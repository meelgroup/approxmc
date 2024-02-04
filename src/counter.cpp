/*
 ApproxMC

 Copyright (c) 2019-2020, Mate Soos and Kuldeep S. Meel. All rights reserved
 Copyright (c) 2009-2018, Mate Soos. All rights reserved.
 Copyright (c) 2015, Supratik Chakraborty, Daniel J. Fremont,
 Kuldeep S. Meel, Sanjit A. Seshia, Moshe Y. Vardi
 Copyright (c) 2014, Supratik Chakraborty, Kuldeep S. Meel, Moshe Y. Vardi

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 */

#include <ctime>
#include <cerrno>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <set>
#include <utility>
#include <fstream>
#include <sys/stat.h>
#include <cstring>
#include <list>
#include <array>
#include <cmath>
#include <complex>

#include "counter.h"
#include "time_mem.h"
#include "GitSHA1.h"
#ifdef CMS_LOCAL_BUILD
#include "arjun.h"
#else
#include <arjun/arjun.h>
#endif

#define verb_print(a, x) \
    do { if (conf.verb >= a) {std::cout << "c " << x << std::endl;} } while (0)

using std::cout;
using std::endl;
using std::map;
using std::make_pair;
using namespace AppMCInt;

bool Counter::solver_add_clause(const vector<Lit>& cl) {
    if (conf.dump_intermediary_cnf) cls_in_solver.push_back(cl);
    return solver->add_clause(cl);
}


bool Counter::solver_add_xor_clause(const vector<Lit>& lits, const bool rhs) {
    if (conf.dump_intermediary_cnf) xors_in_solver.push_back(make_pair(lits, rhs));
    return solver->add_xor_clause(lits, rhs);
}

bool Counter::solver_add_xor_clause(const vector<uint32_t>& vars, const bool rhs) {
    vector<Lit> lits(vars.size());
    for(size_t i = 0; i < vars.size(); i++) lits[i] = Lit(vars[i], false);
    if (conf.dump_intermediary_cnf) xors_in_solver.push_back(make_pair(lits, rhs));
    return solver->add_xor_clause(vars, rhs);
}

Hash Counter::add_hash(uint32_t hash_index, SparseData& sparse_data)
{
    const string random_bits = gen_rnd_bits(conf.sampling_set.size(), hash_index, sparse_data);

    vector<uint32_t> vars;
    for (uint32_t j = 0; j < conf.sampling_set.size(); j++) {
        if (random_bits[j] == '1') vars.push_back(conf.sampling_set[j]);
    }

    solver->new_var();
    const uint32_t act_var = solver->nVars()-1;
    const bool rhs = gen_rhs();
    vars.push_back(act_var);
    solver_add_xor_clause(vars, rhs);
    if (conf.verb_cls) print_xor(vars, rhs);

    return Hash {act_var, vars, rhs};
}

void Counter::ban_one(const uint32_t act_var, const vector<lbool>& model)
{
    vector<Lit> lits;
    lits.push_back(Lit(act_var, false));
    for (const uint32_t var: conf.sampling_set) lits.push_back(Lit(var, model[var] == l_True));
    solver_add_clause(lits);
}

///adding banning clauses for repeating solutions
uint64_t Counter::add_glob_banning_cls(
    const HashesModels* hm
    , const uint32_t act_var
    , const uint32_t num_hashes)
{
    uint64_t repeat = 0;
    uint64_t checked = 0;

    if (hm != nullptr) {
        assert(act_var != std::numeric_limits<uint32_t>::max());
        assert(num_hashes != std::numeric_limits<uint32_t>::max());

        for (uint32_t i = 0; i < hm->glob_model.size(); i++) {
            const SavedModel& sm = hm->glob_model[i];
            //Model was generated with 'sm.hash_num' active
            //We will have 'num_hashes' hashes active

            if (sm.hash_num >= num_hashes) {
                ban_one(act_var, sm.model);
                repeat++;
            } else {
                //Model has to fit all hashes
                checked++;
                bool ok = true;
                for(const auto& h: hm->hashes) {
                    //This hash is number: h.first
                    //Only has to match hashes below current need
                    //note that "h.first" is numbered from 0, so this is a "<" not "<="
                    if (h.first < num_hashes) {
                        ok &= check_model_against_hash(h.second, sm.model);
                        if (!ok) break;
                    }
                }
                if (ok) {
                    //cout << "Found repeat model, had to check " << checked << " hashes" << endl;
                    ban_one(act_var, sm.model);
                    repeat++;
                }
            }
        }
    }
    if (conf.verb) {
        cout << "c [appmc] repeat solutions: " << std::setw(6) << repeat
        << " checked: " << std::setw(6) << checked;
        if (hm) cout << " out of: " << std::setw(6) << hm->glob_model.size();
        cout << endl;
    }
    return repeat;
}

void Counter::dump_cnf_from_solver(const vector<Lit>& assumps, const uint32_t iter, const lbool result) {
    std::string result_str;
    if (result == l_True) result_str = "SAT";
    else if (result == l_False) result_str = "UNSAT";
    else assert(false && "Should not be called with unknown!");

    std::stringstream ss;
    ss << "cnfdump" << "-res-" << result_str << "-iter-" << iter
        << "-active-xors-" << assumps.size() << "-out-" << cnf_dump_no++ << ".cnf";

    std::ofstream f;
    f.open(ss.str(), std::ios::out);
    f << "p cnf " << solver->nVars()+1 << " " << cls_in_solver.size()+xors_in_solver.size()+assumps.size() << endl;
    for(const auto& cl: cls_in_solver) f << cl << " 0" << endl;
    f << "c XORs below" << endl;
    for(const auto& x: xors_in_solver) {
        if (x.first.empty() && x.second == false) continue; // empty && false == tautology
        f << "x ";
        for(uint32_t i = 0; i < x.first.size(); i++) {
            Lit l = x.first[i];
            if (i == 0) l ^= !x.second;
            f << l << " ";
        }
        f << "0" << endl;
    }
    f << "c assumptions below" << endl;
    for(const auto& l: assumps) f << l << " 0" << endl;
    f.close();
}

SolNum Counter::bounded_sol_count(
        uint32_t max_solutions,
        const vector<Lit>* assumps,
        const uint32_t hash_cnt,
        const uint32_t iter,
        HashesModels* hm
) {
    verb_print(1, "[appmc] "
        "[ " << std::setw(7) << std::setprecision(2) << std::fixed << (cpuTimeTotal()-start_time) << " ]"
        << " bounded_sol_count looking for " << std::setw(4) << max_solutions << " solutions"
        << " -- hashes active: " << hash_cnt);

    //Set up things for adding clauses that can later be removed
    vector<Lit> new_assumps;
    if (assumps) {
        assert(assumps->size() == hash_cnt);
        new_assumps = *assumps;
    } else assert(hash_cnt == 0);
    solver->new_var();
    const uint32_t sol_ban_var = solver->nVars()-1;
    new_assumps.push_back(Lit(sol_ban_var, true));

    if (conf.simplify >= 2) {
        verb_print(2, "[appmc] inter-simplifying");
        double my_time = cpuTime();
        solver->simplify(&new_assumps);
        total_inter_simp_time += cpuTime() - my_time;
        verb_print(1, "[appmc] inter-simp finished, total simp time: " << total_inter_simp_time);
    }

    cnf_dump_no = 0;
    const uint64_t repeat = add_glob_banning_cls(hm, sol_ban_var, hash_cnt);
    uint64_t solutions = repeat;
    double last_found_time = cpuTimeTotal();
    vector<vector<lbool>> models;
    while (solutions < max_solutions) {
        lbool ret = solver->solve(&new_assumps, !conf.force_sol_extension);
        assert(ret == l_False || ret == l_True);
        if ((conf.dump_intermediary_cnf >= 2 && ret == l_True) ||
            (conf.dump_intermediary_cnf >= 1 && ret == l_False)) {
            dump_cnf_from_solver(new_assumps, iter, ret);
        }

        if (conf.verb >= 2) {
            cout << "c [appmc] bounded_sol_count ret: " << std::setw(7) << ret;
            if (ret == l_True) cout << " sol no.  " << std::setw(3) << solutions;
            else cout << " No more. " << std::setw(3) << "";
            cout << " T: "
            << std::setw(7) << std::setprecision(2) << std::fixed
            << (cpuTimeTotal()-start_time)
            << " -- hashes act: " << hash_cnt
            << " -- T since last: "
            << std::setw(7) << std::setprecision(2) << std::fixed
            << (cpuTimeTotal()-last_found_time) << endl;
            if (conf.verb >= 4) solver->print_stats();
            last_found_time = cpuTimeTotal();
        }
        if (ret != l_True) break;

        //Add solution to set
        solutions++;
        const vector<lbool> model = solver->get_model();
        check_model(model, hm, hash_cnt);
        models.push_back(model);

        //ban solution
        vector<Lit> lits;
        lits.push_back(Lit(sol_ban_var, false));
        for (const uint32_t var: conf.sampling_set) {
            assert(solver->get_model()[var] != l_Undef);
            lits.push_back(Lit(var, solver->get_model()[var] == l_True));
        }
        if (conf.verb_cls) {
            cout << "c [appmc] Adding banning clause: " << lits << endl;
        }
        solver_add_clause(lits);
    }

    //Save global models
    if (hm && conf.reuse_models) {
        for (const auto& model: models) {
            hm->glob_model.push_back(SavedModel {model, hash_cnt});
        }
    }

    //Remove solution banning
    vector<Lit> cl_that_removes;
    cl_that_removes.push_back(Lit(sol_ban_var, false));
    solver_add_clause(cl_that_removes);

    return SolNum(solutions, repeat);
}

void Counter::print_final_count_stats(ApproxMC::SolCount sol_count)
{
    if (sol_count.hashCount == 0 && sol_count.cellSolCount == 0) {
        cout << "c [appmc] Formula was UNSAT " << endl;
    }

    if (conf.verb > 2) solver->print_stats();
}

ApproxMC::SolCount Counter::solve() {
    orig_num_vars = solver->nVars();
    start_time = cpuTimeTotal();

    open_logfile();
    rnd_engine.seed(conf.seed);

    ApproxMC::SolCount sol_count = count();
    print_final_count_stats(sol_count);

    if (conf.verb) {
        cout << "c [appmc] ApproxMC T: "
        << (cpuTimeTotal() - start_time) << " s"
        << endl;
    }
    return sol_count;
}

vector<Lit> Counter::set_num_hashes(
    uint32_t num_wanted,
    map<uint64_t, Hash>& hashes,
    SparseData& sparse_data
) {
    vector<Lit> assumps;
    for(uint32_t i = 0; i < num_wanted; i++) {
        if (hashes.find(i) != hashes.end()) {
            assumps.push_back(Lit(hashes[i].act_var, true));
        } else {
            Hash h = add_hash(i, sparse_data);
            assumps.push_back(Lit(h.act_var, true));
            hashes[i] = h;
        }
    }
    assert(num_wanted == assumps.size());

    return assumps;
}

void Counter::simplify()
{
    if (conf.verb >= 1) {
        cout << "c [appmc] simplifying" << endl;
    }

    solver->set_sls(1);
    solver->set_intree_probe(1);
    solver->set_full_bve_iter_ratio(conf.var_elim_ratio);
    solver->set_full_bve(1);
    solver->set_scc(1);

    solver->simplify();

    solver->set_sls(0);
    solver->set_full_bve(0);
}

//Set up probabilities, threshold and measurements
void Counter::set_up_probs_threshold_measurements(
    uint32_t& measurements, SparseData& sparse_data)
{
    int best_match = -1;
    bool using_sparse = false;
    double thresh_factor;

    if (conf.sparse) {
        best_match = find_best_sparse_match();
    }

    if (best_match != -1) {
        sparse_data = SparseData(best_match);
        thresh_factor = 1.1;
        using_sparse = true;
    } else {
        thresh_factor = 1.0;
    }

    threshold = int(
        1 +
        thresh_factor*
        9.84*
        (1.0+(1.0/conf.epsilon))*
        (1.0+(1.0/conf.epsilon))*
        (1.0+(conf.epsilon/(1.0+conf.epsilon)))
    );

    verb_print(1, "[appmc] threshold set to " << threshold << " sparse: " << (int)using_sparse);
    measurements = (int)std::ceil(std::log2(3.0/conf.delta)*17);
    for (int count = 0; count < 256; count++) {
        if (constants.iterationConfidences[count] >= 1 - conf.delta) {
            measurements = count*2+1;
            break;
        }
    }
}

bool Counter::find_one_solution()
{
    auto ret = solver->solve();
    return ret == l_True;
}

ApproxMC::SolCount Counter::count()
{
    const int64_t hash_cnt = conf.start_iter;

    SparseData sparse_data(-1);
    HashesModels hm;
    uint32_t measurements;
    set_up_probs_threshold_measurements(measurements, sparse_data);

    verb_print(1, "[appmc] Starting at hash count: " << hash_cnt);

    int64_t prev_measure = hash_cnt;
    num_hash_list.clear();
    num_count_list.clear();

    //See Algorithm 1 in paper "Algorithmic Improvements in Approximate Counting
    //for Probabilistic Inference: From Linear to Logarithmic SAT Calls"
    //https://www.ijcai.org/Proceedings/16/Papers/503.pdf
    for (uint32_t j = 0; j < measurements; j++) {
        one_measurement_count(prev_measure, j, sparse_data, &hm);
        if (prev_measure == 0) {
            // Exact count, no need to measure multiple times.
            break;
        }
        sparse_data.next_index = 0;
        if (conf.simplify >= 1 && j+1 < measurements) simplify();
        hm.clear();
    }
    assert(!num_hash_list.empty() && "UNSAT should not be possible");

    return calc_est_count();
}

ApproxMC::SolCount Counter::calc_est_count()
{
    ApproxMC::SolCount ret_count;
    if (num_hash_list.empty() || num_count_list.empty()) return ret_count;

    const auto min_hash = find_min(num_hash_list);
    auto cnt_it = num_count_list.begin();
    for (auto hash_it = num_hash_list.begin()
        ; hash_it != num_hash_list.end() && cnt_it != num_count_list.end()
        ; hash_it++, cnt_it++
    ) {
        if ((*hash_it) - min_hash > 10) {
            cout << "Internal ERROR: Something is VERY fishy, the difference between each count must"
                " never be this large. Please report this bug to the maintainers" << endl;
            exit(-1);
        }
        *cnt_it *= pow(2, (*hash_it) - min_hash);
    }
    ret_count.valid = true;
    ret_count.cellSolCount = find_median(num_count_list);
    ret_count.hashCount = min_hash;

    return ret_count;
}

int Counter::find_best_sparse_match()
{
    for(int i = 0; i < (int)constants.index_var_maps.size(); i++) {
        if (constants.index_var_maps[i].vars_to_inclusive >= conf.sampling_set.size()) {
            if (conf.verb) {
                cout << "c [sparse] Using match: " << i
                << " sampling set size: " << conf.sampling_set.size()
                << " prev end inclusive is: " << (i == 0 ? -1 : (int)constants.index_var_maps[i-1].vars_to_inclusive)
                << " this end inclusive is: " << constants.index_var_maps[i].vars_to_inclusive
                << " next end inclusive is: " << ((i+1 < (int)constants.index_var_maps.size()) ? ((int)constants.index_var_maps[i+1].vars_to_inclusive) : -1)
                << " sampl size: " << conf.sampling_set.size()
                << endl;
            }

            return i;
        }
    }

    cout << "c [sparse] No match. Using default 0.5" << endl;
    return -1;
}

//See Algorithm 2+3 in paper "Algorithmic Improvements in Approximate Counting
//for Probabilistic Inference: From Linear to Logarithmic SAT Calls"
//https://www.ijcai.org/Proceedings/16/Papers/503.pdf
void Counter::one_measurement_count(
    int64_t& prev_measure,
    const uint32_t iter,
    SparseData sparse_data,
    HashesModels* hm)
{
    if (conf.sampling_set.empty()) {
        num_hash_list.push_back(0);
        num_count_list.push_back(1);
        return;
    }

    //Tells the number of solutions found at hash number N
    //sols_for_hash[N] tells the number of solutions found when N hashes were added
    map<uint64_t,int64_t> sols_for_hash;

    //threshold_sols[hash_num]==1 tells us that at hash_num number of hashes
    //there were found to be FULL threshold number of solutions
    //threshold_sols[hash_num]==0 tells that there were less than threshold
    //number of solutions.
    //if it's not set, we have no clue.
    map<uint64_t,bool> threshold_sols;
    int64_t total_max_xors = conf.sampling_set.size();
    int64_t num_explored = 0;
    int64_t lower_fib = 0;
    int64_t upper_fib = total_max_xors;
    int64_t hash_cnt = prev_measure;
    int64_t hash_prev = hash_cnt;

    //We are doing a galloping search here (see our IJCAI-16 paper for more details).
    //lowerFib is referred to as loIndex and upperFib is referred to as hiIndex
    //The key idea is that we first do an exponential search and then do binary search
    //This is implemented by using two sentinels: lowerFib and upperFib. The correct answer
    // is always between lowFib and upperFib. We do exponential search until upperFib < lowerFib/2
    // Once upperFib < lowerFib/2; we do a binary search.
    while (num_explored < total_max_xors) {
        uint64_t cur_hash_cnt = hash_cnt;
        const vector<Lit> assumps = set_num_hashes(hash_cnt, hm->hashes, sparse_data);

        verb_print(1, "[appmc] "
            "[ " << std::setw(7) << std::setprecision(2) << std::fixed << (cpuTimeTotal()-start_time) << " ]"
            << " round: " << std::setw(2) << iter
            << " hashes: " << std::setw(6) << hash_cnt);
        double my_time = cpuTime();
        SolNum sols = bounded_sol_count(
            threshold + 1, //max no. solutions
            &assumps, //assumptions to use
            hash_cnt,
            iter,
            hm
        );
        const uint64_t num_sols = std::min<uint64_t>(sols.solutions, threshold + 1);
        assert(num_sols <= threshold + 1);
        bool found_full = (num_sols == threshold + 1);
        write_log(
            false, //not sampling
            iter, hash_cnt, found_full, num_sols, sols.repeated,
            cpuTime() - my_time
        );

        if (num_sols < threshold + 1) {
            num_explored = lower_fib + total_max_xors - hash_cnt;

            //one less hash count had threshold solutions
            //this one has less than threshold
            //so this is the real deal!
            if (hash_cnt == 0 ||
                    (threshold_sols.find(hash_cnt-1) != threshold_sols.end()
                    && threshold_sols[hash_cnt-1] == 1)
            ) {
                num_hash_list.push_back(hash_cnt);
                num_count_list.push_back(num_sols);
                prev_measure = hash_cnt;
                return;
            }

            threshold_sols[hash_cnt] = 0;
            sols_for_hash[hash_cnt] = num_sols;
            if (iter > 0 &&
                std::abs(hash_cnt - prev_measure) <= 2
            ) {
                //Doing linear, this is a re-count
                upper_fib = hash_cnt;
                hash_cnt--;
            } else {
                if (hash_prev > hash_cnt) hash_prev = 0;
                upper_fib = hash_cnt;
                if (hash_prev > lower_fib) lower_fib = hash_prev;
                hash_cnt = (upper_fib+lower_fib)/2;
            }
        } else {
            assert(num_sols == threshold + 1);
            num_explored = hash_cnt + total_max_xors - upper_fib;

            //success record for +1 hashcount exists and is 0
            //so one-above hashcount was below threshold, this is above
            //we have a winner -- the one above!
            if (threshold_sols.find(hash_cnt+1) != threshold_sols.end()
                && threshold_sols[hash_cnt+1] == 0
            ) {
                num_hash_list.push_back(hash_cnt+1);
                num_count_list.push_back(sols_for_hash[hash_cnt+1]);
                prev_measure = hash_cnt+1;
                return;
            }

            threshold_sols[hash_cnt] = 1;
            sols_for_hash[hash_cnt] = threshold+1;
            if (iter > 0
                && std::abs(hash_cnt - prev_measure) < 2
            ) {
                //Doing linear, this is a re-count
                lower_fib = hash_cnt;
                hash_cnt++;
            } else if (lower_fib + (hash_cnt-lower_fib)*2 >= upper_fib-1) {

                // Whenever the above condition is satisfied, we are in binary search mode
                lower_fib = hash_cnt;
                hash_cnt = (lower_fib+upper_fib)/2;
            } else {
                // We are in exponential search mode.
                const auto old_hash_cnt = hash_cnt;
                hash_cnt = lower_fib + (hash_cnt-lower_fib)*2;
                if (old_hash_cnt == hash_cnt) hash_cnt++;
            }
        }
        hash_prev = cur_hash_cnt;
    }
}
bool Counter::gen_rhs()
{
    std::uniform_int_distribution<uint32_t> dist{0, 1};
    bool rhs = dist(rnd_engine);
    return rhs;
}

string Counter::gen_rnd_bits(
    const uint32_t size,
    // The name of parameter was changed to indicate that this is the index of hash function
    const uint32_t hash_index,
    SparseData& sparse_data)
{
    string random_bits;
    std::uniform_int_distribution<uint32_t> dist{0, 999};
    uint32_t cutoff = 500;
    if (conf.sparse && sparse_data.table_no != -1) {
        //Do we need to update the probability?
        const auto& table = constants.index_var_maps[sparse_data.table_no];
        const auto next_var_index = table.index_var_map[sparse_data.next_index];
        if (hash_index >= next_var_index) {
            sparse_data.sparseprob = constants.probval[sparse_data.next_index];
            sparse_data.next_index = std::min<uint32_t>(
                sparse_data.next_index+1, table.index_var_map.size()-1);
        }
        assert(sparse_data.sparseprob <= 0.5);
        cutoff = std::ceil(1000.0*sparse_data.sparseprob);
        if (conf.verb > 3) {
            cout << "c [sparse] cutoff: " << cutoff
            << " table: " << sparse_data.table_no
            << " lookup index: " << sparse_data.next_index
            << " hash index: " << hash_index
            << endl;
        }
    }

    while (random_bits.size() < size) {
        bool val = dist(rnd_engine) < cutoff;
        random_bits += '0' + val;
    }
    assert(random_bits.size() >= size);
    return random_bits;
}

void Counter::print_xor(const vector<uint32_t>& vars, const uint32_t rhs)
{
    cout << "c [appmc] Added XOR ";
    for (size_t i = 0; i < vars.size(); i++) {
        cout << vars[i]+1;
        if (i < vars.size()-1) {
            cout << " + ";
        }
    }
    cout << " = " << (rhs ? "True" : "False") << endl;
}

template<class T> inline T Counter::find_median(const vector<T>& nums) {
    assert(!nums.empty());
    auto tmp = nums;

    std::sort(tmp.begin(), tmp.end());
    size_t med_index = tmp.size() / 2;
    if (med_index >= tmp.size()) return tmp[tmp.size() - 1];
    return tmp[med_index];
}

template<class T> inline T Counter::find_min(const vector<T>& nums) {
    T min = std::numeric_limits<T>::max();
    for (const auto a: nums) {
        if (a < min) min = a;
    }
    return min;
}

string scalmc_version_info()
{
    std::stringstream ss;
    ss << "c ApproxMC SHA revision " << AppMCInt::get_version_sha1() << endl;
    ss << "c ApproxMC version " << AppMCInt::get_version_tag() << endl;
    ss << "c ApproxMC compilation env " << AppMCInt::get_compilation_env() << endl;
    #ifdef __GNUC__
    ss << "c ApproxMC compiled with gcc version " << __VERSION__ << endl;
    #else
    ss << "c ApproxMC compiled with non-gcc compiler" << endl;
    #endif

    return ss.str();
}

string Counter::get_version_info() const
{
    string ret = ::scalmc_version_info();
    ret += solver->get_text_version_info();

    return ret;
}

void Counter::open_logfile()
{
    if (!conf.logfilename.empty()) {
        logfile.open(conf.logfilename.c_str());
        if (!logfile.is_open()) {
            cout << "[appmc] Cannot open Counter log file '" << conf.logfilename
                 << "' for writing." << endl;
            exit(1);
        }

        logfile << std::left
        << std::setw(5) << "sampl"
        << " " << std::setw(4) << "iter"
        << " " << std::setw(4) << "hash"
        << " " << std::setw(4) << "full"
        << " " << std::setw(4) << "sols"
        << " " << std::setw(4) << "rep"
        << " " << std::setw(7) << "T"
        << " " << std::setw(7) << "total T"
        << endl;

    }
}

void Counter::write_log(
    bool sampling,
    int iter,
    uint32_t hash_cnt,
    int found_full,
    uint32_t num_sols,
    uint32_t repeat_sols,
    double used_time
) {
    if (!conf.logfilename.empty()) {
        logfile
        << std::left
        << std::setw(5) << (int)sampling
        << " " << std::setw(4) << iter
        << " " << std::setw(4) << hash_cnt
        << " " << std::setw(4) << found_full
        << " " << std::setw(4) << num_sols
        << " " << std::setw(4) << repeat_sols
        << " " << std::setw(7) << std::fixed << std::setprecision(2) << used_time
        << " " << std::setw(7) << std::fixed << std::setprecision(2) << (cpuTimeTotal() - start_time)
        << endl;
    }
}

void Counter::check_model(
    const vector<lbool>& model,
    const HashesModels* const hm,
    const uint32_t hash_cnt
) {
    for(uint32_t var: conf.sampling_set) assert(model[var] != l_Undef);
    if (conf.debug) {
        assert(conf.force_sol_extension);
        assert(conf.dump_intermediary_cnf);
        for(const auto& cl: cls_in_solver) {
            bool sat = false;
            for(const auto& l: cl) {
                assert(model[l.var()] != l_Undef);
                if ((model[l.var()] == l_True && !l.sign()) ||
                    (model[l.var()] == l_False && l.sign())) {sat = true; break;}
            }
            assert(sat);
        }
        for(const auto& x: xors_in_solver) {
            bool sat = !x.second;
            for(const auto& l: x.first) {
                assert(model[l.var()] != l_Undef);
                sat ^= ((model[l.var()]^l.sign()) == l_True);
            }
            assert(sat);
        }
    }

    if (!hm) return;

    for(const auto& h: hm->hashes) {
        //This hash is number: h.first
        //Only has to match hashes at & below
        //Notice that "h.first" is numbered from 0, so it's a "<" not "<="
        if (h.first < hash_cnt) {
            //cout << "Checking model against hash " << h.first << endl;
            assert(check_model_against_hash(h.second, model));
        }
    }
}

bool Counter::check_model_against_hash(const Hash& h, const vector<lbool>& model) {
    bool rhs = false;
    for (auto const& var: h.hash_vars) {
        assert(model[var] != l_Undef);
        rhs ^= model[var] == l_True;
    }
    return rhs == h.rhs;
}
