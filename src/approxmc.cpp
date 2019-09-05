/*
 ApproxMC

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
#include <cstring>
#include <errno.h>
#include <algorithm>
#include <string.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <set>
#include <fstream>
#include <sys/stat.h>
#include <string.h>
#include <list>
#include <array>
#include <cmath>
#include <complex>

#include "approxmc.h"
#include "time_mem.h"
#include "cryptominisat5/cryptominisat.h"
#include "cryptominisat5/solvertypesmini.h"
#include "GitSHA1.h"

using std::cout;
using std::cerr;
using std::endl;
using std::list;
using std::map;

void AppMC::add_hash(
    uint32_t num_xor_cls,
    vector<Lit>& assumps,
    uint32_t total_num_hashes)
{
    const string randomBits =
        gen_rnd_bits(conf.sampling_set.size() * num_xor_cls, total_num_hashes);

    bool rhs;
    vector<uint32_t> vars;

    for (uint32_t i = 0; i < num_xor_cls; i++) {
        //new activation variable
        solver->new_var();
        uint32_t act_var = solver->nVars()-1;
        assumps.push_back(Lit(act_var, true));

        vars.clear();
        vars.push_back(act_var);
        rhs = gen_rhs();

        for (uint32_t j = 0; j < conf.sampling_set.size(); j++) {
            if (randomBits.at(conf.sampling_set.size() * i + j) == '1') {
                vars.push_back(conf.sampling_set[j]);
            }
        }
        solver->add_xor_clause(vars, rhs);
        if (conf.verb_appmc_cls) {
            print_xor(vars, rhs);
        }
    }
}

///adding banning clauses for repeating solutions
void AppMC::add_glob_banning_cls(
    const vector<vector<lbool>>* glob_model
    , const uint32_t act_var)
{
    if (glob_model == NULL)
        return;

    vector<Lit> lits;
    for (uint32_t i = 0; i < glob_model->size(); i++)
    {
        lits.clear();
        lits.push_back(Lit(act_var, false));
        for (const uint32_t var: conf.sampling_set) {
            lits.push_back(Lit(var, (*glob_model)[i][var] == l_True));
        }
        solver->add_clause(lits);
    }
}

int64_t AppMC::bounded_sol_count(
        uint32_t maxSolutions,
        const vector<Lit>* assumps,
        const uint32_t hashCount,
        uint32_t minSolutions,
        vector<vector<lbool>>* glob_model,
        vector<string>* out_solutions
) {
    cout << "[appmc] "
    "[ " << std::setw(7) << std::setprecision(2) << std::fixed
    << (cpuTimeTotal()-startTime)
    << " ]"
    << " bounded_sol_count looking for " << std::setw(4) << maxSolutions << " solutions"
    << " -- hashes active: " << hashCount << endl;

    //Will we need to extend the solution?
    bool must_extend = true;
    if (out_solutions != NULL) {
        must_extend = false;
    } else {
        must_extend = !conf.only_indep_samples;
    }

    //Set up things for adding clauses that can later be removed
    vector<Lit> new_assumps;
    if (assumps) {
        new_assumps = *assumps;
    }
    solver->new_var();
    const uint32_t sol_ban_var = solver->nVars()-1;
    new_assumps.push_back(Lit(sol_ban_var, true));
    add_glob_banning_cls(glob_model, sol_ban_var);

    uint64_t solutions = 0;
    double last_found_time = cpuTimeTotal();
    vector<vector<lbool>> models;
    while (solutions < maxSolutions) {
        lbool ret = solver->solve(&new_assumps, must_extend);
        assert(ret == l_False || ret == l_True);

        if (conf.verb >= 2) {
            cout << "[appmc] bounded_sol_count ret: " << std::setw(7) << ret;
            if (ret == l_True) {
                cout << " sol no.  " << std::setw(3) << solutions;
            } else {
                cout << " No more. " << std::setw(3) << "";
            }
            cout << " T: "
            << std::setw(7) << std::setprecision(2) << std::fixed
            << (cpuTimeTotal()-startTime)
            << " -- hashes act: " << hashCount
            << " -- T since last: "
            << std::setw(7) << std::setprecision(2) << std::fixed
            << (cpuTimeTotal()-last_found_time)
            << endl;
            if (conf.verb >= 3) {
                solver->print_stats();
            }
            last_found_time = cpuTimeTotal();
        }

        if (ret != l_True) {
            break;
        }

        //Add solution to set
        solutions++;
        vector<lbool> model = solver->get_model();
        models.push_back(model);
        if (out_solutions) {
            out_solutions->push_back(get_solution_str(model));
        }

        //ban solution
        vector<Lit> lits;
        lits.push_back(Lit(sol_ban_var, false));
        for (const uint32_t var: conf.sampling_set) {
            assert(solver->get_model()[var] != l_Undef);
            lits.push_back(Lit(var, solver->get_model()[var] == l_True));
        }
        if (conf.verb_appmc_cls) {
            cout << "[appmc] Adding banning clause: " << lits << endl;
        }
        solver->add_clause(lits);
    }

    if (solutions < maxSolutions) {
        //Sampling -- output a random sample of N solutions
        if (solutions >= minSolutions && samples_out != NULL) {
            assert(minSolutions > 0);
            vector<size_t> modelIndices;
            for (uint32_t i = 0; i < models.size(); i++) {
                modelIndices.push_back(i);
            }
            std::shuffle(modelIndices.begin(), modelIndices.end(), randomEngine);

            for (uint32_t i = 0; i < sols_to_return(solutions); i++) {
                const auto& model = models.at(modelIndices.at(i));
                (*samples_out) << get_solution_str(model) << endl << std::flush;
            }
        }

        //Save global models
        if (glob_model) {
            for (const auto& model: models) {
                glob_model->push_back(model);
            }
        }
    }

    //Remove solution banning
    vector<Lit> cl_that_removes;
    cl_that_removes.push_back(Lit(sol_ban_var, false));
    solver->add_clause(cl_that_removes);

    return solutions;
}

int AppMC::solve(AppMCConfig _conf)
{
    conf = _conf;
    orig_num_vars = solver->nVars();
    startTime = cpuTimeTotal();
    sampling = false;

    openLogFile();
    randomEngine.seed(conf.seed);
    if (conf.samples == 0) {
        cout << "[appmc] Using start iteration " << conf.startiter << endl;

        SATCount solCount;
        count(solCount);

        cout << "[appmc] FINISHED AppMC T: "
        << (cpuTimeTotal() - startTime) << " s"
        << endl;

        if (solCount.hashCount == 0 && solCount.cellSolCount == 0) {
            cout << "[appmc] Formula was UNSAT " << endl;
        }

        if (conf.verb > 2) {
            solver->print_stats();
        }

        solCount.print_num_solutions();
    } else {
        if (conf.startiter > conf.sampling_set.size()) {
            cerr << "ERROR: Manually-specified startiter for gen_n_samples"
                 "is larger than the size of the independent set.\n" << endl;
            exit(-1);
        }

        /* Compute threshold via formula from TACAS-15 paper */
        threshold_appmcgen = ceil(4.03 * (1 + (1/conf.kappa)) * (1 + (1/conf.kappa)));

        //No startiter, we have to figure it out
        if (conf.startiter == 0) {
            SATCount solCount;
            std::ostream* backup = samples_out;
            samples_out = NULL;

            count(solCount);
            cout << "[appmc] finished counting solutions in "
            << (cpuTimeTotal() - startTime) << " s" << endl;

            if (solCount.hashCount == 0 && solCount.cellSolCount == 0) {
                cout << "[appmc] The input formula is unsatisfiable." << endl;
                return correctReturnValue(l_False);
            }

            if (conf.verb) {
                solver->print_stats();
            }
            solCount.print_num_solutions();

            double si = round(solCount.hashCount + log2(solCount.cellSolCount)
                + log2(1.8) - log2(threshold_appmcgen)) - 2;
            if (si > 0) {
                conf.startiter = si;
            } else {
                conf.startiter = 0;   /* Indicate ideal sampling case */
            }
            samples_out = backup;
        } else {
            cout << "Using manually-specified startiter for sample generation" << endl;
        }
        generate_samples();
    }

    return correctReturnValue(l_True);
}

void AppMC::set_num_hashes(
    uint32_t num_wanted,
    map<uint64_t,Lit>& hashVars,
    vector<Lit>& assumps
) {
    if (num_wanted < assumps.size()) {
        uint64_t numberToRemove = assumps.size()- num_wanted;
        for (uint64_t i = 0; i<numberToRemove; i++) {
            assumps.pop_back();
        }
    } else {
        if (num_wanted > assumps.size() && assumps.size() < hashVars.size()) {
            for (uint32_t i = assumps.size(); i< hashVars.size() && i < num_wanted; i++) {
                assumps.push_back(hashVars[i]);
            }
        }
        if (num_wanted > hashVars.size()) {
            add_hash(num_wanted-hashVars.size(), assumps, num_wanted);
            for (uint64_t i = hashVars.size(); i < num_wanted; i++) {
                hashVars[i] = assumps[i];
            }
        }
}
}

void AppMC::count(SATCount& ret_count)
{
    ret_count.clear();
    int64_t hashCount = conf.startiter;

    cout << "[appmc] Starting up, initial measurement" << endl;
    if (hashCount == 0) {
        cout << "[appmc] Checking if there are at least threshold+1 solutions..." << endl;
        double myTime = cpuTime();
        int64_t currentNumSolutions = bounded_sol_count(
            conf.threshold+1, //max solutions
            NULL, // no assumptions
            hashCount
        );
        write_log(0, 0,
                  currentNumSolutions == (conf.threshold + 1),
                  currentNumSolutions, 0, cpuTime() - myTime);

        //Din't find at least threshold+1
        if (currentNumSolutions <= conf.threshold) {
            cout << "[appmc] Did not find at least threshold+1 ("
            << conf.threshold << ") we found only " << currentNumSolutions
            << ", i.e. we got exact count" << endl;

            ret_count.cellSolCount = currentNumSolutions;
            ret_count.hashCount = 0;
            return;
        }
        hashCount++;
    }
    cout << "[appmc] Starting at hash count: " << hashCount << endl;

    vector<uint64_t> numHashList;
    vector<int64_t> numCountList;
    int64_t mPrev = hashCount;
    if (false) {
        uint32_t threshold_prev = conf.threshold;
        conf.threshold = 1;
        one_measurement_count(
            numHashList
            , numCountList
            , mPrev
            , -1
        );
        numHashList.clear();
        numCountList.clear();
        conf.threshold = threshold_prev;
        mPrev -= log2(conf.threshold);
        mPrev = std::max<int64_t>(0, mPrev);
        mPrev++;
    }

    //See Algorithm 1 in paper "Algorithmic Improvements in Approximate Counting
    //for Probabilistic Inference: From Linear to Logarithmic SAT Calls"
    //https://www.ijcai.org/Proceedings/16/Papers/503.pdf
    for (uint32_t j = 0; j < conf.measurements; j++) {
        one_measurement_count(
            numHashList
            , numCountList
            , mPrev
            , j
        );
        solver->simplify();
    }
    assert(numHashList.size() > 0 && "UNSAT should not be possible");

    //Median
    auto minHash = findMin(numHashList);
    auto cnt_it = numCountList.begin();
    for (auto hash_it = numHashList.begin()
        ; hash_it != numHashList.end() && cnt_it != numCountList.end()
        ; hash_it++, cnt_it++
    ) {
        *cnt_it *= pow(2, (*hash_it) - minHash);
    }
    int medSolCount = findMedian(numCountList);
    ret_count.cellSolCount = medSolCount;
    ret_count.hashCount = minHash;
    return;
}

//See Algorithm 2+3 in paper "Algorithmic Improvements in Approximate Counting
//for Probabilistic Inference: From Linear to Logarithmic SAT Calls"
//https://www.ijcai.org/Proceedings/16/Papers/503.pdf
void AppMC::one_measurement_count(
    vector<uint64_t>& numHashList,
    vector<int64_t>& numCountList,
    int64_t& mPrev,
    const int iter
)
{
    vector<Lit> assumps;
    map<uint64_t, Lit> hashVars;

    //Tells the number of solutions found at hash number N
    //sols_for_hash[N] tells the number of solutions found when N hashes were added
    map<uint64_t,int64_t> sols_for_hash;

    //threshold_sols[hash_num]==1 tells us that at hash_num number of hashes
    //there were found to be FULL threshold number of solutions
    //threshold_sols[hash_num]==0 tells that there were less than threshold
    //number of solutions.
    //if it's not set, we have no clue.
    map<uint64_t,bool> threshold_sols;

    vector<vector<lbool>> glob_model; //global table storing models
    int64_t repeat = 0; //repeatable solutions

    int64_t total_max_xors = conf.sampling_set.size();
    int64_t numExplored = 0;
    int64_t lowerFib = 0;
    int64_t upperFib = total_max_xors;

    int64_t hashCount = mPrev;
    int64_t hashPrev = hashCount;
    while (numExplored < total_max_xors) {
        uint64_t cur_hash_count = hashCount;
        set_num_hashes(hashCount, hashVars, assumps);
        assert(conf.threshold + 1 >= repeat);

        cout << "[appmc] hashes active: " << std::setw(6) << hashCount << endl;
        double myTime = cpuTime();
        int64_t num_sols = bounded_sol_count(
            conf.threshold + 1 - repeat, //max no. solutions
            &assumps, //assumptions to use
            hashCount,
            1, //ignored
            &glob_model //all solutions put here
        );
        assert(num_sols <= conf.threshold + 1 - repeat);
        bool found_full = (num_sols == conf.threshold + 1 - repeat);
        write_log(iter, hashCount, found_full, num_sols + repeat, repeat,
            cpuTime() - myTime
        );

        if (num_sols < conf.threshold + 1 - repeat) {
            numExplored = lowerFib + total_max_xors - hashCount;

            //one less hash count had threshold solutions
            //this one has less than threshold
            //so this is the real deal!
            if (threshold_sols.find(hashCount-1) != threshold_sols.end()
                && threshold_sols[hashCount-1] == 1
            ) {
                numHashList.push_back(hashCount);
                numCountList.push_back(num_sols+repeat);
                mPrev = hashCount;
                return;
            }

            threshold_sols[hashCount] = 0;
            sols_for_hash[hashCount] = repeat + num_sols;
            if (iter > 0 &&
                std::abs(hashCount - mPrev) <= 2
            ) {
                //Doing linear, this is a re-count
                upperFib = hashCount;
                hashCount--;
            } else {
                if (hashPrev > hashCount) {
                    hashPrev = 0;
                }
                upperFib = hashCount;
                if (hashPrev > lowerFib) {
                    lowerFib = hashPrev;
                }

                //Fast hit
                if (false) {
                    //Trying to hit the right place in case
                    //we got some solutions here -- calculate the right place
                    int64_t diff_delta = 0;
                    if (repeat+num_sols > 0) {
                        diff_delta = log2(conf.threshold/(repeat+num_sols));
                        if (diff_delta == 0){
                            diff_delta = 1;
                        }
                        hashCount -= diff_delta;
                    } else {
                        hashCount = (upperFib+lowerFib)/2;
                    }
                } else {
                    //Slow hit
                    hashCount = (upperFib+lowerFib)/2;
                }
            }
            repeat += num_sols;
        } else {
            assert(num_sols == conf.threshold + 1 - repeat);
            numExplored = hashCount + total_max_xors - upperFib;

            //success record for +1 hashcount exists and is 0
            //so one-above hashcount was below threshold, this is above
            //we have a winner -- the one above!
            if (threshold_sols.find(hashCount+1) != threshold_sols.end()
                && threshold_sols[hashCount+1] == 0
            ) {
                numHashList.push_back(hashCount+1);
                numCountList.push_back(sols_for_hash[hashCount+1]);
                mPrev = hashCount+1;
                return;
            }

            threshold_sols[hashCount] = 1;
            sols_for_hash[hashCount] = conf.threshold+1;
            if (iter > 0
                && std::abs(hashCount - mPrev) < 2
            ) {
                //Doing linear, this is a re-count
                lowerFib = hashCount;
                hashCount++;
            } else if (lowerFib + (hashCount-lowerFib)*2 >= upperFib-1) {
                lowerFib = hashCount;
                hashCount = (lowerFib+upperFib)/2;
            } else {
                hashCount = lowerFib + (hashCount-lowerFib)*2;
            }
        }
        hashPrev = cur_hash_count;
    }
}

void AppMC::generate_samples()
{
    sampling = true;
    assert(samples_out != NULL);
    double genStartTime = cpuTimeTotal();

    hiThresh = ceil(1 + (1.4142136 * (1 + conf.kappa) * threshold_appmcgen));
    loThresh = floor(threshold_appmcgen / (1.4142136 * (1 + conf.kappa)));
    const uint32_t samplesPerCall = sols_to_return(conf.samples);
    const uint32_t callsNeeded =
        conf.samples / samplesPerCall + (bool)(conf.samples % samplesPerCall);
    cout << "Samples requested: " << conf.samples << endl;
    cout << "samples per XOR set:" << samplesPerCall << endl;
    cout << "-> calls needed: " << callsNeeded << endl;

    //TODO WARNING what is this 14???????????????????
    uint32_t callsPerLoop = std::min(solver->nVars() / 14, callsNeeded);
    callsPerLoop = std::max(callsPerLoop, 1U);
    cout << "callsPerLoop:" << callsPerLoop << endl;

    cout << "[appmc] starting sample generation."
    << " loThresh: " << loThresh
    << ", hiThresh: " << hiThresh
    << ", startiter: " << conf.startiter << endl;

    uint32_t samples = 0;
    if (conf.startiter > 0) {
        uint32_t lastSuccessfulHashOffset = 0;
        while(samples < conf.samples) {
            samples += gen_n_samples(
                callsPerLoop,
                &lastSuccessfulHashOffset);
        }
    } else {
        /* Ideal sampling case; enumerate all solutions */
        vector<string> out_solutions;
        const uint32_t count = bounded_sol_count(
            std::numeric_limits<uint32_t>::max() //max no. solutions
            , NULL //assumps is empty
            , 0 //number of hashes (information only)
            , 1 //min num. solutions
            , NULL //gobal model (would be banned)
            , &out_solutions
        );
        assert(count > 0);

        std::uniform_int_distribution<unsigned> uid {0, count-1};
        for (uint32_t i = 0; i < conf.samples; ++i) {
            vector<string>::iterator it = out_solutions.begin();
            for (uint32_t j = uid(randomEngine); j > 0; --j)    // TODO improve hack
            {
                ++it;
            }
            samples++;
            (*samples_out) << *it << endl << std::flush;
        }
    }

    cout
    << "[appmc]"
    << " Time to sample: " << cpuTimeTotal() - genStartTime << " s"
    << " -- Time count+samples: " << cpuTimeTotal() << " s"
    << endl;

    cout << "[appmc] Samples generated: " << samples << endl;
}

uint32_t AppMC::gen_n_samples(
    const uint32_t num_calls
    , uint32_t* lastSuccessfulHashOffset)
{
    uint32_t num_samples = 0;
    uint32_t i = 0;
    while(i < num_calls) {
        uint32_t hashOffsets[3];
        hashOffsets[0] = *lastSuccessfulHashOffset;

        //Specific values
        if (hashOffsets[0] == 0) { // Starting at q-2; go to q-1 then q
            hashOffsets[1] = 1;
            hashOffsets[2] = 2;
        }
        if (hashOffsets[0] == 2) { // Starting at q; go to q-1 then q-2
            hashOffsets[1] = 1;
            hashOffsets[2] = 0;
        }

        vector<Lit> assumps;
        map<uint64_t,Lit> hashVars;
        bool ok;
        for (uint32_t j = 0; j < 3; j++) {
            uint32_t currentHashOffset = hashOffsets[j];
            uint32_t currentHashCount = currentHashOffset + conf.startiter;
            set_num_hashes(currentHashCount, hashVars, assumps);

            double myTime = cpuTime();
            const uint64_t solutionCount = bounded_sol_count(
                hiThresh // max num solutions
                , &assumps //assumptions to use
                , loThresh //min number of solutions (samples not output otherwise)
            );
            ok = (solutionCount < hiThresh && solutionCount >= loThresh);
            write_log(i, currentHashCount, solutionCount == hiThresh,
                      solutionCount, 0, cpuTime()-myTime);

            if (ok) {
                num_samples += sols_to_return(conf.samples);
                *lastSuccessfulHashOffset = currentHashOffset;
                break;
            }
            // Number of solutions too small or too large

            // At q-1, and need to pick next hash count
            if (j == 0 && currentHashOffset == 1) {
                if (solutionCount < loThresh) {
                    // Go to q-2; next will be q
                    hashOffsets[1] = 0;
                    hashOffsets[2] = 2;
                } else {
                    // Go to q; next will be q-2
                    hashOffsets[1] = 2;
                    hashOffsets[2] = 0;
                }
            }
        }

        if (ok) {
            i++;
        }
        solver->simplify();
    }
    return num_samples;
}

////////////////////
//Helper functions
////////////////////

std::string AppMC::get_solution_str(const vector<lbool>& model)
{
    assert(samples_out != NULL);

    std::stringstream  solution;
    if (conf.only_indep_samples) {
        for (uint32_t j = 0; j < conf.sampling_set.size(); j++) {
            uint32_t var = conf.sampling_set[j];
            assert(model[var] != l_Undef);
            solution << ((model[var] != l_True) ? "-":"") << var + 1 << " ";
        }
    } else {
        for(uint32_t var = 0; var < orig_num_vars; var++) {
            assert(model[var] != l_Undef);
            solution << ((model[var] != l_True) ? "-":"") << var + 1 << " ";
        }
    }
    solution << "0";
    return solution.str();
}

bool AppMC::gen_rhs()
{
    std::uniform_int_distribution<uint32_t> dist{0, 1};
    bool rhs = dist(randomEngine);
    //cout << "rnd rhs:" << (int)rhs << endl;
    return rhs;
}

string AppMC::gen_rnd_bits(const uint32_t size, const uint32_t num_hashes)
{
    string randomBits;
    std::uniform_int_distribution<uint32_t> dist{0, 1000};
    uint32_t cutoff = 500;
    if (conf.sparse && num_hashes > 132) {
        double probability = 13.46*std::log(num_hashes)/num_hashes;
        assert(probability < 0.5);
        cutoff = std::ceil(1000.0*probability);
        cout << "[appmc] sparse hashing used, cutoff: " << cutoff << endl;
    }

    while (randomBits.size() < size) {
        bool val = dist(randomEngine) < cutoff;
        randomBits += '0' + val;
    }
    assert(randomBits.size() >= size);

    //cout << "rnd bits: " << randomBits << endl;
    return randomBits;
}

void AppMC::print_xor(const vector<uint32_t>& vars, const uint32_t rhs)
{
    cout << "[appmc] Added XOR ";
    for (size_t i = 0; i < vars.size(); i++) {
        cout << vars[i]+1;
        if (i < vars.size()-1) {
            cout << " + ";
        }
    }
    cout << " = " << (rhs ? "True" : "False") << endl;
}

template<class T>
inline T AppMC::findMedian(vector<T>& numList)
{
    std::sort(numList.begin(), numList.end());
    size_t medIndex = (numList.size() + 1) / 2;
    size_t at = 0;
    if (medIndex >= numList.size()) {
        at += numList.size() - 1;
        return numList[at];
    }
    at += medIndex;
    return numList[at];
}

template<class T>
inline T AppMC::findMin(vector<T>& numList)
{
    T min = std::numeric_limits<T>::max();
    for (const auto a: numList) {
        if (a < min) {
            min = a;
        }
    }
    return min;
}

int AppMC::correctReturnValue(const lbool ret) const
{
    int retval = -1;
    if (ret == l_True) {
        retval = 10;
    } else if (ret == l_False) {
        retval = 20;
    } else if (ret == l_Undef) {
        retval = 15;
    } else {
        std::cerr << "Something is very wrong, output is neither "
        "l_Undef, nor l_False, nor l_True" << endl;
        exit(-1);
    }

    return retval;
}

void printVersionInfoAppMC()
{
    cout << "c AppMC SHA revision " << ::get_version_sha1() << endl;
    cout << "c AppMC version " << ::get_version_tag() << endl;
    cout << "c AppMC compilation env " << ::get_compilation_env() << endl;
    #ifdef __GNUC__
    cout << "c AppMC compiled with gcc version " << __VERSION__ << endl;
    #else
    cout << "c AppMC compiled with non-gcc compiler" << endl;
    #endif
}

void AppMC::printVersionInfo() const
{
    ::printVersionInfoAppMC();
    cout << solver->get_text_version_info();
}

/* Number of solutions to return from one invocation of gen_n_samples. */
uint32_t AppMC::sols_to_return(uint32_t numSolutions)
{
    if (conf.startiter == 0)   // TODO improve hack for ideal sampling case?
        return numSolutions;
    else if (conf.multisample)
        return loThresh;
    else
        return 1;
}

void AppMC::openLogFile()
{
    if (!conf.logfilename.empty()) {
        logfile.open(conf.logfilename.c_str());
        if (!logfile.is_open()) {
            cout << "[appmc] Cannot open AppMC log file '" << conf.logfilename
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

void AppMC::write_log(
    int iter,
    uint32_t hashCount,
    int found_full,
    uint32_t num_sols,
    uint32_t repeat_sols,
    double used_time
)
{
    if (!conf.logfilename.empty()) {
        logfile
        << std::left
        << std::setw(5) << sampling
        << " " << std::setw(4) << iter
        << " " << std::setw(4) << hashCount
        << " " << std::setw(4) << found_full
        << " " << std::setw(4) << num_sols
        << " " << std::setw(4) << repeat_sols
        << " " << std::setw(7) << std::fixed << std::setprecision(2) << used_time
        << " " << std::setw(7) << std::fixed << std::setprecision(2) << (cpuTimeTotal() - startTime)
        << endl;
    }
}
