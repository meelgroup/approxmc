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

#include "appmcconfig.h"
#include <memory>
#include <random>
#include <map>
#include <utility>
#include <cstdint>
#include "approxmc.h"
#include "appmc_constants.h"

using std::string;
using std::vector;
using std::map;
using std::pair;
using namespace CMSat;
using std::unique_ptr;

namespace AppMCInt {

struct SavedModel
{
    SavedModel (const vector<lbool>& _model, uint32_t _hash_num):
        model(_model), hash_num(_hash_num) {}
    const vector<lbool> model;
    uint32_t hash_num;
};

struct Hash {
    Hash () = default;
    Hash(const uint32_t _act_var, const vector<uint32_t>& _hash_vars, const bool _rhs):
        act_var(_act_var), hash_vars(_hash_vars), rhs(_rhs) {}

    uint32_t act_var;
    vector<uint32_t> hash_vars;
    bool rhs;
};

struct HashesModels {
    map<uint64_t, Hash> hashes;
    vector<SavedModel> glob_model; //global table storing models

    void clear() {
        hashes.clear();
        vector<SavedModel> clean_glob_model;
        for(auto const& m: glob_model) {
            if (m.hash_num == 0) clean_glob_model.push_back(m);
        }
        std::swap(glob_model, clean_glob_model);
    }
};

struct SolNum {
    SolNum (uint64_t _solutions, uint64_t _repeated):
        solutions(_solutions), repeated(_repeated) {}
    uint64_t solutions = 0;
    uint64_t repeated = 0;
};

struct SparseData {
    explicit SparseData(int _table_no) : table_no(_table_no) {}

    uint32_t next_index = 0;
    double sparseprob = 0.5;
    int table_no = -1;
};

class Counter {
public:
    Counter(Config& _conf, const unique_ptr<FieldGen>& _fg) : fg(_fg->dup()), conf(_conf) {}
    ApproxMC::SolCount solve();
    string gen_rnd_bits(const uint32_t size,
                        const uint32_t numhashes, SparseData& sparse_data);
    string binary(const uint32_t x, const uint32_t length);
    bool find_one_solution();
    bool gen_rhs();
    uint32_t threshold_appmcgen;
    unique_ptr<SATSolver> solver = nullptr;
    ApproxMC::SolCount calc_est_count();
    const Constants constants;
    bool solver_add_clause(const vector<Lit>& cl);
    bool solver_add_xor_clause(const vector<uint32_t>& vars, const bool rhs);
    bool solver_add_xor_clause(const vector<Lit>& lits, const bool rhs);

private:
    unique_ptr<FieldGen> fg;
    Config& conf;
    ApproxMC::SolCount count();
    void add_appmc_options();
    Hash add_hash(uint32_t total_num_hashes, SparseData& sparse_data);
    SolNum bounded_sol_count(
        uint32_t max_sols,
        const vector<Lit>* assumps,
        const uint32_t hash_cnt,
        const uint32_t iter,
        HashesModels* hm = nullptr
    );
    vector<Lit> set_num_hashes(
        uint32_t num_wanted,
        map<uint64_t, Hash>& hashes,
        SparseData& sparse_data
    );
    void simplify();

    ////////////////
    //Helper functions
    ////////////////
    void dump_cnf_from_solver(const vector<Lit>& assumps, const uint32_t iter, const lbool result);
    void print_xor(const vector<uint32_t>& vars, const uint32_t rhs);
    void one_measurement_count(
        int64_t& prev_measure,
        const unsigned iter,
        SparseData sparse_data,
        HashesModels* hm
    );
    void appmc7_one_measurement_count(
        int64_t& prev_measure,
        const unsigned iter,
        SparseData sparse_data,
        HashesModels* hm
    );
    void call_after_parse();
    void ban_one(const uint32_t act_var, const vector<lbool>& model);
    void check_model(
        const vector<lbool>& model,
        const HashesModels* const hm,
        const uint32_t hash_count
    );
    bool check_model_against_hash(const Hash& h, const vector<lbool>& model);
    uint64_t add_glob_banning_cls(
        const HashesModels* glob_model = nullptr
        , const uint32_t act_var = std::numeric_limits<uint32_t>::max()
        , const uint32_t num_hashes = std::numeric_limits<uint32_t>::max()
    );

    int find_best_sparse_match();
    void set_up_probs_threshold_measurements(uint32_t& measurements, SparseData& sparse_data);
    double calc_error_bound(uint32_t t, double p);

    //Data so we can output temporary count when catching the signal
    vector<uint64_t> num_hash_list;
    vector<double> num_count_list;
    template<class T> T find_median(const vector<T>& nums);
    template<class T> T find_min(const vector<T>& nums);

    ////////////////
    // internal data
    ////////////////
    double start_time;
    std::mt19937 rnd_engine;
    uint32_t orig_num_vars;
    double total_inter_simp_time = 0;
    uint32_t threshold; //precision, it's computed
    uint32_t cnf_dump_no = 0;
    vector<vector<Lit>> cls_in_solver; // needed for accurate dumping
    vector<pair<vector<Lit>, bool>> xors_in_solver; // needed for accurate dumping
    double alpha = -1;
    double beta = -1;
};

}
