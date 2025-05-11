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


#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include <cryptominisat5/cryptominisat.h>
namespace ApproxMC {

#ifdef _WIN32
class __declspec(dllexport) SolCount
#else
class SolCount
#endif
{
    public:
    void clear() {
        SolCount tmp;
        *this = tmp;
    }
    bool valid = false;
    uint32_t hashCount = 0;
    uint32_t cellSolCount = 0;
};

struct AppMCPrivateData;
#ifdef _WIN32
class __declspec(dllexport) AppMC
#else
class AppMC
#endif
{
public:
    AppMC(const std::unique_ptr<CMSat::FieldGen>& _fg);
    ~AppMC();
    ApproxMC::SolCount count();
    bool find_one_solution();

    // Sampling set
    void set_sampl_vars(const std::vector<uint32_t>& vars);
    void set_opt_sampl_vars(const std::vector<uint32_t>& vars);
    bool get_sampl_vars_set() const;
    bool get_opt_sampl_vars_set() const { return false; }
    const std::vector<uint32_t>& get_sampl_vars() const;
    void set_multiplier_weight(const std::unique_ptr<CMSat::Field>& weight);
    const std::unique_ptr<CMSat::Field>& get_multiplier_weight() const;
    void set_weighted(const bool weighted);
    void set_projected(const bool projected);
    void set_lit_weight(const CMSat::Lit& lit, const std::unique_ptr<CMSat::Field>& weight);

    // Adding constraints
    void new_var();
    void new_vars(uint32_t num);
    uint32_t nVars();
    bool add_clause(const std::vector<CMSat::Lit>& lits);
    bool add_red_clause(const std::vector<CMSat::Lit>& lits);
    bool add_xor_clause(const std::vector<CMSat::Lit>& lits, bool rhs);
    bool add_xor_clause(const std::vector<uint32_t>& vars, bool rhs);

    // Information about approxmc
    static std::string get_version_sha1();
    void print_stats(const double start_time);

    //Main options
    void set_verbosity(uint32_t verb);
    void set_seed(uint32_t seed);
    void set_epsilon(double epsilon);
    void set_delta(double delta);
    CMSat::SATSolver* get_solver();

    //Misc options -- do NOT to change unless you know what you are doing!
    void set_start_iter(uint32_t start_iter);
    void set_verb_cls(uint32_t verb_cls);
    void set_var_elim_ratio(double var_elim_ratio);
    void set_detach_xors(uint32_t detach_xors);
    void set_reuse_models(uint32_t reuse_models);
    void set_sparse(uint32_t sparse);
    void set_simplify(uint32_t simplify);
    void set_dump_intermediary_cnf(const int dump_intermediary_cnf);
    void set_debug(int debug);
    void set_force_sol_extension(int val);

    //Querying default values
    const std::vector<uint32_t>& get_sampling_set() const;
    double get_epsilon();
    uint32_t get_seed();
    double get_delta();
    uint32_t get_simplify();
    double get_var_elim_ratio();
    uint32_t get_sparse();
    bool get_reuse_models();

private:
    ////////////////////////////
    // Do not bother with this, it's private
    ////////////////////////////
    AppMCPrivateData* data;
};

}
