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


#ifndef APPROXMC_H__
#define APPROXMC_H__

#include <cstdint>
#include <string>
#include <vector>
#ifdef CMS_LOCAL_BUILD
#include "cryptominisat.h"
#else
#include <cryptominisat5/cryptominisat.h>
#endif
namespace ApproxMC {

#ifdef _WIN32
struct __declspec(dllexport) SolCount
#else
struct SolCount
#endif
{
    void clear()
    {
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
    AppMC();
    ~AppMC();
    void set_projection_set(const std::vector<uint32_t>& vars);
    ApproxMC::SolCount count();
    bool find_one_solution();

    // Adding constraints
    void new_var();
    void new_vars(uint32_t num);
    uint32_t nVars();
    bool add_clause(const std::vector<CMSat::Lit>& lits);
    bool add_xor_clause(const std::vector<uint32_t>& vars, bool rhs);
    bool add_bnn_clause(
        const std::vector<CMSat::Lit>& lits,
        signed cutoff,
        CMSat::Lit out = CMSat::lit_Undef);

    // Information about approxmc
    std::string get_version_info();
    void print_stats(const double start_time);

    //Main options
    void set_up_log(std::string log_file_name);
    void set_verbosity(uint32_t verb);
    void set_detach_warning();
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
    void set_dump_intermediary_cnf(const bool dump_intermediary_cnf);

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

#endif
