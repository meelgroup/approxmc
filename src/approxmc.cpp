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

#include "approxmc.h"
#include "counter.h"
#include "appmc_constants.h"
#include "config.h"
#include <iostream>

using std::cout;
using std::endl;
using namespace AppMCInt;

#if defined _WIN32
    #define DLL_PUBLIC __declspec(dllexport)
#else
    #define DLL_PUBLIC __attribute__ ((visibility ("default")))
    #define DLL_LOCAL  __attribute__ ((visibility ("hidden")))
#endif

namespace ApproxMC {
    struct AppMCPrivateData {
        AppMCPrivateData(): counter(conf) {}
        Config conf;
        Counter counter;
        bool sampl_vars_declared = false;
    };
}

using namespace ApproxMC;

DLL_PUBLIC AppMC::AppMC()
{
    data = new AppMCPrivateData;
    data->counter.solver = new SATSolver();
    data->counter.solver->set_up_for_scalmc();
    data->counter.solver->set_allow_otf_gauss();
}

DLL_PUBLIC AppMC::~AppMC()
{
    delete data->counter.solver;
    delete data;
}

// Helper function, used only in this unit
void setup_sampling_vars(AppMCPrivateData* data)
{
    if (data->conf.verb) {
        cout << "c o [appmc] Sampling set size: " << data->conf.sampl_vars.size() << endl;
        if (data->conf.sampl_vars.size() > 100) {
            cout
            << "c o [appmc] Sampling var set contains over 100 variables, not displaying"
            << endl;
        } else {
            cout << "c o [appmc] Sampling set: ";
            for (auto v: data->conf.sampl_vars) {
                cout << v+1 << " ";
            }
            cout << "0" << endl;
        }
    }

    data->counter.solver->set_sampl_vars(data->conf.sampl_vars);
}

DLL_PUBLIC string AppMC::get_version_info()
{
    return data->counter.get_version_info();
}

DLL_PUBLIC void AppMC::set_up_log(string log_file_name)
{
    data->conf.logfilename = log_file_name;
}

DLL_PUBLIC void AppMC::set_verbosity(uint32_t verb)
{
    data->conf.verb = verb;
    data->counter.solver->set_verbosity(std::max<int>(0, (int)data->conf.verb-2));
}

DLL_PUBLIC void AppMC::set_seed(uint32_t seed)
{
    data->conf.seed = seed;
}

DLL_PUBLIC void AppMC::set_epsilon(double epsilon)
{
    data->conf.epsilon = epsilon;
}

DLL_PUBLIC void AppMC::set_delta(double delta)
{
    data->conf.delta = delta;
}

DLL_PUBLIC void AppMC::set_debug(int debug) { data->conf.debug = debug; }
DLL_PUBLIC void AppMC::set_force_sol_extension(int val) {
    data->conf.force_sol_extension = val;
}

DLL_PUBLIC void AppMC::set_start_iter(uint32_t start_iter)
{
    data->conf.start_iter = start_iter;
}

DLL_PUBLIC void AppMC::set_verb_cls(uint32_t verb_cls)
{
    data->conf.verb_cls = verb_cls;
}

DLL_PUBLIC void AppMC::set_simplify(uint32_t simplify)
{
    data->conf.simplify = simplify;
    if (!simplify) {
        data->counter.solver->set_no_bve();
        data->counter.solver->set_no_bva();
        data->counter.solver->set_scc(0);
        data->counter.solver->set_simplify(0);
    }
}

DLL_PUBLIC void AppMC::set_var_elim_ratio(double var_elim_ratio)
{
    data->conf.var_elim_ratio = var_elim_ratio;
}

DLL_PUBLIC void AppMC::set_reuse_models(uint32_t reuse_models)
{
    data->conf.reuse_models = reuse_models;
}

DLL_PUBLIC void AppMC::set_sparse(uint32_t sparse)
{
    data->conf.sparse = sparse;
}

DLL_PUBLIC double AppMC::get_epsilon()
{
    return data->conf.epsilon;
}

DLL_PUBLIC double AppMC::get_delta()
{
    return data->conf.delta;
}

DLL_PUBLIC uint32_t AppMC::get_simplify()
{
    return data->conf.simplify;
}

DLL_PUBLIC double AppMC::get_var_elim_ratio()
{
    return data->conf.var_elim_ratio;
}

DLL_PUBLIC uint32_t AppMC::get_sparse()
{
    return data->conf.sparse;
}

DLL_PUBLIC uint32_t AppMC::get_seed()
{
    return data->conf.seed;
}

DLL_PUBLIC bool AppMC::get_reuse_models()
{
    return data->conf.reuse_models;
}

DLL_PUBLIC bool AppMC::find_one_solution()
{
    return data->counter.find_one_solution();
}

DLL_PUBLIC ApproxMC::SolCount AppMC::count()
{
    if (!data->sampl_vars_declared) {
        cout << "ERROR: Sampling set was not declared!" << endl;
        exit(-1);
    }
    if (data->conf.verb > 2) {
        cout << "c o [appmc] using seed: " << data->conf.seed << endl;
    }

    if (data->conf.epsilon < 0.0) {
        cout << "ERROR: invalid epsilon" << endl;
        exit(-1);
    }

    if (data->conf.delta <= 0.0 || data->conf.delta > 1.0) {
        cout << "ERROR: invalid delta: " << data->conf.delta << endl;
        exit(-1);
    }

    setup_sampling_vars(data);
    SolCount sol_count = data->counter.solve();
    return sol_count;
}

DLL_PUBLIC void AppMC::set_sampl_vars(const vector<uint32_t>& vars)
{
    data->sampl_vars_declared = true;
    data->conf.sampl_vars_set = true;
    for(const auto& v: vars) {
        if (v >= data->counter.solver->nVars()) {
            std::cout << "ERROR: function set_projection_set() called with variable that is larger than the number of variables inside the solver. Exiting." << endl;
            assert(false);
            exit(-1);
        }
    }
    data->conf.sampl_vars = vars;
}

DLL_PUBLIC const std::vector<uint32_t>& AppMC::get_sampl_vars() const  {
    return data->conf.sampl_vars;
}


DLL_PUBLIC uint32_t AppMC::nVars() {
    return data->counter.solver->nVars();
}

DLL_PUBLIC void AppMC::new_vars(uint32_t num)
{
    data->counter.solver->new_vars(num);
}

DLL_PUBLIC void AppMC::new_var()
{
    data->counter.solver->new_var();
}

DLL_PUBLIC bool AppMC::add_red_clause(const vector<CMSat::Lit>& lits)
{
    return data->counter.solver->add_red_clause(lits);
}

DLL_PUBLIC bool AppMC::add_clause(const vector<CMSat::Lit>& lits)
{
    return data->counter.solver_add_clause(lits);
}

DLL_PUBLIC bool AppMC::add_xor_clause(const vector<Lit>& lits, bool rhs)
{
    return data->counter.solver_add_xor_clause(lits, rhs);
}

DLL_PUBLIC bool AppMC::add_xor_clause(const vector<uint32_t>& vars, bool rhs)
{
    return data->counter.solver_add_xor_clause(vars, rhs);
}

DLL_PUBLIC CMSat::SATSolver* AppMC::get_solver()
{
    return data->counter.solver;
}

DLL_PUBLIC const std::vector<uint32_t>& AppMC::get_sampling_set() const
{
    return data->conf.sampl_vars;
}

DLL_PUBLIC void AppMC::set_dump_intermediary_cnf(const int dump_intermediary_cnf)
{
    data->conf.dump_intermediary_cnf = dump_intermediary_cnf;
}

DLL_PUBLIC void AppMC::print_stats(const double start_time)
{
    data->counter.solver->set_verbosity(1);
    data->counter.solver->print_stats(start_time);
    data->counter.solver->set_verbosity(0);
    if (data->conf.verb > 2) {
        data->counter.solver->set_verbosity(data->conf.verb);
    }
}

DLL_PUBLIC bool AppMC::get_sampl_vars_set() const {
    return data->conf.sampl_vars_set;
}

 DLL_PUBLIC void AppMC::set_multiplier_weight(const mpq_class& weight) {
     data->conf.multiplier_weight = weight;
 }

 DLL_PUBLIC const mpq_class& AppMC::get_multiplier_weight() const {
     return data->conf.multiplier_weight;
 }

DLL_PUBLIC void AppMC::set_weighted(const bool weighted) {
    if (weighted) {
        cout << "ERROR: Weighted ApproxMC not supported" << endl;
        exit(-1);
    }
}

DLL_PUBLIC void AppMC::set_projected(const bool) {
}

DLL_PUBLIC void AppMC::set_lit_weight(const Lit&, const mpq_class&) {
    cout << "ERROR: Weighted ApproxMC is not supported" << endl;
    exit(-1);
}

 DLL_PUBLIC void AppMC::set_opt_sampl_vars(const std::vector<uint32_t>&) {
     // Not interesting for AppMC
 }
