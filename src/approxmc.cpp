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
#include "constants.h"
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
        Counter counter;
        Config conf;
    };
}

using namespace ApproxMC;

DLL_PUBLIC AppMC::AppMC()
{
    data = new AppMCPrivateData;
    data->counter.solver = new SATSolver();
    data->counter.solver->set_up_for_scalmc();
    data->counter.solver->set_allow_otf_gauss();
    data->counter.solver->set_xor_detach(data->conf.cms_detach_xor);
}

DLL_PUBLIC AppMC::~AppMC()
{
    delete data->counter.solver;
    delete data;
}

// Helper function, used only in this unit
void setup_sampling_vars(AppMCPrivateData* data)
{
    if (data->conf.sampling_set.empty()) {
        if (data->conf.verb) {
            cout
            << "c [appmc] WARNING! Sampling set was not declared! We will be **VERY** slow"
            << endl;
        }
        for (size_t i = 0; i < data->counter.solver->nVars(); i++) {
            data->conf.sampling_set.push_back(i);
        }
    }

    if (data->conf.verb) {
        cout << "c [appmc] Sampling set size: " << data->conf.sampling_set.size() << endl;
        if (data->conf.sampling_set.size() > 100) {
            cout
            << "c [appmc] Sampling var set contains over 100 variables, not displaying"
            << endl;
        } else {
            cout << "c [appmc] Sampling set: ";
            for (auto v: data->conf.sampling_set) {
                cout << v+1 << " ";
            }
            cout << "0" << endl;
        }
    }

    data->counter.solver->set_sampling_vars(&(data->conf.sampling_set));
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
    if (verb > 2) {
        data->counter.solver->set_verbosity(data->conf.verb-2);
    }
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

DLL_PUBLIC void AppMC::set_detach_xors(uint32_t detach_xors)
{
    data->conf.cms_detach_xor = detach_xors;
    data->counter.solver->set_xor_detach(data->conf.cms_detach_xor);
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
    if (data->conf.verb > 2) {
        cout << "c [appmc] using seed: " << data->conf.seed << endl;
    }

    if (data->conf.epsilon < 0.0) {
        cout << "[appmc] ERROR: invalid epsilon" << endl;
        exit(-1);
    }

    if (data->conf.delta <= 0.0 || data->conf.delta > 1.0) {
        cout << "[appmc] ERROR: invalid delta" << endl;
        exit(-1);
    }

    setup_sampling_vars(data);
    SolCount sol_count = data->counter.solve(data->conf);
    return sol_count;
}

DLL_PUBLIC void AppMC::set_projection_set(const vector<uint32_t>& vars)
{
    for(const auto& v: vars) {
        if (v >= data->counter.solver->nVars()) {
            std::cout << "ERROR: function set_projection_set() called with variable that is larger than the number of variables inside the solver. Exiting." << endl;
            assert(false);
            exit(-1);
        }
    }
    data->conf.sampling_set = vars;
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

DLL_PUBLIC bool AppMC::add_clause(const vector<CMSat::Lit>& lits)
{
    return data->counter.solver->add_clause(lits);
}

DLL_PUBLIC bool AppMC::add_xor_clause(const vector<uint32_t>& vars, bool rhs)
{
    return data->counter.solver->add_xor_clause(vars, rhs);
}

DLL_PUBLIC bool AppMC::add_bnn_clause(
            const std::vector<CMSat::Lit>& lits,
            signed cutoff,
            Lit out)
{
    return data->counter.solver->add_bnn_clause(lits, cutoff, out);
}

DLL_PUBLIC void AppMC::set_detach_warning()
{
    data->counter.solver->set_verbosity_detach_warning(true);
}

DLL_PUBLIC CMSat::SATSolver* AppMC::get_solver()
{
    return data->counter.solver;
}


DLL_PUBLIC const std::vector<uint32_t>& AppMC::get_sampling_set() const
{
    return data->conf.sampling_set;
}

DLL_PUBLIC void AppMC::set_dump_intermediary_cnf(const bool dump_intermediary_cnf)
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
