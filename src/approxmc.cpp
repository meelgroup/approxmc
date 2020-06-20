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

#include "counter.h"
#include "constants.h"
#include "approxmc/approxmc.h"
#include "config.h"
#include <iostream>

using std::cout;
using std::endl;

namespace ApproxMC {
    struct AppMCPrivateData {
        Counter counter;
        Config conf;
    };
}

using namespace ApproxMC;

AppMC::AppMC()
{
    data = new AppMCPrivateData;
    data->counter.solver = new SATSolver();
    data->counter.solver->set_up_for_scalmc();
    data->counter.solver->set_allow_otf_gauss();
    data->counter.solver->set_xor_detach(data->conf.cms_detach_xor);
}

AppMC::~AppMC()
{
    delete data->counter.solver;
    delete data;
}

void setup_sampling_vars(AppMCPrivateData* data)
{
    if (data->conf.sampling_set.empty()) {
        cout
        << "c [appmc] WARNING! Sampling set was not declared with 'c ind var1 [var2 var3 ..] 0'"
        " notation in the CNF." << endl
        << "c [appmc] we may work substantially worse!" << endl;
        for (size_t i = 0; i < data->counter.solver->nVars(); i++) {
            data->conf.sampling_set.push_back(i);
        }
    }
    cout << "c [appmc] Sampling set size: " << data->conf.sampling_set.size() << endl;
    if (data->conf.sampling_set.size() > 100) {
        cout
        << "c [appmc] Sampling var set contains over 100 variables, not displaying"
        << endl;
    } else {
        cout << "c [appmc] Sampling set: ";
        for (auto v: data->conf.sampling_set) {
            cout << v+1 << ", ";
        }
        cout << endl;
    }
    data->counter.solver->set_sampling_vars(&(data->conf.sampling_set));
}

string AppMC::get_version_info()
{
    return data->counter.get_version_info();
}

void AppMC::set_up_log(string log_file_name)
{
    data->conf.logfilename = log_file_name;
}

void AppMC::set_verbosity(uint32_t verb)
{
    data->conf.verb = verb;
    if (verb > 2) {
        data->counter.solver->set_verbosity(data->conf.verb-2);
    }
}

ApproxMC::SolCount AppMC::count()
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

    if (data->conf.startiter > data->conf.sampling_set.size()) {
        cout << "c [appmc] ERROR: Manually-specified start_iter"
             "is larger than the size of the sampling set.\n" << endl;
        exit(-1);
    }

    SolCount sol_count = data->counter.solve(data->conf);
    return sol_count;
}

void AppMC::set_sampling_set(const vector<uint32_t>& vars)
{
    data->conf.sampling_set = vars;
}


uint32_t AppMC::nVars() {
    return data->counter.solver->nVars();
}


void AppMC::new_vars(uint32_t num)
{
    data->counter.solver->new_vars(num);
}

void AppMC::new_var()
{
    data->counter.solver->new_var();
}

void AppMC::add_clause(const vector<CMSat::Lit>& lits)
{
    data->counter.solver->add_clause(lits);
}

void AppMC::add_xor_clause(const vector<uint32_t>& vars, bool rhs)
{
    data->counter.solver->add_xor_clause(vars, rhs);
}
