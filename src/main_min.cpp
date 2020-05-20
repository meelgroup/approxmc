/*
 ApproxMC and Unigen

 Copyright (c) 2009-2018, Mate Soos. All rights reserved.
 Copyright (c) 2014, Supratik Chakraborty, Kuldeep S. Meel, Moshe Y. Vardi
 Copyright (c) 2015, Supratik Chakraborty, Daniel J. Fremont,
 Kuldeep S. Meel, Sanjit A. Seshia, Moshe Y. Vardi

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

#include <string>
#include <vector>
using std::string;
using std::vector;

#include "approxmcconfig.h"
#include "time_mem.h"
#include "approxmc.h"
#include <cryptominisat5/cryptominisat.h>
#include "cryptominisat5/dimacsparser.h"
#include "cryptominisat5/streambuffer.h"

using namespace CMSat;
using std::cout;
using std::cerr;
using std::endl;
AppMC* appmc = NULL;

SolverConf satconf;
AppMCConfig conf;

#if defined _WIN32
    #define DLL_PUBLIC __declspec(dllexport)
#else
    #define DLL_PUBLIC __attribute__ ((visibility ("default")))
    #define DLL_LOCAL  __attribute__ ((visibility ("hidden")))
#endif

void set_sampling_vars()
{
    if (conf.sampling_set.empty()) {
        cout
        << "c [appmc] WARNING! Sampling set was not declared with 'vp var1 [var2 var3 ..] 0'"
        " notation in the PCNF." << endl
        << "c [appmc] we may work substantially worse!" << endl;
        for (size_t i = 0; i < appmc->solver->nVars(); i++) {
            conf.sampling_set.push_back(i);
        }
    }
    cout << "c [appmc] Sampling set size: " << conf.sampling_set.size() << endl;
    cout << "c [appmc] Sampling set: ";
    for (auto v: conf.sampling_set) {
        cout << v+1 << ", ";
    }
    cout << endl;
    appmc->solver->set_sampling_set(&conf.sampling_set);
}

int solve(const char* cnf)
{
    appmc = new AppMC;
    appmc->printVersionInfo();
    appmc->solver = new SATSolver;
    cout << "c [appmc] using seed: " << conf.seed << endl;
    conf.logfilename = "log.txt";

    DimacsParser<StreamBuffer<const char*, CH>> parser(appmc->solver, NULL, 2);
    if (!parser.parse_DIMACS(cnf, false)) {
        exit(-1);
    }
    conf.sampling_set.swap(parser.sampling_set);

    set_sampling_vars();
    appmc->solve(conf);
    return 0;
}

extern "C" {

DLL_PUBLIC int csolve(const char *input) {
  return solve(input);
}

}
