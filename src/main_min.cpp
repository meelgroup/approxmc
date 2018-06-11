/*
 ScalMC and ScalGen

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

#include "scalmcconfig.h"
#include "time_mem.h"
#include "scalmc.h"
#include <cryptominisat5/cryptominisat.h>
#include <cryptominisat5/solverconf.h>
#include "cryptominisat5/dimacsparser.h"
#include "cryptominisat5/streambuffer.h"

using namespace CMSat;
using std::cout;
using std::cerr;
using std::endl;
ScalMC* scalmc = NULL;

SolverConf satconf;
ScalMCConfig conf;

#if defined _WIN32
    #define DLL_PUBLIC __declspec(dllexport)
#else
    #define DLL_PUBLIC __attribute__ ((visibility ("default")))
    #define DLL_LOCAL  __attribute__ ((visibility ("hidden")))
#endif

void set_indep_vars()
{
    if (conf.independent_vars.empty()) {
        cout
        << "[scalmc] WARNING! No independent vars were set using 'c ind var1 [var2 var3 ..] 0'"
        "notation in the CNF." << endl
        << "[scalmc] we may work substantially worse!" << endl;
        for (size_t i = 0; i < scalmc->solver->nVars(); i++) {
            conf.independent_vars.push_back(i);
        }
    }
    cout << "[scalmc] Num independent vars: " << conf.independent_vars.size() << endl;
    cout << "[scalmc] Independent vars: ";
    for (auto v: conf.independent_vars) {
        cout << v+1 << ", ";
    }
    cout << endl;
    scalmc->solver->set_independent_vars(&conf.independent_vars);
}

int solve(const char* cnf)
{
    scalmc = new ScalMC;
    scalmc->printVersionInfo();
    scalmc->solver = new SATSolver;
    cout << "[scalmc] using seed: " << conf.seed << endl;
    conf.logfilename = "log.txt";

    DimacsParser<StreamBuffer<const char*, CH>> parser(scalmc->solver, NULL, 2);
    if (!parser.parse_DIMACS(cnf, false)) {
        exit(-1);
    }
    conf.independent_vars.swap(parser.independent_vars);

    set_indep_vars();
    scalmc->solve(conf);
    return 0;
}

extern "C" {

DLL_PUBLIC int csolve(const char *input) {
  return solve(input);
}

}
