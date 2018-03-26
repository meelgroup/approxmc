/*
 ScalMC and ScalMC

 Copyright (c) 2009-2015, Mate Soos. All rights reserved.
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


#ifndef ScalMC_H_
#define ScalMC_H_

//#include "main.h"
#include <boost/program_options.hpp>
#include <fstream>
#include <random>
#include <map>
#include <cstdint>
#include <cryptominisat5/cryptominisat.h>
#include <cryptominisat5/solverconf.h>

using std::string;
using std::vector;

#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace CMSat;

struct SATCount {
    void clear()
    {
        SATCount tmp;
        *this = tmp;
    }
    uint32_t hashCount = 0;
    uint32_t cellSolCount = 0;
};

class ScalMC {
public:
    ScalMC(int _argc, char** _argv):
        scalmc_options("ScalMC options")
        , argc(_argc)
        , argv(_argv)
    {
        //must_interrupt.store(false, std::memory_order_relaxed);
        solver = new SATSolver;
    }

    int solve();
    void add_supported_options();

    po::options_description scalmc_options = po::options_description("ScalMC options");
    po::options_description help_options;
    po::variables_map vm;
    po::positional_options_description p;

private:
    SolverConf conf;
    bool count(SATCount& count);
    void add_scalmc_options();
    bool ScalScalMC(SATCount& count);
    bool AddHash(uint32_t num_xor_cls, vector<Lit>& assumps);
    void SetHash(uint32_t clausNum, std::map<uint64_t,Lit>& hashVars, vector<Lit>& assumps);

    void printVersionInfo() const;
    int correctReturnValue(const lbool ret) const;

    int64_t BoundedSATCount(uint32_t maxSolutions, const vector<Lit>& assumps);
    lbool BoundedSAT(
        uint32_t maxSolutions, uint32_t minSolutions
        , vector<Lit>& assumptions
        , std::map<std::string, uint32_t>& solutionMap
        , uint32_t* solutionCount
    );
    string GenerateRandomBits(uint32_t size);

    void readInAFile(SATSolver* solver2, const string& filename);
    void readInStandardInput(SATSolver* solver2);

    //config
    std::string logfile;

    double startTime;
    std::map< std::string, std::vector<uint32_t>> globalSolutionMap;
    bool openLogFile();
    //std::atomic<bool> must_interrupt;
    void call_after_parse();

    uint32_t start_iter = 0;
    uint32_t pivot = 72; //precision
    uint32_t tScalMC = 9; //confidence of 0.81
    double   loopTimeout = 2500;
    int      unset_vars = 0;
    std::ofstream cusp_logf;
    std::mt19937 randomEngine;
    SATSolver* solver;
    vector<uint32_t> independent_vars;
    unsigned verb = 1;
    double total_runtime; //runTime
    int what_to_break = 0;

    int argc;
    char** argv;
};


#endif //ScalMC_H_
