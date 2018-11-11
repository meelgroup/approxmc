/*
 ApproxMC and ScalGen

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


#ifndef AppMC_H_
#define AppMC_H_

#include "approxmcconfig.h"
#include <fstream>
#include <random>
#include <map>
#include <cstdint>
#include <cryptominisat5/cryptominisat.h>

using std::string;
using std::vector;
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

class AppMC {
public:
    AppMC()
    {
    }

    ~AppMC()
    {
    }

    int solve(AppMCConfig _conf);


    uint32_t ScalGen();
    string GenerateRandomBits(const uint32_t size, const uint32_t numhashes);
    string binary(const uint32_t x, const uint32_t length);
    uint32_t SolutionsToReturn(uint32_t numSolutions);
    void generate_samples();
    bool gen_rhs();
    uint32_t ScalGen(
        uint32_t samples
        , uint32_t sampleCounter
        , std::map<string, uint32_t>& solutionMap
        , uint32_t* lastSuccessfulHashOffset
        , double timeReference
    );
    int ScalGenCall(
        uint32_t samples
        , uint32_t sampleCounter
        , std::map<string, uint32_t>& solutionMap
        , uint32_t* lastSuccessfulHashOffset
        , double timeReference
    );
    uint32_t loThresh;
    uint32_t hiThresh;
    uint32_t threshold_scalgen;
    SATSolver* solver = NULL;
    void printVersionInfo() const;

private:
    AppMCConfig conf;
    bool count(SATCount& count);
    void add_scalmc_options();
    bool ScalAppMC(SATCount& count);
    bool add_hash(uint32_t num_xor_cls, vector<Lit>& assumps, uint32_t total_num_hashes);
    void SetHash(uint32_t clausNum, std::map<uint64_t,Lit>& hashVars, vector<Lit>& assumps);

    int correctReturnValue(const lbool ret) const;
    void output_samples();

    void add_solution_to_map(
        const vector<lbool>& model
        , std::map<std::string, uint32_t>* solutionMap) const;

    int64_t bounded_sol_count(
        uint32_t maxSolutions,
        const vector<Lit>& assumps,
        const uint32_t hashCount,
        std::map<std::string, uint32_t>* solutionMap = NULL,
        uint32_t minSolutions = 1
    );

    void readInAFile(SATSolver* solver2, const string& filename);
    void readInStandardInput(SATSolver* solver2);


    double startTime;
    std::map< std::string, std::vector<uint32_t>> globalSolutionMap;
    void openLogFile();
    void call_after_parse();

    std::ofstream logfile;
    std::mt19937 randomEngine;
    double total_runtime; //runTime

    int argc;
    char** argv;
};


#endif //AppMC_H_
