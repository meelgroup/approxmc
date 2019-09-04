/*
 ApproxMC and gen_n_samples

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
using std::map;
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
    string gen_rnd_bits(const uint32_t size, const uint32_t numhashes);
    string binary(const uint32_t x, const uint32_t length);
    uint32_t sols_to_return(uint32_t numSolutions);
    void generate_samples();
    bool gen_rhs();
    uint32_t gen_n_samples(
        const uint32_t samples
        , uint32_t* lastSuccessfulHashOffset
    );
    uint32_t loThresh;
    uint32_t hiThresh;
    uint32_t threshold_appmcgen;
    SATSolver* solver = NULL;
    void printVersionInfo() const;
    void set_samples_file(std::ostream* os);

private:
    AppMCConfig conf;
    void count(SATCount& count);
    void add_appmc_options();
    bool ScalAppMC(SATCount& count);
    void add_hash(uint32_t num_xor_cls,
                  vector<Lit>& assumps,
                  uint32_t total_num_hashes);
    void set_num_hashes(
        uint32_t num_wanted,
        map<uint64_t,Lit>& hashVars,
        vector<Lit>& assumps
    );

    //Helper functions
    template<class T> T findMedian(vector<T>& numList);
    template<class T> T findMin(vector<T>& numList);
    void print_xor(const vector<uint32_t>& vars, const uint32_t rhs);
    int correctReturnValue(const lbool ret) const;
    std::string get_solution_str(const vector<lbool>& model);
    void one_measurement_count(
        vector<uint64_t>& numHashList,
        vector<int64_t>& numCountList,
        uint64_t& hashPrev,
        uint64_t& mPrev,
        uint64_t& hashCount,
        const uint32_t iter
    );

    void add_glob_banning_cls(
        const vector<vector<lbool>>* glob_model
        , const uint32_t act_var);

    int64_t bounded_sol_count(
        uint32_t maxSolutions,
        const vector<Lit>* assumps,
        const uint32_t hashCount,
        uint32_t minSolutions = 1,
        std::vector<vector<lbool>>* glob_model = NULL,
        vector<string>* out_solutions = NULL
    );

    void readInAFile(SATSolver* solver2, const string& filename);
    void readInStandardInput(SATSolver* solver2);


    double startTime;
    std::ostream* samples_out = NULL;
    void openLogFile();
    void call_after_parse();

    std::ofstream logfile;
    std::mt19937 randomEngine;
    double total_runtime; //runTime
    uint32_t orig_num_vars;

    int argc;
    char** argv;
};

inline void AppMC::set_samples_file(std::ostream* out)
{
    samples_out = out;
}


#endif //AppMC_H_
