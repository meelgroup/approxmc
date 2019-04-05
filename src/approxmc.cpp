/*
 ApproxMC

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

#include <ctime>
#include <cstring>
#include <errno.h>
#include <algorithm>
#include <string.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <set>
#include <fstream>
#include <sys/stat.h>
#include <string.h>
#include <list>
#include <array>
#include <cmath>
#include <complex>

#include "approxmc.h"
#include "time_mem.h"
#include "cryptominisat5/cryptominisat.h"
#include "cryptominisat5/solvertypesmini.h"
#include "GitSHA1.h"

using std::cout;
using std::cerr;
using std::endl;
using std::list;
using std::map;

void print_xor(const vector<uint32_t>& vars, const uint32_t rhs)
{
    cout << "[appmc] Added XOR ";
    for (size_t i = 0; i < vars.size(); i++) {
        cout << vars[i]+1;
        if (i < vars.size()-1) {
            cout << " + ";
        }
    }
    cout << " = " << (rhs ? "True" : "False") << endl;
}

void AppMC::openLogFile()
{
    if (!conf.logfilename.empty()) {
        logfile.open(conf.logfilename.c_str());
        if (!logfile.is_open()) {
            cout << "[appmc] Cannot open AppMC log file '" << conf.logfilename
                 << "' for writing." << endl;
            exit(1);
        }
    }
}

template<class T>
inline T findMedian(vector<T>& numList)
{
    std::sort(numList.begin(), numList.end());
    size_t medIndex = (numList.size() + 1) / 2;
    size_t at = 0;
    if (medIndex >= numList.size()) {
        at += numList.size() - 1;
        return numList[at];
    }
    at += medIndex;
    return numList[at];
}

template<class T>
inline T findMin(vector<T>& numList)
{
    T min = std::numeric_limits<T>::max();
    for (const auto a: numList) {
        if (a < min) {
            min = a;
        }
    }
    return min;
}

bool AppMC::add_hash(uint32_t num_xor_cls, vector<Lit>& assumps, uint32_t total_num_hashes)
{
    const string randomBits =
        GenerateRandomBits(conf.sampling_set.size() * num_xor_cls, total_num_hashes);

    bool rhs;
    vector<uint32_t> vars;

    for (uint32_t i = 0; i < num_xor_cls; i++) {
        //new activation variable
        solver->new_var();
        uint32_t act_var = solver->nVars()-1;
        assumps.push_back(Lit(act_var, true));

        vars.clear();
        vars.push_back(act_var);
        rhs = gen_rhs();

        for (uint32_t j = 0; j < conf.sampling_set.size(); j++) {
            if (randomBits.at(conf.sampling_set.size() * i + j) == '1') {
                vars.push_back(conf.sampling_set[j]);
            }
        }
        solver->add_xor_clause(vars, rhs);
        if (conf.verb_appmc_cls) {
            print_xor(vars, rhs);
        }
    }
    return true;
}

int64_t AppMC::bounded_sol_count(
        uint32_t maxSolutions,
        const vector<Lit>& assumps,
        const uint32_t hashCount
) {
    cout << "[appmc] "
    "[ " << std::setw(7) << std::setprecision(2) << std::fixed
    << (cpuTimeTotal()-total_runtime)
    << " ]"
    << " bounded_sol_count looking for " << std::setw(4) << maxSolutions << " solutions"
    << " -- hashes active: " << hashCount << endl;

    //Set up things for adding clauses that can later be removed
    vector<lbool> model;
    vector<Lit> new_assumps(assumps);
    solver->new_var();
    uint32_t act_var = solver->nVars()-1;
    new_assumps.push_back(Lit(act_var, true));
    if (hashCount > 2) {
        solver->simplify(&new_assumps);
    }

    uint64_t solutions = 0;
    lbool ret;
    double last_found_time = cpuTimeTotal();
    while (solutions < maxSolutions) {
        ret = solver->solve(&new_assumps);
        assert(ret == l_False || ret == l_True);

        if (conf.verb >=2 ) {
            cout << "[appmc] bounded_sol_count ret: " << std::setw(7) << ret;
            if (ret == l_True) {
                cout << " sol no.  " << std::setw(3) << solutions;
            } else {
                cout << " No more. " << std::setw(3) << "";
            }
            cout << " T: "
            << std::setw(7) << std::setprecision(2) << std::fixed << (cpuTimeTotal()-total_runtime)
            << " -- hashes act: " << hashCount
            << " -- T since last: "
            << std::setw(7) << std::setprecision(2) << std::fixed << (cpuTimeTotal()-last_found_time)
            << endl;
            last_found_time = cpuTimeTotal();
        }

        if (ret != l_True) {
            break;
        }
        model = solver->get_model();

        if (solutions < maxSolutions) {
            vector<Lit> lits;
            lits.push_back(Lit(act_var, false));
            for (const uint32_t var: conf.sampling_set) {
                if (solver->get_model()[var] != l_Undef) {
                    lits.push_back(Lit(var, solver->get_model()[var] == l_True));
                } else {
                    assert(false);
                }
            }
            if (conf.verb_appmc_cls) {
                cout << "[appmc] Adding banning clause: " << lits << endl;
            }
            solver->add_clause(lits);
        }
        solutions++;
    }

    //Remove clauses added
    vector<Lit> cl_that_removes;
    cl_that_removes.push_back(Lit(act_var, false));
    solver->add_clause(cl_that_removes);

    assert(ret != l_Undef);
    return solutions;
}

bool AppMC::gen_rhs()
{
    std::uniform_int_distribution<uint32_t> dist{0, 1};
    bool rhs = dist(randomEngine);
    //cout << "rnd rhs:" << (int)rhs << endl;
    return rhs;
}

string AppMC::GenerateRandomBits(const uint32_t size, const uint32_t num_hashes)
{
    string randomBits;
    std::uniform_int_distribution<uint32_t> dist{0, 1000};
    uint32_t cutoff = 500;
    if (conf.sparse && num_hashes > 132) {
        double probability = 13.46*std::log(num_hashes)/num_hashes;
        assert(probability < 0.5);
        cutoff = std::ceil(1000.0*probability);
        cout << "[appmc] sparse hashing used, cutoff: " << cutoff << endl;
    }

    while (randomBits.size() < size) {
        bool val = dist(randomEngine) < cutoff;
        randomBits += '0' + val;
    }
    assert(randomBits.size() >= size);

    //cout << "rnd bits: " << randomBits << endl;
    return randomBits;
}

int AppMC::solve(AppMCConfig _conf)
{
    conf = _conf;

    openLogFile();
    randomEngine.seed(conf.seed);
    total_runtime = cpuTimeTotal();
    cout << "[appmc] Using start iteration " << conf.start_iter << endl;

    SATCount solCount;
    bool finished = count(solCount);
    assert(finished);

    cout << "[appmc] FINISHED AppMC T: " << (cpuTimeTotal() - startTime) << " s" << endl;
    if (solCount.hashCount == 0 && solCount.cellSolCount == 0) {
        cout << "[appmc] Formula was UNSAT " << endl;
    }

    if (conf.verb > 2) {
        solver->print_stats();
    }

    cout << "[appmc] Number of solutions is: "
    << solCount.cellSolCount
     << " x 2^" << solCount.hashCount << endl;

    return correctReturnValue(l_True);
}

void AppMC::SetHash(uint32_t clausNum, std::map<uint64_t,Lit>& hashVars, vector<Lit>& assumps)
{
    if (clausNum < assumps.size()) {
        uint64_t numberToRemove = assumps.size()- clausNum;
        for (uint64_t i = 0; i<numberToRemove; i++) {
            assumps.pop_back();
        }
    } else {
        if (clausNum > assumps.size() && assumps.size() < hashVars.size()) {
            for (uint32_t i = assumps.size(); i< hashVars.size() && i < clausNum; i++) {
                assumps.push_back(hashVars[i]);
            }
        }
        if (clausNum > hashVars.size()) {
            add_hash(clausNum-hashVars.size(), assumps, clausNum);
            for (uint64_t i = hashVars.size(); i < clausNum; i++) {
                hashVars[i] = assumps[i];
            }
        }
    }
}

bool AppMC::count(SATCount& count)
{
    count.clear();
    vector<uint64_t> numHashList;
    vector<int64_t> numCountList;
    vector<Lit> assumps;

    uint64_t hashCount = conf.start_iter;
    uint64_t hashPrev = 0;
    uint64_t mPrev = 0;

    double myTime = cpuTimeTotal();
    cout << "[appmc] Starting up, initial measurement" << endl;
    if (hashCount == 0) {
        int64_t currentNumSolutions = bounded_sol_count(conf.threshold+1, assumps, count.hashCount);
        if (!conf.logfilename.empty()) {
            logfile << "appmc:"
            <<"0:0:"
            << std::fixed << std::setprecision(2) << (cpuTimeTotal() - myTime) << ":"
            << (int)(currentNumSolutions == (conf.threshold + 1)) << ":"
            << currentNumSolutions << endl;
        }

        //Din't find at least threshold+1
        if (currentNumSolutions <= conf.threshold) {
            cout << "[appmc] Did not find at least threshold+1 (" << conf.threshold << ") we found only " << currentNumSolutions << ", exiting AppMC" << endl;
            count.cellSolCount = currentNumSolutions;
            count.hashCount = 0;
            return true;
        }
        hashCount++;
    }

    for (uint32_t j = 0; j < conf.measurements; j++) {
        map<uint64_t,int64_t> countRecord;
        map<uint64_t,uint32_t> succRecord;
        map<uint64_t,Lit> hashVars; //map assumption var to XOR hash


        //Note, the rank of a random NxN matrix is not N of course. It has an expected
        //rank that is of course lower than N. So we need to shoot higher.
        //https://math.stackexchange.com/questions/324150/expected-rank-of-a-random-binary-matrix
        //Apparently this question is analyzed in Kolchin's book Random Graphs in sect. 3.2.
        //Thanks to Yash Pote to digging this one out. Very helpful.
        uint64_t total_max_xors = std::ceil((double)conf.sampling_set.size()*1.2)+5;

        uint64_t numExplored = 0;
        uint64_t lowerFib = 0;
        uint64_t upperFib = total_max_xors;

        while (numExplored < total_max_xors) {
            cout << "[appmc] Explored: " << std::setw(4) << numExplored
                 << " ind set size: " << std::setw(6) << conf.sampling_set.size() << endl;
            myTime = cpuTimeTotal();
            uint64_t swapVar = hashCount;
            SetHash(hashCount,hashVars,assumps);
            cout << "[appmc] hashes active: " << std::setw(6) << hashCount << endl;
            int64_t currentNumSolutions = bounded_sol_count(conf.threshold + 1, assumps, hashCount);

            //cout << currentNumSolutions << ", " << threshold << endl;
            if (!conf.logfilename.empty()) {
                logfile << "appmc:"
                << j << ":" << hashCount << ":"
                << std::fixed << std::setprecision(2) << (cpuTimeTotal() - myTime) << ":"
                << (int)(currentNumSolutions == (conf.threshold + 1)) << ":"
                << currentNumSolutions << endl;
            }

            if (currentNumSolutions <= conf.threshold) {
                numExplored = lowerFib + total_max_xors - hashCount;

                //check success record if it exists
                if (succRecord.find(hashCount-1) != succRecord.end()
                    && succRecord[hashCount-1] == 1
                ) {
                    numHashList.push_back(hashCount);
                    numCountList.push_back(currentNumSolutions);
                    mPrev = hashCount;
                    //less than threshold solutions
                    break;
                }

                //No success record
                succRecord[hashCount] = 0;
                countRecord[hashCount] = currentNumSolutions;
                if (std::abs<int64_t>((int64_t)hashCount - (int64_t)mPrev) <= 2
                    && mPrev != 0
                ) {
                    upperFib = hashCount;
                    hashCount--;
                } else {
                    if (hashPrev > hashCount) {
                        hashPrev = 0;
                    }
                    upperFib = hashCount;
                    if (hashPrev > lowerFib) {
                        lowerFib = hashPrev;
                    }
                    hashCount = (upperFib+lowerFib)/2;
                }
            } else {
                assert(currentNumSolutions == conf.threshold+1);
                numExplored = hashCount + total_max_xors - upperFib;

                //Check if success record for +1 hashcount exists and is 0
                if (succRecord.find(hashCount+1) != succRecord.end()
                    && succRecord[hashCount+1] == 0
                ) {
                    numHashList.push_back(hashCount+1);
                    numCountList.push_back(countRecord[hashCount+1]);
                    mPrev = hashCount+1;
                    break;
                }

                //No success record of hashCount+1 or it's not 0
                succRecord[hashCount] = 1;
                if (std::abs<int64_t>((int64_t)hashCount - (int64_t)mPrev) < 2
                    && mPrev!=0
                ) {
                    lowerFib = hashCount;
                    hashCount ++;
                } else if (lowerFib + (hashCount - lowerFib)*2 >= upperFib-1) {
                    lowerFib = hashCount;
                    hashCount = (lowerFib+upperFib)/2;
                } else {
                    //cout << "hashPrev: " << hashPrev << " hashCount: " << hashCount << endl;
                    hashCount = lowerFib + (hashCount -lowerFib)*2;
                }
            }
            hashPrev = swapVar;
        }
        assumps.clear();
        hashCount = mPrev;
    }
    if (numHashList.size() == 0) {
        //UNSAT
        return true;
    }

    auto minHash = findMin(numHashList);
    auto cnt_it = numCountList.begin();
    for (auto hash_it = numHashList.begin()
        ; hash_it != numHashList.end() && cnt_it != numCountList.end()
        ; hash_it++, cnt_it++
    ) {
        *cnt_it *= pow(2, (*hash_it) - minHash);
    }
    int medSolCount = findMedian(numCountList);

    count.cellSolCount = medSolCount;
    count.hashCount = minHash;
    return true;
}

///////////
// Useful helper functions
///////////

void printVersionInfoAppMC()
{
    cout << "c AppMC SHA revision " << ::get_version_sha1() << endl;
    cout << "c AppMC version " << ::get_version_tag() << endl;
    cout << "c AppMC compilation env " << ::get_compilation_env() << endl;
    #ifdef __GNUC__
    cout << "c AppMC compiled with gcc version " << __VERSION__ << endl;
    #else
    cout << "c AppMC compiled with non-gcc compiler" << endl;
    #endif
}

void AppMC::printVersionInfo() const
{
    ::printVersionInfoAppMC();
    cout << solver->get_text_version_info();
}

int AppMC::correctReturnValue(const lbool ret) const
{
    int retval = -1;
    if (ret == l_True) {
        retval = 10;
    } else if (ret == l_False) {
        retval = 20;
    } else if (ret == l_Undef) {
        retval = 15;
    } else {
        std::cerr << "Something is very wrong, output is neither l_Undef, nor l_False, nor l_True" << endl;
        exit(-1);
    }

    return retval;
}
