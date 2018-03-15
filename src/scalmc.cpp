/*
 CUSP and ScalMC

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

#if defined(__GNUC__) && defined(__linux__)

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <fenv.h>
#endif

#include <stdio.h>
#include <ctime>
#include <cstring>
#include <errno.h>
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
#include <complex>

#include "scalmc.h"
#include "time_mem.h"
#include "cryptominisat5/cryptominisat.h"
#include "cryptominisat5/dimacsparser.h"
#include "cryptominisat5/streambuffer.h"
#include "cryptominisat5/solvertypesmini.h"
#include "signalcode.h"

using std::cout;
using std::cerr;
using std::endl;
using boost::lexical_cast;
using std::list;
using std::map;

string binary(unsigned x, uint32_t length)
{
    uint32_t logSize = (x == 0 ? 1 : log2(x) + 1);
    string s;
    do {
        s.push_back('0' + (x & 1));
    } while (x >>= 1);

    for (uint32_t i = logSize; i < length; i++) {
        s.push_back('0');
    }
    std::reverse(s.begin(), s.end());

    return s;
}

string ScalMC::GenerateRandomBits(uint32_t size)
{
    string randomBits;
    std::uniform_int_distribution<unsigned> uid {0, 2147483647U};
    uint32_t i = 0;
    while (i < size) {
        i += 31;
        randomBits += binary(uid(randomEngine), 31);
    }
    return randomBits;
}

void ScalMC::add_scalmc_options()
{
    scalmc_options.add_options()
    ("help,h", "Prints help")
    ("seed,s", po::value< int >(), "Seed")
    ("pivotAC", po::value(&pivot)->default_value(pivot)
        , "Number of solutions to check for")
    ("mode", po::value(&searchMode)->default_value(searchMode)
        ,"Seach mode. ApproxMX = 0, ScalMC = 1")
    ("tApproxMC", po::value(&tApproxMC)->default_value(tApproxMC)
        , "Number of measurements")
    ("start", po::value(&start_iter)->default_value(start_iter),
         "")
    ("looptout", po::value(&loopTimeout)->default_value(loopTimeout)
        , "Timeout for one measurement, consisting of finding pivotAC solutions")
    ("cuspLogFile", po::value(&cuspLogFile)->default_value(cuspLogFile),
         "Log of SCALMC iterations")
    ("unset", po::value(&unset_vars)->default_value(unset_vars),
         "Try to ask the solver to unset some independent variables, thereby"
         "finding more than one solution at a time")
    ("input", po::value< vector<string> >(), "file(s) to read")
    ("verb,v", po::value(&verb)->default_value(verb), "verbosity")
    ;

    help_options.add(scalmc_options);
    //help_options_complicated.add(scalmc_options);
}

void ScalMC::add_supported_options()
{
    add_scalmc_options();
    p.add("input", 1);

    try {
        po::store(po::command_line_parser(argc, argv).options(help_options).positional(p).run(), vm);
        if (vm.count("help"))
        {
            cout
            << "Approximate counter" << endl;

            cout
            << "scalmc [options] inputfile [drat-trim-file]" << endl << endl;

            cout << help_options << endl;
            std::exit(0);
        }
        po::notify(vm);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::unknown_option> >& c
    ) {
        cerr
        << "ERROR: Some option you gave was wrong. Please give '--help' to get help" << endl
        << "       Unkown option: " << c.what() << endl;
        std::exit(-1);
    } catch (boost::bad_any_cast &e) {
        std::cerr
        << "ERROR! You probably gave a wrong argument type" << endl
        << "       Bad cast: " << e.what()
        << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::invalid_option_value> > what
    ) {
        cerr
        << "ERROR: Invalid value '" << what.what() << "'" << endl
        << "       given to option '" << what.get_option_name() << "'"
        << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::multiple_occurrences> > what
    ) {
        cerr
        << "ERROR: " << what.what() << " of option '"
        << what.get_option_name() << "'"
        << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::required_option> > what
    ) {
        cerr
        << "ERROR: You forgot to give a required option '"
        << what.get_option_name() << "'"
        << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::too_many_positional_options_error> > what
    ) {
        cerr
        << "ERROR: You gave too many positional arguments. Only at most two can be given:" << endl
        << "       the 1st the CNF file input, and optinally, the 2nd the DRAT file output" << endl
        << "    OR (pre-processing)  1st for the input CNF, 2nd for the simplified CNF" << endl
        << "    OR (post-processing) 1st for the solution file" << endl
        ;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::ambiguous_option> > what
    ) {
        cerr
        << "ERROR: The option you gave was not fully written and matches" << endl
        << "       more than one option. Please give the full option name." << endl
        << "       The option you gave: '" << what.get_option_name() << "'" <<endl
        << "       The alternatives are: ";
        for(size_t i = 0; i < what.alternatives().size(); i++) {
            cout << what.alternatives()[i];
            if (i+1 < what.alternatives().size()) {
                cout << ", ";
            }
        }
        cout << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::invalid_command_line_syntax> > what
    ) {
        cerr
        << "ERROR: The option you gave is missing the argument or the" << endl
        << "       argument is given with space between the equal sign." << endl
        << "       detailed error message: " << what.what() << endl
        ;
        std::exit(-1);
    }
}

void print_xor(const vector<uint32_t>& vars, const uint32_t rhs)
{
    cout << "Added XOR ";
    for (size_t i = 0; i < vars.size(); i++) {
        cout << vars[i]+1;
        if (i < vars.size()-1) {
            cout << " + ";
        }
    }
    cout << " = " << (rhs ? "True" : "False") << endl;
}

bool ScalMC::openLogFile()
{
    cusp_logf.open(cuspLogFile.c_str());
    if (!cusp_logf.is_open()) {
        cout << "Cannot open ScalMC log file '" << cuspLogFile
             << "' for writing." << endl;
        exit(1);
    }
    return true;
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

bool ScalMC::AddHash(uint32_t num_xor_cls, vector<Lit>& assumps)
{
    string randomBits = GenerateRandomBits((independent_vars.size() + 1) * num_xor_cls);
    bool rhs = true;
    vector<uint32_t> vars;

    for (uint32_t i = 0; i < num_xor_cls; i++) {
        //new activation variable
        solver->new_var();
        uint32_t act_var = solver->nVars()-1;
        assumps.push_back(Lit(act_var, true));

        vars.clear();
        vars.push_back(act_var);
        rhs = (randomBits[(independent_vars.size() + 1) * i] == '1');

        for (uint32_t j = 0; j < independent_vars.size(); j++) {
            if (randomBits[(independent_vars.size() + 1) * i + j+1] == '1') {
                vars.push_back(independent_vars[j]);
            }
        }
        solver->add_xor_clause(vars, rhs);
        if (verb >= 3) {
            print_xor(vars, rhs);
        }
    }
    return true;
}

int64_t ScalMC::BoundedSATCount(uint32_t maxSolutions, const vector<Lit>& assumps)
{
    cout << "BoundedSATCount looking for " << maxSolutions << " solutions"
    << " T:" << std::setprecision(2) << std::fixed << (cpuTime()-myTime) << endl;

    //Set up things for adding clauses that can later be removed
    solver->new_var();
    uint32_t act_var = solver->nVars()-1;
    vector<Lit> new_assumps(assumps);
    new_assumps.push_back(Lit(act_var, true));

    double start_time = cpuTime();
    uint64_t solutions = 0;
    lbool ret;
    while (solutions < maxSolutions) {
        //solver->set_max_confl(10*1000*1000);
        double this_iter_timeout = loopTimeout-(cpuTime()-start_time);
        solver->set_timeout_all_calls(this_iter_timeout);
        ret = solver->solve(&new_assumps);
        if (verb >= 3) {
            cout << "Found one" << endl;
        }
        if (ret != l_True)
            break;

        size_t num_undef = 0;
        if (solutions < maxSolutions) {
            vector<Lit> lits;
            lits.push_back(Lit(act_var, false));
            for (const uint32_t var: independent_vars) {
                if (solver->get_model()[var] != l_Undef) {
                    lits.push_back(Lit(var, solver->get_model()[var] == l_True));
                } else {
                    num_undef++;
                }
            }
            solver->add_clause(lits);
        }
        if (num_undef) {
            cout << "WOW Num undef:" << num_undef << endl;
        }

        //Try not to be crazy, 2**30 solutions is enough
        if (num_undef <= 30) {
            solutions += 1U << num_undef;
        } else {
            solutions += 1U << 30;
            cout << "WARNING, in this cut there are > 2**30 solutions indicated by the solver!" << endl;
        }
    }
    if (solutions > maxSolutions) {
        solutions = maxSolutions;
    }

    //Remove clauses added
    vector<Lit> cl_that_removes;
    cl_that_removes.push_back(Lit(act_var, false));
    solver->add_clause(cl_that_removes);

    //Timeout
//     if (ret == l_Undef) {
//         must_interrupt.store(false, std::memory_order_relaxed);
//         return -1;
//     }
    return solutions;
}

bool ScalMC::ApproxMC(SATCount& count)
{
    count.clear();
    int64_t currentNumSolutions = 0;
    vector<uint64_t> numHashList;
    vector<int64_t> numCountList;
    vector<Lit> assumps;
    for (uint32_t j = 0; j < tApproxMC; j++) {
        uint64_t hashCount;
        uint32_t repeatTry = 0;
        for (hashCount = 0; hashCount < solver->nVars(); hashCount++) {
            cout << "-> Hash Count " << hashCount << endl;
            double myTime = cpuTimeTotal();
            currentNumSolutions = BoundedSATCount(pivot + 1, assumps);

            //cout << currentNumSolutions << ", " << pivot << endl;
            cusp_logf << "ApproxMC:" << searchMode << ":"
                      << j << ":" << hashCount << ":"
                      << std::fixed << std::setprecision(2) << (cpuTimeTotal() - myTime) << ":"
                      << (int)(currentNumSolutions == (pivot + 1)) << ":"
                      << currentNumSolutions << endl;
            //Timeout!
            if (currentNumSolutions < 0) {
                //Remove all hashes
                assumps.clear();

                if (repeatTry < 2) {    /* Retry up to twice more */
                    assert(hashCount > 0);
                    AddHash(hashCount, assumps); //add new set of hashes
                    solver->simplify(&assumps);
                    hashCount --;
                    repeatTry += 1;
                    cout << "Timeout, try again -- " << repeatTry << endl;
                } else {
                    //this set of hashes does not work, go up
                    AddHash(hashCount + 1, assumps);
                    solver->simplify(&assumps);
                    cout << "Timeout, moving up" << endl;
                }
                continue;
            }

            if (currentNumSolutions < pivot + 1) {
                //less than pivot solutions
                break;
            }

            //Found all solutions needed
            AddHash(1, assumps);
        }
        assumps.clear();
        numHashList.push_back(hashCount);
        numCountList.push_back(currentNumSolutions);
        solver->simplify(&assumps);
    }
    if (numHashList.size() == 0) {
        //UNSAT
        return true;
    }

    auto minHash = findMin(numHashList);
    auto hash_it = numHashList.begin();
    auto cnt_it = numCountList.begin();
    for (; hash_it != numHashList.end() && cnt_it != numCountList.end()
            ; hash_it++, cnt_it++
        ) {
        *cnt_it *= pow(2, (*hash_it) - minHash);
    }
    int medSolCount = findMedian(numCountList);

    count.cellSolCount = medSolCount;
    count.hashCount = minHash;
    return true;
}

void ScalMC::readInAFile(SATSolver* solver2, const string& filename)
{
    solver2->add_sql_tag("filename", filename);
    #ifndef USE_ZLIB
    FILE * in = fopen(filename.c_str(), "rb");
    DimacsParser<StreamBuffer<FILE*, FN> > parser(solver, NULL, 2);
    #else
    gzFile in = gzopen(filename.c_str(), "rb");
    DimacsParser<StreamBuffer<gzFile, GZ> > parser(solver, NULL, 2);
    #endif

    if (in == NULL) {
        std::cerr
        << "ERROR! Could not open file '"
        << filename
        << "' for reading: " << strerror(errno) << endl;

        std::exit(1);
    }

    if (!parser.parse_DIMACS(in, false)) {
        exit(-1);
    }

    independent_vars.swap(parser.independent_vars);
    call_after_parse();

    #ifndef USE_ZLIB
        fclose(in);
    #else
        gzclose(in);
    #endif
}

void ScalMC::readInStandardInput(SATSolver* solver2)
{
    cout
    << "c Reading from standard input... Use '-h' or '--help' for help."
    << endl;

    #ifndef USE_ZLIB
    FILE * in = stdin;
    #else
    gzFile in = gzdopen(0, "rb"); //opens stdin, which is 0
    #endif

    if (in == NULL) {
        std::cerr << "ERROR! Could not open standard input for reading" << endl;
        std::exit(1);
    }

    #ifndef USE_ZLIB
    DimacsParser<StreamBuffer<FILE*, FN> > parser(solver2, NULL, 2);
    #else
    DimacsParser<StreamBuffer<gzFile, GZ> > parser(solver2, NULL, 2);
    #endif

    if (!parser.parse_DIMACS(in, false)) {
        exit(-1);
    }

    #ifdef USE_ZLIB
        gzclose(in);
    #endif
}

int ScalMC::solve()
{
    myTime = cpuTime();
    //set seed
    if (vm.count("seed") == 0) {
        cerr << "ERROR: You must provide a seed value with the '-s NUM' option" << endl;
        exit(-1);
    }
    unsigned int seed = vm["seed"].as<int>();
    randomEngine.seed(seed);

    openLogFile();
    startTime = cpuTimeTotal();

    //solver = new SATSolver(&must_interrupt);
    solver = new SATSolver();
    solverToInterrupt = solver;

    /*conf.reconfigure_at = 0;
    conf.reconfigure_val = 15;*/
    solver->set_allow_otf_gauss();
    if (verb > 2) {
        solver->set_verbosity(verb-2);
    }
    CMSat::GaussConf gconf;
    gconf.max_matrix_rows = 3000;
    gconf.decision_until = 3000;
    gconf.max_num_matrixes = 1;
    gconf.min_matrix_rows = 5;
    gconf.autodisable = false;
    gconf.only_nth_gauss_save = 10;
    solver->set_gauss_config(gconf);

    if (unset_vars) {
        solver->set_greedy_undef();
    }
    printVersionInfo();
    vector<string> inp = vm["input"].as<vector<string> >();
    if (inp.size() >= 1) {
        if (inp.size() > 1) {
            cout << "ERROR: can only parse in one file" << endl;
        }
        readInAFile(solver, inp[0].c_str());
    } else {
        readInStandardInput(solver);
    }

    if (start_iter > independent_vars.size()) {
        cout << "ERROR: Manually-specified start_iter"
             "is larger than the size of the independent set.\n" << endl;
        return -1;
    }

    SATCount solCount;
    cout << "Using start iteration " << start_iter << endl;

    bool finished = false;
    if (searchMode == 0) {
        finished = ApproxMC(solCount);
    } else {
        finished = ScalApproxMC(solCount);
    }
    cout << "ApproxMC finished in " << (cpuTimeTotal() - startTime) << " s" << endl;
    if (!finished) {
        cout << " (TIMED OUT)" << endl;
        return 0;
    }

    if (solCount.hashCount == 0 && solCount.cellSolCount == 0) {
        cout << "The input formula is unsatisfiable." << endl;
        return correctReturnValue(l_False);
    }

    if (verb) {
        solver->print_stats();
    }

    cout << "Number of solutions is: " << solCount.cellSolCount
         << " x 2^" << solCount.hashCount << endl;

    return correctReturnValue(l_True);
}

int main(int argc, char** argv)
{
    #if defined(__GNUC__) && defined(__linux__)
    feenableexcept(FE_INVALID   |
                   FE_DIVBYZERO |
                   FE_OVERFLOW
                  );
    #endif

    ScalMC main(argc, argv);
    main.add_supported_options();
    return main.solve();
}

void ScalMC::call_after_parse()
{
    if (independent_vars.empty()) {
        cout
        << "c WARNING! No independent vars were set using 'c ind var1 [var2 var3 ..] 0'"
        "notation in the CNF." << endl
        << " c ScalMC may work substantially worse!" << endl;
        for (size_t i = 0; i < solver->nVars(); i++) {
            independent_vars.push_back(i);
        }
    }
    solver->set_independent_vars(&independent_vars);
}

//For ScalApproxMC only
void ScalMC::SetHash(uint32_t clausNum, std::map<uint64_t,Lit>& hashVars, vector<Lit>& assumps)
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
            AddHash(clausNum-hashVars.size(),assumps);
            for (uint64_t i = hashVars.size(); i < clausNum; i++) {
                hashVars[i] = assumps[i];
            }
        }
    }
}

//For ScalApproxMC only
bool ScalMC::ScalApproxMC(SATCount& count)
{
    count.clear();
    vector<uint64_t> numHashList;
    vector<int64_t> numCountList;
    vector<Lit> assumps;

    uint64_t hashCount = start_iter;
    uint64_t hashPrev = 0;
    uint64_t mPrev = 0;

    double myTime = cpuTimeTotal();
    cout << "ScalApproxMC: Starting up, initial measurement" << endl;
    if (hashCount == 0) {
        int64_t currentNumSolutions = BoundedSATCount(pivot+1,assumps);
        cusp_logf << "ApproxMC:"<< searchMode<<":"<<"0:0:"
                  << std::fixed << std::setprecision(2) << (cpuTimeTotal() - myTime) << ":"
                  << (int)(currentNumSolutions == (pivot + 1)) << ":"
                  << currentNumSolutions << endl;

        //Din't find at least pivot+1
        if (currentNumSolutions <= pivot) {
            cout << "Did not find at least pivot+1 (" << pivot << ") we found only " << currentNumSolutions << ", exiting ScalApproxMC" << endl;
            count.cellSolCount = currentNumSolutions;
            count.hashCount = 0;
            return true;
        }
        hashCount++;
    }

    for (uint32_t j = 0; j < tApproxMC; j++) {
        map<uint64_t,int64_t> countRecord;
        map<uint64_t,uint32_t> succRecord;
        map<uint64_t,Lit> hashVars; //map assumption var to XOR hash

        uint32_t repeatTry = 0;
        uint64_t numExplored = 1;
        uint64_t lowerFib = 0, upperFib = independent_vars.size();

        while (numExplored < independent_vars.size()) {
            cout << "Num Explored: " << numExplored
                 << " ind set size: " << independent_vars.size() << endl;
            myTime = cpuTimeTotal();
            uint64_t swapVar = hashCount;
            SetHash(hashCount,hashVars,assumps);
            cout << "Number of XOR hashes active: " << hashCount << endl;
            int64_t currentNumSolutions = BoundedSATCount(pivot + 1, assumps);

            //cout << currentNumSolutions << ", " << pivot << endl;
            cusp_logf << "ApproxMC:" << searchMode<<":"
                      << j << ":" << hashCount << ":"
                      << std::fixed << std::setprecision(2) << (cpuTimeTotal() - myTime) << ":"
                      << (int)(currentNumSolutions == (pivot + 1)) << ":"
                      << currentNumSolutions << endl;
            //Timeout!
            if (currentNumSolutions < 0) {
                //Remove all hashes
                assumps.clear();
                hashVars.clear();

                if (repeatTry < 2) {    /* Retry up to twice more */
                    assert(hashCount > 0);
                    SetHash(hashCount,hashVars,assumps);
                    solver->simplify(&assumps);
                    hashCount --;
                    repeatTry += 1;
                    cout << "Timeout, try again -- " << repeatTry << endl;
                } else {
                    //this set of hashes does not work, go up
                    SetHash(hashCount + 1, hashVars, assumps);
                    solver->simplify(&assumps);
                    cout << "Timeout, moving up" << endl;
                }
                hashCount = swapVar;
                continue;
            }

            if (currentNumSolutions < pivot + 1) {
                numExplored = lowerFib+independent_vars.size()-hashCount;
                if (succRecord.find(hashCount-1) != succRecord.end()
                    && succRecord[hashCount-1] == 1
                ) {
                    numHashList.push_back(hashCount);
                    numCountList.push_back(currentNumSolutions);
                    mPrev = hashCount;
                    //less than pivot solutions
                    break;
                }
                succRecord[hashCount] = 0;
                countRecord[hashCount] = currentNumSolutions;
                if (std::abs<int64_t>((int64_t)hashCount - (int64_t)mPrev) <= 2 && mPrev != 0) {
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
                assert(currentNumSolutions == pivot+1);

                numExplored = hashCount + independent_vars.size()-upperFib;
                if (succRecord.find(hashCount+1) != succRecord.end()
                    && succRecord[hashCount+1] == 0
                ) {
                    numHashList.push_back(hashCount+1);
                    numCountList.push_back(countRecord[hashCount+1]);
                    mPrev = hashCount+1;
                    break;
                }
                succRecord[hashCount] = 1;
                if (std::abs<int64_t>((int64_t)hashCount - (int64_t)mPrev) < 2 && mPrev!=0) {
                    lowerFib = hashCount;
                    hashCount ++;
                } else if (lowerFib + (hashCount - lowerFib)*2 >= upperFib-1) {
                    lowerFib = hashCount;
                    hashCount = (lowerFib+upperFib)/2;
                } else {
                    //printf("hashPrev:%d hashCount:%d\n",hashPrev, hashCount);
                    hashCount = lowerFib + (hashCount -lowerFib)*2;
                }
            }
            hashPrev = swapVar;
        }
        assumps.clear();
        solver->simplify(&assumps);
        hashCount =mPrev;
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

void ScalMC::printVersionInfo() const
{
    cout << "c CryptoMiniSat version " << solver->get_version() << endl;
    cout << "c CryptoMiniSat SHA revision " << solver->get_version_sha1() << endl;
    cout << "c CryptoMiniSat compilation env " << solver->get_compilation_env() << endl;
    #ifdef __GNUC__
    cout << "c compiled with gcc version " << __VERSION__ << endl;
    #else
    cout << "c compiled with non-gcc compiler" << endl;
    #endif
}

int ScalMC::correctReturnValue(const lbool ret) const
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
