/***************************************************************************************[Solver.cc]
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2009, Niklas Sorensson
Copyright (c) 2009-2012, Mate Soos
Copyright (c) 2014-16, Kuldeep S. Meel
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 **************************************************************************************************/



#include <ctime>
#include <cstring>
#include <errno.h>
#include <string.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <omp.h>
#include <map>
#include <set>
#include "cmsat/constants.h"
#include <signal.h>
#include <time.h>
#include<list>
#include <cmath>

#include <signal.h>
#include "cmsat/Main.h"
#include "cmsat/time_mem.h"
#include "cmsat/constants.h"
#include "cmsat/DimacsParser.h"

#if defined(__linux__)
#include <fpu_control.h>
#endif

#include "cmsat/Main.h"

using namespace CMSat;

Main::Main(int _argc, char** _argv) :
numThreads(1)
, grouping(false)
, debugLib(false)
, debugNewVar(false)
, printResult(true)
, max_nr_of_solutions(1)
, fileNamePresent(false)
, twoFileNamesPresent(false)
, argc(_argc)
, argv(_argv) {
}


double findMean(std::list<int> numList) {
    double sum = 0;
    for (std::list<int>::iterator it = numList.begin(); it != numList.end(); it++) {
        sum += *it;
    }
    return (sum * 1.0 / numList.size());
}

double findMedian(std::list<int> numList) {
    numList.sort();
    int medIndex = int((numList.size()) / 2);
    std::list<int>::iterator it = numList.begin();
    if (medIndex >= (int) numList.size()) {
        std::advance(it, numList.size() - 1);
        return double(*it);
    }
    std::advance(it, medIndex);
    return double(*it);
}

int findMin(std::list<int> numList) {
    int min = INT_MAX;
    for (std::list<int>::iterator it = numList.begin(); it != numList.end(); it++) {
        if ((*it) < min) {
            min = *it;
        }
    }
    return min;
}
std::map<uint32_t, Solver*> solversToInterrupt;
std::set<uint32_t> finished;
timer_t mytimer;
bool need_clean_exit;
bool timerSetFirstTime;
void start_timer(int num) {
    struct itimerspec value;
    value.it_value.tv_sec = num; //waits for n seconds before sending timer signal
    value.it_value.tv_nsec = 0;
    value.it_interval.tv_sec = 0; //exipire once
    value.it_interval.tv_nsec = 0;
    if (timerSetFirstTime){
    timer_create(CLOCK_REALTIME, NULL, &mytimer);
  //  timer_delete(mytimer);
    }
    timerSetFirstTime = false;
    timer_settime(mytimer, 0, &value, NULL);
}
void SIGINT_handler(int) {
#pragma omp critical
  exit(1);
}
void SIGALARM_handler(int) {
#pragma omp critical
    {
        Solver& solver = *solversToInterrupt.begin()->second;
        printf("\n");
        std::cerr << "*** INTERRUPTED ***" << std::endl;
        if (solver.conf.needToDumpLearnts || solver.conf.needToDumpOrig || need_clean_exit) {
            solver.needToInterrupt = true;
            std::cerr << "*** Please wait. We need to interrupt cleanly" << std::endl;
            std::cerr << "*** This means we might need to finish some calculations" << std::endl;
        } else {
            if (solver.conf.verbosity >= 1) solver.printStats();
            exit(1);
        }
    }
}
void Main::readInAFile(const std::string& filename, Solver& solver) {
#pragma omp single
    if (solver.conf.verbosity >= 1) {
        std::cout << "c Reading file '" << filename << "'" << std::endl;
    }
#ifdef DISABLE_ZLIB
    FILE * in = fopen(filename.c_str(), "rb");
#else
    gzFile in = gzopen(filename.c_str(), "rb");
#endif // DISABLE_ZLIB

#pragma omp single
    if (in == NULL) {
        std::cout << "ERROR! Could not open file '" << filename << "' for reading" << std::endl;
        exit(1);
    }

    DimacsParser parser(&solver, debugLib, debugNewVar, grouping);
    parser.parse_DIMACS(in);

#ifdef DISABLE_ZLIB
    fclose(in);
#else
    gzclose(in);
#endif // DISABLE_ZLIB
}

void Main::readInStandardInput(Solver& solver) {
    if (solver.conf.verbosity >= 1) {
        std::cout << "c Reading from standard input... Use '-h' or '--help' for help." << std::endl;
    }
#ifdef DISABLE_ZLIB
    FILE * in = stdin;
#else
    gzFile in = gzdopen(fileno(stdin), "rb");
#endif // DISABLE_ZLIB

    if (in == NULL) {
        std::cout << "ERROR! Could not open standard input for reading" << std::endl;
        exit(1);
    }

    DimacsParser parser(&solver, debugLib, debugNewVar, grouping);
    parser.parse_DIMACS(in);

#ifndef DISABLE_ZLIB
    gzclose(in);
#endif // DISABLE_ZLIB
}

void Main::parseInAllFiles(Solver& solver) {
    double myTime = cpuTime();

    //First read normal extra files
    if ((debugLib || debugNewVar) && filesToRead.size() > 0) {
        std::cout << "debugNewVar and debugLib must both be OFF to parse in extra files" << std::endl;
        exit(-1);
    }
    for (uint32_t i = 0; i < filesToRead.size(); i++) {
        readInAFile(filesToRead[i].c_str(), solver);
    }

    //Then read the main file or standard input
    if (!fileNamePresent) {
        readInStandardInput(solver);
    } else {
        string filename = argv[(twoFileNamesPresent ? argc - 2 : argc - 1)];
        readInAFile(filename, solver);
    }

    if (solver.conf.verbosity >= 1) {
        std::cout << "c Parsing time: "
                << std::fixed << std::setw(5) << std::setprecision(2) << (cpuTime() - myTime)
                << " s" << std::endl;
    }
}

void Main::printUsage(char** argv) {
#ifdef DISABLE_ZLIB
    printf("USAGE: %s [options] <input-file> <result-output-file>\n\n  where input is plain DIMACS.\n\n", argv[0]);
#else
    printf("USAGE: %s [options] <input-file> <result-output-file>\n\n  where input may be either in plain or gzipped DIMACS.\n\n", argv[0]);
#endif // DISABLE_ZLIB
    printf("OPTIONS:\n\n");
    printf(" --delta = <float> provide value of 1-confidence\n");   
    printf(" --epsilon = <float> tolerance\n");
    printf(" --searchMode = {0,1}. Sets search mode to linear (0) or Logarithmic (1) (default=1)\n");
    printf("  --verbosity      = {0,1,2}\n");
/*
    printf("  --polarity-mode  = {true,false,rnd,auto} [default: auto]. Selects the default\n");
    printf("                     polarity mode. Auto is the Jeroslow&Wang method\n");
    //printf("  -decay         = <num> [ 0 - 1 ]\n");
    printf("  --rnd-freq       = <num> [ 0 - 1 ]\n");
   
    printf("  --randomize      = <seed> [0 - 2^32-1] Sets random seed, used for picking\n");
    printf("                     decision variables (default = 0)\n");
    printf("  --restrict       = <num> [1 - varnum] when picking random variables to branch\n");
    printf("                     on, pick one that in the 'num' most active vars useful\n");
    printf("                     for cryptographic problems, where the question is the key,\n");
    printf("                     which is usually small (e.g. 80 bits)\n");
    printf("  --gaussuntil     = <num> Depth until which Gaussian elimination is active.\n");
    printf("                     Giving 0 switches off Gaussian elimination\n");
    printf("  --restarts       = <num> [1 - 2^32-1] No more than the given number of\n");
    printf("                     restarts will be performed during search\n");
    printf("  --nonormxorfind    Don't find and collect >2-long xor-clauses from\n");
    printf("                     regular clauses\n");
    printf("  --nobinxorfind     Don't find and collect 2-long xor-clauses from\n");
    printf("                     regular clauses\n");
    printf("  --noregbxorfind    Don't regularly find and collect 2-long xor-clauses\n");
    printf("                     from regular clauses\n");
    printf("  --doextendedscc    Do strongly conn. comp. finding using non-exist. bins\n");
    printf("  --noconglomerate   Don't conglomerate 2 xor clauses when one var is dependent\n");
    printf("  --nosimplify       Don't do regular simplification rounds\n");
    printf("  --greedyunbound    Greedily unbound variables that are not needed for SAT\n");
    printf("  --debuglib         Solve at specific 'c Solver::solve()' points in the CNF\n");
    printf("                     file. Used to debug file generated by Solver's\n");
    printf("                     needLibraryCNFFile() function\n");
    printf("  --debugnewvar      Add new vars at specific 'c Solver::newVar()' points in \n");
    printf("                     the CNF file. Used to debug file generated by Solver's\n");
    printf("                     needLibraryCNFFile() function.\n");
    printf("  --novarreplace     Don't perform variable replacement. Needed for programmable\n");
    printf("                     solver feature\n");
    printf("  --restart        = {auto, static, dynamic}   Which kind of restart strategy to\n");
    printf("                     follow. Default is auto\n");
    printf("  --dumplearnts    = <filename> If interrupted or reached restart limit, dump\n");
    printf("                     the learnt clauses to the specified file. Maximum size of\n");
    printf("                     dumped clauses can be specified with next option.\n");
    printf("  --maxdumplearnts = [0 - 2^32-1] When dumping the learnts to file, what\n");
    printf("                     should be maximum length of the clause dumped. Useful\n");
    printf("                     to make the resulting file smaller. Default is 2^32-1\n");
    printf("                     note: 2-long XOR-s are always dumped.\n");
    printf("  --dumporig       = <filename> If interrupted or reached restart limit, dump\n");
    printf("                     the original problem instance, simplified to the\n");
    printf("                     current point.\n");
    printf("  --alsoread       = <filename> Also read this file in\n");
    printf("                     Can be used to re-read dumped learnts, for example\n");
    printf("  --maxsolutions     Search for given amount of solutions\n");
    printf("                     Can only be used in single-threaded more (\"--threads=1\")\n");
    printf("  --pavgbranch       Print average branch depth\n");
    printf("  --nofailedlit      Don't search for failed literals, and don't search for lits\n");
    printf("                     propagated both by 'varX' and '-varX'\n");
    printf("  --noheuleprocess   Don't try to minimise XORs by XOR-ing them together.\n");
    printf("                     Algo. as per global/local substitution in Heule's thesis\n");
    printf("  --nosatelite       Don't do clause subsumption, clause strengthening and\n");
    printf("                     variable elimination (implies -novarelim and -nosubsume1).\n");
    printf("  --noxorsubs        Don't try to subsume xor-clauses.\n");
    printf("  --nosolprint       Don't print the satisfying assignment if the solution\n");
    printf("                     is SAT\n");
    printf("  --novarelim        Don't perform variable elimination as per Een and Biere\n");
    printf("  --nosubsume1       Don't perform clause contraction through resolution\n");
#ifdef USE_GAUSS
    printf("  --nomatrixfind     Don't find distinct matrixes. Put all xors into one\n");
    printf("                     big matrix\n");
    printf("  --noordercol       Don't order variables in the columns of Gaussian\n");
    printf("                     elimination. Effectively disables iterative reduction\n");
    printf("                     of the matrix\n");
    printf("  --noiterreduce     Don't reduce iteratively the matrix that is updated\n");
    printf("  --maxmatrixrows    [0 - 2^32-1] Set maximum no. of rows for gaussian matrix.\n");
    printf("                     Too large matrixes should bee discarded for\n");
    printf("                     reasons of efficiency. Default: %d\n", gaussconfig.maxMatrixRows);
    printf("  --minmatrixrows  = [0 - 2^32-1] Set minimum no. of rows for gaussian matrix.\n");
    printf("                     Normally, too small matrixes are discarded for\n");
    printf("                     reasons of efficiency. Default: %d\n", gaussconfig.minMatrixRows);
    printf("  --savematrix     = [0 - 2^32-1] Save matrix every Nth decision level.\n");
    printf("                     Default: %d\n", gaussconfig.only_nth_gauss_save);
    printf("  --maxnummatrixes = [0 - 2^32-1] Maximum number of matrixes to treat.\n");
    printf("                     Default: %d\n", gaussconfig.maxNumMatrixes);
#endif //USE_GAUSS
    //printf("  --addoldlearnts  = Readd old learnts for failed variable searching.\n");
    //printf("                     These learnts are usually deleted, but may help\n");
    printf("  --nohyperbinres    Don't add binary clauses when doing failed lit probing.\n");
    printf("  --noremovebins     Don't remove useless binary clauses\n");
    printf("  --noremlbins       Don't remove useless learnt binary clauses\n");
    printf("  --nosubswithbins   Don't subsume with binary clauses\n");
    printf("  --nosubswithnbins  Don't subsume with non-existent binary clauses\n");
    printf("  --noclausevivif    Don't do perform clause vivification\n");
    printf("  --nosortwatched    Don't sort watches according to size: bin, tri, etc.\n");
    printf("  --nolfminim        Don't do on-the-fly self-subsuming resolution\n");
    printf("                     (called 'strong minimisation' in PrecoSat)\n");
    printf("  --nocalcreach      Don't calculate reachability and interfere with\n");
    printf("                     variable decisions accordingly\n");
    printf("  --nobxor           Don't find equivalent lits during failed lit search\n");
    printf("  --norecotfssr      Don't perform recursive/transitive OTF self-\n");
    printf("                     subsuming resolution\n");
    printf("  --nocacheotfssr    Don't cache 1-level equeue. Less memory used, but\n");
    printf("                     disables trans OTFSSR, adv. clause vivifier, etc.\n");
    printf("  --nootfsubsume     Don't do on-the-fly subsumption after conf. gen.\n");
#ifdef ENABLE_UNWIND_GLUE
    printf("  --maxgluedel       Automatically delete clauses over max glue. See '--maxglue'\n");
    printf("  --maxglue        = [0 - 2^%d-1] default: %d. Glue value above which we\n", MAX_GLUE_BITS, conf.maxGlue);
#endif //ENABLE_UNWIND_GLUE
    printf("                     throw the clause away on backtrack.\n");
    printf("  --threads        = Num threads (default is 1)\n");
    printf("  --plain            Get rid of all simplification algorithms\n");
    printf("  --maxconfl       = [0..2^63-1] Maximum number of conflicts to do\n");
    printf("  --maxtime        = [0..] Maximum number of seconds to run after which we exit cleanly\n");
    printf("  --switchoffsubs  = Number of variables after which to switch off subsumption and all related algorithms. Saves time. Default: %ld\n", conf.switch_off_subsumer_max_vars);
*/
    printf("\n");
}

const char* Main::hasPrefix(const char* str, const char* prefix) {
    int len = strlen(prefix);
    if (strncmp(str, prefix, len) == 0)
        return str + len;
    else
        return NULL;
}

void Main::parseCommandLine() {
    const char* value;
    char tmpFilename[201];
    tmpFilename[0] = '\0';
    uint32_t unparsedOptions = 0;
    bool needTwoFileNames = false;
    conf.verbosity = 0;
    need_clean_exit = false;

    for (int i = 0; i < argc; i++) {
        if ((value = hasPrefix(argv[i], "--polarity-mode="))) {
            if (strcmp(value, "true") == 0)
                conf.polarity_mode = polarity_true;
            else if (strcmp(value, "false") == 0)
                conf.polarity_mode = polarity_false;
            else if (strcmp(value, "rnd") == 0)
                conf.polarity_mode = polarity_rnd;
            else if (strcmp(value, "auto") == 0)
                conf.polarity_mode = polarity_auto;
            else {
                printf("ERROR! unknown polarity-mode %s\n", value);
                exit(0);
            }

        } else if ((value = hasPrefix(argv[i], "--rnd-freq="))) {
            double rnd;
            if (sscanf(value, "%lf", &rnd) <= 0 || rnd < 0 || rnd > 1) {
                printf("ERROR! illegal rnRSE ERROR!d-freq constant %s\n", value);
                exit(0);
            }
            conf.random_var_freq = rnd;

            /*} else if ((value = hasPrefix(argv[i], "--decay="))) {
                double decay;
                if (sscanf(value, "%lf", &decay) <= 0 || decay <= 0 || decay > 1) {
                    printf("ERROR! illegal decay constant %s\n", value);
                    exit(0);
                }
                conf.var_decay = 1 / decay;*/
        } else if ((value = hasPrefix(argv[i], "--epsilon="))) {
            float epsilon;

            if (sscanf(value, "%f", &epsilon) < 0 || epsilon > 1) {
                printf("ERROR! Illegal epsilon %s\n", value);
            }

            conf.epsilon = epsilon;
       }else if ((value = hasPrefix(argv[i], "--delta="))) {
            float delta;

            if (sscanf(value, "%f", &delta) < 0 || delta > 1) {
                printf("ERROR! Illegal delta %d\n", delta);
            }
            conf.delta = delta;

         }else if ((value = hasPrefix(argv[i], "--searchMode="))){
            int tmp;
            if (sscanf(value, "%d", &tmp) < 0 || tmp < 0) {
                   std::cout << "ERROR! illegal search mode: " << tmp << std::endl;
                  exit(0);
            }
            if (tmp ==0){
              conf.approxMCMode = LINEAR;
            }

       }else if ((value = hasPrefix(argv[i], "--startIteration="))) {
            int t;
            if (sscanf(value, "%d", &t) < 0) {
                printf("ERROR! Illegal startIteration %d\n", t);
            }
            conf.startIteration = t;
        } else if ((value = hasPrefix(argv[i], "--verbosity="))) {
            int verbosity = (int) strtol(value, NULL, 10);
            if (verbosity == EINVAL || verbosity == ERANGE) {
                printf("ERROR! illegal verbosity level %s\n", value);
                exit(0);
            }
            conf.verbosity = verbosity;
        } else if ((value = hasPrefix(argv[i], "--randomize="))) {
            int seed;
            if (sscanf(value, "%d", &seed) < 0) {
                printf("ERROR! illegal seed %s\n", value);
                exit(0);
            }
            conf.origSeed = seed;
        } else if ((value = hasPrefix(argv[i], "--restrict="))) {
            int branchTo;
            if (sscanf(value, "%d", &branchTo) < 0 || branchTo < 1) {
                printf("ERROR! illegal restricted pick branch number %d\n", branchTo);
                exit(0);
            }
            conf.restrictPickBranch = branchTo;
        } else if ((value = hasPrefix(argv[i], "--gaussuntil="))) {
            int until;
            if (sscanf(value, "%d", &until) < 0) {
                printf("ERROR! until %s\n", value);
                exit(0);
            }
            gaussconfig.decision_until = until;
        } else if ((value = hasPrefix(argv[i], "--restarts="))) {
            int maxrest;
            if (sscanf(value, "%d", &maxrest) < 0 || maxrest == 0) {
                printf("ERROR! illegal maximum restart number %d\n", maxrest);
                exit(0);
            }
            conf.maxRestarts = maxrest;
        } else if ((value = hasPrefix(argv[i], "--dumplearnts="))) {
            if (sscanf(value, "%200s", tmpFilename) < 0 || strlen(tmpFilename) == 0) {
                printf("ERROR! wrong filename '%s'\n", tmpFilename);
                exit(0);
            }
            conf.learntsFilename.assign(tmpFilename);
            conf.needToDumpLearnts = true;
        } else if ((value = hasPrefix(argv[i], "--dumporig="))) {
            if (sscanf(value, "%200s", tmpFilename) < 0 || strlen(tmpFilename) == 0) {
                printf("ERROR! wrong filename '%s'\n", tmpFilename);
                exit(0);
            }
            conf.origFilename.assign(tmpFilename);
            conf.needToDumpOrig = true;
        } else if ((value = hasPrefix(argv[i], "--logFile="))) {
            if (sscanf(value, "%200s", tmpFilename) < 0 || strlen(tmpFilename) == 0) {
                printf("ERROR! wrong fileName '%s'\n", tmpFilename);
                exit(0);
            }
            conf.shouldLog = true;
            conf.logFilename.assign(tmpFilename);
        } else if ((value = hasPrefix(argv[i], "--alsoread="))) {
            if (sscanf(value, "%400s", tmpFilename) < 0 || strlen(tmpFilename) == 0) {
                printf("ERROR! wrong filename '%s'\n", tmpFilename);
                exit(0);
            }
            filesToRead.push_back(tmpFilename);
        } else if ((value = hasPrefix(argv[i], "--maxdumplearnts="))) {
            if (!conf.needToDumpLearnts) {
                printf("ERROR! -dumplearnts=<filename> must be first activated before issuing -maxdumplearnts=<size>\n");
                exit(0);
            }
            int tmp;
            if (sscanf(value, "%d", &tmp) < 0 || tmp < 0) {
                std::cout << "ERROR! wrong maximum dumped learnt clause size is illegal: " << tmp << std::endl;
                exit(0);
            }
            conf.maxDumpLearntsSize = (uint32_t) tmp;
        } else if ((value = hasPrefix(argv[i], "--maxsolutions="))) {
            int tmp;
            if (sscanf(value, "%d", &tmp) < 0 || tmp < 0) {
                std::cout << "ERROR! wrong maximum number of solutions is illegal: " << tmp << std::endl;
                exit(0);
            }
            max_nr_of_solutions = (uint32_t) tmp;

        } else if ((value = hasPrefix(argv[i], "--pavgbranch"))) {
            conf.doPrintAvgBranch = true;
        } else if ((value = hasPrefix(argv[i], "--greedyunbound"))) {
            conf.greedyUnbound = true;
        } else if ((value = hasPrefix(argv[i], "--nonormxorfind"))) {
            conf.doFindXors = false;
        } else if ((value = hasPrefix(argv[i], "--nobinxorfind"))) {
            conf.doFindEqLits = false;
        } else if ((value = hasPrefix(argv[i], "--noregbxorfind"))) {
            conf.doRegFindEqLits = false;
        } else if ((value = hasPrefix(argv[i], "--doextendedscc"))) {
            conf.doExtendedSCC = true;
        } else if ((value = hasPrefix(argv[i], "--noconglomerate"))) {
            conf.doConglXors = false;
        } else if ((value = hasPrefix(argv[i], "--nosimplify"))) {
            conf.doSchedSimp = false;
        } else if ((value = hasPrefix(argv[i], "--debuglib"))) {
            debugLib = true;
        } else if ((value = hasPrefix(argv[i], "--debugnewvar"))) {
            debugNewVar = true;
        } else if ((value = hasPrefix(argv[i], "--novarreplace"))) {
            conf.doReplace = false;
        } else if ((value = hasPrefix(argv[i], "--nofailedlit"))) {
            conf.doFailedLit = false;
        } else if ((value = hasPrefix(argv[i], "--nodisablegauss"))) {
            gaussconfig.dontDisable = true;
        } else if ((value = hasPrefix(argv[i], "--maxnummatrixes="))) {
            int maxNumMatrixes;
            if (sscanf(value, "%d", &maxNumMatrixes) < 0) {
                printf("ERROR! maxnummatrixes: %s\n", value);
                exit(0);
            }
            gaussconfig.maxNumMatrixes = maxNumMatrixes;
        } else if ((value = hasPrefix(argv[i], "--noheuleprocess"))) {
            conf.doHeuleProcess = false;
        } else if ((value = hasPrefix(argv[i], "--nosatelite"))) {
            conf.doSatELite = false;
        } else if ((value = hasPrefix(argv[i], "--noxorsubs"))) {
            conf.doXorSubsumption = false;
        } else if ((value = hasPrefix(argv[i], "--nohyperbinres"))) {
            conf.doHyperBinRes = false;
        } else if ((value = hasPrefix(argv[i], "--novarelim"))) {
            conf.doVarElim = false;
        } else if ((value = hasPrefix(argv[i], "--nosubsume1"))) {
            conf.doSubsume1 = false;
        } else if ((value = hasPrefix(argv[i], "--nomatrixfind"))) {
            gaussconfig.noMatrixFind = true;
        } else if ((value = hasPrefix(argv[i], "--noiterreduce"))) {
            gaussconfig.iterativeReduce = false;
        } else if ((value = hasPrefix(argv[i], "--noordercol"))) {
            gaussconfig.orderCols = false;
        } else if ((value = hasPrefix(argv[i], "--maxmatrixrows="))) {
            int rows;
            if (sscanf(value, "%d", &rows) < 0 || rows < 0) {
                printf("ERROR! maxmatrixrows: %s\n", value);
                exit(0);
            }
            gaussconfig.maxMatrixRows = (uint32_t) rows;
        } else if ((value = hasPrefix(argv[i], "--minmatrixrows="))) {
            int rows;
            if (sscanf(value, "%d", &rows) < 0 || rows < 0) {
                printf("ERROR! minmatrixrows: %s\n", value);
                exit(0);
            }
            gaussconfig.minMatrixRows = rows;
        } else if ((value = hasPrefix(argv[i], "--savematrix"))) {
            int every;
            if (sscanf(value, "%d", &every) < 0) {
                printf("ERROR! savematrix: %s\n", value);
                exit(0);
            }
            gaussconfig.only_nth_gauss_save = every;
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "-help") == 0 || strcmp(argv[i], "--help") == 0) {
            printUsage(argv);
            exit(0);
        } else if ((value = hasPrefix(argv[i], "--restart="))) {
            if (strcmp(value, "auto") == 0)
                conf.fixRestartType = auto_restart;
            else if (strcmp(value, "static") == 0)
                conf.fixRestartType = static_restart;
            else if (strcmp(value, "dynamic") == 0)
                conf.fixRestartType = dynamic_restart;
            else {
                printf("ERROR! unknown restart type %s\n", value);
                exit(0);
            }
        } else if ((value = hasPrefix(argv[i], "--nosolprint"))) {
            printResult = false;
            //} else if ((value = hasPrefix(argv[i], "--addoldlearnts"))) {
            //    conf.readdOldLearnts = true;
        } else if ((value = hasPrefix(argv[i], "--nohyperbinres"))) {
            conf.doHyperBinRes = false;
        } else if ((value = hasPrefix(argv[i], "--noremovebins"))) {
            conf.doRemUselessBins = false;
        } else if ((value = hasPrefix(argv[i], "--nosubswithnbins"))) {
            conf.doSubsWNonExistBins = false;
        } else if ((value = hasPrefix(argv[i], "--nosubswithbins"))) {
            conf.doSubsWBins = false;
        } else if ((value = hasPrefix(argv[i], "--noclausevivif"))) {
            conf.doClausVivif = false;
        } else if ((value = hasPrefix(argv[i], "--nosortwatched"))) {
            conf.doSortWatched = false;
        } else if ((value = hasPrefix(argv[i], "--nolfminim"))) {
            conf.doMinimLearntMore = false;
        } else if ((value = hasPrefix(argv[i], "--nocalcreach"))) {
            conf.doCalcReach = false;
        } else if ((value = hasPrefix(argv[i], "--norecotfssr"))) {
            conf.doMinimLMoreRecur = false;
        } else if ((value = hasPrefix(argv[i], "--nocacheotfssr"))) {
            conf.doCacheOTFSSRSet = false;
            conf.doCacheOTFSSR = false;
        } else if ((value = hasPrefix(argv[i], "--nootfsubsume"))) {
            conf.doOTFSubsume = false;
        } else if ((value = hasPrefix(argv[i], "--noremlbins"))) {
            conf.doRemUselessLBins = false;
        } else if ((value = hasPrefix(argv[i], "--maxconfl="))) {
            int maxconfl = 0;
            if (sscanf(value, "%d", &maxconfl) < 0 || maxconfl < 2) {
                printf("ERROR! max confl: %s\n", value);
                exit(-1);
            }
            conf.maxConfl = maxconfl;
        } else if ((value = hasPrefix(argv[i], "--maxTotalTime="))) {
            int maxtime = 0;
            if (sscanf(value, "%d", &maxtime) < 0 || maxtime < 2) {
                printf("ERROR! max time is too small: %s\n", value);
                exit(-1);
            }
            conf.totalTimeout = maxtime;
        } else if ((value = hasPrefix(argv[i], "--maxLoopTime="))) {
            int maxtime = 0;
            if (sscanf(value, "%d", &maxtime) < 0 || maxtime < 2) {
                printf("ERROR! max time is too small: %s\n", value);
                exit(-1);
            }
            conf.loopTimeout = maxtime;
        }else if ((value = hasPrefix(argv[i], "--switchoffsubs="))) {
            long vars = 0;
            if (sscanf(value, "%ld", &vars) < 0 || vars < 2) {
                printf("ERROR! max time is too small: %s\n", value);
                exit(-1);
            }
            conf.switch_off_subsumer_max_vars = vars;
        } else if ((value = hasPrefix(argv[i], "--plain"))) {
            conf.isPlain = true;
            conf.doOTFSubsume = false;
            conf.doFindXors = false;
            conf.doFindEqLits = false;
            conf.doRegFindEqLits = false;
            conf.doExtendedSCC = false;
            conf.doConglXors = false;
            conf.doSchedSimp = false;
            conf.doReplace = false;
            conf.doFailedLit = false;
            conf.doHeuleProcess = false;
            conf.doSatELite = false;
            conf.doXorSubsumption = false;
            conf.doVarElim = false;
            //nomatrixfind
            gaussconfig.orderCols = false;
            gaussconfig.iterativeReduce = false;
            conf.doHyperBinRes = false;
            conf.doRemUselessBins = false;
            conf.doRemUselessLBins = false;
            conf.doSubsWBins = false;
            conf.doSubsWNonExistBins = false;
            conf.doClausVivif = false;
            conf.doCalcReach = false;
            conf.doBXor = false;
            conf.doMinimLMoreRecur = false;
            conf.doMinimLearntMore = false;
            conf.doCacheOTFSSR = false;
        } else if ((value = hasPrefix(argv[i], "--nobxor"))) {
            conf.doBXor = false;
#ifdef ENABLE_UNWIND_GLUE
        } else if ((value = hasPrefix(argv[i], "--maxglue="))) {
            int glue = 0;
            if (sscanf(value, "%d", &glue) < 0 || glue < 2) {
                printf("ERROR! maxGlue: %s\n", value);
                exit(0);
            }
            if (glue >= (1 << MAX_GLUE_BITS) - 1) {
                std::cout << "Due to memory-packing limitations, max glue cannot be more than "
                        << ((1 << MAX_GLUE_BITS) - 2) << std::endl;
                exit(-1);
            }
            conf.maxGlue = (uint32_t) glue;
        } else if ((value = hasPrefix(argv[i], "--maxgluedel"))) {
            conf.doMaxGlueDel = true;
#endif //ENABLE_UNWIND_GLUE
        } else if ((value = hasPrefix(argv[i], "--threads="))) {
            numThreads = 0;
            if (sscanf(value, "%d", &numThreads) < 0 || numThreads < 1) {
                printf("ERROR! numThreads: %s\n", value);
                exit(0);
            }
        } else if (strncmp(argv[i], "-", 1) == 0 || strncmp(argv[i], "--", 2) == 0) {
            printf("ERROR! unknown flag %s\n", argv[i]);
            exit(0);
        } else {
            //std::std::cout << "argc:" << argc << " i:" << i << ", value:" << argv[i] << std::endl;
            unparsedOptions++;
            if (unparsedOptions == 2) {
                if (!(argc <= i + 2)) {
                    std::cout << "You must give the input file as either:" << std::endl;
                    std::cout << " -- last option if you want the output to the console" << std::endl;
                    std::cout << " -- or one before the last option" << std::endl;
                    std::cout << "It appears that you did neither. Maybe you forgot the '--' from an option?" << std::endl;
                    exit(-1);
                }
                fileNamePresent = true;
                if (argc == i + 2) needTwoFileNames = true;
            }
            if (unparsedOptions == 3) {
                if (!(argc <= i + 1)) {
                    std::cout << "You must give the output file as the last option. Exiting" << std::endl;
                    exit(-1);
                }
                twoFileNamesPresent = true;
            }
            if (unparsedOptions == 4) {
                std::cout << "You gave more than two filenames as parameters." << std::endl;
                std::cout << "The first one is interpreted as the input, the second is the output." << std::endl;
                std::cout << "However, the third one I cannot do anything with. EXITING" << std::endl;
                exit(-1);
            }
        }
    }
    if (conf.verbosity >= 1) {
        if (twoFileNamesPresent) {
            std::cout << "c Outputting solution to file: " << argv[argc - 1] << std::endl;
        } else {
            std::cout << "c Outputting solution to console" << std::endl;
        }
    }


    if (unparsedOptions == 2 && needTwoFileNames == true) {
        std::cout << "Command line wrong. You probably frogot to add " << std::endl
                << "the '--'  in front of one of the options, or you started" << std::endl
                << "your output file with a hyphen ('-'). Exiting." << std::endl;
        exit(-1);
    }
    if (!debugLib) conf.libraryUsage = false;
}

FILE* Main::openOutputFile() {
    FILE* res = NULL;
    if (twoFileNamesPresent) {
        char* filename = argv[argc - 1];
        res = fopen(filename, "wb");
        if (res == NULL) {
            int backup_errno = errno;
            printf("Cannot open %s for writing. Problem: %s", filename, strerror(backup_errno));
            exit(1);
        }
    }

    return res;
}

FILE* Main::openLogFile() {
    FILE* res = NULL;
    if (!conf.shouldLog) {
        return res;
    }
    res = fopen(conf.logFilename.c_str(), "wb");
    if (res == NULL) {
        int backup_errno = errno;
        printf("Cannot open %s for writing. Problem: %s\n", conf.logFilename.c_str(), strerror(backup_errno));
        exit(1);
    }
    return res;
}

void Main::setDoublePrecision(const uint32_t verbosity) {
#if defined(_FPU_EXTENDED) && defined(_FPU_DOUBLE)
    fpu_control_t oldcw, newcw;
    _FPU_GETCW(oldcw);
    newcw = (oldcw & ~_FPU_EXTENDED) | _FPU_DOUBLE;
    _FPU_SETCW(newcw);
#pragma omp single
    if (verbosity >= 1) {
        printf("c WARNING: for repeatability, setting FPU to use double precision\n");
    }
#endif
}

void Main::printVersionInfo(const uint32_t verbosity) {
#pragma omp single
    if (verbosity >= 1) {
        printf("c This is ApproxMC %s\n", VERSION);
#ifdef __GNUC__
        printf("c compiled with gcc version %s\n", __VERSION__);
#else
        printf("c compiled with non-gcc compiler\n");
#endif
    }
}

int Main::correctReturnValue(const lbool ret) const {
    int retval = -1;
    if (ret == l_True) retval = 10;
    else if (ret == l_False) retval = 20;
    else if (ret == l_Undef) retval = 15;
    else {
        std::cerr << "Something is very wrong, output is neither l_Undef, nor l_False, nor l_True" << std::endl;
        exit(-1);
    }

#ifdef NDEBUG
    // (faster than "return", which will invoke the destructor for 'Solver')
    exit(retval);
#endif
    return retval;
}

std::string binary(int x, uint32_t length) {
    uint32_t logSize = log2(x) + 1;
    std::string s;
    do {
        s.push_back('0' + (x & 1));
    } while (x >>= 1);
    for (uint32_t i = logSize; i < (uint32_t) length; i++) {
        s.push_back('0');
    }
    std::reverse(s.begin(), s.end());

    return s;

}

bool Main::GenerateRandomBits(string &randomBits, uint32_t size) {
    std::uniform_int_distribution<int> uid{0, 2147483647};
    uint32_t i = 0;
    while (i < size) {
        i += 31;
        randomBits += binary(uid(rd), 31);
    }
    return true;
}

int Main::GenerateRandomNum(int maxRange) {
    std::uniform_int_distribution<int> uid(0, maxRange);
    int value = -1;
    while (value < 0 || value > maxRange) {
        value = uid(rd);
    }
    return value;
}

bool Main::AddHash(uint32_t clausNum, std::map<int,Lit> &hashAssumptionVars,Solver &s, vec<Lit> &assumptions){
    string randomBits;
    GenerateRandomBits(randomBits, (s.independentSet.size() + 1) * clausNum);
    bool xorEqualFalse = false;
    Var activationVar;
    vec<Lit> lits;
    for (uint32_t i = 0; i < clausNum; i++) {
        lits.clear();
        activationVar = s.newVar();
        assumptions.push(Lit(activationVar, true));
       
        hashAssumptionVars[assumptions.size()-1] = Lit(activationVar,true);
        lits.push(Lit(activationVar, false));
        xorEqualFalse = (randomBits[(s.independentSet.size() + 1) * i] == 1);

        for (uint32_t j = 0; j < s.independentSet.size(); j++) {

            if (randomBits[(s.independentSet.size() + 1) * i + j] == '1') {
                lits.push(Lit(s.independentSet[j], true));
            }
        }
        s.addXorClause(lits, xorEqualFalse);
    }
    return true;
}


bool Main::SetHash(uint32_t clausNum, std::map<int,Lit> &hashAssumptionVars, Solver &s, vec<Lit> &assumptions){
    if (clausNum < assumptions.size())
    {
      uint32_t numberToRemove = assumptions.size()- clausNum;
      for (uint32_t i = 0; i<numberToRemove; i++){
        assumptions.pop();
      }
    }
    
    else{
      if (clausNum > assumptions.size() && assumptions.size() < hashAssumptionVars.size())
    {
        for (uint32_t i = 0; i<hashAssumptionVars.size()-assumptions.size(); i++)
        {
            assumptions.push(hashAssumptionVars[i]);
        }
    }
    if (clausNum > hashAssumptionVars.size())
    {
      return AddHash(clausNum-hashAssumptionVars.size(),hashAssumptionVars,s,assumptions);
    }
    }
   return true; 
}
void Main::printResultFunc(Solver &S, vec<lbool> solutionModel, const lbool ret, FILE* res) {
    if (res != NULL && printResult) {
        if (ret == l_True) {
            fprintf(res, "v ");
            for (Var var = 0; var != S.nOrigVars(); var++)
                if (solutionModel[var] != l_Undef)
                    fprintf(res, "%s%d ", (S.model[var] == l_True) ? "" : "-", var + 1);
            fprintf(res, "0\n");
            fflush(res);
        }
    } else {

        if (ret == l_True && printResult) {
            std::stringstream toPrint;
            toPrint << "v ";
            for (Var var = 0; var != S.nOrigVars(); var++)
                if (solutionModel[var] != l_Undef)
                    toPrint << ((solutionModel[var] == l_True) ? "" : "-") << var + 1 << " ";
            toPrint << "0" << std::endl;
            std::cout << toPrint.str();
        }
    }
}

int32_t Main::BoundedSATCount(uint32_t maxSolutions, Solver &solver, vec<Lit> &assumptions) {
    unsigned long current_nr_of_solutions = 0;
    lbool ret = l_True;
    Var activationVar = solver.newVar();
    vec<Lit> allSATAssumptions;
    if (!assumptions.empty()) {
        assumptions.copyTo(allSATAssumptions);
    }
    allSATAssumptions.push(Lit(activationVar, true));
    //signal(SIGALRM, SIGALARM_handler);
     start_timer(conf.loopTimeout);
    while (current_nr_of_solutions < maxSolutions && ret == l_True) {

        ret = solver.solve(allSATAssumptions);
        current_nr_of_solutions++;
        if (ret == l_True && current_nr_of_solutions < maxSolutions) {
            vec<Lit> lits;
            lits.push(Lit(activationVar, false));
            for (uint32_t j = 0; j < solver.independentSet.size(); j++) {
                Var var = solver.independentSet[j];
                if (solver.model[var] != l_Undef) {
                    lits.push(Lit(var, (solver.model[var] == l_True) ? true : false));
                }
            }
            solver.addClause(lits);
        }
    }
    vec<Lit> cls_that_removes;
    cls_that_removes.push(Lit(activationVar, false));
    solver.addClause(cls_that_removes);
    if (ret == l_Undef){
        solver.needToInterrupt = false;
        return -1*current_nr_of_solutions;
    }
    return current_nr_of_solutions;
}


SATCount Main::ApproxMC(Solver &solver, FILE* resLog) {
    int32_t currentNumSolutions = 0;
    uint32_t  hashCount = 0, swapVar = 0, mPrev = 0, hashPrev = 0;
    std::list<int> numHashList, numCountList, medianComputeList;
    vec<Lit> assumptions;
    SATCount solCount;
    solCount.cellSolCount = 0;
    solCount.hashCount = 0;
    double elapsedTime = 0;
    int repeatTry = 0,minHash=0,medSolCount=0;
    map<int,int> succRecord, countRecord;
    std::map<int,Lit> hashVars;
    uint32_t numExplored = 1, lowerFib = 0, upperFib = solver.independentSet.size();
    currentNumSolutions = BoundedSATCount(conf.pivotApproxMC+1, solver, assumptions);
    if (currentNumSolutions <= conf.pivotApproxMC){
        solCount.cellSolCount = currentNumSolutions;
        solCount.hashCount = 0;
        return solCount;
    }
    hashCount = 1;
    for (uint32_t j = 0; j < conf.tApproxMC; j++) {
        Solver solver(conf, gaussconfig);
    solversToInterrupt.clear();
    solversToInterrupt[0] = &solver;
    need_clean_exit = true;
    setDoublePrecision(conf.verbosity);
    parseInAllFiles(solver);
 
        hashVars.clear();
        succRecord.clear();
        countRecord.clear();
        hashPrev = 0;
        lowerFib = 0;
        if (conf.approxMCMode == LINEAR){
          hashCount = 1;
        }
        numExplored = 1; 
        upperFib = solver.independentSet.size();
        while (numExplored < solver.independentSet.size())
        {
            if (conf.approxMCMode == LINEAR){
                numExplored = hashCount;
            }
            double currentTime = totalTime();
            elapsedTime = currentTime-startTime;
            if (elapsedTime > conf.totalTimeout - 3000){
                break;
            }
            swapVar = hashCount;
            double myTime = totalTime();
            SetHash(hashCount, hashVars, solver,assumptions);
            //BoundedSATCount returns -1 for timeout, 0 for UNSAT
            currentNumSolutions = BoundedSATCount(conf.pivotApproxMC + 1, solver, assumptions);
            myTime = totalTime() - myTime;
            //printf("%f\n", myTime);
            //printf("%d %d\n",currentNumSolutions,conf.pivotApproxMC);
            if (conf.shouldLog) {
                fprintf(resLog, "ApproxMC:%d:%d:%d:%f:%d:%d\n",conf.approxMCMode, j, hashCount, myTime,
                    (currentNumSolutions == (int32_t)(conf.pivotApproxMC + 1)),currentNumSolutions);
                fflush(resLog);
            }
            if (currentNumSolutions <= -1){
                assumptions.clear();
                hashVars.clear();
                if (repeatTry < 2){     /* Retry up to twice more */
                    AddHash(hashCount, hashVars, solver,assumptions);
                    repeatTry += 1;
                } else{
                    repeatTry = 0;
                    hashCount++;
                }
            }
            else if (currentNumSolutions <= conf.pivotApproxMC){
                if (conf.approxMCMode == LINEAR){
                  numHashList.push_back(hashCount);
                  numCountList.push_back(currentNumSolutions);
                  break;
                }
                numExplored = lowerFib+solver.independentSet.size()-hashCount;
                if (succRecord.find(hashCount-1) != succRecord.end()){
                  if (succRecord[hashCount-1] == 1){
                    numHashList.push_back(hashCount);
                    numCountList.push_back(currentNumSolutions);
                    mPrev = hashCount;
                    break;
                  }
                }
                //printf("mPrev:%d\n",mPrev);
                succRecord[hashCount] = 0;
                countRecord[hashCount] = currentNumSolutions;
                if (std::abs((int)hashCount -(int)mPrev) <= 2 && mPrev != 0){
                    upperFib = hashCount;
                    hashCount --;
                }
                else{
                    if (hashPrev > hashCount){
                      hashPrev = 0;
                    }
                    upperFib = hashCount;
                    lowerFib = hashPrev;
                    hashCount = (upperFib+lowerFib)/2;
                }
            }
            else if (currentNumSolutions == conf.pivotApproxMC + 1) {
                if (conf.approxMCMode == LINEAR){
                    assumptions.clear();
                    hashVars.clear();
                    hashCount ++;
                    continue;
                  }
                numExplored = hashCount + solver.independentSet.size()-upperFib;
                if (succRecord.find(hashCount+1) != succRecord.end()){
                    if (succRecord[hashCount+1] == 0){
                    numHashList.push_back(hashCount+1);
                    numCountList.push_back(countRecord[hashCount+1]);
                    mPrev = hashCount+1;
                    break;
                  }
                }
                succRecord[hashCount] = 1;
                if (std::abs((int)hashCount - (int)mPrev) < 2 && mPrev!=0){
                    lowerFib = hashCount;
                    hashCount ++;
                }
                else if (lowerFib + (hashCount - lowerFib)*2 >= upperFib-1){
                    lowerFib = hashCount;
                    hashCount = (lowerFib+upperFib)/2;
                }else{
                    //printf("hashPrev:%d hashCount:%d\n",hashPrev, hashCount);
                    hashCount = lowerFib + (hashCount -lowerFib)*2;
                }
            }
          
          hashPrev = swapVar;
        }
	minHash = findMin(numHashList);
	medianComputeList.clear();
	for (std::list<int>::iterator it1 = numHashList.begin(), it2 = numCountList.begin();
	     it1 != numHashList.end() && it2 != numCountList.end(); it1++, it2++) {
	  	medianComputeList.push_back((*it2)* pow(2, (*it1) - minHash));
	}
	medSolCount = findMedian(medianComputeList);
	if (j%2 == 0){
	printf(" With confidence %f, Esimate of the number of solutions is: %d * 2^%d\n",
                confidence[j+1], medSolCount, minHash);
	}
	assumptions.clear();
        if (elapsedTime > conf.totalTimeout - 3000){
            break;
        }
        hashCount = mPrev;
    }
    if (numHashList.size() == 0){
        return solCount;
    }

    solCount.cellSolCount = medSolCount;
    solCount.hashCount = minHash;
    return solCount;
}

bool Main::printSolutions(Solver& solver, FILE* res){
    for (map< std::string, uint32_t>:: iterator it = solutionMap.begin();
                                    it != solutionMap.end(); it++)
            fprintf(res, "%s:%d\n ", it->first.c_str(),it->second);
            
    fflush(res);
     
    return true;
}
int Main::singleThreadSolve() {

   int val; 
   std::string line;
   char * pch;

   std::ifstream probmapfile;
   if (conf.approxMCMode == LINEAR)
     {
       probmapfile.open("./ProbMapFile_40.txt");
     }
   else
  {
     probmapfile.open("./ProbMapFile_36.txt");

  }
   if (probmapfile.is_open()){
     while (getline(probmapfile, line)) {
       pch = strtok(strdup(line.c_str()), ":");
       val = std::atoi(pch);
       pch = strtok(NULL, ":");
       if (std::atof(pch) > (1 - conf.delta))
	        break;
     }
   }
   probmapfile.close();
   if (val == 0) {
     std::cerr << "Probmapfile failed, delta is " << conf.delta << std::endl;
     exit(-1);
   }
   conf.tApproxMC = val;
   std::cout<<"t ApproxMC:"<<conf.tApproxMC<<std::endl;
     confidence = (float *) malloc(sizeof(float)*(2+conf.tApproxMC));
   if(confidence == NULL){  
     printf("Out of memory, could not allocate confidence list\n");  
     exit(-1);  
   }  
   
   
   if (conf.approxMCMode == LINEAR)
     {
       probmapfile.open("./ProbMapFile_40.txt");
     }
   else{
     probmapfile.open("./ProbMapFile_36.txt");
   }
   if (probmapfile.is_open()){
     int iter = 0;
     while (getline(probmapfile, line) && iter < conf.tApproxMC) {
       pch = strtok(strdup(line.c_str()), ":");
       iter = std::atoi(pch);
       pch = strtok(NULL, ":");
       confidence[iter] = std::atof(pch);
       //  iter++;
     }
   }
   probmapfile.close();
  if (conf.approxMCMode == LINEAR){
    conf.pivotApproxMC = 2*ceil(4.94*(1+(1/conf.epsilon))*(1+(1/conf.epsilon)));
  }
  else
  {
    conf.pivotApproxMC = int(1 + 9.84*(1+(1/conf.epsilon))*(1+(1/conf.epsilon))*(1+(conf.epsilon/(1+conf.epsilon))));
  }
    sampleCounter = 0;
    startTime = totalTime();
    timerSetFirstTime = true;
    Solver solver(conf, gaussconfig);
    solversToInterrupt.clear();
    solversToInterrupt[0] = &solver;
    need_clean_exit = true;
    printVersionInfo(conf.verbosity);
    setDoublePrecision(conf.verbosity);
    parseInAllFiles(solver);
    FILE* res = openOutputFile();
    FILE* resLog = openLogFile();
    lbool ret = l_True;
    if (conf.startIteration > solver.independentSet.size()){
        conf.startIteration = 0;
        string filename = argv[(twoFileNamesPresent ? argc - 2 : argc - 1)];
        std::cerr<<filename<<std::endl;
    }

    SATCount solCount;
    solCount = ApproxMC(solver, resLog);
    std::cout<<"With confidence at least:"<<1-conf.delta<<"  "<<"Approximate count is :"<<solCount.cellSolCount<<"x2^"<<solCount.hashCount<<std::endl;
    if (conf.needToDumpOrig) {
        if (ret != l_False) {
            solver.addAllXorAsNorm();
        }
        if (ret == l_False && conf.origFilename == "stdout") {
            std::cout << "p cnf 0 1" << std::endl;
            std::cout << "0";
        } else if (ret == l_True && conf.origFilename == "stdout") {
            std::cout << "p cnf " << solver.model.size() << " " << solver.model.size() << std::endl;
            for (uint32_t i = 0; i < solver.model.size(); i++) {
                std::cout << (solver.model[i] == l_True ? "" : "-") << i + 1 << " 0" << std::endl;
            }
        } else {
            if (!solver.dumpOrigClauses(conf.origFilename)) {
                std::cout << "Error: Cannot open file '" << conf.origFilename << "' to write learnt clauses!" << std::endl;
                exit(-1);
            }
            if (conf.verbosity >= 1)
                std::cout << "c Simplified original clauses dumped to file '"
                    << conf.origFilename << "'" << std::endl;
        }
    }
    if (ret == l_Undef && conf.verbosity >= 1) {
        std::cout << "c Not finished running -- signal caught or maximum restart reached" << std::endl;
    }
    if (conf.verbosity >= 1) solver.printStats();

    // printResultFunc(solver, ret, res, current_nr_of_solutions == 1);

    return correctReturnValue(ret);
}

/**
@brief For correctly and gracefully exiting

It can happen that the user requests a dump of the learnt clauses. In this case,
the program must wait until it gets to a state where the learnt clauses are in
a correct state, then dump these and quit normally. This interrupt hander
is used to achieve this
 */

int main(int argc, char** argv) {
    Main main(argc, argv);
    main.parseCommandLine();
    signal(SIGINT, SIGINT_handler);
    signal(SIGALRM, SIGALARM_handler);
    try{
        return main.singleThreadSolve();

    }

    catch(std::bad_alloc) {
        std::cerr << "Memory manager cannot handle the load. Sorry. Exiting." << std::endl;
        exit(-1);
    }

    catch(std::out_of_range oor) {
        std::cerr << oor.what() << std::endl;
        exit(-1);
    }

    catch(CMSat::DimacsParseError dpe) {
        std::cerr << "PARSE ERROR!" << dpe.what() << std::endl;
        exit(3);
    }
    return 0;
}
