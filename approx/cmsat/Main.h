/***************************************************************************************[Solver.cc]
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2009, Niklas Sorensson
Copyright (c) 2009-2012, Mate Soos

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

#ifndef MAIN_H
#define MAIN_H

#include<iostream>
#include <sstream>
#include<random>
#include <string>
#include <vector>
#include <map>
#ifndef DISABLE_ZLIB
#include <zlib.h>
#endif // DISABLE_ZLIB
#include <chrono>
#include <ctime>
#include "cmsat/Solver.h"
#include "cmsat/SharedData.h"
#include "cmsat/DimacsParser.h"
namespace CMSat {

    using std::string;

    struct SATCount {
        uint32_t hashCount;
        uint32_t cellSolCount;
    };

    class Main {
    public:
        Main(int argc, char** argv);

        void parseCommandLine();

        int singleThreadSolve();
        int oneThreadSolve();
        int multiThreadSolve();

        int numThreads;

    private:

        void printUsage(char** argv);
        const char* hasPrefix(const char* str, const char* prefix);
        void printResultFunc(const Solver& S, const lbool ret, FILE* res, const bool firstSolution);

        //File reading
        void readInAFile(const std::string& filename, Solver& solver);
        void readInStandardInput(Solver& solver);
        void parseInAllFiles(Solver& solver);
        FILE* openOutputFile();
        FILE* openLogFile();

        int singleThreadUniGenCall(uint32_t samples, FILE* res, FILE* resLog);
        void setDoublePrecision(const uint32_t verbosity);
        void printVersionInfo(const uint32_t verbosity);
        int correctReturnValue(const lbool ret) const;

        SolverConf conf;
        GaussConf gaussconfig;

        float *confidence;

        bool grouping;
        bool debugLib;
        bool debugNewVar;
        bool printResult;
        uint32_t max_nr_of_solutions;
        bool fileNamePresent;
        bool twoFileNamesPresent;
        std::vector<std::string> filesToRead;

        SharedData sharedData;

        int argc;
        char** argv;
         SATCount ApproxMC(Solver &solver, FILE* resLog);
        bool AddHash(uint32_t clausNum, std::map<int,Lit>&hashAssumptionVars,Solver &s, vec<Lit> &assumptions);
        bool SetHash(uint32_t clausNum, std::map<int,Lit>&hashAssumptionVars, Solver &s, vec<Lit> &assumptions);
        int32_t BoundedSATCount(uint32_t maxSolutions, Solver &solver, vec<Lit> &assumptions);
         bool GenerateRandomBits(string &randomBits, uint32_t size);
        int GenerateRandomNum(int maxRange);
        void printResultFunc(Solver &S, vec<lbool> solutionModel, const lbool ret, FILE* res);
        vec<Lit> lits; ///<To reduce temporary creation overhead    
        vector<vec<lbool> > modelsSet;
        vec<lbool> model;
        vec<Lit> assumptions;
        vec<Lit> allSATAssumptions;
        time_t  startTime;
        uint32_t sampleCounter;
        std::map< std::string, uint32_t> solutionMap;  
        bool printSolutions(Solver &s, FILE* res);
        std::random_device rd;
        void SeedEngine(std::mt19937 &randomEngine);
    };

}

#endif //MAIN_H
