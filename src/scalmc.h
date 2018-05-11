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


    uint32_t ScalGen();
    string GenerateRandomBits(const uint32_t size);
    string binary(const uint32_t x, const uint32_t length);
    uint32_t SolutionsToReturn(uint32_t numSolutions);
    void generate_samples();
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
    std::string sampleFilename;
    uint32_t startIterationUG = 0;
    uint32_t multisample = 1;
    uint32_t loThresh;
    uint32_t hiThresh;
    bool sparse = false;
    double kappa = 0.638;      /* Corresponds to epsilon=16 */
    uint32_t samples = 0;
    uint32_t callsPerSolver = 0;
    uint32_t pivotUniGen; //TODO rename scalgen
    int scalmc = 1;

private:
    SolverConf conf;
    bool count(SATCount& count);
    void add_scalmc_options();
    bool ScalScalMC(SATCount& count);
    bool AddHash(uint32_t num_xor_cls, vector<Lit>& assumps);
    void SetHash(uint32_t clausNum, std::map<uint64_t,Lit>& hashVars, vector<Lit>& assumps);

    void printVersionInfo() const;
    int correctReturnValue(const lbool ret) const;

    int64_t BoundedSATCount(
        uint32_t maxSolutions,
        uint32_t minSolutions,
        const vector<Lit>& assumps,
        const uint32_t hashCount,
        std::map<std::string, uint32_t>* solutionMap = NULL
    );

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
    std::ofstream cusp_logf;
    std::mt19937 randomEngine;
    int learn_type = 0;
    int maple = 1;
    uint32_t num_threads = 1;
    SATSolver* solver;
    vector<uint32_t> independent_vars;
    unsigned verb = 1;
    unsigned verb_scalmc_cls = 0;
    double total_runtime; //runTime
    int what_to_break = 0;
    int dosimp = 1;

    int argc;
    char** argv;
};

std::array<double,256> iterationConfidences = {{
    0.64, 0.704512, 0.7491026944, 0.783348347699,
    0.81096404252, 0.833869604432, 0.853220223135, 0.869779929746,
    0.884087516258, 0.896540839559, 0.907443973174, 0.917035558936,
    0.92550684748, 0.933013712405, 0.939684956024, 0.945628233538,
    0.950934391703, 0.955680718969, 0.9599334282, 0.963749585636,
    0.967178632062, 0.970263598173, 0.973042086923, 0.975547075737,
    0.977807577642, 0.979849190627, 0.98169455749, 0.98336375333,
    0.984874614022, 0.98624301618, 0.987483116952, 0.988607560325,
    0.989627655353, 0.990553530695, 0.991394269076, 0.992158024647,
    0.99285212571, 0.993483164877, 0.994057078393, 0.994579216081,
    0.995054403148, 0.995486994909, 0.995880925327, 0.996239750132,
    0.996566685198, 0.99686464074, 0.99713625183, 0.997383905666,
    0.997609765964, 0.997815794807, 0.998003772226, 0.998175313784,
    0.998331886363, 0.998474822355, 0.998605332446, 0.998724517108,
    0.998833376973, 0.998932822179, 0.999023680805, 0.999106706494,
    0.999182585332, 0.999251942072, 0.999315345767, 0.999373314859,
    0.999426321798, 0.999474797213, 0.99951913371, 0.999559689296,
    0.9995967905, 0.999630735199, 0.99966179518, 0.999690218475,
    0.999716231474, 0.999740040852, 0.999761835313, 0.999781787183,
    0.999800053855, 0.999816779104, 0.999832094286, 0.999846119426,
    0.999858964211, 0.999870728891, 0.999881505109, 0.999891376644,
    0.999900420098, 0.999908705519, 0.999916296969, 0.99992325304,
    0.999929627331, 0.999935468875, 0.999940822533, 0.999945729354,
    0.999950226904, 0.999954349561, 0.999958128793, 0.999961593401,
    0.999964769754, 0.999967681991, 0.999970352216, 0.999972800666,
    0.999975045876, 0.999977104817, 0.999978993036, 0.999980724771,
    0.999982313065, 0.999983769866, 0.999985106122, 0.99998633186,
    0.99998745627, 0.999988487774, 0.999989434086, 0.999990302279,
    0.999991098834, 0.99999182969, 0.999992500293, 0.999993115634,
    0.999993680288, 0.999994198449, 0.999994673963, 0.999995110355,
    0.999995510858, 0.999995878437, 0.999996215809, 0.999996525467,
    0.999996809697, 0.999997070596, 0.999997310086, 0.999997529931,
    0.999997731749, 0.999997917023, 0.999998087116, 0.999998243274,
    0.999998386645, 0.999998518279, 0.99999863914, 0.999998750113,
    0.999998852009, 0.999998945574, 0.999999031492, 0.999999110388,
    0.99999918284, 0.999999249374, 0.999999310476, 0.999999366591,
    0.999999418127, 0.999999465459, 0.99999950893, 0.999999548857,
    0.99999958553, 0.999999619214, 0.999999650153, 0.999999678573,
    0.999999704678, 0.999999728658, 0.999999750687, 0.999999770922,
    0.999999789512, 0.999999806589, 0.999999822278, 0.999999836692,
    0.999999849933, 0.999999862099, 0.999999873277, 0.999999883546,
    0.999999892982, 0.999999901651, 0.999999909617, 0.999999916936,
    0.999999923661, 0.999999929841, 0.999999935519, 0.999999940737,
    0.999999945532, 0.999999949938, 0.999999953987, 0.999999957708,
    0.999999961128, 0.99999996427, 0.999999967158, 0.999999969812,
    0.999999972252, 0.999999974493, 0.999999976554, 0.999999978447,
    0.999999980188, 0.999999981788, 0.999999983258, 0.999999984609,
    0.999999985851, 0.999999986993, 0.999999988043, 0.999999989007,
    0.999999989894, 0.999999990709, 0.999999991458, 0.999999992147,
    0.99999999278, 0.999999993362, 0.999999993897, 0.999999994389,
    0.999999994841, 0.999999995257, 0.999999995639, 0.99999999599,
    0.999999996313, 0.99999999661, 0.999999996883, 0.999999997134,
    0.999999997364, 0.999999997577, 0.999999997772, 0.999999997951,
    0.999999998116, 0.999999998267, 0.999999998407, 0.999999998535,
    0.999999998653, 0.999999998761, 0.999999998861, 0.999999998952,
    0.999999999036, 0.999999999114, 0.999999999185, 0.999999999251,
    0.999999999311, 0.999999999366, 0.999999999417, 0.999999999464,
    0.999999999507, 0.999999999547, 0.999999999583, 0.999999999616,
    0.999999999647, 0.999999999676, 0.999999999702, 0.999999999726,
    0.999999999748, 0.999999999768, 0.999999999786, 0.999999999804,
    0.999999999819, 0.999999999834, 0.999999999847, 0.999999999859,
    0.999999999871, 0.999999999881, 0.999999999891, 0.999999999899,
    0.999999999907, 0.999999999915, 0.999999999922, 0.999999999928,
    0.999999999934, 0.999999999939, 0.999999999944, 0.999999999948
    }};


#endif //ScalMC_H_
