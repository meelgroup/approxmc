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

#include <boost/program_options.hpp>
using boost::lexical_cast;
namespace po = boost::program_options;
using std::string;
using std::vector;

#if defined(__GNUC__) && defined(__linux__)
#include <fenv.h>
#endif

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

AppMCConfig conf;
po::options_description appmc_options = po::options_description("ApproxMC options");
po::options_description help_options;
po::variables_map vm;
po::positional_options_description p;

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

void add_appmc_options()
{
    std::ostringstream my_epsilon;
    std::ostringstream my_delta;
    my_epsilon << std::setprecision(8) << conf.epsilon;
    my_delta << std::setprecision(8) << conf.delta;


    appmc_options.add_options()
    ("help,h", "Prints help")
    ("version", "Print version info")
    ("input", po::value< vector<string> >(), "file(s) to read")
    ("verb,v", po::value(&conf.verb)->default_value(conf.verb), "verbosity")
    ("seed,s", po::value(&conf.seed)->default_value(conf.seed), "Seed")
    ("epsilon", po::value(&conf.epsilon)->default_value(conf.epsilon, my_epsilon.str())
        , "epsilon parameter as per PAC guarantees")
    ("delta", po::value(&conf.delta)->default_value(conf.delta, my_delta.str())
        , "delta parameter as per PAC guarantees; 1-delta is the confidence")
    ("start", po::value(&conf.start_iter)->default_value(conf.start_iter),
         "Start at this many XORs")
    ("log", po::value(&conf.logfilename)->default_value(conf.logfilename),
         "Logs of ApproxMC execution")
    ("th", po::value(&conf.num_threads)->default_value(conf.num_threads),
         "How many solving threads to use per solver call")
    ("vcl", po::value(&conf.verb_appmc_cls)->default_value(conf.verb_appmc_cls)
        ,"Print banning clause + xor clauses. Highly verbose.")
    ("sparse", po::value(&conf.sparse)->default_value(conf.sparse)
        , "Generate sparse XORs when possible")
    ("startiter", po::value(&conf.startiter)->default_value(conf.startiter)
        , "If positive, use instead of startiter computed by AppMC")
    ;

    help_options.add(appmc_options);
}

void add_supported_options(int argc, char** argv)
{
    add_appmc_options();
    p.add("input", 1);

    try {
        po::store(po::command_line_parser(argc, argv).options(help_options).positional(p).run(), vm);
        if (vm.count("help"))
        {
            cout
            << "Probably Approximate counter" << endl;

            cout
            << "approxmc [options] inputfile" << endl << endl;

            cout << help_options << endl;
            std::exit(0);
        }

        if (vm.count("version")) {
            appmc->printVersionInfo();
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
        boost::exception_detail::error_info_injector<po::invalid_option_value> >& what
    ) {
        cerr
        << "ERROR: Invalid value '" << what.what() << "'" << endl
        << "       given to option '" << what.get_option_name() << "'"
        << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::multiple_occurrences> >& what
    ) {
        cerr
        << "ERROR: " << what.what() << " of option '"
        << what.get_option_name() << "'"
        << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::required_option> >& what
    ) {
        cerr
        << "ERROR: You forgot to give a required option '"
        << what.get_option_name() << "'"
        << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::too_many_positional_options_error> >& what
    ) {
        cerr
        << "ERROR: You gave too many positional arguments. Only the input CNF can be given as a positional option." << endl;
        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::ambiguous_option> >& what
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
        boost::exception_detail::error_info_injector<po::invalid_command_line_syntax> >& what
    ) {
        cerr
        << "ERROR: The option you gave is missing the argument or the" << endl
        << "       argument is given with space between the equal sign." << endl
        << "       detailed error message: " << what.what() << endl
        ;
        std::exit(-1);
    }
}

void readInAFile(SATSolver* solver2, const string& filename)
{
    solver2->add_sql_tag("filename", filename);
    #ifndef USE_ZLIB
    FILE * in = fopen(filename.c_str(), "rb");
    DimacsParser<StreamBuffer<FILE*, FN> > parser(solver, NULL, 2);
    #else
    gzFile in = gzopen(filename.c_str(), "rb");
    DimacsParser<StreamBuffer<gzFile, GZ> > parser(appmc->solver, NULL, 2);
    #endif

    if (in == NULL) {
        std::cerr
        << "ERROR! Could not open file '"
        << filename
        << "' for reading: " << strerror(errno) << endl;

        std::exit(-1);
    }

    if (!parser.parse_DIMACS(in, false)) {
        exit(-1);
    }

    conf.sampling_set.swap(parser.sampling_vars);

    #ifndef USE_ZLIB
        fclose(in);
    #else
        gzclose(in);
    #endif
}

void readInStandardInput(SATSolver* solver2)
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

    conf.sampling_set.swap(parser.sampling_vars);

    #ifdef USE_ZLIB
        gzclose(in);
    #endif
}

void set_sampling_vars()
{
    if (conf.sampling_set.empty()) {
        cout
        << "[appmc] WARNING! Sampling set was not declared with 'c ind var1 [var2 var3 ..] 0'"
        " notation in the CNF." << endl
        << "[appmc] we may work substantially worse!" << endl;
        for (size_t i = 0; i < appmc->solver->nVars(); i++) {
            conf.sampling_set.push_back(i);
        }
    }
    cout << "[appmc] Sampling set size: " << conf.sampling_set.size() << endl;
    if (conf.sampling_set.size() > 100) {
        cout
        << "[appmc] Sampling var set contains over 100 variables, not displaying"
        << endl;
    } else {
        cout << "[appmc] Sampling set: ";
        for (auto v: conf.sampling_set) {
            cout << v+1 << ", ";
        }
        cout << endl;
    }
    appmc->solver->set_sampling_vars(&conf.sampling_set);
}

int main(int argc, char** argv)
{
    #if defined(__GNUC__) && defined(__linux__)
    feenableexcept(FE_INVALID   |
                   FE_DIVBYZERO |
                   FE_OVERFLOW
                  );
    #endif

    appmc = new AppMC;
    add_supported_options(argc, argv);
    appmc->printVersionInfo();

    cout << "[appmc] using seed: " << conf.seed << endl;

    if (vm.count("log") == 0) {
        if (vm.count("input") != 0) {
            conf.logfilename = vm["input"].as<vector<string> >()[0] + ".log";
            cout << "[appmc] Logfile name not given, assumed to be " << conf.logfilename << endl;
        } else {
            std::cerr << "[appmc] ERROR: You must provide the logfile name" << endl;
            exit(-1);
        }
    }

    //startTime = cpuTimeTotal();
    appmc->solver = new SATSolver();
    appmc->solver->set_up_for_scalmc();

    if (conf.verb > 2) {
        appmc->solver->set_verbosity(conf.verb-2);
    }
    appmc->solver->set_allow_otf_gauss();

    if (conf.num_threads > 1) {
        appmc->solver->set_num_threads(conf.num_threads);
    }

    if (conf.epsilon < 0.0) {
        cout << "[appmc] ERROR: invalid epsilon" << endl;
        exit(-1);
    }
    conf.threshold = int(1 + 9.84*(1+(1/conf.epsilon))*(1+(1/conf.epsilon))*(1+(conf.epsilon/(1+conf.epsilon))));

    if (conf.delta <= 0.0 || conf.delta > 1.0) {
        cout << "[appmc] ERROR: invalid delta" << endl;
        exit(-1);
    }
    conf.measurements = (int)std::ceil(std::log2(3.0/conf.delta)*17);
    for (int count = 0; count < 256; count++) {
        if(iterationConfidences[count] >= 1 - conf.delta){
            conf.measurements = count*2+1;
            break;
        }
    }

    //parsing the input
    if (vm.count("input") != 0) {
        vector<string> inp = vm["input"].as<vector<string> >();
        if (inp.size() > 1) {
            cout << "[appmc] ERROR: can only parse in one file" << endl;
        }
        readInAFile(appmc->solver, inp[0].c_str());
    } else {
        readInStandardInput(appmc->solver);
    }
    set_sampling_vars();

    if (conf.start_iter > conf.sampling_set.size()) {
        cout << "[appmc] ERROR: Manually-specified start_iter"
             "is larger than the size of the sampling set.\n" << endl;
        exit(-1);
    }

    return appmc->solve(conf);
}
