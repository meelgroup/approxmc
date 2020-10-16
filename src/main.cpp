/*
 ApproxMC

 Copyright (c) 2019-2020, Mate Soos and Kuldeep S. Meel. All rights reserved
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
#include <signal.h>
#include <gmp.h>

#include "approxmc.h"
#include "time_mem.h"
#include <cryptominisat5/dimacsparser.h>
#include <cryptominisat5/streambuffer.h>
#include <arjun/arjun.h>

using namespace CMSat;
using std::cout;
using std::cerr;
using std::endl;
ApproxMC::AppMC* appmc = NULL;
ArjunNS::Arjun* arjun = NULL;

po::options_description main_options = po::options_description("Main options");
po::options_description improvement_options = po::options_description("Improvement options");
po::options_description misc_options = po::options_description("Misc options");
po::options_description help_options;
po::variables_map vm;
po::positional_options_description p;
uint32_t verbosity = 1;
uint32_t seed;
double epsilon;
double delta;
string logfilename;
uint32_t start_iter = 0;
uint32_t verb_cls = 0;
uint32_t simplify;
double var_elim_ratio;
uint32_t detach_xors = 1;
uint32_t reuse_models = 1;
uint32_t force_sol_extension = 0;
uint32_t sparse = 0;
vector<uint32_t> sampling_vars;
int ignore_sampl_set = 0;
int do_arjun = 1;
int debug_arjun = 0;
int arjun_incidence_sort;

void add_appmc_options()
{
    ApproxMC::AppMC tmp;
    epsilon = tmp.get_epsilon();
    delta = tmp.get_delta();
    simplify = tmp.get_simplify();
    var_elim_ratio = tmp.get_var_elim_ratio();
    sparse = tmp.get_sparse();
    seed = tmp.get_seed();

    ArjunNS::Arjun tmpa;
    arjun_incidence_sort = tmpa.get_incidence_sort();

    std::ostringstream my_epsilon;
    std::ostringstream my_delta;
    std::ostringstream my_var_elim_ratio;

    my_epsilon << std::setprecision(8) << epsilon;
    my_delta << std::setprecision(8) << delta;
    my_var_elim_ratio << std::setprecision(8) << var_elim_ratio;
 
    main_options.add_options()
    ("help,h", "Prints help")
    ("input", po::value< vector<string> >(), "file(s) to read")
    ("verb,v", po::value(&verbosity)->default_value(1), "verbosity")
    ("seed,s", po::value(&seed)->default_value(seed), "Seed")
    ("version", "Print version info")

    ("epsilon", po::value(&epsilon)->default_value(epsilon, my_epsilon.str())
        , "epsilon parameter as per PAC guarantees")
    ("delta", po::value(&delta)->default_value(delta, my_delta.str())
        , "delta parameter as per PAC guarantees; 1-delta is the confidence")
    ("log", po::value(&logfilename),
         "Logs of ApproxMC execution")
    ("ignore", po::value(&ignore_sampl_set)->default_value(ignore_sampl_set)
        , "Ignore given sampling set and recompute it with Arjun")
    ("arjun", po::value(&do_arjun)->default_value(do_arjun)
        , "Use arjun to minimize sampling set")
    ("arjuninc", po::value(&arjun_incidence_sort)->default_value(arjun_incidence_sort)
        , "Select incidence sorting. Probe-based is 3. Simple incidence-based is 1. Component-to-other-component based is 5. Random is 5")
    ("debugarjun", po::value(&debug_arjun)->default_value(debug_arjun)
        , "Use CNF from Arjun, but use sampling set from CNF")
    ;

    improvement_options.add_options()
    ("sparse", po::value(&sparse)->default_value(sparse)
        , "0 = (default) Do not use sparse method. 1 = Generate sparse XORs when possible.")
    ("detachxor", po::value(&detach_xors)->default_value(detach_xors)
        , "Detach XORs in CMS")
    ("reusemodels", po::value(&reuse_models)->default_value(reuse_models)
        , "Reuse models while counting solutions")
    ("forcesolextension", po::value(&force_sol_extension)->default_value(force_sol_extension)
        , "Use trick of not extending solutions in the SAT solver to full solution")
    ;

    misc_options.add_options()
    ("verbcls", po::value(&verb_cls)->default_value(verb_cls)
        ,"Print banning clause + xor clauses. Highly verbose.")
    ("simplify", po::value(&simplify)->default_value(simplify)
        , "Simplify agressiveness")
    ("velimratio", po::value(&var_elim_ratio)->default_value(var_elim_ratio)
        , "Variable elimination ratio for each simplify run")
    ;

    help_options.add(main_options);
    help_options.add(improvement_options);
    help_options.add(misc_options);
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
            cout << appmc->get_version_info();
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

// void SIGINT_handler(int)
// {
//     if (!appmc) {
//         return;
//     }
//     SolCount solCount = appmc->calc_est_count();
//     if (!solCount.valid) {
//         cout << "c did not manage to get a single measurement, we have no estimate of the count" << endl;
//         exit(-1);
//     }
//     cout << "c Below count is NOT FULLY APPROXIMIATE due to early-abort!" << endl;
//     solCount.print_num_solutions();
//     exit(-1);
// }

template<class T>
void read_in_file(const string& filename, T* myreader)
{
    #ifndef USE_ZLIB
    FILE * in = fopen(filename.c_str(), "rb");
    DimacsParser<StreamBuffer<FILE*, FN>, T> parser(myreader, NULL, verbosity);
    #else
    gzFile in = gzopen(filename.c_str(), "rb");
    DimacsParser<StreamBuffer<gzFile, GZ>, T> parser(myreader, NULL, verbosity);
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

    sampling_vars = parser.sampling_vars;

    #ifndef USE_ZLIB
    fclose(in);
    #else
    gzclose(in);
    #endif
}

void print_indep_set(const vector<uint32_t>& indep_set, uint32_t orig_sampling_set_size)
{
    cout << "vp ";
    for(const uint32_t s: indep_set) {
        cout << s+1 << " ";
    }
    cout << "0" << endl;

    cout << "c set size: " << std::setw(8)
    << indep_set.size()
    << " fraction of original: "
    <<  std::setw(6) << std::setprecision(4)
    << (double)indep_set.size()/(double)orig_sampling_set_size
    << endl << std::flush;
}

template<class T>
void read_stdin(T* myreader)
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
    DimacsParser<StreamBuffer<FILE*, FN>, T> parser(myreader, NULL, verbosity);
    #else
    DimacsParser<StreamBuffer<gzFile, GZ>, T> parser(myreader, NULL, verbosity);
    #endif

    if (!parser.parse_DIMACS(in, false)) {
        exit(-1);
    }

    sampling_vars = parser.sampling_vars;

    #ifdef USE_ZLIB
    gzclose(in);
    #endif
}

void print_num_solutions(uint32_t cellSolCount, uint32_t hashCount)
{
    cout << "c [appmc] Number of solutions is: "
    << cellSolCount << "*2**" << hashCount << endl;

    mpz_t num_sols;
    mpz_init (num_sols);
    mpz_ui_pow_ui(num_sols, 2, hashCount);
    mpz_mul_ui(num_sols, num_sols, cellSolCount);

    cout << "s mc " << std::flush;
    mpz_out_str(0, 10, num_sols);
    cout << endl;
    mpz_clear(num_sols);

}

void get_cnf_from_arjun()
{
    bool ret = true;
    const uint32_t orig_num_vars = arjun->get_orig_num_vars();
    appmc->new_vars(orig_num_vars);
    arjun->start_getting_small_clauses(
        std::numeric_limits<uint32_t>::max(),
        std::numeric_limits<uint32_t>::max(),
        false);
    vector<Lit> clause;
    while (ret) {
        ret = arjun->get_next_small_clause(clause);
        if (!ret) {
            break;
        }

        bool ok = true;
        for(auto l: clause) {
            if (l.var() >= orig_num_vars) {
                ok = false;
                break;
            }
        }

        if (ok) {
            appmc->add_clause(clause);
        }

    }
    arjun->end_getting_small_clauses();
}

template<class T>
void read_input_cnf(T* reader)
{
    //Init Arjun, read in file, get minimal indep set
    if (vm.count("input") != 0) {
        vector<string> inp = vm["input"].as<vector<string> >();
        if (inp.size() > 1) {
            cout << "[appmc] ERROR: you must only give one CNF as input" << endl;
            exit(-1);
        }

        read_in_file(inp[0].c_str(), reader);
    } else {
        read_stdin(reader);
    }
}

void minimize_sampling_set()
{
    uint32_t orig_sampling_set_size;
    if (sampling_vars.empty() || ignore_sampl_set) {
        orig_sampling_set_size = arjun->start_with_clean_sampling_set();
    } else {
        orig_sampling_set_size = arjun->set_starting_sampling_set(sampling_vars);
    }

    sampling_vars = arjun->get_indep_set();
    print_indep_set(sampling_vars , orig_sampling_set_size);
}

void set_approxmc_options()
{
    //Main options
    appmc->set_verbosity(verbosity);
    if (verbosity > 2) {
        appmc->set_detach_warning();
    }
    appmc->set_seed(seed);
    appmc->set_epsilon(epsilon);
    appmc->set_delta(delta);

    //Improvement options
    appmc->set_detach_xors(detach_xors);
    appmc->set_reuse_models(reuse_models);
    appmc->set_force_sol_extension(force_sol_extension);
    appmc->set_sparse(sparse);

    //Misc options
    appmc->set_start_iter(start_iter);
    appmc->set_verb_cls(verb_cls);
    appmc->set_simplify(simplify);
    appmc->set_var_elim_ratio(var_elim_ratio);

    if (logfilename != "") {
        appmc->set_up_log(logfilename);
        cout << "c [appmc] Logfile set " << logfilename << endl;
    }
}

void print_orig_sampling_vars(const vector<uint32_t>& orig_sampling_vars)
{
    if (!orig_sampling_vars.empty()) {
        cout << "Original sampling vars: ";
        for(auto v: orig_sampling_vars) {
            cout << v << " ";
        }
        cout << endl;
        cout << "Orig sampling vars size: " << orig_sampling_vars.size() << endl;
    } else {
        cout << "No original sampling vars given" << endl;
        cout << "Orig sampling vars size: " << arjun->nVars() << endl;
    }
}

int main(int argc, char** argv)
{
    #if defined(__GNUC__) && defined(__linux__)
    feenableexcept(FE_INVALID   |
                   FE_DIVBYZERO |
                   FE_OVERFLOW
                  );
    #endif
//     signal(SIGALRM, SIGINT_handler);
//     signal(SIGTERM, SIGINT_handler);
    double start_time = cpuTime();

    //Reconstruct the command line so we can emit it later if needed
    string command_line;
    for(int i = 0; i < argc; i++) {
        command_line += string(argv[i]);
        if (i+1 < argc) {
            command_line += " ";
        }
    }

    appmc = new ApproxMC::AppMC;
    add_supported_options(argc, argv);
    if (verbosity) {
        cout << appmc->get_version_info();
        cout << "c executed with command line: " << command_line << endl;
    }

    set_approxmc_options();

    if (do_arjun) {
        //Arjun-based minimization
        arjun = new ArjunNS::Arjun;
        arjun->set_seed(seed);
        arjun->set_verbosity(verbosity);
        arjun->set_incidence_sort(arjun_incidence_sort);
        read_input_cnf(arjun);
        print_orig_sampling_vars(sampling_vars);
        auto old_sampling_vars = sampling_vars;

        minimize_sampling_set();
        get_cnf_from_arjun();
        if (debug_arjun) {
            sampling_vars = old_sampling_vars;
        }
        delete arjun;
    } else {
        read_input_cnf(appmc);
        if (ignore_sampl_set) {
            sampling_vars.clear();
            for(uint32_t i = 0; i < appmc->nVars(); i++) {
                sampling_vars.push_back(i);
            }
        }
    }
    appmc->set_projection_set(sampling_vars);

    //Count with ApproxMC
    auto sol_count = appmc->count();
    appmc->print_stats(start_time);

    print_num_solutions(sol_count.cellSolCount, sol_count.hashCount);
    delete appmc;
}
