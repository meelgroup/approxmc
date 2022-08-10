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

#include "time_mem.h"
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
po::options_description arjun_options = po::options_description("Arjun options");
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
int dump_intermediary_cnf = 0;

//Arjun
vector<uint32_t> sampling_vars;
bool sampling_vars_found = false;
int ignore_sampl_set = 0;
int do_arjun = 1;
int debug_arjun = 0;
int arjun_incidence_sort;
int cont_recomp_indep_set = 0;
int with_e = 0;
int do_empty_occ = 1;
int arjun_irreg = 1;
int arjun_emtpy = 1;
int arjun_mirror_emtpy = 1;

void add_appmc_options()
{
    ApproxMC::AppMC tmp;
    epsilon = tmp.get_epsilon();
    delta = tmp.get_delta();
    simplify = tmp.get_simplify();
    var_elim_ratio = tmp.get_var_elim_ratio();
    sparse = tmp.get_sparse();
    seed = tmp.get_seed();
    cont_recomp_indep_set = tmp.get_cont_recomp_indep_set();

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
        , "Tolerance parameter, i.e. what level of error we expect of the result? CAV-2020 paper default is 0.80")
    ("delta", po::value(&delta)->default_value(delta, my_delta.str())
        , "Confidence parameter, i.e. how sure are we of the result? CAV-2020 paper default is 0.20")
    ("log", po::value(&logfilename),
         "Logs of ApproxMC execution")
    ("ignore", po::value(&ignore_sampl_set)->default_value(ignore_sampl_set)
        , "Ignore given sampling set and recompute it with Arjun")
    ;

    ArjunNS::Arjun tmpa;
    arjun_incidence_sort = tmpa.get_incidence_sort();

    arjun_options.add_options()
    ("arjun", po::value(&do_arjun)->default_value(do_arjun)
        , "Use arjun to minimize sampling set")
    ("arjuninc", po::value(&arjun_incidence_sort)->default_value(arjun_incidence_sort)
        , "Select incidence sorting. Probe-based is 3. Simple incidence-based is 1. Component-to-other-component based is 5. Random is 5")
    ("arjundebug", po::value(&debug_arjun)->default_value(debug_arjun)
        , "Use CNF from Arjun, but use sampling set from CNF")
    ("arjuncontrecomp", po::value(&cont_recomp_indep_set)->default_value(cont_recomp_indep_set)
        , "Continiously, at every XOR addition, recompute the independent set through Arjun")
    ("arjunirreg", po::value(&arjun_irreg)->default_value(arjun_irreg)
        , "Arjun should use irregular gates")
    ("arjunempty", po::value(&arjun_emtpy)->default_value(arjun_emtpy)
        , "Arjun should use check for empty occurrence variables")
    ("arjunmirror", po::value(&arjun_mirror_emtpy)->default_value(arjun_mirror_emtpy)
        , "Arjun should use mirror variables to improve empty detection")
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
    ("withe", po::value(&with_e)->default_value(with_e)
        , "Eliminate variables and simplify CNF as well")
    ("emptyocc", po::value(&do_empty_occ)->default_value(do_empty_occ)
        , "Don't count over empty occ variables, instead just multiply")
    ;


    misc_options.add_options()
    ("verbcls", po::value(&verb_cls)->default_value(verb_cls)
        ,"Print banning clause + xor clauses. Highly verbose.")
    ("simplify", po::value(&simplify)->default_value(simplify)
        , "Simplify agressiveness")
    ("velimratio", po::value(&var_elim_ratio)->default_value(var_elim_ratio)
        , "Variable elimination ratio for each simplify run")
    ("dumpintercnf", po::value(&dump_intermediary_cnf)->default_value(dump_intermediary_cnf)
        , "Dump intermediary CNFs during solving into files cnf_dump-X.cnf. If set to 1 only UNSAT is dumped, if set to 2, all are dumped");
    ;

    help_options.add(main_options);
    help_options.add(improvement_options);
    help_options.add(arjun_options);
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
    sampling_vars_found = parser.sampling_vars_found;

    #ifndef USE_ZLIB
    fclose(in);
    #else
    gzclose(in);
    #endif
}

inline double stats_line_percent(double num, double total)
{
    if (total == 0) {
        return 0;
    } else {
        return num/total*100.0;
    }
}

void print_final_indep_set(
    const vector<uint32_t>& indep_set, uint32_t orig_sampling_set_size, const vector<uint32_t>& empty_occs)
{
    cout << "c ind ";
    for(const uint32_t s: indep_set) cout << s+1 << " ";
    cout << "0" << endl;

    cout
    << "c [arjun] final set size:      " << std::setw(7) << indep_set.size()
    << " percent of original: "
    <<  std::setw(6) << std::setprecision(4)
    << stats_line_percent(indep_set.size(), orig_sampling_set_size)
    << " %" << endl
    << "c [arjun] of which empty occs: " << std::setw(7) << empty_occs.size()
    << " percent of original: "
    <<  std::setw(6) << std::setprecision(4)
    << stats_line_percent(empty_occs.size(), orig_sampling_set_size)
    << " %" << endl;
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
    sampling_vars_found = parser.sampling_vars_found;

    #ifdef USE_ZLIB
    gzclose(in);
    #endif
}

void print_num_solutions(uint32_t cellSolCount, uint32_t hashCount)
{
    cout << "c [appmc] Number of solutions is: "
    << cellSolCount << "*2**" << hashCount << endl;
    if (cellSolCount == 0) {
        cout << "s UNSATISFIABLE" << endl;
    } else {
        cout << "s SATISFIABLE" << endl;
    }

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

    vector<Lit> lits;
    for(const auto& bnn: arjun->get_bnns()) {
        if (bnn) {
            lits.clear();
            lits.insert(lits.end(), bnn->begin(), bnn->end());
            appmc->add_bnn_clause(lits, bnn->cutoff, bnn->out);
        }
    }
}

// void get_cnf_and_sampl_from_arjun_fully_simplified_renumbered()
// {
//     auto cnf_sampl = arjun->get_fully_simplified_renumbered_cnf(sampling_vars, arjun->get_orig_num_vars());
//
//     //Get num vars
//     uint32_t max_var = 0;
//     for(const auto& cl: cnf_sampl.first) {
//         for(const auto& l: cl) {
//             if (l.var()+1 > max_var) {
//                 max_var = l.var()+1;
//             }
//         }
//     }
//     appmc->new_vars(max_var);
//
//     //Add clauses
//     for(const auto& cl: cnf_sampl.first) {
//         appmc->add_clause(cl);
//     }
//     sampling_vars = cnf_sampl.second;
// }

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

uint32_t set_up_sampling_set()
{
    uint32_t orig_sampling_set_size;
    if (!sampling_vars_found || ignore_sampl_set) {
        orig_sampling_set_size = arjun->start_with_clean_sampling_set();
    } else {
        orig_sampling_set_size = arjun->set_starting_sampling_set(sampling_vars);
    }

    return orig_sampling_set_size;
}

void set_approxmc_options()
{
    //Main options
    appmc->set_verbosity(verbosity);
    if (verbosity >= 2) {
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
    appmc->set_dump_intermediary_cnf(dump_intermediary_cnf);

    //Arjun options
    appmc->set_cont_recomp_indep_set(cont_recomp_indep_set);

    if (logfilename != "") {
        appmc->set_up_log(logfilename);
        cout << "c [appmc] Logfile set " << logfilename << endl;
    }
}

template<class T>
void print_orig_sampling_vars(const vector<uint32_t>& orig_sampling_vars, T* ptr)
{
    if (!orig_sampling_vars.empty()) {
        cout << "c Original sampling vars: ";
        for(auto v: orig_sampling_vars) {
            cout << v << " ";
        }
        cout << endl;
        cout << "c [appmc] Orig sampling vars size: " << orig_sampling_vars.size() << endl;
    } else {
        cout << "c [appmc] No original sampling vars given" << endl;
        cout << "c [appmc] Orig sampling vars size: " << ptr->nVars() << endl;
    }
}

void transfer_unit_clauses_from_arjun()
{
    vector<Lit> cl(1);
    auto units = arjun->get_zero_assigned_lits();
    for(const auto& unit: units) {
        if (unit.var() < appmc->nVars()) {
            cl[0] = unit;
            appmc->add_clause(cl);
        }
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

    uint32_t offset_count_by_2_pow = 0;
    if (do_arjun) {
        //Arjun-based minimization
        arjun = new ArjunNS::Arjun;
        arjun->set_seed(seed);
        arjun->set_verbosity(verbosity);
        arjun->set_incidence_sort(arjun_incidence_sort);
        arjun->set_simp(simplify);
        arjun->set_irreg_gate_based(arjun_irreg);
        arjun->set_empty_occs_based(arjun_emtpy);
        arjun->set_mirror_empty(arjun_mirror_emtpy);

        if (verbosity) cout << "c Arjun SHA revision " <<  arjun->get_version_info() << endl;

        read_input_cnf(arjun);
        print_orig_sampling_vars(sampling_vars, arjun);
        auto old_sampling_vars = sampling_vars;
        uint32_t orig_sampling_set_size = set_up_sampling_set();
        sampling_vars = arjun->get_indep_set();
        vector<uint32_t> empty_occ_sampl_vars;
        if (do_empty_occ) empty_occ_sampl_vars = arjun->get_empty_occ_sampl_vars();
        print_final_indep_set(sampling_vars , orig_sampling_set_size, empty_occ_sampl_vars);
        if (!with_e) {
            if (do_empty_occ) {
                std::set<uint32_t> sampl_vars_set;
                sampl_vars_set.insert(sampling_vars.begin(), sampling_vars.end());
                for(auto const& v: empty_occ_sampl_vars) {
                    assert(sampl_vars_set.find(v) != sampl_vars_set.end()); // this is guaranteed by arjun
                    sampl_vars_set.erase(v);
                }
                offset_count_by_2_pow = empty_occ_sampl_vars.size();
                sampling_vars.clear();
                sampling_vars.insert(sampling_vars.end(), sampl_vars_set.begin(), sampl_vars_set.end());
            }
            get_cnf_from_arjun();
            transfer_unit_clauses_from_arjun();
        } else {
            auto ret = arjun->get_fully_simplified_renumbered_cnf(sampling_vars, empty_occ_sampl_vars, arjun->nVars(), false);
            //const std::pair<vector<vector<Lit>>, uint32_t>& cnf, const vector<uint32_t>& sampl_set, const uint32_t multiply = 0
            const auto& cls = std::get<0>(ret).first;
            const auto& max_var = std::get<0>(ret).second;
            appmc->new_vars(max_var);
            for(const auto& cl: cls) appmc->add_clause(cl);
            sampling_vars = std::get<1>(ret);
            offset_count_by_2_pow = std::get<2>(ret);
        }
        if (debug_arjun) {
            sampling_vars = old_sampling_vars;
            offset_count_by_2_pow = 0;
        }
        delete arjun;
    } else {
        read_input_cnf(appmc);
        if (ignore_sampl_set || !sampling_vars_found) {
            sampling_vars.clear();
            for(uint32_t i = 0; i < appmc->nVars(); i++) sampling_vars.push_back(i);
        }
        print_final_indep_set(sampling_vars , 0, vector<uint32_t>());
    }

    ApproxMC::SolCount sol_count;
    if (!sampling_vars.empty()) {
        appmc->set_projection_set(sampling_vars);
        sol_count = appmc->count();
        appmc->print_stats(start_time);
    } else {
        bool ret = appmc->find_one_solution();
        sol_count.hashCount = 0;
        if (ret) sol_count.cellSolCount = 1;
        else sol_count.cellSolCount = 0;
    }
    cout << "c [appmc+arjun] Total time: " << (cpuTime() - start_time) << endl;
    print_num_solutions(sol_count.cellSolCount, sol_count.hashCount+offset_count_by_2_pow);

    delete appmc;
}
