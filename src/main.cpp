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

#include <string>
#include <vector>
#if defined(__GNUC__) && defined(__linux__)
#include <cfenv>
#endif
#include <cstdint>
#include <set>
#include <gmp.h>

#include "time_mem.h"
#include "approxmc.h"
#include "time_mem.h"
#include <cryptominisat5/dimacsparser.h>
#include <cryptominisat5/streambuffer.h>
#include <arjun/arjun.h>
#include "src/argparse.hpp"

using namespace CMSat;
using std::cout;
using std::cerr;
using std::endl;
using std::set;
using std::string;
using std::vector;
ApproxMC::AppMC* appmc = nullptr;
argparse::ArgumentParser program = argparse::ArgumentParser("approxmc");

uint32_t verb = 1;
uint32_t seed;
double epsilon;
double delta;
string logfilename;
uint32_t start_iter = 0;
uint32_t verb_cls = 0;
uint32_t simplify;
double var_elim_ratio;
uint32_t reuse_models = 1;
uint32_t force_sol_extension = 0;
uint32_t sparse = 0;
int dump_intermediary_cnf = 0;
int arjun_gates = 0;
bool debug = false;

//Arjun
ArjunNS::SimpConf simp_conf;
int with_e = 0;
bool all_indep = false;
bool do_arjun = true;

#define myopt(name, var, fun, hhelp) \
    program.add_argument(name) \
        .action([&](const auto& a) {var = std::fun(a.c_str());}) \
        .default_value(var) \
        .help(hhelp)
#define myopt2(name1, name2, var, fun, hhelp) \
    program.add_argument(name1, name2) \
        .action([&](const auto& a) {var = std::fun(a.c_str());}) \
        .default_value(var) \
        .help(hhelp)

void add_appmc_options()
{
    ApproxMC::AppMC tmp;
    epsilon = tmp.get_epsilon();
    delta = tmp.get_delta();
    simplify = tmp.get_simplify();
    var_elim_ratio = tmp.get_var_elim_ratio();
    sparse = tmp.get_sparse();
    seed = tmp.get_seed();


    myopt2("-v", "--verb", verb, atoi, "Verbosity");
    myopt2("-s", "--seed", seed, atoi, "Seed");
    myopt2("-e", "--epsilon", epsilon, stod,
            "Tolerance parameter, i.e. how close is the count from the correct count? "
            "Count output is within bounds of (exact_count/(1+e)) < count < (exact_count*(1+e)). "
            "So e=0.8 means we'll output at most 180%% of exact count and at least 55%% of exact count. "
            "Lower value means more precise.");
    myopt2("-d", "--delta", delta, stod, "Confidence parameter, i.e. how sure are we of the result? "
            "(1-d) = probability the count is within range as per epsilon parameter. "
            "So d=0.2 means we are 80%% sure the count is within range as specified by epsilon. "
            "The lower, the higher confidence we have in the count.");

    /* arjun_options.add_options() */
    myopt("--arjun", do_arjun, atoi, "Use arjun to minimize sampling set");

    /* improvement_options.add_options() */
    myopt("--sparse", sparse, atoi,
            "0 = (default) Do not use sparse method. 1 = Generate sparse XORs when possible.");
    myopt("--reusemodels", reuse_models, atoi, "Reuse models while counting solutions");
    myopt("--forcesolextension", force_sol_extension, atoi,
            "Use trick of not extending solutions in the SAT solver to full solution");
    myopt("--withe", with_e, atoi, "Eliminate variables and simplify CNF as well");
    myopt("--eiter1", simp_conf.iter1, atoi, "Num iters of E on 1st round");
    myopt("--eiter2", simp_conf.iter2, atoi, "Num iters of E on 1st round");
    myopt("--evivif", simp_conf.oracle_vivify, atoi, "E vivif");
    myopt("--esparsif", simp_conf.oracle_sparsify, atoi, "E sparsify");
    myopt("--egetreds", simp_conf.oracle_vivify_get_learnts, atoi, "Get redundant from E");

    /* misc_options.add_options() */
    myopt("--verbcls", verb_cls, atoi, "Print banning clause + xor clauses. Highly verbose.");
    myopt("--simplify", simplify, atoi, "Simplify aggressiveness");
    myopt("--velimratio", var_elim_ratio, stod, "Variable elimination ratio for each simplify run");
    myopt("--dumpintercnf", dump_intermediary_cnf, atoi,
            "Dump intermediary CNFs during solving into files cnf_dump-X.cnf. If set to 1 only UNSAT is dumped, if set to 2, all are dumped");
    myopt("--log", logfilename, string, "Put logs of ApproxMC execution to this file");
    myopt("--debug", debug, atoi, "Turn on more heavy internal debugging");

    program.add_argument("inputfile").remaining().help("input CNF");
}

void add_supported_options(int argc, char** argv) {
    add_appmc_options();
    try {
        program.parse_args(argc, argv);
        if (program.is_used("--help")) {
            cout << "Probilistic Approcimate Counter" << endl << endl
            << "approxmc [options] inputfile" << endl;
            cout << program << endl;
            exit(0);
        }
    }
    catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        exit(-1);
    }

    if (program["version"] == true) {
        cout << appmc->get_version_info();
        exit(0);
    }
}

template<class T> void read_in_file(const string& filename, T* myreader)
{
    #ifndef USE_ZLIB
    FILE * in = fopen(filename.c_str(), "rb");
    DimacsParser<StreamBuffer<FILE*, FN>, T> parser(myreader, nullptr, verb);
    #else
    gzFile in = gzopen(filename.c_str(), "rb");
    DimacsParser<StreamBuffer<gzFile, GZ>, T> parser(myreader, nullptr, verb);
    #endif

    if (in == nullptr) {
        std::cerr
        << "ERROR! Could not open file '"
        << filename
        << "' for reading: " << strerror(errno) << endl;

        std::exit(-1);
    }

    if (!parser.parse_DIMACS(in, true)) exit(-1);

    #ifndef USE_ZLIB
    fclose(in);
    #else
    gzclose(in);
    #endif
}

inline double stats_line_percent(double num, double total)
{
    if (total == 0) { return 0;
    } else { return num/total*100.0;
    }
}

void print_final_indep_set(const vector<uint32_t>& indep_set, uint32_t orig_sampling_set_size)
{
    cout << "c ind ";
    for(const uint32_t s: indep_set) cout << s+1 << " ";
    cout << "0" << endl;

    cout
    << "c [arjun] final set size:      " << std::setw(7) << indep_set.size()
    << " percent of original: "
    <<  std::setw(6) << std::setprecision(4)
    << stats_line_percent(indep_set.size(), orig_sampling_set_size)
    << " %" << endl;
}

template<class T> void read_stdin(T* myreader) {
    cout << "c Reading from standard input... Use '-h' or '--help' for help." << endl;

    #ifndef USE_ZLIB
    FILE * in = stdin;
    #else
    gzFile in = gzdopen(0, "rb"); //opens stdin, which is 0
    #endif

    if (in == nullptr) {
        std::cerr << "ERROR! Could not open standard input for reading" << endl;
        std::exit(1);
    }

    #ifndef USE_ZLIB
    DimacsParser<StreamBuffer<FILE*, FN>, T> parser(myreader, nullptr, verb);
    #else
    DimacsParser<StreamBuffer<gzFile, GZ>, T> parser(myreader, nullptr, verb);
    #endif

    if (!parser.parse_DIMACS(in, false)) exit(-1);

    #ifdef USE_ZLIB
    gzclose(in);
    #endif
}

void print_num_solutions(uint32_t cell_sol_cnt, uint32_t hash_count, const mpq_class& mult)
{
    cout << "c [appmc] Number of solutions is: "
    << cell_sol_cnt << "*2**" << hash_count << "*" << mult << endl;
    if (cell_sol_cnt == 0) cout << "s UNSATISFIABLE" << endl;
    else cout << "s SATISFIABLE" << endl;

    mpz_class num_sols(2);
    mpz_pow_ui(num_sols.get_mpz_t(), num_sols.get_mpz_t(), hash_count);
    num_sols *= cell_sol_cnt;
    mpq_class final = mult * num_sols;

    cout << "s mc " << num_sols << endl;
}

void set_approxmc_options()
{
    //Main options
    appmc->set_verbosity(verb);
    appmc->set_seed(seed);
    appmc->set_epsilon(epsilon);
    appmc->set_delta(delta);

    //Improvement options
    appmc->set_reuse_models(reuse_models);
    appmc->set_sparse(sparse);

    //Misc options
    appmc->set_start_iter(start_iter);
    appmc->set_verb_cls(verb_cls);
    appmc->set_simplify(simplify);
    appmc->set_var_elim_ratio(var_elim_ratio);
    appmc->set_dump_intermediary_cnf(dump_intermediary_cnf);
    appmc->set_force_sol_extension(force_sol_extension);
    if (debug) {
        appmc->set_force_sol_extension(1);
        appmc->set_debug(1);
        appmc->set_dump_intermediary_cnf(std::max(dump_intermediary_cnf, 1));
    }

    if (!logfilename.empty()) {
        appmc->set_up_log(logfilename);
        cout << "c [appmc] Logfile set " << logfilename << endl;
    }
}

template<class T> void parse_file(const std::string& filename, T* reader) {
  #ifndef USE_ZLIB
  FILE * in = fopen(filename.c_str(), "rb");
  DimacsParser<StreamBuffer<FILE*, CMSat::FN>, T> parser(reader, nullptr, 0);
  #else
  gzFile in = gzopen(filename.c_str(), "rb");
  DimacsParser<StreamBuffer<gzFile, CMSat::GZ>, T> parser(reader, nullptr, 0);
  #endif
  if (in == nullptr) {
      std::cout << "ERROR! Could not open file '" << filename
      << "' for reading: " << strerror(errno) << endl;
      std::exit(-1);
  }
  if (!parser.parse_DIMACS(in, true)) exit(-1);
  #ifndef USE_ZLIB
  fclose(in);
  #else
  gzclose(in);
  #endif

  if (!reader->get_sampl_vars_set()) {
    all_indep = true;
    vector<uint32_t> tmp;
    for(uint32_t i = 0; i < reader->nVars(); i++) tmp.push_back(i);
    reader->set_sampl_vars(tmp); // will automatically set the opt_sampl_vars
  } else {
    // Check if CNF has all vars as indep. Then its's all_indep
    set<uint32_t> tmp;
    for(auto const& s: reader->get_sampl_vars()) {
      if (s >= reader->nVars()) {
        cout << "ERROR: Sampling var " << s+1 << " is larger than number of vars in formula: "
          << reader->nVars() << endl;
        exit(-1);
      }
      tmp.insert(s);
    }
    if (tmp.size() == reader->nVars()) all_indep = true;
    if (!reader->get_opt_sampl_vars_set()) {
      reader->set_opt_sampl_vars(reader->get_sampl_vars());
    }
  }
}

int main(int argc, char** argv)
{
    #if defined(__GNUC__) && defined(__linux__)
    feenableexcept(FE_INVALID   | FE_DIVBYZERO | FE_OVERFLOW);
    #endif
    double start_time = cpuTime();

    //Reconstruct the command line so we can emit it later if needed
    string command_line;
    for(int i = 0; i < argc; i++) {
        command_line += string(argv[i]);
        if (i+1 < argc) command_line += " ";
    }

    appmc = new ApproxMC::AppMC;
    add_supported_options(argc, argv);
    if (verb) {
        cout << appmc->get_version_info();
        cout << "c executed with command line: " << command_line << endl;
    }
    set_approxmc_options();

    ArjunNS::SimplifiedCNF cnf;
    auto files = program.get<std::vector<std::string>>("inputfile");
    string fname = files[0];
    if (do_arjun) {
        parse_file(fname, &cnf);
        const auto orig_sampl_vars = cnf.sampl_vars;
        double my_time = cpuTime();
        ArjunNS::Arjun arjun;
        arjun.set_verb(verb);
        arjun.set_or_gate_based(arjun_gates);
        arjun.set_xor_gates_based(arjun_gates);
        arjun.set_ite_gate_based(arjun_gates);
        arjun.set_irreg_gate_based(arjun_gates);
        arjun.only_backbone(cnf);
        arjun.only_run_minimize_indep(cnf);
        bool do_extend_indep = false;
        bool do_unate = false;
        all_indep = false;
        bool do_bce = false;
        int sbva_steps = 1000;
        int sbva_cls_cutoff = 4;
        int sbva_lits_cutoff = 5;
        int sbva_tiebreak = 1;
        if (with_e) {
            arjun.elim_to_file(cnf, all_indep, do_extend_indep, do_bce, do_unate, simp_conf, sbva_steps, sbva_cls_cutoff, sbva_lits_cutoff, sbva_tiebreak);
        }
        appmc->new_vars(cnf.nVars());
        appmc->set_sampl_vars(cnf.sampl_vars);
        for(const auto& c: cnf.clauses) appmc->add_clause(c);
        for(const auto& c: cnf.red_clauses) appmc->add_red_clause(c);
        appmc->set_multiplier_weight(cnf.multiplier_weight);
        print_final_indep_set(cnf.sampl_vars, orig_sampl_vars.size());
        cout << "c [arjun] Arjun finished. T: " << (cpuTime() - my_time) << endl;
    } else {
        parse_file(fname, appmc);
        print_final_indep_set(appmc->get_sampl_vars(), appmc->get_sampl_vars().size());
    }

    ApproxMC::SolCount sol_count;
    sol_count = appmc->count();
    appmc->print_stats(start_time);
    cout << "c [appmc+arjun] Total time: " << (cpuTime() - start_time) << endl;
    print_num_solutions(sol_count.cellSolCount, sol_count.hashCount, appmc->get_multiplier_weight());

    delete appmc;
}
