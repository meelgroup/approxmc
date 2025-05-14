[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![build](https://github.com/meelgroup/approxmc/workflows/build/badge.svg)

# ApproxMC6: Approximate Model Counter
ApproxMCv6 is a state-of-the-art approximate model counter using
[Arjun](https://github.com/meelgroup/arjun) and
[CryptoMiniSat](https://github.com/msoos/cryptominisat) to give probabilistic
approximate model counts to problems of size and complexity that were not
possible before.

This work is the culmination of work by a number of people, including but not
limited to, Mate Soos, Jiong Yang, Stephan Gocht, Yash Pote, and Kuldeep S.
Meel. Publications: published [in
AAAI-19](https://www.cs.toronto.edu/~meel/Papers/aaai19-sm.pdf), [in
CAV2020](https://www.cs.toronto.edu/~meel/Papers/cav20-sgm.pdf), and [in
CAV2023](https://arxiv.org/pdf/2305.09247).

ApproxMC handles CNF formulas and performs approximate counting.
1. If you are interested in exact model counting, visit our exact counter
   [Ganak](http://github.com/meelgroup/ganak)
2. If you are instead interested in DNF formulas, visit our approximate DNF
   counter [Pepin](https://github.com/meelgroup/pepin).

## Installation
We recommend using a prebuilt binary from our [release
page](https://github.com/meelgroup/approxmc/releases) which contains binaries
for many different platforms. The second best thing to use is Nix. Simply [install
nix](https://nixos.org/download/) and then:
```shell
git clone https://github.com/meelgroup/approxmc
cd approxmc
nix-shell
```

Then you will have `approxmc` binary available and ready to use.

## Providing a Projection Set
For some applications, one is not interested in solutions over all the
variables and instead interested in counting the number of unique solutions to
a subset of variables, called sampling set (also called a "projection set").
ApproxMC allows you to specify the sampling set using the following modified
version of DIMACS format:

```plain
$ cat myfile.cnf
c p show 1 3 4 6 7 8 10 0
p cnf 500 1
3 4 0
```
Above, using the `c p show` line, we declare that only variables 1, 3, 4, 6, 7,
8 and 10 form part of the sampling set out of the CNF's 500 variables
`1,2...500`. This line must end with a 0. The solution that ApproxMC will be
giving is essentially answering the question: how many different combination of
settings to this variables are there that satisfy this problem? Naturally, if
your sampling set only contains 7 variables, then the maximum number of
solutions can only be at most 2^7 = 128. This is true even if your CNF has
thousands of variables.

## Running ApproxMC
In our case, the maximum number of solutions could at most be 2^7=128, but our
CNF should be restricting this. Let's see:

```plain
$ approxmc --seed 5 myfile.cnf
c ApproxMC version 3.0
[...]
c CryptoMiniSat SHA revision [...]
c Using code from 'When Boolean Satisfiability Meets Gauss-E. in a Simplex Way'
[...]
[appmc] using seed: 5
[appmc] Sampling set size: 7
[appmc] Sampling set: 1, 3, 4, 6, 7, 8, 10,
[appmc] Using start iteration 0
[appmc] [    0.00 ] bounded_sol_count looking for   73 solutions -- hashes active: 0
[appmc] [    0.01 ] bounded_sol_count looking for   73 solutions -- hashes active: 1
[appmc] [    0.01 ] bounded_sol_count looking for   73 solutions -- hashes active: 0
[...]
[appmc] FINISHED ApproxMC T: 0.04 s
c [appmc] Number of solutions is: 48*2**1
s mc 96
```
ApproxMC reports that we have approximately `96 (=48*2)` solutions to the CNF's
independent support. This is because for variables 3 and 4 we have banned the
`false,false` solution, so out of their 4 possible settings, one is banned.
Therefore, we have `2^5 * (4-1) = 96` solutions.

## Guarantees
ApproxMC provides so-called "PAC", or Probably Approximately Correct,
guarantees. In less fancy words, the system guarantees that the solution found
is within a certain tolerance (called "epsilon") with a certain probability
(called "delta"). The default tolerance and probability, i.e. epsilon and delta
values, are set to 0.8 and 0.2, respectively. Both values are configurable.

## How to use the Python interface
Install using pip:
```bash
pip install pyapproxmc
```

Then you can use it as:
```python
import pyapproxmc
c = pyapproxmc.Counter()
c.add_clause([1,2,3])
c.add_clause([3,20])
count = c.count()
print("Approximate count is: %d*2**%d" % (count[0], count[1]))
```

The above will print that `Approximate count is: 11*2**16`. Since the largest
variable in the clauses was 20, the system contained 2\*\*20 (i.e. 1048576)
potential models. However, some of these models were prohibited by the two
clauses, and so only approximately 11*2\*\*16 (i.e. 720896) models remained.

If you want to count over a projection set, you need to call
`count(projection_set)`, for example:
```python
import pyapproxmc
c = pyapproxmc.Counter()
c.add_clause([1,2,3])
c.add_clause([3,20])
count = c.count(range(1,10))
print("Approximate count is: %d*2**%d" % (count[0], count[1]))
```

This now prints `Approximate count is: 7*2**6`, which corresponds to the
approximate count of models, projected over variables 1..10.

### Library usage
The system can be used as a library:

```c++
#include <approxmc/approxmc.h>
#include <vector>
#include <cassert>

using std::vector;
using namespace ApproxMC;
using namespace CMSat;

int main() {
    AppMC appmc;
    appmc.new_vars(10);

    vector<Lit> clause;

    //add "-3 4 0"
    clause.clear();
    clause.push_back(Lit(2, true));
    clause.push_back(Lit(3, false));
    appmc.add_clause(clause);

    //add "3 -4 0"
    clause.clear();
    clause.push_back(Lit(2, false));
    clause.push_back(Lit(3, true));
    appmc.add_clause(clause);

    SolCount c = appmc.count();
    uint32_t cnt = std::pow(2, c.hashCount)*c.cellSolCount;
    assert(cnt == std::pow(2, 9));

    return 0;
}
```

## ApproxMC5: Sparse-XOR based Approximate Model Counter
Note: this is beta version release, not recommended for general use. We are
currently working on a tight integration of sparse XORs into ApproxMC based on
our [LICS-20](http://www.cs.toronto.edu/~meel/Papers/lics20-ma.pdf) paper. You
can turn on the sparse XORs using the flag `--sparse 1` but beware as reported in
LICS-20 paper, this may slow down solving in some cases. It is likely to give a
significant speedup if the number of solutions is very large.


## Issues, questions, bugs, etc.
Please click on "issues" at the top and [create a new
issue](https://github.com/meelgroup/mis/issues/new). All issues are responded
to promptly.

## How to Cite
If you use ApproxMC, please cite the following papers:
[AAAI-19](https://www.cs.toronto.edu/~meel/Papers/aaai19-sm.pdf), [in
CAV2020](https://www.cs.toronto.edu/~meel/Papers/cav20-sgm.pdf), and [in
CAV2023](https://arxiv.org/pdf/2305.09247).
[CAV20](https://dblp.uni-trier.de/rec/conf/cav/SoosGM20.html?view=bibtex),
[AAAI19](https://www.cs.toronto.edu/~meel/bib/SM19.bib) and
[IJCAI16](https://www.cs.toronto.edu/~meel/bib/CMV16.bib). If you use sparse
XORs, please also cite the
[LICS20](https://www.cs.toronto.edu/~meel/publications/AM20.bib) paper.
ApproxMC builds on a series of papers on hashing-based approach: [Related
Publications](https://www.cs.toronto.edu/~meel/publications.html)

The benchmarks used in our evaluation can be found [here](https://zenodo.org/records/10449477).
