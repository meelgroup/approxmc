[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![build](https://github.com/meelgroup/approxmc/workflows/build/badge.svg)

# ApproxMC7: Approximate Model Counter
ApproxMCv7 is a state-of-the-art approximate model counter using
[Arjun](https://github.com/meelgroup/arjun) and
[CryptoMiniSat](https://github.com/msoos/cryptominisat) to give probabilistic
approximate model counts to problems of size and complexity that were not
possible before.

This work is the culmination of work by a number of people. See publication list at the end of this README for more details.

ApproxMC handles CNF formulas and performs approximate counting.
1. If you are interested in exact model counting, visit our exact counter
   [Ganak](http://github.com/meelgroup/ganak)
2. If you need to count a weighted CNF formula, you need to preprocess your CNF
   using [our tool](https://github.com/meelgroup/weighted-to-unweighted) to
   convert it to an unweighted CNF formula. Then you can use ApproxMC to count it.
3. If you are interested in DNF formulas, visit our approximate DNF
   counter [Pepin](https://github.com/meelgroup/pepin).

Notice that for some formula families, Ganak is faster than ApproxMC. It
can be worthwhile to try both tools on your instances.

## Installation
You can try out ApproxMC [from your browser](https://www.msoos.org/approxmc/).

It is strongly recommended to not build, but to use the precompiled
binaries as in our [release](https://github.com/meelgroup/approxmc/releases).
The second best thing to use is Nix. Simply [install
nix](https://nixos.org/download/) and then:
```shell
nix shell github:meelgroup/approxmc#approxmc
```

Then you will have `approxmc` binary available and ready to use.

### Building statically

To build a static binary, you first need to build GMP with position-independent code
enabled. GMP's hand-optimised assembly is normally compiled without `-fPIC` (fine for
a native static binary), but `-fPIC` is required whenever the static `.a` is linked
into a shared object — for example, the Python extension (`.so`). Passing `--with-pic`
to GMP's `configure` makes both uses work:

```shell
wget https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz
tar xf gmp-6.3.0.tar.xz
cd gmp-6.3.0
./configure --enable-static --enable-cxx --enable-shared --with-pic
make -j8
sudo make install
cd ..
```

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

The legacy directive `c ind 1 3 4 6 7 8 10 0` declares the same thing.
ApproxMC accepts both forms but they cannot appear together in the same file;
prefer `c p show` in new files.

## Running ApproxMC
In our case, the maximum number of solutions could at most be 2^7=128, but our
CNF should be restricting this. Let's see:

```plain
$ approxmc --seed 5 myfile.cnf
c o CMS SHA1: ...
c o Arjun SHA1: ...
c o ApproxMC SHA1: ...
c o Using code from 'When Boolean Satisfiability Meets Gauss-E. in a Simplex Way'
[...]
c o [arjun-simp] Removed set       : 0 new size: 7
c o [arjun-simp] get-empties removed: 5 perc: 71.43 total empties now: 5 T: 0.00
[...]
c o [appmc] Sampling set size: 2
c o [appmc] Sampling set: 3 4 0
c o [appmc] threshold set to 72 sparse: 0
c o [appmc] [    0.00 ] bounded_sol_count looking for   73 solutions -- hashes active: 0
[...]
c o [appmc] Counted without XORs, i.e. we got exact count
c o [appmc] ApproxMC T: 0.00 s
c [appmc] Number of solutions is: 3*2**0*32
s SATISFIABLE
s mc 96
```
ApproxMC reports `96` solutions on the final `s mc` line. The
`Number of solutions is: 3*2**0*32` line above it is in the format
`<cell>*2**<hash>*<mult>`, where `<mult>` is a multiplier produced by Arjun's
preprocessing pass (enabled by default). In this run Arjun shrank the
7-variable sampling set down to just `{3, 4}` — the only ones the clause
constrains — so the 5 removed variables contribute `2^5 = 32` as the
multiplier. The remaining two variables have one banned setting
(`false, false`) out of four, hence `cell = 3`. The final count is
`3 * 2^0 * 32 = 96`. Passing `--arjun 0` would skip this minimization.

### Common command-line options

Run `approxmc --help` for the full list. The most commonly used flags:

| Flag | Default | Description |
|------|---------|-------------|
| `-s, --seed <N>` | 1 | Random seed |
| `-e, --epsilon <F>` | 0.8 | Tolerance: count is in `[exact/(1+e), exact*(1+e)]` |
| `-d, --delta <F>` | 0.2 | Confidence: count is within tolerance with prob. `1-d` |
| `-v, --verb <N>` | 1 | Verbosity (0 = quiet) |
| `--sparse <0/1>` | 0 | Generate sparse XOR constraints (see below) |
| `--version` |  | Print version info and exit |

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

Or build from source (requires GMP: `apt-get install libgmp-dev` / `brew install gmp`):
```bash
git clone --recurse-submodules https://github.com/meelgroup/approxmc
cd approxmc
python -m venv venv
venv/bin/pip install .
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

Note: `count()` may only be called once per `Counter` instance.

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
approximate count of models, projected over variables 1..9.

See [python/README.md](python/README.md) for the full Python API documentation.

### Library usage
The system can be used as a library:

```c++
#include <approxmc/approxmc.h>
#include <arjun/arjun.h>
#include <cryptominisat5/solvertypesmini.h>
#include <cassert>
#include <cmath>
#include <memory>
#include <vector>

using std::vector;
using namespace ApproxMC;
using namespace CMSat;

int main() {
    // AppMC takes a FieldGen used by the multiplier-weight machinery.
    // For unweighted counting, use Arjun's mpq-backed generator.
    std::unique_ptr<FieldGen> fg = std::make_unique<ArjunNS::FGenMpq>();
    AppMC appmc(fg);
    appmc.new_vars(10);

    vector<Lit> clause;

    //add "-3 4 0"
    clause = {Lit(2, true), Lit(3, false)};
    appmc.add_clause(clause);

    //add "3 -4 0"
    clause = {Lit(2, false), Lit(3, true)};
    appmc.add_clause(clause);

    // The sampling set must be set explicitly (no default).
    // Here we count over all 10 variables.
    vector<uint32_t> sampl;
    for (uint32_t i = 0; i < 10; i++) sampl.push_back(i);
    appmc.set_sampl_vars(sampl);

    SolCount c = appmc.count();
    uint64_t cnt = (uint64_t)std::pow(2, c.hashCount) * c.cellSolCount;
    assert(cnt == (uint64_t)std::pow(2, 9));

    return 0;
}
```

Link against `libapproxmc`, `libarjun`, and `libcryptominisat5`. The
`<approxmc/approxmc.h>` header is the only public ApproxMC header you need;
`<arjun/arjun.h>` is required for the `FGenMpq` field generator passed to the
`AppMC` constructor.

## Sparse XOR Counting
You can turn on the sparse XORs using the flag `--sparse 1` but beware as reported in
[LICS-20](http://www.cs.toronto.edu/~meel/Papers/lics20-ma.pdf) paper,
this may slow down solving in some cases. It is likely to give a
significant speedup if the number of solutions is very large.

## How to Cite
If you use ApproxMC, please cite the following papers:
[AAAI-25](https://ojs.aaai.org/index.php/AAAI/article/view/33231)
, [CAV-23](https://arxiv.org/pdf/2305.09247)
, [AAAI-19](https://www.cs.toronto.edu/~meel/Papers/aaai19-sm.pdf)
, [CAV-20](https://www.cs.toronto.edu/~meel/Papers/cav20-sgm.pdf)
, [CAV-20](https://dblp.uni-trier.de/rec/conf/cav/SoosGM20.html?view=bibtex)
, [LICS-20](https://www.cs.toronto.edu/~meel/publications/AM20.bib)
, [AAAI-19](https://www.cs.toronto.edu/~meel/publications/SM19.bib)
, [IJCAI-16](https://www.cs.toronto.edu/~meel/publications/CMV16.bib)
, [CP-13](https://www.cs.toronto.edu/~meel/publications/CMV13b.bib)

Additional related publications can be found [here](https://www.cs.toronto.edu/~meel/publications.html)

## Benchmarks Used
The benchmarks used in our evaluation can be found [here](https://zenodo.org/records/10449477)
