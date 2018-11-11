# ApproxMC Approximate Model Counter

ApproxMC is a state-of-the-art approximate model counter utilizing an improved version of CryptoMiniSat to give approximate model counts to problems of size and complexity that were not possible before. This work is by Kuldeep Meel and Mate Soos, as published in AAAI-19. A large part of the work is in CryptoMiniSat, here: https://github.com/msoos/cryptominisat


## How to Build
To build on Linux, you will need the following:
```
sudo apt-get install build-essential cmake
sudo apt-get install zlib1g-dev libboost-program-options-dev libm4ri-dev
```

Then, build CryptoMiniSat and ApproxMC:
```
git clone https://github.com/msoos/cryptominisat
cd cryptominisat
mkdir build && cd build
cmake cmake -DUSE_GAUSS=ON ..
make
sudo make install
cd ../..
git clone https://github.com/meelgroup/approxmc/
cd approxmc
mkdir build && cd build
cmake ..
make
sudo make install
```

## How to Use

First, you must translate your problem to CNF. ApproxMC takes in a modified version of DIMACS where you can specify your sampling set like this:

```
$ cat myfile.cnf
c ind 1 3 4 6 7 8 10 0
p cnf 500 1
3 4 0
```
Above, using the `c ind` line, we declare that only variables 1, 3, 4, 6, 7, 8 and 10 form part of the sampling set out of the CNF's 500 variables `1,2...500`. This line must end with a 0. The solution that ApproxMC will be giving is essentially answering the question: how many different combination of settings to this variables are there that satisfy this problem? Naturally, if your sampling set only contains 7 variables, then the maximum number of solutions can only be at most 2^7 = 128. This is true even if your CNF has thousands of variables.

In our case, the maximum number of solutions could at most be 2^7=128, but our CNF should be restricting this. Let's see:

```
$ approxmc --seed 5 myfile.cnf
c ApproxMC version 2.5
[...]
c CryptoMiniSat SHA revision 17a1aed4956848404e33d514eef257ca1ed2382b
c CMS is MIT licensed
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
[appmc] Number of solutions is: 48 x 2^1
```
ApproxMC reports that have approximately `96 (=48*2)` solutions to the CNF's independent support. This is because for variables 3 and 4 we have banned the `false,false` solution, so out of their 4 possible settings, one is banned. Therefore, we have `2^5 * (4-1) = 96` solutions.

## Issues, Bugs, Wishes, etc

In case you encounter a bug, or a problem solving is very slow, etc. please create a [GitHub issue](https://github.com/meelgroup/approxmc/issues). We will deal with *every* case, and will resolve them all. The authors of this code have resolved over 500 GitHub issues on other software they release. All issues will be dealt with promptly.

## Old Versions

The old version, 2.0 is available under the branch "ver2". Please check out the releases for the 2.x versions under GitHub [releases](https://github.com/meelgroup/approxmc/releases). Please read the README of the old release to know how to compile the code. Old releases should easily compile.
