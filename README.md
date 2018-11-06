# ScalMC, an state-of-the-art approximate Model Counter

ScalMC is a state-of-the-art approximate model counter that uses an improved version of CryptoMiniSat to give approximate model counts to problems of size and complexity that were not possible before. This work is by Kuldeep Meel and Mate Soos, as published in AAAI-19. A large part of the work is in CryptoMiniSat, here: https://github.com/msoos/cryptominisat


## How to build
To build on Linux, you will need the following:
```
sudo apt-get install build-essential cmake
sudo apt-get install zlib1g-dev libboost-program-options-dev libm4ri-dev
```

Then, build CryptoMiniSat and ScalMC:
```
git clone https://github.com/msoos/cryptominisat
cd cryptominisat
mkdir build && cd build
cmake cmake -DUSE_GAUSS=ON ..
make
sudo make install
cd ../..
git clone https://github.com/meelgroup/scalmc/
cd scalmc
mkdir build && cd build
cmake ..
make
sudo make install
```

## How to use

First, you must translate your problem to CNF. ScalMC takes in a modified version of DIMACS where you can specify your independent support (essentially, independent variables) like this:

```
$ cat myfile.cnf
c ind 1 3 4 0
p cnf 5 2
1 2 0
3 4 0
4 5 0
```
Where we declare that only variables 1, 3, and 4 form part of the independent support with the `c ind` line. This line must end with a 0. The solution that ScalMC will be giving is essentially answering the question: how many different combination of settings to this variables are there that satisfy this problem? Naturally, if you independent support only contains 5 variables, then the maximum number of solutions can only be at most 2^5 = 32. This is true even if your CNF has thousands of variables.

In our case, the maximum number of solutions would be 2^3 = 8, but our CNF should be restricting this. Let's see:

```
4 scalmc --seed 5 myfile.cnf
c ScalMC SHA revision ea21bfaaa97cf2aa6d7864083cf9597848202f39
[...]
c CryptoMiniSat SHA revision 17a1aed4956848404e33d514eef257ca1ed2382b
[scalmc] using seed: 5
[scalmc] Num independent vars: 17
[scalmc] Independent vars: 10, 13, 15, 16, 25, 28, 39, 41, 43, 45, 5, 53, 6, 69, 78, 9, 93, 
[scalmc] Using start iteration 0
[scalmc] [    0.00 ] bounded_sol_count looking for   73 solutions -- hashes active: 0
[scalmc] [    0.01 ] bounded_sol_count looking for   73 solutions -- hashes active: 1
[scalmc] [    0.01 ] bounded_sol_count looking for   73 solutions -- hashes active: 2
[...]
[scalmc] FINISHED ScalMC T: 0.22 s
[scalmc] Number of solutions is: 65 x 2^8
```

For the same seed you will always get the same results. You can check all the options with -h.
