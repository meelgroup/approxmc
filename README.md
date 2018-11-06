# ScalMC -- an state-of-the-art approximate Model Counter

ScalMC is a state-of-the-art approximate model counter that uses an improved version of CryptoMiniSat to give approximate model counts to problems of size and complexity that were not possible before.

To build perform the following:
```
git clone https://github.com/msoos/cryptominisat
cd cryptominisat
mkdir build && cd build
cmake cmake -DUSE_GAUSS=ON ..
sudo make install
cd ../..
git clone https://github.com/meelgroup/scalmc/
cd scalmc
mkdir build && cd build
cmake ..
make
```

Then run ScalMC like:

```
./scalmc --seed 5 myfile.cnf
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
