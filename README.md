# ApproxMC

ApproxMC2 is an approximate model counter for CNF formulas based on our IJCAI-16  paper: http://www.comp.nus.edu.sg/~meel/Papers/ijcai16_counting.pdf

The version 1 can be invoked by setting searchMode to 0 (Use "-h" for
more information when running approxmc binary).
The version 1 was based on our CP-13 paper. 


## Licensing

ApproxMC is released under MIT License. For more details, please see the file `LICENSE`.

## Building

You will need:

* automake
* autotools
* autoconfig
* make
* gcc

Then:

```
git clone https://github.com/meelgroup/approxmc
git checkout ver2
cd approxmc
cd approx
make -f Makefile.cvs
./configure
make
```

## Running

```
./approxmc --epsilon=0.8 --delta=0.2 --gaussuntil=400 myfile.cnf

```
The file ProbMapFile.txt should be in the same directory from which
approxmc is being invoked. 
Both epsilon and delta should be floats between 0 and 1.
Epsilon represents the accuracy of the result with respect to the
exact answer. (1-delta) is the confidence.

ApproxMC, run on a given CNF file, returns an approximate #SAT result, within epsilon percent of the exact count, with a confidence of (1-delta).


## Contact
Kuldeep Meel ([kuldeep@rice.edu](mailto:kuldeep@rice.edu))

