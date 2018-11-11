# ApproxMC #

### What is it? ###
ApproxMC2 is an approximate model counter for CNF formulas based on our IJCAI-16  paper: http://www.comp.nus.edu.sg/~meel/Papers/ijcai16_counting.pdf

The version 1 can be invoked by setting searchMode to 0 (Use "-h" for
more information when running approxmc binary).
The version 1 was based on our CP-13 paper. 

For more details on our hashing-based approach to sampling and counting, please visit: http://www.kuldeepmeel.com

### The Latest Version ###
The latest version can be downloaded from [https://bitbucket.org/kuldeepmeel/approxmc](https://bitbucket.org/kuldeepmeel/approxmc).


### Licensing ###
ApproxMC is released under MIT License. For more details, please see the file `LICENSE`. 

## USAGE ##

```
#!shell

 ./approxmc --epsilon=0.8 --delta=0.2 --gaussuntil=400 <input file>

```
The file ProbMapFile.txt should be in the same directory from which
approxmc is being invoked. 
Both epsilon and delta should be floats between 0 and 1.
Epsilon represents the accuracy of the result with respect to the
exact answer. (1-delta) is the confidence.

ApproxMC, run on a given CNF file, returns an approximate #SAT result, within epsilon percent of the exact count, with a confidence of (1-delta).


### Questions/Feedback/Comments ###
Report issues using bitbucket. Please do not email me. 

### Contact ###


  1. Kuldeep Meel ([kuldeep@rice.edu](mailto:kuldeep@rice.edu))


Enjoy!
