# pyapproxmc: bindings to the ApproxMC model counter

This directory provides Python bindings to ApproxMC on the C++ level,
i.e. when importing pyapproxmc, the ApproxMC counter becomes part of the
Python process itself.


## Installing

```
pip install pyapproxmc
```

## Compiling
If you don't want to use the pip package, you can compile it:

```
apt-get install python-dev
cd python
git clone https://github.com/msoos/cryptominisat
git clone https://github.com/meelgroup/arjun
cd ..
python -m build

You will then find the files under "dist/".
```

## Usage

```
import pyapproxmc
c = pyapproxmc.Counter()
c.add_clause([1,2,3])
c.add_clause([3,20])
count = c.count()
print("Approximate count is: %d*2**%d" % (count[0], count[1]))
```

The above will print that `Approximate count is: 88*2**13`. Since the largest variable in the clauses was 20, the system contained 2**20 (i.e. 1048576) potential models. However, some of these models were prohibited by the two clauses, and so only approximately 88*2**13 (i.e. 720896) models remained.

If you want to count over a projection set, you need to call `count(projection_set)`, for example:

```
import pyapproxmc
c = pyapproxmc.Counter()
c.add_clause([1,2,3])
c.add_clause([3,20])
count = c.count(range(1,10))
print("Approximate count is: %d*2**%d" % (count[0], count[1]))
```

This now prints `Approximate count is: 56*2**3`, which corresponds to the approximate count of models, projected over variables 1..10.

## Counter Object

You can give the following arguments to `Counter`:
* `seed` -- sets the random seed
* `verbosity` -- sets the verbosity of the system (default = 0)
* `epsilon` -- Tolerance parameter, i.e. sets how approximate the returned count is. Default = 0.8
* `delta` -- Confidence parameter, i.e. sets how probabilistically correct the returned count is. Default = 0.20

