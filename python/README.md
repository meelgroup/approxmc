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

The `pyapproxmc` module has one object, `Counter` that has two functions
`count` and `add_clause`.

The function ``add_clause()`` takes an iterable list of literals such as
``[1, 2]`` which represents the truth ``1 or 2 = True``. For example,
``add_clause([1])`` sets variable ``1`` to ``True``.

The function `count()` counts the number of solutions to the system of constraints
that have been added with `add_clause()`:

```
>>> from pyapproxmc import Counter
>>> s = Counter()
>>> s.add_clause([1, 2])
>>> cells, hashes = s.count()
>>> print("There are approx ", cells*2**hashes, " solutions")
There are 55 solutions, approximately
```

The return value is a tuple of cells and hashes. Which gives how many solutions
there are, probabilistically approximately

You can give the following arguments to `Counter`:
* `seed` -- sets the random seed
* `verbosity` -- sets the verbosity of the system (default = 0)
* `epsilon` -- Tolerance parameter, i.e. sets how approximate the returned count is. Default = 0.8
* `delta` -- Confidence parameter, i.e. sets how probabilistically correct the returned count is. Default = 0.20

