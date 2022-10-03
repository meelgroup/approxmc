# pyapproxmc: bindings to the ApproxMC model counter

This directory provides Python bindings to CryptoMiniSat on the C++ level,
i.e. when importing pycryptosat, the CryptoMiniSat solver becomes part of the
Python process itself.

## Compiling
The pyapproxmc python package compiles separately from ApproxMC, the binary.

In order to compile, install the python developer tools:

```
apt-get install python-dev
```

Then:

```
cd python
git clone https://github.com/msoos/cryptominisat
git clone https://github.com/meelgroup/arjun
cd ..
python -m build

You will then find the files under "dist/".
```

## Usage

The ``pyapproxmc`` module has one object, ``Counter`` that has two functions
``count`` and ``add_clause``.

The funcion ``add_clause()`` takes an iterable list of literals such as
``[1, 2]`` which represents the truth ``1 or 2 = True``. For example,
``add_clause([1])`` sets variable ``1`` to ``True``.

The function ``count()`` solves the system of equations that have been added
with ``add_clause()``:

```
   >>> from pyapproxmc import Counter
   >>> s = Counter()
   >>> s.add_clause([1, 2])
   >>> cells, hashes = s.count()
   >>> print "There are ", cells*2**hashes, " solutions, approximately"
   There are 55 solutions, approximately
```

The return value is a tuple of cells and hashes. Which gives how many solutions
there are, probabilistically approximately

You can give the following arguments to `Counter`:
* `seed` -- sets the random seed
* `verbosity` -- sets the verbosity of the system (default = 0)
* `epsilon` -- Tolerance parameter, i.e. sets how approximate the returned count is. Default = 0.8
* `delta` -- Confidence parameter, i.e. sets how probabilistically correct the returned count is. Default = 0.20

