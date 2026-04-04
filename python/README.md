# pyapproxmc: Python bindings to ApproxMC

Python bindings to [ApproxMC](https://github.com/meelgroup/approxmc), an
approximate model counter with PAC (Probably Approximately Correct) guarantees.
ApproxMC counts the number of satisfying assignments of CNF formulas.

## Installing

```bash
pip install pyapproxmc
```

## Building from source

```bash
# Install build tools
pip install scikit-build-core build

# Clone the repository
git clone --recurse-submodules https://github.com/meelgroup/approxmc
cd approxmc

# Build and install into a virtual environment
python -m venv .venv && source .venv/bin/activate
pip install .
```

The build system uses scikit-build-core and CMake. All C++ dependencies
(CryptoMiniSat, Arjun, and their sub-dependencies) are fetched and built
automatically by CMake during the `pip install` step. The only system library
required is [GMP](https://gmplib.org/):

- **Linux**: `sudo apt-get install libgmp-dev`  (or `dnf install gmp-devel`)
- **macOS**: `brew install gmp`

## Usage

```python
import pyapproxmc

c = pyapproxmc.Counter()
c.add_clause([1, 2, 3])
c.add_clause([3, 20])
count = c.count()
print("Approximate count is: %d*2**%d" % (count[0], count[1]))
```

`count()` may only be called **once** per `Counter` instance; create a new
`Counter` if you need to count a different formula.

The above prints `Approximate count is: 11*2**16`. Since the largest variable
in the clauses is 20, the formula has at most 2\*\*20 models; the two clauses
restrict this to approximately 11\*2\*\*16 ≈ 720896 models.

### Counting over a projection set

Pass a projection set (sampling set) to `count()` to count models projected
onto a subset of variables:

```python
import pyapproxmc

c = pyapproxmc.Counter()
c.add_clause([1, 2, 3])
c.add_clause([3, 20])
count = c.count(range(1, 10))
print("Approximate count is: %d*2**%d" % (count[0], count[1]))
```

This prints `Approximate count is: 7*2**6`, the approximate count projected
over variables 1–9.

### Adding clauses from arrays

For performance-critical code, multiple clauses can be added in one call via
`add_clauses()` with a flat, zero-terminated `array.array`. This uses the
buffer protocol and avoids per-element Python overhead:

```python
import pyapproxmc
from array import array

c = pyapproxmc.Counter()
# Clauses [1,2,3] and [3,20], each terminated by 0
c.add_clauses(array('i', [1, 2, 3, 0, 3, 20, 0]))
count = c.count()
```

Note: `add_clause()` (singular) also accepts any iterable including arrays,
but iterates element by element and does not use the buffer protocol.

## Counter constructor parameters

| Parameter   | Default | Description |
|-------------|---------|-------------|
| `seed`      | 1       | Random seed for reproducibility |
| `verbosity` | 0       | Output verbosity (0 = silent) |
| `epsilon`   | 0.8     | Tolerance: how approximate the count is |
| `delta`     | 0.2     | Confidence: probability the count is within tolerance |

Example:

```python
c = pyapproxmc.Counter(seed=42, epsilon=0.5, delta=0.1)
```

## Version

```python
import pyapproxmc
print(pyapproxmc.__version__)
```
