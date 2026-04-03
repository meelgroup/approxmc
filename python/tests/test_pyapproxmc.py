#!/usr/bin/env python3
import sys
from array import array
from pathlib import Path

import pytest

import pyapproxmc
from pyapproxmc import Counter

# Seed-pinned tests below are intentionally locked to specific seeds; they
# serve as regression guards.  If the counting algorithm changes the expected
# values must be updated deliberately.


# ---------------------------------------------------------------------------
# Module-level attributes
# ---------------------------------------------------------------------------

def test_version_string():
    assert isinstance(pyapproxmc.__version__, str)
    assert isinstance(pyapproxmc.VERSION, str)
    assert pyapproxmc.__version__ == pyapproxmc.VERSION
    parts = pyapproxmc.__version__.split(".")
    assert len(parts) == 3, f"Expected 3 version parts, got: {pyapproxmc.__version__}"
    assert all(p.isdigit() for p in parts), f"Non-numeric version part: {pyapproxmc.__version__}"


# ---------------------------------------------------------------------------
# Constructor validation
# ---------------------------------------------------------------------------

def test_invalid_epsilon_zero():
    with pytest.raises(ValueError, match="epsilon"):
        Counter(epsilon=0)


def test_invalid_epsilon_negative():
    with pytest.raises(ValueError, match="epsilon"):
        Counter(epsilon=-1.0)


def test_invalid_delta_negative():
    with pytest.raises(ValueError, match="delta"):
        Counter(delta=-0.1)


def test_invalid_delta_one():
    with pytest.raises(ValueError, match="delta"):
        Counter(delta=1.0)


def test_invalid_delta_greater_than_one():
    with pytest.raises(ValueError, match="delta"):
        Counter(delta=2.0)


def test_invalid_verbosity():
    with pytest.raises(ValueError, match="verbosity"):
        Counter(verbosity=-1)


# ---------------------------------------------------------------------------
# add_clause input validation
# ---------------------------------------------------------------------------

def test_add_clause_zero_literal():
    c = Counter(seed=1)
    with pytest.raises(ValueError, match="non-zero"):
        c.add_clause([0])


def test_add_clause_zero_literal_in_middle():
    c = Counter(seed=1)
    with pytest.raises(ValueError, match="non-zero"):
        c.add_clause([1, 0, 2])


def test_add_clause_wrong_type():
    c = Counter(seed=1)
    with pytest.raises(TypeError):
        c.add_clause(["a"])


def test_add_clause_non_iterable():
    c = Counter(seed=1)
    with pytest.raises(TypeError):
        c.add_clause(42)


# ---------------------------------------------------------------------------
# count() called twice must raise; re-init must reset that flag
# ---------------------------------------------------------------------------

def test_count_twice_raises():
    c = Counter(seed=1)
    c.add_clause([1, 2])
    c.count()
    with pytest.raises(ValueError, match="once"):
        c.count()


def test_reinit_resets_count_called():
    c = Counter(seed=1)
    c.add_clause([1, 2])
    c.count()
    # Re-initialise via __init__ — count_called must be cleared
    Counter.__init__(c, seed=1)
    c.add_clause([1, 2])
    result = c.count()
    assert result is not None


# ---------------------------------------------------------------------------
# UNSAT formula / trivially SAT / empty projection
# ---------------------------------------------------------------------------

def test_unsat():
    c = Counter(seed=1)
    c.add_clause([1])
    c.add_clause([-1])
    cell_count, _ = c.count()
    assert cell_count == 0


def test_empty_formula():
    # No clauses added: trivially satisfiable over an empty variable set.
    c = Counter(seed=1)
    cell_count, hash_count = c.count()
    assert cell_count == 1
    assert hash_count == 0


def test_empty_projection():
    # Projecting onto zero variables: exactly one assignment (the empty one).
    c = Counter(seed=1)
    c.add_clause([1, 2])
    cell_count, hash_count = c.count([])
    assert cell_count == 1
    assert hash_count == 0


# ---------------------------------------------------------------------------
# Functional correctness (deterministic with fixed seed)
# ---------------------------------------------------------------------------

def test_minimal():
    counter = Counter(seed=2157, epsilon=0.8, delta=0.2)
    counter.add_clause(list(range(1, 100)))

    significand, exponent = counter.count()
    assert significand * 2**exponent == 512 * 2**90


def test_sampling_set():
    counter = Counter(seed=2157, epsilon=0.8, delta=0.2)
    counter.add_clause(range(1, 100))
    assert counter.count(list(range(1, 50))) == (64, 43)


def test_projection_keyword_argument():
    # count() must also accept 'projection' as a keyword argument
    counter = Counter(seed=2157, epsilon=0.8, delta=0.2)
    counter.add_clause(range(1, 100))
    assert counter.count(projection=list(range(1, 50))) == (64, 43)


def test_real_example():
    counter = Counter(seed=120, epsilon=0.8, delta=0.2)

    cnf_file = Path(__file__).parent / "test_1.cnf"
    with open(cnf_file) as test_cnf:
        # Skip sampling set and header lines
        lines = test_cnf.readlines()[2:]

        for line in lines:
            literals = [int(i) for i in line.split()[:-1]]
            counter.add_clause(literals)

    assert counter.count(list(range(1, 21))) == (64, 14)


# ---------------------------------------------------------------------------
# add_clauses — array (buffer) path
# ---------------------------------------------------------------------------

def test_add_clauses_array_minimal():
    counter = Counter(seed=2157, epsilon=0.8, delta=0.2)
    clauses = array('i', list(range(1, 100)) + [0])
    counter.add_clauses(clauses)

    significand, exponent = counter.count()
    assert significand * 2**exponent == 512 * 2**90


def test_add_clauses_array_real_example():
    counter = Counter(seed=120, epsilon=0.8, delta=0.2)
    flat = []

    cnf_file = Path(__file__).parent / "test_1.cnf"
    with open(cnf_file) as test_cnf:
        lines = test_cnf.readlines()[2:]
        for line in lines:
            flat += [int(i) for i in line.split()]

    counter.add_clauses(array('i', flat))
    significand, exponent = counter.count(list(range(1, 21)))
    assert significand * 2**exponent == 64 * 2**14


def test_add_clauses_array_last_not_terminated():
    c = Counter(seed=1)
    with pytest.raises(ValueError, match="zero"):
        c.add_clauses(array('i', [1, 2, 3]))  # missing trailing 0


def test_add_clauses_array_long():
    # 'l' (signed long) format path in _add_clauses_from_buffer
    counter = Counter(seed=2157, epsilon=0.8, delta=0.2)
    counter.add_clauses(array('l', list(range(1, 100)) + [0]))
    significand, exponent = counter.count()
    assert significand * 2**exponent == 512 * 2**90


def test_add_clauses_array_longlong():
    # 'q' (signed long long) format path in _add_clauses_from_buffer
    counter = Counter(seed=2157, epsilon=0.8, delta=0.2)
    counter.add_clauses(array('q', list(range(1, 100)) + [0]))
    significand, exponent = counter.count()
    assert significand * 2**exponent == 512 * 2**90


def test_add_clauses_array_unsupported_format():
    c = Counter(seed=1)
    with pytest.raises(ValueError, match="format"):
        c.add_clauses(array('f', [1.0, 0.0]))



# ---------------------------------------------------------------------------
# add_clauses — iterable-of-iterables (list of lists) path
# ---------------------------------------------------------------------------

def test_add_clauses_list_of_lists():
    # Same formula as test_minimal but via add_clauses([[...]])
    counter = Counter(seed=2157, epsilon=0.8, delta=0.2)
    counter.add_clauses([list(range(1, 100))])
    significand, exponent = counter.count()
    assert significand * 2**exponent == 512 * 2**90


def test_add_clauses_list_of_lists_multiple():
    counter = Counter(seed=2157, epsilon=0.8, delta=0.2)
    counter.add_clauses([list(range(1, 50)), list(range(50, 100))])
    cell_count, _ = counter.count()
    assert cell_count > 0  # just check it runs without error


def test_add_clauses_list_wrong_type():
    c = Counter(seed=1)
    with pytest.raises(TypeError):
        c.add_clauses([["a", "b"]])


if __name__ == '__main__':
    ret = pytest.main([__file__, '-v'] + sys.argv[1:])
    raise SystemExit(ret)
