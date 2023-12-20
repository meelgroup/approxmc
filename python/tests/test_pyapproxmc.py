import sys

import pytest

from pyapproxmc import Counter

def test_minimal():
    counter = Counter(seed=2157, epsilon=0.8, delta=0.2)
    counter.add_clause(list(range(1,100)))

    significand, exponent = counter.count()
    print(f'count: {significand} * 2**{exponent}')
    assert significand * 2**exponent == 512 * 2**90

def test_sampling_set():
    counter = Counter(seed=2157, epsilon=0.8, delta=0.2)
    counter.add_clause(range(1,100))
    assert counter.count(list(range(1,50))) == (64, 43)

def test_real_example():
    counter = Counter(seed=120, epsilon=0.8, delta=0.2)

    with open("test_1.cnf") as test_cnf:
        # Pop sampling set and metadata lines
        lines = test_cnf.readlines()[2:]

        # Add clauses to counter
        for line in lines:
            literals = [int(i) for i in line.split()[:-1]]
            counter.add_clause(literals)

    assert counter.count(list(range(1,21))) == (64,14)


if __name__ == '__main__':
    pytest.main([__file__, '-v'] + sys.argv[1:])
