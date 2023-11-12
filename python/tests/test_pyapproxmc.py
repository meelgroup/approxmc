from pyapproxmc import Counter

def minimal_test():
    counter = Counter(seed=2157, epsilon=0.8, delta=0.2)
    counter.add_clause(list(range(1,100)))
    assert counter.count() == (512, 90)

def sampling_set_test():
    counter = Counter(seed=2157, epsilon=0.8, delta=0.2)
    counter.add_clause(range(1,100))
    assert counter.count(list(range(1,50))) == (64, 43)

def real_example_test():
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
    minimal_test()
    sampling_set_test()
    real_example_test()
