#!/bin/python3

import pyapproxmc
import sys

fname = sys.argv[1]
print("using file: ", fname)
solver = pyapproxmc.Counter(verbosity=1)

num_cls = 0
p_found = False
p = []
with open(fname, "r") as f:
    for line in f:
        line = line.strip()
        if len(line) == 0: continue
        if line[0] == "c":
            line = line.split(" ")
            if len(line) < 3: continue
            if line[1] != "ind": continue
            p_found = True
            for l in line[2:]:
                l = int(l)
                if l != 0:
                    p.append(l)

            continue
        if line[0] == "p": continue
        # print(line)
        clause = []
        for l in line.split(" "):
            # print(l)
            l = int(l)
            if l == 0: break
            clause.append(l)
        solver.add_clause(clause)
        num_cls+=1
print("num cls: ", num_cls)
if not p_found:
    print("ERRROR: no projection found!!")
    exit(-1)
print("projection: " , p)

ret = solver.count(projection=p)
print("ret:", ret)
