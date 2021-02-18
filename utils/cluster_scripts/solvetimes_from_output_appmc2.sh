#!/bin/bash

xzgrep "at least" *.out.xz | sed "s/no_w.cnf.gz.*/no_w.cnf.gz/" > solved.csv
ls *.out.xz | sed "s/.cnf.*/.cnf/" > all_files.csv
rm solveTimes.csv
for a in `cat solved.csv`; do
    time=`xzgrep "User time" "$a.timeout.xz" | awk '{print $4}'`
    echo "$time $a" >> solveTimes.csv
done
