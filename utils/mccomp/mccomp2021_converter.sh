#!/bin/bash

file=$1
mc=`grep "^c t " $file`
echo "c o found header: $mc"

if [[ "$mc" = "c t mc" ]]; then
    echo "c o this is a regular model counting file"
    grep -v "^c" $file |  ./approxmc > sol.txt
    sed -E "s/^(.)/c o \1/" sol.txt

    sat=`grep "^s .*SATIS" sol.txt`
    count=`grep "^s mc" sol.txt | awk '{print $3}'`
    log_10_count=`echo "scale=15; l($count)/l(10)" | bc -l `

    echo $sat
    echo "c s type mc"
    echo "c s log10-estimate $log_10_count"
    echo "c s approx arb int $count"

elif [[ "$mc" = "c t pmc" ]]; then
    echo "c o this is a projected model counting file"
    grep "c p show" $file | sed -E "s/c p show (.*)/c ind \1 0/" > ind.txt
    grep -v "^c" $file >> clean_file.txt
    cat clean_file.txt ind.txt | ./approxmc  > sol.txt
    sed -E "s/^(.)/c o \1/" sol.txt

    sat=`grep "^s .*SATIS" sol.txt`
    count=`grep "^s mc" sol.txt | awk '{print $3}'`
    log_10_count=`echo "scale=15; l($count)/l(10)" | bc -l `

    echo $sat
    echo "c s type pmc"
    echo "c s log10-estimate $log_10_count"
    echo "c s approx arb int $count"

elif [[ "$mc" = "c t wmc" ]]; then
    echo "c o this is a weighted model counting file"
    echo "ERROR: we cannot deal with this, exiting"
    exit -1
else
    echo "ERROR: Header not found, this is not an MCComp 2021 file!"
    exit -1
fi
