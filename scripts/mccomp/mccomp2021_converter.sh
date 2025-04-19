#!/bin/bash

file=$1
mc=`grep "^c t " $file`
echo "c o found header: $mc"

solfile=$(mktemp)
indfile=$(mktemp)
cleanfile=$(mktemp)
echo "c o solfile: $solfile  indfile: $indfile  cleanfile: $cleanfile"


if [[ "$mc" = "c t mc" ]]; then
    echo "c o this is a regular model counting file"
    grep -v "^c" $file |  ./approxmc > $solfile
    sed -E "s/^(.)/c o \1/" $solfile

    sat=`grep "^s .*SATIS" $solfile`
    count=`grep "^s mc" $solfile | awk '{print $3}'`
    log_10_count=`echo "scale=15; l($count)/l(10)" | bc -l `

    echo $sat
    echo "c s type mc"
    echo "c s log10-estimate $log_10_count"
    echo "c s approx arb int $count"

elif [[ "$mc" = "c t pmc" ]]; then
    echo "c o this is a projected model counting file"
    grep "c p show" $file | sed -E "s/c p show (.*)/c ind \1 0/" > $indfile
    grep -v "^c" $file > $cleanfile
    cat $cleanfile $indfile | ./approxmc  > $solfile
    sed -E "s/^(.)/c o \1/" $solfile

    sat=`grep "^s .*SATIS" $solfile`
    count=`grep "^s mc" $solfile | awk '{print $3}'`
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
