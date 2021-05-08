#!/bin/bash

file=$1
mc=`grep "^c t " $file`
echo "c o found header: $mc"

if [ -z ${TMPDIR+x} ];
then
    TMPDIR="."
    echo "c o TMP dir not defined, setting it to $TMPDIR"
else
    echo "c o TMP dir defined: $TMPDIR"
fi

rm -f "$TMPDIR/ind.txt" "$TMPDIR/sol.txt" "$TMPDIR/clean_file.txt"


if [[ "$mc" = "c t mc" ]]; then
    echo "c o this is a regular model counting file"
    grep -v "^c" $file |  ./approxmc > "$TMPDIR/sol.txt"
    sed -E "s/^(.)/c o \1/" "$TMPDIR/sol.txt"

    sat=`grep "^s .*SATIS" "$TMPDIR/sol.txt"`
    count=`grep "^s mc" "$TMPDIR/sol.txt" | awk '{print $3}'`
    log_10_count=`echo "scale=15; l($count)/l(10)" | bc -l `

    echo $sat
    echo "c s type mc"
    echo "c s log10-estimate $log_10_count"
    echo "c s approx arb int $count"

elif [[ "$mc" = "c t pmc" ]]; then
    echo "c o this is a projected model counting file"
    grep "c p show" $file | sed -E "s/c p show (.*)/c ind \1 0/" > "$TMPDIR/ind.txt"
    grep -v "^c" $file > "$TMPDIR/clean_file.txt"
    cat "$TMPDIR/clean_file.txt" "$TMPDIR/ind.txt" | ./approxmc  > "$TMPDIR/sol.txt"
    sed -E "s/^(.)/c o \1/" "$TMPDIR/sol.txt"

    sat=`grep "^s .*SATIS" "$TMPDIR/sol.txt"`
    count=`grep "^s mc" "$TMPDIR/sol.txt" | awk '{print $3}'`
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

rm -f "$TMPDIR/ind.txt" "$TMPDIR/sol.txt" "$TMPDIR/clean_file.txt"
