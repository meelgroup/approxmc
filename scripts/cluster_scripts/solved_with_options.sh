#!/bin/bash

check="$1"
FILE="20180321_110706599_p_cnf_320_1120.cnf"
FILE="001-80-12-sc2014.cnf"
FILE="01A-1.cnf.gz.no_w.cnf"
FILE="01A-1.cnf.gz.no_w.cnf"

rm out
for d in `ls | grep $check`; do
    echo $d;
    PAR2=`cat $d/PAR2score`
    name=`xzgrep "Command being" $d/${FILE}.gz.timeout.xz`
    numsolved=`wc -l $d/solved.csv | awk '{print $1}'`
    numALL=`wc -l $d/allFiles.csv  | awk '{print $1}'`
    #c Arjun SHA revision 8885476de181e6efe9ab3038ded86ed0bcd9cfec
    revArj=`xzgrep -i "Arjun.*revision" $d/${FILE}.gz.out.xz | awk '{print $5}' | cut -c1-7`
    revCMS=`xzgrep -i "CMS.*revision" $d/${FILE}.gz.out.xz | awk '{print $5}' | cut -c1-7`
    revApp=`xzgrep -i "AppMC.*revision" $d/${FILE}.gz.out.xz | awk '{print $5}' | cut -c1-7`
    echo "$PAR2 $d $numsolved $numALL $revApp $revCMS $revArj $name" >> out
done
sed "s/${FILE}.*//" out | sed "s/Command.*time.*approxmc//"  | sed "s/.solved.csv//" |                 sed "s/ *Command being timed.*approxmc//" | sed "s/\t/ /g" | sed "s/\t/ /g" | sort -n
