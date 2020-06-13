xzgrep -i "terminated" *.timeout.xz |  sed -e "s/.gz[^ ]*//" | awk '{print $1 " " $5}' > signals.csv
ls -- *.out.xz > allFiles_xz.csv
ls -- *.out.xz | sed "s/.gz.*/.gz/" > allFiles.csv
xzgrep "FINISHED" *.out.xz | sed -e "s/.gz[^ ]*//" | awk '{print $5 " " $1}' > solveTimes.csv
awk '{print $2}' solveTimes.csv > solved.csv
grep -v -f solved.csv allFiles.csv | sed "s/.gz.*/.gz/" > unsolved.csv
cat unsolved.csv | awk '{print "5000.00 " $1}' >> solveTimes.csv
awk '{print $2 " " $1}' solveTimes.csv | sort > solveTimes_rev.csv
awk '{if ($1=="5000.00") {x+=10000} else {x += $1};} END {printf "%d\n", x}' solveTimes.csv > PAR2score
