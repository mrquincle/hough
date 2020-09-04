#!/bin/sh


mkdir -p ../build/results

for f in ~/data/lines/*.data.txt; do

	echo $f

	../build/DetectLineModule $f

	./get_assignment.sh output.txt > ../build/assignment.txt

	fbase=${f##*/}
	fresult=${fbase%.*.*}
	frfile=../build/results/$fresult.result.txt
	echo "Print to $frfile"
	~/workspace/clustering/build/bin/clustering ../build/assignment.txt > $frfile

done
