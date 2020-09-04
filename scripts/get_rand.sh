#!/bin/sh

opath=${1:? "Usage: $0 <output path>"}

rm -f $opath/purity.txt
rm -f $opath/rand.txt
rm -f $opath/adjusted_rand.txt
rm -f $opath/mirkin.txt
rm -f $opath/hubert.txt

for f in ../build/results/*; do
	< "$f" grep '^Purity' | cut -f2 -d':' | tr -d ' ' >> $opath/purity.txt
	< "$f" grep '^Rand Index' | cut -f2 -d':' | tr -d ' ' >> $opath/rand.txt
	< "$f" grep '^Adjusted Rand Index' | cut -f2 -d':' | tr -d ' ' >> $opath/adjusted_rand.txt
	< "$f" grep '^Mirkin Index' | cut -f2 -d':' | tr -d ' ' >> $opath/mirkin.txt
	< "$f" grep '^Hubert Index' | cut -f2 -d':' | tr -d ' ' >> $opath/hubert.txt
done

wc -l $opath/purity.txt
wc -l $opath/rand.txt
wc -l $opath/adjusted_rand.txt
wc -l $opath/mirkin.txt
wc -l $opath/hubert.txt
