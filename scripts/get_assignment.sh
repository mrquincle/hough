#!/bin/sh

ifile=${1:? "Usage: $0 <input file>"}

cat $ifile | sort -n | uniq -c | sort -nr | sort -u -n -k 2 | cut -f2,5
