#!/bin/bash

# Make filelist for each chromosome. Note that the order is important! The order is actually alphabetical so we can simply sort by name.

#for i in `seq 1 10`; do
for i in 2 3 4 5 6 7 8 9 10 11 12 14 15; do 
	ls *Chr${i}.mpileup | sort > Chr${i}_filelist.txt
done

