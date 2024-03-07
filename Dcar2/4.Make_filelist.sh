#!/bin/bash

# Make filelist for each chromosome. Note that the order is important! The order is actually alphabetical so we can simply sort by name.

for i in `seq 1 15`; do
	ls *Chr${i}.mpileup | sort > Chr${i}_filelist.txt
done

