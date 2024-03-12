#!/usr/bin/python
# coding=utf-8

# This script extracts consecutive outlier windows to define candidate regions


# Requires:
# 1: filein with list of scaffolds containing oultier windows
# 2: filein with stats for outliers
# 3: window size

# Example:
#python ../5_Extract_overlapping_windows_v2.py list_scaffolds_with_hmm_outliers.txt Overlapping_Outliers_ite_500_pop1.pop3_pop2.pop6_pop4.pop5_pop1.pop4_pop1.pop6_pop2.pop3_pop2.pop4_pop3.pop5_pop5.pop6.txt 500


import sys

# This function finds continuous series within a list
def group(L):
    first = last = L[0]
    for n in L[1:]:
        if n - int(sys.argv[3]) == last: # Part of the group, bump the end  
            last = n
        else: # Not part of the group, yield current group and start a new
            yield first, last
            first = last = n
    yield first, last # Yield the last group


filein1=open(sys.argv[1],'r')
filein2=open(sys.argv[2],'r')

for scaff in filein1: 
	scaff=scaff.rstrip("\n")
	print scaff

	filein2=open(sys.argv[2],'r')
	x=[]
	for line in filein2:
		line=line.rstrip("\n")
		l=line.split("\t")
		if l[0] == scaff:
			x.append(int(l[1]))
	print x		
	ovlap_list=list(group(x))


# ---> Creates a file containing the significant window spans
	fileout=open(scaff+"_signif_window.txt",'w')
	for l in ovlap_list:
		for i in l:
			fileout.write(str(i)+"\t")
		fileout.write("\n")

	filein2.close()
	fileout.close()

# ---> Creates a list of positions between the widow spans foud above (to use in subsequent script)
	filein3=open(scaff+"_signif_window.txt",'r')
	r=[]
	for line in filein3:
		line=line.rstrip("\n")
		l=line.split("\t")
		r.extend(range(int(l[0]),int(l[1])))

	fileout=open(scaff+"_overlapping_signif_windows.txt",'w')
	for i in r:
		fileout.write(str(i)+"\t")
	fileout.write("\n")

	filein3.close()
	fileout.close()

