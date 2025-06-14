#!/usr/bin/python
# coding=utf-8
from __future__ import division


# This script makes a table with frequencies of High / Low protein classes and population coordinates
# It is used to plot pie charts on a map


# Requires:
# 1:A table with population coordinates
# 2:A fasta file with one class of protein (e.g. High)
# 3:A second fasta file with alternate class of proteins (e.g. Low) 
# For semplicity, we use haplotype sequences separate in high and low classes

#How to run:
#python 1_prepare_freqs.py \
#coordinates.txt \
#Dsyl_CEN_haplos_highClass.fasta \
#Dsyl_CEN_haplos_lowClass.fasta  \
#Protein_freqs.txt


import sys
import re


# make dictionary for popualtion names
pop_dict={}
filein1=open(sys.argv[1],'r')
for line in filein1:
	line=line.rstrip('\n')
	line=line.split('\t')
	name=line[0]
	pop_dict[name]=line[1::]
filein1.close()
#print(pop_dict)	
#print(pop_dict.keys())


# make dictionary for High class by assigning N of haplotypes in each pop
dHigh={}
filein1=open(sys.argv[2],'r')
for line in filein1:
	if line.startswith('>'):
		#print line
		line=line.lstrip('>')
		pop=re.split('_[0-9]+.+Seq',line)[0]
		#print pop
		if pop in dHigh.keys():
			dHigh[pop]=dHigh[pop]+1
		else:
			dHigh[pop]=1
#print(dHigh)

# make dictionary for Low class by assigning N of haplotypes in each pop		
dLow={}
filein1=open(sys.argv[3],'r')
for line in filein1:
	if line.startswith('>'):
		#print line
		line=line.lstrip('>')
		pop=re.split('_[0-9]+.+Seq',line)[0]
		#print pop
		if pop in dLow.keys():
			dLow[pop]=dLow[pop]+1
		else:
			dLow[pop]=1
#print(dLow)

# make table with haplotype prequencies
fileout=open(sys.argv[4],'w')
fileout.write('Population\tLatitude\tLongitude\tAltitude\tNr_Total\tNr_High\tNr_Low\tFreqHigh\tFreqLow\n')
for pop in pop_dict.keys():
#	print(pop)
	if pop in dHigh.keys() | dLow.keys():
		if not pop in dHigh.keys():
			NrH=0
		else:	
			NrH=dHigh[pop]	
	
		if not pop in dLow.keys():
			NrL=0
		else:	
			NrL=dLow[pop]	
	
		Tot=NrH+NrL
		FreqH=float(NrH/Tot)
		FreqL=float(NrL/Tot)
		fileout.write(pop+'\t'+'\t'.join(pop_dict[pop][1::])+'\t'+str(Tot)+'\t'+str(NrH)+'\t'+str(NrL)+'\t'+str(FreqH)+'\t'+str(FreqL)+'\n')
		
		
			
		
