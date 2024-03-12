#!/usr/bin/python
# coding=utf-8

# This script looks if the broken spams listed in the *significant_windows   files are broken because of lack 
# of coverage of in-between windows, or beacause differentiation has really a drop in-between peaks

# REQUIRES: 
#1: file with list of outliers with switch to state 'list_scaffolds_with_hmm_outliers.txt '
#2: file with all fst values after filtering for comverage of window. Tipically 'combined_pairs_filter_250bp_sorted.fst'
#3: step size used for fst sliding windows
#4: output file

import sys




filein1=open(sys.argv[1],'r')
fileout1=open(sys.argv[4],'w')

step=int(sys.argv[3])

# creates list of scaffolds
SCAFFS=[]
for line in filein1:
	print line
	line=line.rstrip("\n")
	SCAFFS.append(line)
print '---> scaffolds with outliers are:'
print SCAFFS

for scaff in SCAFFS:
	print '************ '+scaff+' ************'
	POS=[]
	filein2=open(sys.argv[2],'r')
	for line in filein2:
		if line.startswith(scaff):
			pos=line.split('\t')[1]
			POS.append(int(pos))
	print '---> covered position in fst file are:'
	print POS
	
	filein3=open(scaff+'_signif_window.txt','r')
	SPANS=[]
	for line in filein3:
		line=line.rstrip("\n")
		SPANS.append(line)			
	print '---> listed spans are:'
	print SPANS
	
	if len(SPANS) == 1:
		span=SPANS[0].split('\t')
		if span[0] != span[1]:
			fileout1.write(scaff+'\t')
			fileout1.write(span[0]+'\t'+span[1]+'\n')
		
	if len(SPANS) > 1:
		x=0
		y='null'
		for n in range(0,len(SPANS)-1):
			print n
			#print SPANS[n]
			#print SPANS[n+1]
			span1=SPANS[n].split('\t')
			span2=SPANS[n+1].split('\t')
			if span1[0] != span1[1]:
#				if span2[0] != span2[1]:
				if x==0:
					fileout1.write(scaff+'\t'+span1[0]+'\t')
					x=1
					y=n
				if span2[0] != span2[1]:
					print '---> consecutive spans >0 are:'
					print 'span1 is '+str(SPANS[n])
					print 'span2 is '+str(SPANS[n+1])
					for w in range(int(span1[1])+step,int(span2[0]),step):
						print 'checking if this w is in fst file: '+str(w)
						if w in POS:
							print 'FOUND '+str(w)
							fileout1.write(span1[1]+'\n')
							x=0
						else:
							print 'I did not find this w in fst file: '+str(w)
				
			
		fileout1.write(span2[1]+'\n')
			
						
		
			