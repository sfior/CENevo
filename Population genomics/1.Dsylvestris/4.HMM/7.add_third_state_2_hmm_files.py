#!/usr/bin/python
# coding=utf-8

# This script adds a '3rd state' to windows that are significant across elevational pairs

import sys

chrs=sys.argv[1].split(",")
pairs=sys.argv[2].split(",")

d={}
filein1=open('candidate.hmm.continous.sweeps.per.chrm.txt','r')
for line in filein1:
	if line.startswith('scaffold'):
		line=line.rstrip('\n')
		l=line.split('\t')
		s=l[0]
		t=l[1:3]
		d[s]=t
print d

for c in chrs:
	print c
	for p in pairs:
		print p
		filein2=open('hmm.chr'+c+'/'+p+'_chr_'+c+'.hmm.txt','r')
		fileout2=open('hmm.chr'+c+'/'+p+'_chr_'+c+'_3rd_state.hmm.txt','w')
		for line in filein2:
			line=line.rstrip('\n')
			l=line.split('\t')
			if l[0] in d.keys():
				#print l
				scaff=l[0]
				#print d[scaff]
				#print scaff
				#print l[1]
				#print d[scaff][0]
				#print d[scaff][1]
				if int(l[1]) >= int(d[scaff][0]):
					if int(l[1]) <= int(d[scaff][1]):
						if l[5] == '1':
							#print 'pippo'
					 		fileout2.write('\t'.join(l[0:5])+'\t'+str(3)+'\n')
						else:
							fileout2.write(line+'\n') 	
					else:
						fileout2.write(line+'\n') 	
				else:
					fileout2.write(line+'\n') 
			else:
				fileout2.write(line+'\n')		 	
		fileout2.close()					 	
					 	
					 	
					 	