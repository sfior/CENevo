
# This scripts takes the combined output of the iterations and adjsuts signs to keep trios consistent

# Requires:
# arg1: Trio file
# arg2: output of collection of the iterations (from script 5_combine_ites_perTrio.sh)
# arg3: name of output

# Example:
#python3 6_adjust_ABBA-BABA_sign.py trios_BBAA.txt   collected_trios_BBAA.txt   collected_trios_adjustedSign_BBAA.txt


import sys
filein1=open(sys.argv[1],'r')
filein2=open(sys.argv[2],'r')
fileout=open(sys.argv[3],'w')

next(filein1) # skip header
trios=[]
for line in filein1:
	l=line.split('\t')
	trios.append(l[0:3]) 
print("This is the list of trios:")
print(trios)
print(len(trios))

line=filein2.readline()
l=line.split('\t')
fileout.write('\t'.join(l[0:3]+l[7::])) 


print (" >>>>>>>>>>   Starting reading the iterations  <<<<<<<<<<<<<" )

n=0
print(" >>>>>>>>>>  Taking trio  <<<<<<<<<<<<<" )
trio=trios[n]
print(trio)
print('\n')

d={}
d[trio[0]]='P1'
d[trio[1]]='P2'
d[trio[2]]='P3'
print('Dictionary is:')
print(d)
print('\n')

for line in filein2:
	line = line.strip()
	if not line:
		if n < len(trios)-1:
			n=n+1
			trio=trios[n]
			print(" >>>>>>>>>>  Taking new trio  <<<<<<<<<<<<<" )
			print(trio)
			print('\n')
			d[trio[0]]='P1'
			d[trio[1]]='P2'
			d[trio[2]]='P3'
			print('New dictionary is:')
			print(d)
			print('\n')
			fileout.write('\n')
			continue  # skip empty or whitespace-only lines
		else:
			print("Reached end of trios")

	l = line.split('\t')

    # Skip lines that are too short
	if len(l) < 3:
		continue
	print('Reading line:')
	print(l)

	observedTrio=[d[l[0]],d[l[1]],d[l[2]]]
	print('observedTrio is:')
	print(observedTrio)
	print('\n')
	
	BBAA=l[7]
	ABBA=l[8]
	BABA=l[9]
	
	fileout.write('\t'.join(trio) + '\t')
	
	if observedTrio == ['P1', 'P2', 'P3']:
		fileout.write('\t'.join([BBAA,ABBA,BABA]) + '\n')
	elif observedTrio == ['P2', 'P3', 'P1']:
		fileout.write('\t'.join([BABA,BBAA,ABBA]) + '\n')
	elif observedTrio == ['P2', 'P1', 'P3']:
		fileout.write('\t'.join([BBAA,BABA,ABBA]) + '\n')
	elif observedTrio == ['P3', 'P1', 'P2']:
		fileout.write('\t'.join([ABBA,BABA,BBAA]) + '\n')
	elif observedTrio == ['P1', 'P3', 'P2']:
		fileout.write('\t'.join([BABA,ABBA,BBAA]) + '\n')
	elif observedTrio == ['P3', 'P2', 'P1']:
		fileout.write('\t'.join([ABBA,BBAA,BABA]) + '\n')

