#!/usr/bin/python

### This script converts amino acid sequences (e.g. from Geneious) to ACTG characters for use in haplotype network software ###
### It outputs a phylip file ###

## Example usage (in Python): aminoacid2ACTG("test.txt", "output")
## Example usage (in Linux command line): python aminoacid2ACTG.py input.txt output
## Coded in Python 3.x

# Import modules
import numpy
import os
import argparse

# Change to working directory
os.chdir(".")

#%%
############################################## MAIN FUNCTION ##############################################

# We define the main function. It requires 2 arguments: an input file (see below for details) and an output file name (phylip format).
def aminoacid2ACTG(input_filename, output_filename):
	
	# Open input file:
	# The input file is a file containing the amino acid sequences with one sequence per line (sequential format), with the string preceded by the unique sample name.
	with open(input_filename, 'r') as f:
		lines_raw = f.read().splitlines()
		line_raw_split = [line.split(": ") for line in lines_raw]
		sequences_temp = [line[1] for line in line_raw_split]       # To remove leading and ending whitespaces.
		sequences_raw = [line.strip() for line in sequences_temp]
		names_raw = [line[0] for line in line_raw_split]
	# For strict Phylip format. This requires sample names to be exactly 10 characters long. If they are shorter, we fill/append with "_" characters. If longer, we trunctate.
#		names_temp = [h[:10] for h in names_raw] 
#		names = [r + "".join(["_"]*(10 - len(r))) for r in names_temp] 
	# For relaxed Phylip format (allow for full length sample name, followed by a space)
		names = [h.replace(" translation","") + " " for h in names_raw]
	# Or if you just want to retain the first X (e.g. 13) characters in the sample name:	
#		names = [h.replace(" translation","")[:13] + " " for h in names_raw]
		
	# We then define the variables
	sequence_length = len(min(sequences_raw, key=len))
	no_samples = len(names)
	
	# We then convert the input list of rows (which correspond to the amino acid sequence of one sample) to a list of columns (corresponding to the exhibited amino acid per site across all samples)
	sites = []
	for s in range(sequence_length):
		site = [l[s] for l in sequences_raw]
		sites.append(site)
	
	# We now substitue codon letters with nucleotide base letters.
	# Here, we've coded to allow up to 6 different amino acids per site. 
	converted_sites = []
	for i in sites:
		if len(set(i)) == 1:
			converted_site = [list(numpy.random.choice(["A","C","T","G"],1,replace=True))[0]]*no_samples
		elif len(set(i)) > 1:   
			converted_site = ["A" if x==list(set(i))[0] else "C" if x==list(set(i))[1] else "T" if x==list(set(i))[min(len(set(i))-1, 2)] else "G" if x==list(set(i))[min(len(set(i))-1, 3)] else "N" if x==list(set(i))[min(len(set(i))-1, 4)] else "?" if x==list(set(i))[min(len(set(i))-1, 5)] else x for x in i]
		converted_sites.append(converted_site)	
	
	# We then need to convert the list of columns (sites) back to a list of rows (samples)
	sequences = []		
	for k in range(no_samples):
		sequence = [c[k] for c in converted_sites]
		sequences.append(sequence)
		final_sequences = [''.join(m) for m in sequences]
	
	# We then concatenate the sample names to their respective (converted) sequences, for the output	
	sample_output_lines = []
	if len(names) == len(final_sequences):
		for v in range(no_samples):
			sample_output_line = [str(names[v]), str(final_sequences[v])]
			sample_output_lines.append(sample_output_line)
			samples = ["".join(j) for j in sample_output_lines]
	
	# We define the header for the output file (Phylip requires this to be "x y" where x is the number of samples and y is the sequence length, separated by whitespace)
	file_header = [" " + str(no_samples) + " " + str(sequence_length)]	
	
	# Finally, we output the results to an output file:
	with open(output_filename, "w") as output:
		header_line = file_header[0] + "\n"
		output.writelines([header_line])		
		for line in samples:
			output.write("%s\n" % line)
	os.rename(output_filename, output_filename + ".phy")

#%%

# We then allow for parsing of arguments from the command line using argparse 
parser = argparse.ArgumentParser(description="Amino acid to ACTG character converter", add_help=True)

parser.add_argument("input_filename", action="store", help="The main input file, which contains the amino acid sequences with one sequence per line, preceded by the unique sample name.")
parser.add_argument("output_filename", action="store", help="The output file name (file will be in Phylip format).")

args = parser.parse_args()	

input_filename = args.input_filename
output_filename = args.output_filename

if __name__ == "__main__":
	aminoacid2ACTG(input_filename, output_filename)
	print("Conversion completed succesfully! Phylip output file has been produced.")
	
#%%