#!/cluster/apps/gdc/python/3.6.1/bin/python3.6

# This script adds flanking regions to a bed file
# Example usage:python addFlankingRegions2BED.py input.bed -f 1000 -o output.bed

# Load necessary modules
#import os
import pandas
import argparse

#os.chdir("/Users/luqman/Desktop")
#os.chdir("./")

# Parse arguments
parser = argparse.ArgumentParser(description="Add flanking regions to bed file", add_help=True)
parser.add_argument("filename", action="store", help="The input bed file.")
parser.add_argument("-f", action="store", dest="length_flanking", default=0, type=int, help="Length of flanking regions.")
parser.add_argument("-o", action="store", dest="output", help="Output file name.")
args = parser.parse_args()
filename = args.filename
length_flanking = args.length_flanking
output = args.output

# We calculate the scaffold size, which will be used to define the upper and lower bound positions for the flanking regions
def scaffold_size(x):
	scaffold_name = str(x)
	first_split = scaffold_name.split("_size")
	second_split = first_split[1].split("_")
	if len(second_split) <= 1:
		scaffold_size = int(second_split[0])
	elif len(second_split) > 1:
		# We add 1 because BED start positions are zero-based and BED end positions are one-based, see: http://bedtools.readthedocs.io/en/latest/content/overview.html
		scaffold_size = int(second_split[2]) - int(second_split[1]) + 1
	return scaffold_size

# We define a function that puts a lower bound (0) on the start position, as positions starts from 0 (BED start position is 0-based, see: http://bedtools.readthedocs.io/en/latest/content/overview.html)
def min_0(x, length_flanking):
	if x > length_flanking:		
		output = x - length_flanking
	elif x <= length_flanking:
		output = 0
	return output

# We define a function that puts an upper bound (len(chr)) on the end position, as position must lie within chromosome
def max_pos(end_pos, chr_length, length_flanking):
	if end_pos < chr_length - length_flanking:
		output = end_pos + length_flanking
	elif end_pos >= chr_length - length_flanking:		
		output = chr_length
	return output

# Check that end position does not exceed scaffold size
def len_test(x):
	if int(x[2]) > int(x["chr_length"]):
		output = 1
	else:
		output = 0
	return output

# Define the main function, which adds flanking regions left and right of start and end positions respectively
def add_flanking_regions(filename, length_flanking, output):	
	bedfile = pandas.read_csv(filename, sep='\t', header=None)
	bedfile[len(bedfile.columns)] = bedfile[0].apply(lambda x: scaffold_size(x))
	bedfile.rename(columns={ bedfile.columns[len(bedfile.columns)-1]: "chr_length" }, inplace=True)
	bedfile[1] = bedfile[1].apply(lambda x: min_0(x,length_flanking))
	bedfile[2] = bedfile.apply(lambda x: max_pos(x[2], x["chr_length"],length_flanking), axis=1)
	bedfile["check"] = bedfile.apply(len_test,axis=1)
	if bedfile["check"].sum() > 0:
		print("Error: some end positions are greater than the length of the scaffold.")	
	bedfile.drop("chr_length", axis=1, inplace=True)
	bedfile.drop("check", axis=1, inplace=True)
	bedfile.to_csv(output, sep='\t', header=False, index=False)

# Run the function
if __name__ == "__main__":
	add_flanking_regions(filename, length_flanking, output)
	
