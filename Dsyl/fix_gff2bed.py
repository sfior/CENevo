#!/cluster/apps/gdc/python/3.6.1/bin/python3.6

# This script corrects inconsistency between the gff position indexing and bed indexing.
# Example usage:python fix_gff2bed.py input.bed -o output.bed

# Load necessary modules
#import os
import pandas
import argparse

#os.chdir("/Users/luqman/Desktop")
#os.chdir("./")

# Parse arguments
parser = argparse.ArgumentParser(description="Fix position indexing of gff2bed produced file", add_help=True)
parser.add_argument("filename", action="store", help="The input bed file.")
parser.add_argument("-o", action="store", dest="output", help="Output file name.")
args = parser.parse_args()
filename = args.filename
output = args.output

# This function corrects inconsistency between the gff position indexing and bed. The gff2bed script converts 1-based, closed [start, end] General Feature Format v3 (GFF3) to sorted, 0-based, half-open [start-1, end) extended BED-formatted data.
# chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0. See: https://genome.ucsc.edu/FAQ/FAQformat.html
# chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
# The problem here is that, for some reason, this gff file is of the format closed [start+1, end+1], i.e. chromStart = 2 and chromEnd = 101. Hence we need to subtract 1 from both chromStart and Chromend in the resulting bed file. 

def fix_gff2bed(filename, output):
	bedfile = pandas.read_csv(filename, sep='\t', header=None)
	bedfile[1] = bedfile[1].apply(lambda x: x - 1)
	bedfile[2] = bedfile[2].apply(lambda x: x - 1)
	bedfile.to_csv(output, sep='\t', header=False, index=False)

# Run the function
if __name__ == "__main__":
	fix_gff2bed(filename, output)
	
