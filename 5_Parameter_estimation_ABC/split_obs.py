#!/cluster/apps/gdc/python/3.6.1/bin/python3.6

### This script splits an .obs file (for use in ABCtoolbox) into multiple smaller ones, to facilitate parallelising across multiple CPUs ###
# Example usage: python split_obs.py filename num_of_splits split_num[optional]

# Import modules and set working directory
import sys
#import os
#os.chdir("/Users/luqman/Desktop")


# Define input variables
obs_filename = sys.argv[1]
nsplit = int(sys.argv[2])

with open(obs_filename) as f:
	obs_file_raw = f.read().splitlines()

# Remove header
obs_file_raw_noHeader = obs_file_raw[1:]
# Header
obs_file_header = obs_file_raw[:1]
# Define split size	
split_size = round ( len(obs_file_raw_noHeader) / nsplit )

# Split the file into nsplit files
# Initiate empty dict
obs_file_splits = {}	
split_start = 0
# First add the first nsplit-1 splits to dict (we don't add all splits here because the file may not be exactly divisable by an integer nsplits)
for split in range(nsplit-1):
	#print(split+1)
	obs_split_temp = obs_file_raw_noHeader[split_start:(split_size+split_start)]
	obs_split = obs_file_header + obs_split_temp
	split_start += split_size
	obs_file_splits["Split" + str(split+1)] = obs_split
# Then we add the final (remainder) split
final_split_size = len(obs_file_raw_noHeader) - ((nsplit - 1) * split_size)
obs_split_final_temp = obs_file_raw_noHeader[(len(obs_file_raw_noHeader) - final_split_size):]
obs_split_final = obs_file_header + obs_split_final_temp
obs_file_splits["Split" + str(nsplit)] = obs_split_final

# Define output file name
output_prefix = obs_filename.split(".")

# Should we want to output all splits
if len(sys.argv) == 3:

	# Then we write each splits into separate files
	for split in range(nsplit):
		#print(split+1)
		output_filename = output_prefix[0] + "_" + str(split+1) + ".obs"
		with open(output_filename, "w") as splitobsfile:	
			for item in obs_file_splits["Split" + str(split+1)]:
				splitobsfile.write("%s\n" % item)


# Alternatively, should we want to output a specific split				
if len(sys.argv) == 4:		

	# Define input variables
	split_num = int(sys.argv[3])	
	
	# Then write out
	output_filename = output_prefix[0] + "_" + str(split_num) + ".obs"
	with open(output_filename, "w") as splitobsfile:	
		for item in obs_file_splits["Split" + str(split_num)]:
			splitobsfile.write("%s\n" % item)		

print("Completed succesfully!")