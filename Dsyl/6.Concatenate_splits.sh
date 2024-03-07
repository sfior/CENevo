#!/bin/bash
### Concatenate split regions into single (Chr 16 file)

head -n1 split_1/*1kbOverlapWindows.obs > Chr16_allScaffolds_FullSS_start1_end9999999_windowsize5000_M2_Q0.obs
for i in `seq 1 15`; do
	tail -n+2 split_${i}/*1kbOverlapWindows.obs >> Chr16_allScaffolds_FullSS_start1_end9999999_windowsize5000_M2_Q0.obs;
	cat split_${i}/*1kbOverlapWindows.scaffoldlist >> Chr16_allScaffolds_FullSS_start1_end9999999_windowsize5000_M2_Q0.scaffoldlist;
done
