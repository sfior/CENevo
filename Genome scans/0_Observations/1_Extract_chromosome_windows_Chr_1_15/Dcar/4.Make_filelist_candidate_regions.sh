#!/bin/bash

# Make filelist for each chromosome. Note that the order is important! The order is actually alphabetical so we can simply sort by name.


ls *FT.mpileup | sort > FT_filelist.txt
ls *TCF1.mpileup | sort > TCF1_filelist.txt

