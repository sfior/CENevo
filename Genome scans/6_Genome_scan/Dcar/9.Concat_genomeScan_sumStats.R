##### This scripts combines the genome scan results with the observed summary statistics calculated per window #####
# Recall that since some observed windows fail during ABC estimation, there may be less estimation windows than observed windows, hence we'll have to match and merge based on scaffold and window names.

### Load libraries
library(stringr)

### Define input variables
# Species (Dsyl or Dcar)
species <- "Dcar"
# Number of chromosomes
nChr <- 16
# Define input files and directories
if (species == "Dsyl") {
  Chr_scaffold_list <- c("Chr1_allScaffolds_FullSS_start1_end7187273_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved", "Chr2_allScaffolds_FullSS_start1_end5149401_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved", "Chr3_allScaffolds_FullSS_start1_end8134089_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved", "Chr4_allScaffolds_FullSS_start1_end5958733_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved", "Chr5_allScaffolds_FullSS_start1_end6199402_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved", "Chr6_allScaffolds_FullSS_start1_end6750011_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved", "Chr7_allScaffolds_FullSS_start1_end5209322_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved", "Chr8_allScaffolds_FullSS_start1_end5355932_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved", "Chr9_allScaffolds_FullSS_start1_end5836137_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved", "Chr10_allScaffolds_FullSS_start1_end5817725_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved", "Chr11_allScaffolds_FullSS_start1_end4529089_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved", "Chr12_allScaffolds_FullSS_start1_end4630051_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved", "Chr13_allScaffolds_FullSS_start1_end7920842_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved", "Chr14_allScaffolds_FullSS_start1_end8187860_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved", "Chr15_allScaffolds_FullSS_start1_end5372557_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved", "Chr16_allScaffolds_FullSS_start1_end9999999_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved")  #RECON2 newSS 1kb overlapping windows
  obs_dir <- "/cluster/project/gdc/people/lhirzi/GenomeScan/Dsyl/"
  genomescan_dir <- "/cluster/project/gdc/people/lhirzi/GenomeScan/Dsyl/GenomeScan_resultsTables/"
} else if (species == "Dcar") {
  Chr_scaffold_list <- c("Chr1_allScaffolds_FullSS_start1_end5522201_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved", "Chr2_allScaffolds_FullSS_start1_end5623217_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved", "Chr3_allScaffolds_FullSS_start1_end5547344_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved", "Chr4_allScaffolds_FullSS_start1_end5248954_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved", "Chr5_allScaffolds_FullSS_start1_end4732141_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved", "Chr6_allScaffolds_FullSS_start1_end4217913_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved", "Chr7_allScaffolds_FullSS_start1_end4006302_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved", "Chr8_allScaffolds_FullSS_start1_end7220078_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved", "Chr9_allScaffolds_FullSS_start1_end7995774_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved", "Chr10_allScaffolds_FullSS_start1_end5278334_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved", "Chr11_allScaffolds_FullSS_start1_end6967360_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved", "Chr12_allScaffolds_FullSS_start1_end5814850_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved", "Chr13_allScaffolds_FullSS_start1_end5930082_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved", "Chr14_allScaffolds_FullSS_start1_end5148409_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved", "Chr15_allScaffolds_FullSS_start1_end4927259_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved", "Chr16_allScaffolds_FullSS_start1_end9999999_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved")  #RECON2 newSS 1kb overlapping windows
  obs_dir <- "/cluster/project/gdc/people/lhirzi/GenomeScan/Dcar/"
  genomescan_dir <- "/cluster/project/gdc/people/lhirzi/GenomeScan/Dcar/GenomeScan_resultsTables/"
}

### Run for over n chromosomes
for (Chr in seq(1,nChr)) {
  
  ### First we import obs scaffold list. "obs" here refers to observed windows.
  window_pos<-read.delim(paste0(obs_dir, "Chr_", Chr, "/", Chr_scaffold_list[Chr],".scaffoldlist"), sep = "S", header = FALSE)
  # Select relevant data
  window_pos<-as.vector(window_pos[,4])
  # Remove leading and trailing (useless) characters
  scaffold_list.temp <- lapply(window_pos, function(x) substring(x, 2, nchar(x)-4))
  scaffold_list <- as.data.frame(t(as.data.frame(scaffold_list.temp)))
  rownames(scaffold_list)<-c(1:length(scaffold_list[,1]))
  # Split into scaffold name and window range
  scaffold_df <- as.data.frame(str_split_fixed(scaffold_list$V1, "_window_", 2))
  scaffold_df.windows <- as.data.frame(str_split_fixed(scaffold_df$V2, "_", 2)) 
  # Append columns to dataframe
  scaffold_df_concat <- as.data.frame(scaffold_df$V1)
  scaffold_df_concat[,2] <- scaffold_df.windows$V1
  scaffold_df_concat[,3] <- scaffold_df.windows$V2
  # Rename columns
  colnames(scaffold_df_concat) <- c("scaffold", "window_start", "window_end")
  
  ### To this obs scaffold list, we append the obs summary statistics. 
  # Import sumstats file
  raw_SS<-read.table(paste0(obs_dir, "Chr_", Chr, "/", Chr_scaffold_list[Chr],".obs"), header = TRUE)
  # Append to scaffold_df_concat
  scaffold_df_concat <- cbind(scaffold_df_concat,raw_SS)
  
  ### Then we merge scaffold_df_concat with the genome scan results.
  # Import genome scan results file
  genomeScan_raw<-read.table(paste0(genomescan_dir, "genomeScan_table.Chr", Chr, ".txt"), header = TRUE)
  # Assign column for which to match the two dfs. Here we'll match and merge on a column of scaffold_name+window_start+window_end
  genomeScan_raw[,14] <- paste(genomeScan_raw[,1], genomeScan_raw[,4], genomeScan_raw[,5], sep = "_")
  colnames(genomeScan_raw)[14] <- "match_column"
  scaffold_df_concat[,71] <- paste(scaffold_df_concat[,1], scaffold_df_concat[,2], scaffold_df_concat[,3], sep = "_")
  colnames(scaffold_df_concat)[71] <- "match_column"
  # Merge the summary statistics df with the genome scan df by the match_column.
  genomeScan_final <- merge(genomeScan_raw, scaffold_df_concat, by.x = "match_column", by.y = "match_column")
  
  ### Clean up merged df
  # Drop columns
  genomeScan_final <- genomeScan_final[,-c(1,15,16,17)]
  # Rename columns
  colnames(genomeScan_final)[1] <- "scaffold"
  colnames(genomeScan_final)[4] <- "window_start"
  colnames(genomeScan_final)[5] <- "window_end"
  colnames(genomeScan_final)[6] <- "-log(1-p_migHL)"
  colnames(genomeScan_final)[7] <- "-log(1-p_migLH)"
  colnames(genomeScan_final)[8] <- "-log(1-p_joint_migHL_migLH)"
  
  ### Write output
  #sapply(genomeScan_final, class)
  genomeScan_final$window_start <- as.numeric(as.character(genomeScan_final$window_start))
  genomeScan_final$window_end <- as.numeric(as.character(genomeScan_final$window_end))
  
  # Sort by chromosome position (first) and windows (second). Should already be sorted.
  #genomeScan_table <- genomeScan_table[with(genomeScan_table, order(Chr_position, window_start)), ]
  
  # Write out table
  write.table(genomeScan_final, file = paste0(genomescan_dir, "genomeScan_wSumStats.Chr",Chr,".txt"), row.names = FALSE)
}
