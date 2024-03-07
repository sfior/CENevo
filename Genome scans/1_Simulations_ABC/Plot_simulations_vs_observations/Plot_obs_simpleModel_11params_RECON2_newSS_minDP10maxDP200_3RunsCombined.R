# Import libraries
library(ggplot2)
library(gridExtra)

# Set working directory
setwd("/Users/luqman/Desktop")

# Import data
#ABC_rej<-read.delim("./concatenated_results_simpleModel_11params_RECON2_newSS_minDP10maxDP200_Dsyl_3RunsCombined_PLS20.txt")  # DSYL
ABC_rej<-read.delim("./concatenated_results_simpleModel_11params_RECON2_newSS_minDP10maxDP200_Dcar_3RunsCombined_PLS20.txt")  # DCAR
ABC_rej<-ABC_rej[c(1:3150000),c(26:ncol(ABC_rej))]
#ABC_obs<-read.delim("./GW_Neutral_5000_1000_anchored_50kb_start1_end7194149_windowsize5000_M2_Q0_simpleModel_11params_RECON2_newSS_minDP10maxDP200_Dsyl_3RunsCombined_PLS20.obs")  # DSYL
ABC_obs<-read.delim("./GW_Neutral_5000_0_anchored_50kb_start1_end6684312_windowsize5000_M2_Q0_simpleModel_11params_RECON2_newSS_minDP10maxDP200_Dcar_3RunsCombined_PLS20.obs")  # DCAR
ABC_obs<-ABC_obs[,c(2:ncol(ABC_obs))]

# Number of PLS compenents (to plot!)
num_PLS <- 20
#num_PLS <- 15

# Ensures an even number of PLS componenets.Since this plotter plots 2D plots, it requires an even number of PLS components to plot. So in case of an odd number, just plot n-1 PLS components.
if (num_PLS %% 2 !=0) {num_PLS <- num_PLS - 1}

# For storing ggplot objects in a list using a loop: https://stackoverflow.com/questions/31993704/storing-ggplot-objects-in-a-list-from-within-loop-in-r
plot_list <- list()
for (i in seq(1, (num_PLS/2))) {
  local({
    i <- i
    plot_list[[i]] <<- ggplot(data = ABC_rej, aes(x=ABC_rej[,2*i-1], y=ABC_rej[,2*i]) ) +
      geom_hex(bins = 35) + scale_fill_gradientn(colours=c("gray85","gray15"),name = "sim count",na.value=NA) +
      geom_hex(data = ABC_obs, bins = 70, aes(x=ABC_obs[,2*i-1], y=ABC_obs[,2*i], alpha=..count..), fill="red") +
      theme_bw() + xlab(paste0("PLS ",2*i-2)) + ylab(paste0("PLS ",2*i-1)) + labs("asd")
  })
}
grid.arrange(grobs = plot_list, ncol=3, top = "Overlap of simulated & observed PLS-transformed summary statistics")


##############################################################################
#MANUAL

# Plot sim vs obs (PLS-transformed)
# plot_PLS_01 <- ggplot(ABC_rej, aes(x=ABC_rej$LinearCombination_0, y=ABC_rej$LinearCombination_1) ) +
#   geom_hex(bins = 100) +
#   geom_point(data = ABC_obs, aes(x=ABC_obs$LinearCombination_0, y=ABC_obs$LinearCombination_1, col = 'red', alpha = 0.01, shape = "."), show.legend = FALSE) + 
#   theme_bw() + xlab("PLS 0") + ylab("PLS 1")
# 
# plot_PLS_23 <- ggplot(ABC_rej, aes(x=ABC_rej$LinearCombination_2, y=ABC_rej$LinearCombination_3) ) +
#   geom_hex(bins = 100) +
#   geom_point(data = ABC_obs, aes(x=ABC_obs$LinearCombination_2, y=ABC_obs$LinearCombination_3, col = 'red', alpha = 0.01, shape = "."), show.legend = FALSE) + 
#   theme_bw() + xlab("PLS 2") + ylab("PLS 3")
# 
# plot_PLS_45 <- ggplot(ABC_rej, aes(x=ABC_rej$LinearCombination_4, y=ABC_rej$LinearCombination_5) ) +
#   geom_hex(bins = 100) +
#   geom_point(data = ABC_obs, aes(x=ABC_obs$LinearCombination_4, y=ABC_obs$LinearCombination_5, col = 'red', alpha = 0.01, shape = "."), show.legend = FALSE) + 
#   theme_bw() + xlab("PLS 4") + ylab("PLS 5")
# 
# plot_PLS_67 <- ggplot(ABC_rej, aes(x=ABC_rej$LinearCombination_6, y=ABC_rej$LinearCombination_7) ) +
#   geom_hex(bins = 100) +
#   geom_point(data = ABC_obs, aes(x=ABC_obs$LinearCombination_6, y=ABC_obs$LinearCombination_7, col = 'red', alpha = 0.01, shape = "."), show.legend = FALSE) + 
#   theme_bw() + xlab("PLS 6") + ylab("PLS 7")
# 
# plot_PLS_89 <- ggplot(ABC_rej, aes(x=ABC_rej$LinearCombination_8, y=ABC_rej$LinearCombination_9) ) +
#   geom_hex(bins = 100) +
#   geom_point(data = ABC_obs, aes(x=ABC_obs$LinearCombination_8, y=ABC_obs$LinearCombination_9, col = 'red', alpha = 0.01, shape = "."), show.legend = FALSE) + 
#   theme_bw() + xlab("PLS 8") + ylab("PLS 9")
# 
# grid.arrange(plot_PLS_01, plot_PLS_23, plot_PLS_45, plot_PLS_67, plot_PLS_89, ncol=3)
# 
# ################################################################################
# 
# plot_PLS_01 <- ggplot(ABC_rej, aes(x=ABC_rej$LinearCombination_0, y=ABC_rej$LinearCombination_1) ) +
#   geom_hex(bins = 100) +
#   geom_point(data = ABC_obs, aes(x=ABC_obs$LinearCombination_0, y=ABC_obs$LinearCombination_1, col = 'red', alpha = 0.01, shape = "."), show.legend = FALSE) + 
#   theme_bw() + xlab("PLS 0") + ylab("PLS 1")
# 
# plot_PLS_23 <- ggplot(ABC_rej, aes(x=ABC_rej$LinearCombination_2, y=ABC_rej$LinearCombination_3) ) +
#   geom_hex(bins = 100) +
#   geom_point(data = ABC_obs, aes(x=ABC_obs$LinearCombination_2, y=ABC_obs$LinearCombination_3, col = 'red', alpha = 0.01, shape = "."), show.legend = FALSE) + 
#   theme_bw() + xlab("PLS 2") + ylab("PLS 3")
# 
# plot_PLS_45 <- ggplot(ABC_rej, aes(x=ABC_rej$LinearCombination_4, y=ABC_rej$LinearCombination_5) ) +
#   geom_hex(bins = 100) +
#   geom_point(data = ABC_obs, aes(x=ABC_obs$LinearCombination_4, y=ABC_obs$LinearCombination_5, col = 'red', alpha = 0.01, shape = "."), show.legend = FALSE) + 
#   theme_bw() + xlab("PLS 4") + ylab("PLS 5")
# 
# plot_PLS_67 <- ggplot(ABC_rej, aes(x=ABC_rej$LinearCombination_6, y=ABC_rej$LinearCombination_7) ) +
#   geom_hex(bins = 100) +
#   geom_point(data = ABC_obs, aes(x=ABC_obs$LinearCombination_6, y=ABC_obs$LinearCombination_7, col = 'red', alpha = 0.01, shape = "."), show.legend = FALSE) + 
#   theme_bw() + xlab("PLS 6") + ylab("PLS 7")
# 
# plot_PLS_89 <- ggplot(ABC_rej, aes(x=ABC_rej$LinearCombination_8, y=ABC_rej$LinearCombination_9) ) +
#   geom_hex(bins = 100) +
#   geom_point(data = ABC_obs, aes(x=ABC_obs$LinearCombination_8, y=ABC_obs$LinearCombination_9, col = 'red', alpha = 0.01, shape = "."), show.legend = FALSE) + 
#   theme_bw() + xlab("PLS 8") + ylab("PLS 9")
# 
# grid.arrange(plot_PLS_01, plot_PLS_23, plot_PLS_45, plot_PLS_67, plot_PLS_89, ncol=3)

