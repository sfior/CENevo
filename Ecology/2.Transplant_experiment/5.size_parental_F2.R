library(ggplot2)
library(reshape2)

# This script makes plots for size of parental plants growing in the transplant experiment, can compares size with F2 plants

## Size of parental plants in multi-year transplant
df_wild<-read.table('size_parentals.txt',sep='\t',header=T)
df_wild$altitude <- factor(df_wild$altitude, levels = c("low", "high"))  
data<-df_wild_FZ[,c("ros_s16","ros_e16","ros_s17","ros_s18","ros_e18","ros_s19","ros_e19","ros_s20","ros_e20","ros_s21","ros_e21","ros_s22","ros_e22","ros_s23","ros20.23","altitude","site_altitude")]
# Reshape the data into long format
melted_data <- melt(data, id.vars = c("altitude", "site_altitude"))
# Make plot:
colors <- c("cyan", "darkorange1", "blue", "red")
melted_data$fill_factor <- paste(melted_data$altitude, melted_data$site_altitude, sep="_")
melted_data$fill_factor <- factor(melted_data$fill_factor, levels = c("high_high_site", "low_high_site", "high_low_site", "low_low_site"))  # this is to have correct colors in plots automatically

p <- ggplot(melted_data, aes(x = variable, y = value, fill = fill_factor, color = fill_factor)) +
  geom_boxplot(color = "black") +  # Add the color argument here
  labs(x = "Life stage", y = "Size") +
  labs(fill = melted_data$fill_factor, color = melted_data$fill_factor) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "right", axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 11), axis.title = element_text(size = 12)) + 
  scale_x_discrete(labels = c('s16','e16','s17','s18','e18','s19','e19','s20','e20','s21','e21','s22','e22','s23','e20-s23'))
p
num_boxes <- length(unique(melted_data$variable))
for (i in seq(1, num_boxes, by = 1)) {
  p <- p + geom_vline(xintercept = i - 0.5, linetype = "dashed", color = "black")
}
p


# Size of F2 plants in multi-year transplant
df_f2<-read.table('size_F2.txt',header=T,sep='\t')
df_f2$CEN <- factor(df_f2$CEN, levels = c("H-L", "H-H"))  
# Back-transform the standardised mean values using 3rd year:
hs<-subset(df_f2,site=='Findeln')
ls<-subset(df_f2,site=='Zeneggen')

m_hs_yr3<-mean(hs$start21_size, na.rm=T)
sd_hs_yr3<-sd(hs$start21_size, na.rm=T)
bktr_hs_yr3<-(hs$size_mean_scaled*sd_hs_yr3) + m_hs_yr3

m_ls_yr3<-mean(ls$start21_size, na.rm=T)
sd_ls_yr3<-sd(ls$start21_size, na.rm=T)
bktr_ls_yr3<-(ls$size_mean_scaled*sd_ls_yr3) + m_ls_yr3


# make combined plot of size for wild and F2
hshg<-subset(df_wild,site=='findeln' & altitude =='high')
wild_ros_hs<-hshg$ros20.23/10 # adjust measuring unit
wild_F2_hs <- data.frame(ros = c(bktr_hs_yr3, wild_ros_hs), pop = c(rep('F2', length(bktr_hs_yr3)), rep('wild', length(wild_ros_hs))))
plot<-ggplot(wild_F2_hs, aes(x = pop, y = ros, fill=pop)) +
  geom_violin(trim=F) + stat_summary(fun=mean, geom="point", shape=16, size=2,color='black') + 
  xlab("") +
  ylab("Plant size") +
  ggtitle("High environment") +
  scale_fill_manual(values = c("grey", "blue")) +
  theme_minimal() +
  coord_flip() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
plot


lslg<-subset(df_wild,site=='zeneggen' & altitude =='low')
wild_ros_hs<-lslg$ros20.23/10
wild_F2_hs <- data.frame(ros = c(bktr_hs_yr3, wild_ros_hs), pop = c(rep('F2', length(bktr_hs_yr3)), rep('wild', length(wild_ros_hs))))

plot<-ggplot(wild_F2_ls, aes(x = pop, y = ros,fill=pop)) +
  geom_violin(trim=F) + stat_summary(fun=mean, geom="point", shape=16, size=2,color='black') + 
  xlab("") +
  ylab("Plant size") +
  ggtitle("Low environment") +
  scale_fill_manual(values = c("grey", "red")) +
  theme_minimal() + 
  coord_flip() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
plot

