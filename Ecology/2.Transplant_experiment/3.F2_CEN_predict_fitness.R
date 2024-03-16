#############Set WD and load packages and data
library(dplyr)
library(aster)
#setwd("/path/to/wd")
df <- read.csv("F2_CEN_fitness_pheno_data_for_aster.csv")
#############Standardize traits, define function
standardize_subset <- function(df, trait) {
  df_ros19 <- df[df[, trait] > 0, ];df_no_ros19 <- df[df[, trait] == 0, ]
  standardize_subset <- function(df_subset) {
    df_subset[, trait] <- scale(df_subset[, trait], center = TRUE, scale = TRUE)
    return(df_subset)
  }
  df_ros19_fin_hh <- df_ros19[df_ros19$site == "Findeln" & df_ros19$CEN_1 == "H-H", ]
  df_ros19_fin_hl <- df_ros19[df_ros19$site == "Findeln" & df_ros19$CEN_1 == "H-L", ]
  df_ros19_zen_hh <- df_ros19[df_ros19$site == "Zeneggen" & df_ros19$CEN_1 == "H-H", ]
  df_ros19_zen_hl <- df_ros19[df_ros19$site == "Zeneggen" & df_ros19$CEN_1 == "H-L", ]
  df_ros19_fin_hh <- standardize_subset(df_ros19_fin_hh);df_ros19_fin_hl <- standardize_subset(df_ros19_fin_hl)
  df_ros19_zen_hh <- standardize_subset(df_ros19_zen_hh);df_ros19_zen_hl <- standardize_subset(df_ros19_zen_hl)
  df_ros19 <- bind_rows(df_ros19_fin_hh, df_ros19_fin_hl, df_ros19_zen_hh, df_ros19_zen_hl)
  df <- bind_rows(df_ros19, df_no_ros19)
  df[, trait][df[, trait] == 0] <- NA;df[, trait] <- as.numeric(df[, trait])
  return(df)
}
#############Run function over focal traits 
traits_list <- c("start19_size", "start20_size", "start21_size", "day19", "day20", "day21")
#Loop through each trait and standardize it
for (trait in traits_list) {
  df <- standardize_subset(df, trait)
}

#############Calculate mean trait values across growing seasons 
cols_to_mean <- c("barcode", "start19_size", "start20_size", "start21_size");df_subset <- df[, cols_to_mean]
df_subset$mean_3cols <- rowMeans(df_subset[, 2:4], na.rm = TRUE)
df_subset$mean_2cols <- rowMeans(df_subset[, c("start19_size", "start20_size")], na.rm = TRUE)
df_subset$mean_2cols[is.na(df_subset$start20_size)] <- rowMeans(df_subset[is.na(df_subset$start20_size), c("start19_size", "start21_size")], na.rm = TRUE)

mean_of_height_df<-df_subset[c("barcode", "mean_3cols")]
colnames(mean_of_height_df)[2]  <- "mean_size"
df <- merge(df, mean_of_height_df,  by="barcode")

cols_to_mean <- c("barcode", "day19", "day20", "day21"); df_subset <- df[, cols_to_mean]
df_subset$mean_3cols <- rowMeans(df_subset[, 2:4], na.rm = TRUE)
df_subset$mean_2cols <- rowMeans(df_subset[, c("day19", "day20")], na.rm = TRUE)
df_subset$mean_2cols[is.na(df_subset$day20)] <- rowMeans(df_subset[is.na(df_subset$day20), c("day19", "day21")], na.rm = TRUE)

mean_of_height_df<-df_subset[c("barcode", "mean_3cols")]
colnames(mean_of_height_df)[2]  <- "mean_day"
df <- merge(df, mean_of_height_df,  by="barcode")

df<-df[!is.na(df$mean_size),];df<-df[!is.na(df$mean_day),];df<-df[!is.na(df$seed_count21),]
df<-df[!is.na(df$seed_count20),];df<-df[!is.na(df$flowered21),]


#############Reshape data and fit aster model
vars <- c("surv_start19", "surv_start20","surv_start21" ,"flowered19", "flowered20","flowered21","seed_count19", "seed_count20", "seed_count21")
redf1 <- reshape(df, varying = list(vars),
                 direction = "long", timevar = "varb", times = as.factor(vars),v.names = "resp")
redf1 <- data.frame(redf1, root = 1)
pred <- c(0, 1, 2, 1, 2, 3, 4, 5, 6)
fam <- c(1, 1, 1, 1, 1, 1, 2, 2, 2)
sapply(fam.default(), as.character)[fam]

all(pred < seq(along = pred)) 
redf1$cluster <- as.factor(redf1$cluster)
redf1$site <- as.factor(redf1$site)

m1 <- aster(resp ~ varb+CEN_1*mean_day*mean_size*site+cluster,pred, fam, varb, id, root, data = redf1)
m2 <- aster(resp ~ varb+CEN_1*mean_size*mean_size+1,pred, fam, varb, id, root, data = redf1)
anova(m2,m1)
summary(m1, info.tol=1e-20)

#############Predict individual fitness
pout6 <- predict(m1)
pout6 <- matrix(pout6, nrow = nrow(m1$x), ncol = ncol(m1$x))
colnames(pout6) <- colnames(m1$x)
aster_output<-cbind(pout6,df)
colnames(aster_output)[26]<-"count1"
colnames(aster_output)[27]<-"count2"
colnames(aster_output)[28]<-"count3"
aster_output$fitness <- aster_output[, grep("seed_", colnames(aster_output))]
aster_output$comb_fit <- apply(aster_output$fitness, 1, "sum")
hist(aster_output$comb_fit)

#############Relativize fitness within transplant site and genotype
findeln <- subset(aster_output, site=="Findeln");zeneggen <- subset(aster_output, site=="Zeneggen")
fin_hh <- subset(findeln, CEN_1=="H-H");fin_hl <- subset(findeln, CEN_1=="H-L")
zen_hh <- subset(zeneggen, CEN_1=="H-H");zen_hl <- subset(zeneggen, CEN_1=="H-L")
fin_hh$rel.fitness <- fin_hh$comb_fit/mean(fin_hh$comb_fit);fin_hl$rel.fitness <- fin_hl$comb_fit/mean(fin_hl$comb_fit)
zen_hh$rel.fitness <- zen_hh$comb_fit/mean(zen_hh$comb_fit);zen_hl$rel.fitness <- zen_hl$comb_fit/mean(zen_hl$comb_fit)
data <- rbind(fin_hh, fin_hl, zen_hh, zen_hl)

#############Save relative fitness in new df for downstream analyses
data_aster_relative_fitness <- data  %>%  dplyr::select(barcode,rel.fitness)

