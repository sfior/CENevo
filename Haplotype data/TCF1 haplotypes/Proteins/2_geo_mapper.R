# This scipt plots frequencies of protein variants as pie charts on a map

library(maps)
library(mapdata)
library(maptools)  #for shapefiles
library(scales)  #for transparency
library(mapplots)
library(raster)
library(RColorBrewer)

# Import data
data<-read.table("Protein_freqs.txt",header=T,sep='\t')

# add jitter to the pie charts, in case they are too close
data[,2][2]<-46.2  #Gibidumsee
data[,3][2]<-8

#dev.new(width=15, height=15)
# Define the colours for the pie chart (number of colours should be equal to number of haplotypes)
my_col<-c(alpha('blue', 0.6), alpha('green', 0.6), alpha('red',0.6), alpha('black',0.6),alpha('magenta',0.6),alpha('yellow',0.6),alpha('orange',0.6))

# Download map of Switzerland
swiss_dem <- getData('alt' , country='CHE', mask=TRUE)
swiss <- getData('GADM' , country="CHE", level=0)
# Or for a high resolution DEM
alps.hiRes <- readRDS("ELEV.GMTED.alps.hiRes.rds")

plot(swiss_dem)
#plot(alps.hiRes, col = brewer.pal(9, "Greys"))

# Crop
e <- extent(7, 8.5, 45.55, 46.55)
#switzerland <- crop(swiss_dem, e)
switzerland <- crop(alps.hiRes, e)

# Plot the map
#par(mar=c(5.1,4.1,4.1,2.1))
plot(switzerland, col=rev(grey(0:100/100)))

# Add the pie charts.
for (n in 1:nrow(data)){
  line<-data[n,]
  add.pie(z=c(line[,5],line[,6],line[,7],line[,8],line[,9],line[,10],line[,11]), x=line[,3], y=line[,2], radius=sqrt(line[,4]/4000), col=my_col, labels='')
}




