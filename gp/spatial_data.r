require("spBayes")
require("fields")
data("NETemp.dat")


z<-  (NETemp.dat$y.92<20)+1
plot(NETemp.dat$UTMX, NETemp.dat$UTMY, col=z)


colSums(NETemp.dat>30)


library("rgdal")
library("raster")
inFile <- "~/Downloads/mammals-critically-endangered/mammals_critically_endangered.tif"
library(gdalUtils) # For gdal_translate()
s <-raster(inFile)

GDALinfo(inFile)

SP27GTIF <- readGDAL(inFile, output.dim=c(200,200))

image(matrix(SP27GTIF@data$band1,200,200))


outFile <- "~/Downloads/mammals-critically-endangered/subset.tif"
ex <- extent(c(686040.1, 689715.9, 238156.3, 241774.2))

gdal_translate(inFile, outFile,
               projwin=c(xmin(ex), ymax(ex), xmax(ex), ymin(ex)))
