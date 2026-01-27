#WRITTEN BY ELKIN TENORIO
#ACCESS @ https://github.com/Elkin01/Parrots_NicheEvolution/
rm(list=ls())

library(raster)
library(terra)
library(ecospat)
library(sf)
library(sp)
library(ade4)
library(rgdal)
library(rgeos)

## Construction of global environmental  space

elev=raster('~/Library/CloudStorage/OneDrive-CityUniversityofNewYork/Documentos/Estadistica R SIG/Capas SIG/Rasters/Altitude_World.tif')
  
bio1=raster('~/Library/CloudStorage/OneDrive-CityUniversityofNewYork/Documentos/Estadistica R SIG/Capas SIG/Rasters/WorldClim/World/wc2.1_30s_bio_1.tif')
bio2=raster('~/Library/CloudStorage/OneDrive-CityUniversityofNewYork/Documentos/Estadistica R SIG/Capas SIG/Rasters/WorldClim/World/wc2.1_30s_bio_2.tif')
bio3=raster('~/Library/CloudStorage/OneDrive-CityUniversityofNewYork/Documentos/Estadistica R SIG/Capas SIG/Rasters/WorldClim/World/wc2.1_30s_bio_3.tif')
bio4=raster('~/Library/CloudStorage/OneDrive-CityUniversityofNewYork/Documentos/Estadistica R SIG/Capas SIG/Rasters/WorldClim/World/wc2.1_30s_bio_4.tif')
bio5=raster('~/Library/CloudStorage/OneDrive-CityUniversityofNewYork/Documentos/Estadistica R SIG/Capas SIG/Rasters/WorldClim/World/wc2.1_30s_bio_5.tif')
bio6=raster('~/Library/CloudStorage/OneDrive-CityUniversityofNewYork/Documentos/Estadistica R SIG/Capas SIG/Rasters/WorldClim/World/wc2.1_30s_bio_6.tif')
bio7=raster('~/Library/CloudStorage/OneDrive-CityUniversityofNewYork/Documentos/Estadistica R SIG/Capas SIG/Rasters/WorldClim/World/wc2.1_30s_bio_7.tif')
bio8=raster('~/Library/CloudStorage/OneDrive-CityUniversityofNewYork/Documentos/Estadistica R SIG/Capas SIG/Rasters/WorldClim/World/wc2.1_30s_bio_8.tif')
bio9=raster('~/Library/CloudStorage/OneDrive-CityUniversityofNewYork/Documentos/Estadistica R SIG/Capas SIG/Rasters/WorldClim/World/wc2.1_30s_bio_9.tif')
bio10=raster('~/Library/CloudStorage/OneDrive-CityUniversityofNewYork/Documentos/Estadistica R SIG/Capas SIG/Rasters/WorldClim/World/wc2.1_30s_bio_10.tif')
bio11=raster('~/Library/CloudStorage/OneDrive-CityUniversityofNewYork/Documentos/Estadistica R SIG/Capas SIG/Rasters/WorldClim/World/wc2.1_30s_bio_11.tif')
bio12=raster('~/Library/CloudStorage/OneDrive-CityUniversityofNewYork/Documentos/Estadistica R SIG/Capas SIG/Rasters/WorldClim/World/wc2.1_30s_bio_12.tif')
bio13=raster('~/Library/CloudStorage/OneDrive-CityUniversityofNewYork/Documentos/Estadistica R SIG/Capas SIG/Rasters/WorldClim/World/wc2.1_30s_bio_13.tif')
bio14=raster('~/Library/CloudStorage/OneDrive-CityUniversityofNewYork/Documentos/Estadistica R SIG/Capas SIG/Rasters/WorldClim/World/wc2.1_30s_bio_14.tif')
bio15=raster('~/Library/CloudStorage/OneDrive-CityUniversityofNewYork/Documentos/Estadistica R SIG/Capas SIG/Rasters/WorldClim/World/wc2.1_30s_bio_15.tif')
bio16=raster('~/Library/CloudStorage/OneDrive-CityUniversityofNewYork/Documentos/Estadistica R SIG/Capas SIG/Rasters/WorldClim/World/wc2.1_30s_bio_16.tif')
bio17=raster('~/Library/CloudStorage/OneDrive-CityUniversityofNewYork/Documentos/Estadistica R SIG/Capas SIG/Rasters/WorldClim/World/wc2.1_30s_bio_17.tif')
bio18=raster('~/Library/CloudStorage/OneDrive-CityUniversityofNewYork/Documentos/Estadistica R SIG/Capas SIG/Rasters/WorldClim/World/wc2.1_30s_bio_18.tif')
bio19=raster('~/Library/CloudStorage/OneDrive-CityUniversityofNewYork/Documentos/Estadistica R SIG/Capas SIG/Rasters/WorldClim/World/wc2.1_30s_bio_19.tif')
  
bio=stack(bio1, bio2, bio3, bio4, bio5, bio6, bio7, bio8, bio9, bio10, bio11, bio12, bio13, bio14, bio15, bio16, bio17, bio18, bio19)

elev2=aggregate(elev, 5, 'mean')

mundo=st_read('~/Library/CloudStorage/OneDrive-CityUniversityofNewYork/Documentos/Estadistica R SIG/Capas SIG/Shapes/Political_World_Borders/TM_WORLD_BORDERS_SIMPL-0.3/TM_WORLD_BORDERS_SIMPL-0.3.shp', 'TM_WORLD_BORDERS_SIMPL-0.3')
mundo <- mundo[1] 
mundo=as_Spatial(mundo)

#Selection of 100000 random points worldwide to build the background environmental space for all species
coordenadas=sampleRandom(elev, 100000, na.rm=TRUE, cells=TRUE, sp=TRUE) 

filas=as.data.frame(coordenadas)
  
#plot(elev)
#plot(coordenadas, add=T)
  
clim=data.frame(filas[,3:4], extract(bio, coordenadas))
clim=na.omit(clim)
colnames(clim)=c('long', 'lat', colnames(clim)[3:21])

#load coordinate for each species

CoorSpp=list.files('~/Dropbox (AMNH)/parrot_morpho_evolution/ebird_data_with_headers/Filtered Ebird records', pattern = '1_')

setwd('~/Dropbox (AMNH)/parrot_morpho_evolution/ebird_data_with_headers/Filtered Ebird records')

occ.sp=data.frame()

for (i in 97:length(CoorSpp))
  
{
  print(i)
  
  spData=read.table(CoorSpp[i], header=T)
  
  spp=spData[,c('SCIENTIFIC.NAME', 'LONGITUDE', 'LATITUDE')]
  
  spp=spp[!duplicated(spp), ]
  
  Coor=spp[,c('LONGITUDE', 'LATITUDE')]
  
  coordinates(Coor)=c('LONGITUDE', 'LATITUDE')
  
  occ.sp_test=data.frame(spp, extract(bio, Coor))
  
  occ.sp_test=na.omit(occ.sp_test)
  colnames(occ.sp_test)=c('species', 'long', 'lat', colnames(occ.sp_test)[4:22])
  
  occ.data <- cbind(occ.sp_test, occ.sp_test[, 1]) # add species names
  
  occ.sp=rbind(occ.sp, occ.data)
  
}


# list of species
sp.list <- unique(occ.sp[, 1])
sp.nbocc <- c()

for (i in 1:length(sp.list)) {
  sp.nbocc <- c(sp.nbocc, length(which(occ.sp[, 1] == sp.list[i])))
}

#sp.nbocc <- c(sp.nbocc, length(which(occ.sp[, 1] == sp.list[1])))

# calculate the nb of occurences per species
sp.list <- sp.list[sp.nbocc > 4] # remove species with less than 5 occurences


# selection of variables to include in the analyses
Xvar <- c(3:21)
nvar <- length(Xvar)

#################################### PCA-ENVIRONMENT ##################################
data <- rbind(occ.sp[, Xvar + 1], clim[, Xvar])
w <- c(rep(0, nrow(occ.sp)), rep(1, nrow(clim)))
pca.cal <- dudi.pca(data, row.w = w, center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)


#################################################
#################################################
# Loop for niche estimation
#################################################
#################################################


pdf('Niche_plots.pdf', width=25, height=10)

NicheResults=data.frame()

for (i in 314:length(CoorSpp))
  
{
  
  print(i)

row.sp1 <- which(occ.sp[, 1] == sp.list[i]) # rows in data corresponding to sp1

scores.clim <- pca.cal$li[(nrow(occ.sp) + 1):nrow(data), ] # scores for global climate

scores.sp1 <- pca.cal$li[row.sp1, ] # scores for sp1

# calculation of occurence density and test of niche equivalency and similarity
# with the default kernel method
z1 <- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.sp1, R = 100)

Niche_0=z1$z.uncor
Niche_0.5=z1$z.uncor

Niche_0[which(values(Niche_0)==0)]<-NA
Niche_0.5[which(values(Niche_0.5)<0.5)]<-NA

r.to.poly_0.0 <- sf::as_Spatial(sf::st_as_sf(stars::st_as_stars(Niche_0), 
                                         as_points = FALSE, merge = TRUE))

r.to.poly_0.5 <- sf::as_Spatial(sf::st_as_sf(stars::st_as_stars(Niche_0.5), 
                                         as_points = FALSE, merge = TRUE))

st_write(st_as_sf(r.to.poly_0.0), "/Users/elkintenorio/Dropbox (AMNH)/parrot_morpho_evolution/ebird_data_with_headers/Filtered Ebird records/Shapefiles P0", layer=paste(sp.list[i], '_0.0'), driver = "ESRI Shapefile")
st_write(st_as_sf(r.to.poly_0.5), "/Users/elkintenorio/Dropbox (AMNH)/parrot_morpho_evolution/ebird_data_with_headers/Filtered Ebird records/Shapefiles P0.5", layer=paste(sp.list[i], '_0.5'), driver = "ESRI Shapefile")

NicheResults[i,1]<-sp.list[i]
NicheResults[i,2]<-length(which(values(z1$z.uncor)>0))/length(values(z1$z.uncor)) #he density of occurrence of the species
NicheResults[i,3]<-length(which(values(z1$z.uncor)>0.5))/length(values(z1$z.uncor)) #he density of occurrence of the species
NicheResults[i,4]<-gCentroid(r.to.poly_0.0)@coords[1]
NicheResults[i,5]<-gCentroid(r.to.poly_0.0)@coords[2]
NicheResults[i,6]<-gCentroid(r.to.poly_0.5)@coords[1]
NicheResults[i,7]<-gCentroid(r.to.poly_0.5)@coords[2]


#length(which(values(z1$z.cor)>0))/length(values(z1$z.cor)) #the occupancy of the environment by the species (density of occurrences divided by the desinty of environment in the study area

par(mfrow=c(1,3))
plot(mundo)
occ.sp2=occ.sp[which(occ.sp[, 1] == sp.list[i]), c('long', 'lat')]
coordinates(occ.sp2)=c('long', 'lat')
plot(occ.sp2, add=T, pch=21, col='red')

EnvNiche_0=length(which(values(z1$z.uncor)>0))/length(values(z1$z.uncor))
EnvNiche_0.5=length(which(values(z1$z.uncor)>0.5))/length(values(z1$z.uncor))

plot(z1$z.uncor, main = paste(sp.list[i], "with default kernel", sep=' ')) #he density of occurrence of the species
text(x=4, y=4.6, paste("EnvNiche_0 =", EnvNiche_0, sep=" "))
text(x=4, y=3.8, paste("EnvNiche_0.5 =", EnvNiche_0.5, sep=" "))

hist(extract(elev, occ.sp2), xlim=c(0, 6500))


}

dev.off()

colnames(NicheResults)=c('Sp','Total niche','Niche over 50%','XCentroid0.0','yCentroid0.0','XCentroid0.5','yCentroid0.5')

write.table(NicheResults, file = 'NicheEstimations.txt', sep = '\t')


library(ggplot2)
library(factoextra)
library(ade4)

summary(pca.cal)

fviz_eig(pca.cal)

fviz_pca_var(pca.cal,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)     # Avoid text overlapping
  
s.corcircle(pca.cal$co)

eig.val <- get_eigenvalue(pca.cal)
eig.val

# Results for Variables - Loads
res.var <- get_pca_var(pca.cal)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 

quartz()
par(oma=c(1,2,1,0))
plot(pca.cal$tab[,1], pca.cal$tab[,2], pch=21, ylab='Principal component 2', xlab='Principal component 1', cex.lab=1.4)
cor(pca.cal$tab[,1], pca.cal$tab[,2])
