#output summary stats for color data, for species, tribes, families, biogeographic regions, continents

#library and useful function
{
library(here)
library(plyr)
library(dplyr)
library(picante)

#for making tables
#library(insight)
#library(gt)

genRandColors <- function(n = 10, plot = F)
{
  library(RColorBrewer)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'div',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  randCols <- sample(col_vector, n)
  if(plot == T){pie(rep(1,n), col=randCols)}
  return(randCols)
}

errorNoNA <- function(x = NA) {plotrix::std.error(x = x, na.rm = T)}

}

#load files and connect
{
#set working directory to parent directory
#data_pointer <- <YOUR PATH>

data_pointer <- "c:/Users/jtm378/AMNH Dropbox/Jon Merwin/AllParrotHomeAnalysis/Scripts/GITHUB_UPLOAD/"
setwd(data_pointer)

fullMorphRaw <- read.csv(file=("data/MorphometricMeasurements_8.8.25.csv"))

#write supplement table
#write.csv(fullMorphRaw, "Table_S1_raw_morphometric_measurements.csv")

#load taxon key
metadata_taxonomy <- read.csv("data/tip_hotswap_7.16.25.csv")

fullMorphRaw$tree_spp <- mapvalues(fullMorphRaw$CATALOG.NUMBER, from = metadata_taxonomy$Number, to = metadata_taxonomy$Species_in_merwin_et_al)
fullMorphRaw$tree_spp[which(is.na(fullMorphRaw$tree_spp))] <- fullMorphRaw$SPECIES[which(is.na(fullMorphRaw$tree_spp))]
fullMorphRaw$tree_spp <- mapvalues(fullMorphRaw$tree_spp, from = fullMorphRaw$CATALOG.NUMBER, to = fullMorphRaw$SPECIES)

#group_by_catalog number and get mean and std errors
#plot relationships between standard errors for each measurement
as_tibble(fullMorphRaw) %>% group_by(CATALOG.NUMBER) %>% summarise_if(is.numeric, plotrix::std.error, na.rm = T) -> morphDatError
plot(log(morphDatError[,c(3:9)]))

as_tibble(fullMorphRaw) %>% group_by(CATALOG.NUMBER) %>% summarise_if(is.numeric, mean, na.rm = T) -> fullMorphDat
fullMorphDat$tree_spp <- mapvalues(fullMorphDat$CATALOG.NUMBER, from = metadata_taxonomy$Number, to = metadata_taxonomy$Species_in_merwin_et_al)
fullMorphDat$tree_spp <- mapvalues(fullMorphDat$tree_spp, from = fullMorphRaw$CATALOG.NUMBER, to = fullMorphRaw$SPECIES)

#output summaries of morphometric measurements
summary(na.omit(fullMorphDat$BILL.LENGTH))
sd(na.omit(fullMorphDat$BILL.LENGTH))
summary(na.omit(fullMorphDat$BILL.WIDTH)); sd(na.omit(fullMorphDat$BILL.WIDTH))
summary(na.omit(fullMorphDat$BILL.HEIGHT)); sd(na.omit(fullMorphDat$BILL.HEIGHT))
summary(na.omit(fullMorphDat$TARSUS.LENGTH)); sd(na.omit(fullMorphDat$TARSUS.LENGTH))
summary(na.omit(fullMorphDat$TAIL)); sd(na.omit(fullMorphDat$TAIL))
}

#comparison of JTM and EL as measurer
{
png("measurement_ratio_by_measurer_comparison.png", res = 150, width = 20, height = 20, units = "in")
measColors <- plyr::mapvalues(x = fullMorphRaw$MEASURER, from = c("JTM", "EL"), to = c("lightblue", "orange"))
par(mar = c(2,2,2,2))
plot(log(fullMorphRaw[,c(5:11)]), col = measColors, cex = 1, pch = 16)
#text(log(fullMorphRaw[,c(6,10)]), labels = fullMorphRaw$CATALOG.NUMBER, cex = .5, pos = 4)
legend("bottomright", legend = c("EL", "JTM"), fill = c("orange", "lightblue"), xpd = T, cex = 1, horiz = TRUE, seg.len=1, bty = 'n')
dev.off()

system("open measurement_ratio_by_measurer_comparison.png")

#examine relationship between measurements and their error
billLengthmeasurementToErrorRatio <- fullMorphDat$BILL.LENGTH/morphDatError$BILL.LENGTH
billHeightmeasurementToErrorRatio <- fullMorphDat$BILL.HEIGHT/morphDatError$BILL.HEIGHT
billWidthmeasurementToErrorRatio <- fullMorphDat$BILL.WIDTH/morphDatError$BILL.WIDTH
TarsusmeasurementToErrorRatio <- fullMorphDat$TARSUS.LENGTH/morphDatError$TARSUS.LENGTH
tailMeasurementToErrorRatio <- fullMorphDat$TAIL/morphDatError$TAIL
primariesMeasurementToErrorRatio <- fullMorphDat$WING.LENGTH.TO.PRIMARIES/morphDatError$WING.LENGTH.TO.PRIMARIES
secondariesMeasurementToErrorRatio <- fullMorphDat$WING.LENGTH.TO.SECONDARIES/morphDatError$WING.LENGTH.TO.SECONDARIES

#make plots of error for each measurement
png(here::here("Figure_output/errorPlots.png"), width = 2000, height = 1000)
par(mfrow = c(1,2))
plot(billLengthmeasurementToErrorRatio[order(billLengthmeasurementToErrorRatio)], pch = 16, cex = 1,
     main = "Measurement/Error")
points(billHeightmeasurementToErrorRatio[order(billHeightmeasurementToErrorRatio)], pch = 16, cex = 1, col = "gold4")
points(billWidthmeasurementToErrorRatio[order(billWidthmeasurementToErrorRatio)], pch = 16, cex = 1, col = "lightblue")
points(TarsusmeasurementToErrorRatio[order(TarsusmeasurementToErrorRatio)], pch = 16, cex = 1, col = "darkgreen")
points(tailMeasurementToErrorRatio[order(tailMeasurementToErrorRatio)], pch = 16, cex = 1, col = "blue")
points(primariesMeasurementToErrorRatio[order(primariesMeasurementToErrorRatio)], pch = 16, cex = 1, col = "red")
points(secondariesMeasurementToErrorRatio[order(secondariesMeasurementToErrorRatio)], pch = 16, cex = 1, col = "green")
errLegCols <- c("gold4", "lightblue", "darkgreen", "blue", "red", "green")
errLegTraits <- c("Bill Height", "Bill Width", "Tarsus", "Tail", "Primaries", "Secondaries")
legend("topleft", legend = errLegTraits, fill = errLegCols, border = "transparent")

plot(log(fullMorphDat$BILL.LENGTH)~log(morphDatError$BILL.LENGTH), pch = 15, cex = 1, ylim = c(0, 10), xlim = c(-5, 5),
     main = "Per-Bird Error vs. Measurement")
points(log(fullMorphDat$BILL.HEIGHT)~log(morphDatError$BILL.HEIGHT), pch = 15, cex = 1, col = "gold4")
points(log(fullMorphDat$BILL.WIDTH)~log(morphDatError$BILL.WIDTH), pch = 15, cex = 1, col = "lightblue")
points(log(fullMorphDat$TARSUS.LENGTH)~log(morphDatError$TARSUS.LENGTH), pch = 15, cex =1, col = "darkgreen")
points(log(fullMorphDat$TAIL)~log(morphDatError$TAIL), pch = 15, cex =1, col = "blue")
points(log(fullMorphDat$WING.LENGTH.TO.PRIMARIES)~log(morphDatError$WING.LENGTH.TO.PRIMARIES), pch = 15, cex =1, col = "red")
graphics.off()
}

#process morph data, calculate kipps, output summary tables.
{
#calculate stats about kipps for good measurements
nobadKipps <- fullMorphDat[fullMorphDat$WING.LENGTH.TO.PRIMARIES > fullMorphDat$WING.LENGTH.TO.SECONDARIES, ]
hist(nobadKipps$WING.LENGTH.TO.PRIMARIES-nobadKipps$WING.LENGTH.TO.SECONDARIES)
summary(nobadKipps$WING.LENGTH.TO.PRIMARIES-nobadKipps$WING.LENGTH.TO.SECONDARIES)
sd(nobadKipps$WING.LENGTH.TO.PRIMARIES-nobadKipps$WING.LENGTH.TO.SECONDARIES)
summary(na.omit(fullMorphDat$TAIL))
head(nobadKipps[order(nobadKipps$WING.LENGTH.TO.PRIMARIES-nobadKipps$WING.LENGTH.TO.SECONDARIES),])

fullMorphDat$KIPPS_RATIO <- NA
fullMorphDat$KIPPS_RATIO <- (fullMorphDat$WING.LENGTH.TO.PRIMARIES-fullMorphDat$WING.LENGTH.TO.SECONDARIES)/fullMorphDat$WING.LENGTH.TO.PRIMARIES

##export summary tables, adding in sex, catalog number, species ID 
#Species <- plyr::mapvalues(x = as.character(fullMorphDat$CATALOG.NUMBER), from = as.character(fullMorphRaw$CATALOG.NUMBER), to = as.character(fullMorphRaw$SPECIES))
Species <- fullMorphDat$tree_spp
Sex <- plyr::mapvalues(x = as.character(fullMorphDat$CATALOG.NUMBER), from = as.character(fullMorphRaw$CATALOG.NUMBER), to = as.character(fullMorphRaw$SEX))
fullMorphDat %>% tibble::add_column(Species, .before = 1) %>% tibble::add_column(Sex, .after = 2) %>% select(-c("MEASUREMENT"))  %>% group_by(Species) %>% arrange(Species) -> summMorphTable_Export
write.csv(summMorphTable_Export, "2_processing_color_and_morphology_data/data/datasets_to_collate/Morph_collate.csv")
#format a nice table for sharing, requires table making packages
#summMorphTable_Export %>% format_table(pretty_names = T) %>% export_table(format = "html")

system("open Table_S2_summarized_morphometric_measurements.csv")

#connecting color and measurement data 
#load cleaned color data 
#this one has both males and females
cleanedColorDat <- read.csv("data/sexes_combined_spp_color_spaces.csv")

#cleanedColorDat$Sex <- plyr::mapvalues(x = as.character(cleanedColorDat$X), from = as.character(fullMorphRaw$CATALOG.NUMBER), to = as.character(fullMorphRaw$SEX))
#write.csv(cleanedColorDat, "Table_S3_summarized_color_measurements.csv")
#color data by specimen with SEX and a column
#cleanedColorDat %>% format_table(pretty_names = T) %>% export_table(format = "html")

#now colspan
#summStatsColorNoNA
head(cleanedColorDat[order(cleanedColorDat$colspan.m, decreasing = T),])
tail(cleanedColorDat[order(cleanedColorDat$colspan.m, decreasing = T),c("spp", "colspan.m")], 10)

#hue disparity
# summStatsColorNoNA
head(cleanedColorDat[order(cleanedColorDat$huedisp.m, decreasing = T),c("spp", "huedisp.m")])
tail(cleanedColorDat[order(cleanedColorDat$huedisp.m, decreasing = T),c("spp", "huedisp.m")], 10)

#saturation (mean.ra)
# summStatsColorNoNA
head(cleanedColorDat[order(cleanedColorDat$mean.ra, decreasing = T),c("spp", "mean.ra")], 10)
tail(cleanedColorDat[order(cleanedColorDat$mean.ra, decreasing = T),c("spp", "mean.ra")], 10)

# summStatsColorNoNA
head(cleanedColorDat[order(cleanedColorDat$a.vol, decreasing = T),c("spp", "a.vol")], 10)
tail(cleanedColorDat[order(cleanedColorDat$a.vol, decreasing = T),c("spp", "a.vol")], 10)
}

# #below commented-out scripts to output certain data for MUSSCRAT analysis
{
# allDat <- read.csv("data/PGLS_comparative_dataset.csv")
# summary(allDat)
# colnames(allDat)
# 
# par(mfrow = c(1,1))
# plot(allDat$XCentroid0.5~allDat$yCentroid0.5, pch = 16, col = as.factor(allDat$Tribe))
# summary(lm(allDat$XCentroid0.0~allDat$yCentroid0.0, pch = 16))
# legend("topleft", legend = levels(as.factor(allDat$Tribe)), fill = 1:length(unique(allDat$Tribe)), border = "transparent", cex = .5, bg = "transparent", box.col = "transparent")
# 
# #output data and tree for MUSSCRAT analysis
# allColDatPCA <- apply(X = allDat[,c(15:27)], MARGIN = 2, FUN = function(X) as.numeric(X))
# rownames(allColDatPCA) <- allDat$TreeNames
# colnames(allColDatPCA)
# 
# #color data PCA
# colPCA <- prcomp(na.exclude(allColDatPCA), center = T, scale. = T)
# dotplot(colPCA$rotation[,1:3], col = c("black", "red", "blue"), pch = 16, cex = 2)
# summary(colPCA)
# 
# pairs(colPCA$x[,1:3])
# 
# dataForBrian <- allDat[,c()]
# 
# colnames(allDat)
# colsForBrian <- c(2, 69, 70, 11:13, 23, 24, 28, 29, 30, 31, 76, 77)
# datForBrian <- allDat[,colsForBrian]
# rownames(datForBrian) <- allDat$TreeNames

# convertNames <- read.csv(file = here::here("Data/MetadataHelpers/newTreeNamesForPubMatching.csv"))
# testNames <- plyr::mapvalues(root$tip.label, from = convertNames$morphTaxNames, to = convertNames$TreeNames)
# 
# root$tip.label <- testNames
# 
# matchedTree <- picante::match.phylo.data(root, datForBrian)
# 
# 
# write.csv(x = datForBrian, file = here::here("MuSSCRat_DataForBrian_1.30.2023"))
# write.tree(phy = matchedTree$phy, file = here::here("MuSSCRat_TreeForBrian_1.30.2023"))
 }






