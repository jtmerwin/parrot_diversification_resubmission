#TO DO CLEAN UP FOR PUB!

#import packages and useful functions...set wd
{
require("geiger")
require("plyr")
require("dplyr")
library("randomcoloR")
library("here")
library("colorspace")
library("pavo")
library("ggplot2")

#useful method for matching later 
lastDropper <- function(string = NULL, seper = " ")
{
  require("stringr")
  lab <- word(string, 1, -2, sep = seper)
  return(lab)
}

#set working directory to parent directory
#data_pointer <- <YOUR PATH>
data_pointer <- "c:/Users/jtm378//AMNH Dropbox/Jon Merwin/AllParrotHomeAnalysis/Scripts/GITHUB_UPLOAD/"
setwd(data_pointer)
#import tree

#taxonomy guide import
metadata_taxonomy <- read.csv("data/Specimen_metadata_pub_newsamples.csv")
tre <- read.nexus("data/parrot_timetree.tre")
rootTre <- root.phylo(tre, "Caracara_cheriway_b5328")

functional_names <- intersect(gsub(lastDropper(tre$tip.label, seper = "_"), pattern = "_", replacement = " "), metadata_taxonomy$Species_in_merwin_et_al) 
broken_names_1 <- setdiff(x= gsub(lastDropper(tre$tip.label, seper = "_"), pattern = "_", replacement = " "), y = metadata_taxonomy$Species_in_merwin_et_al) 
broken_names_2 <- setdiff(y = gsub(lastDropper(tre$tip.label, seper = "_"), pattern = "_", replacement = " "), x = metadata_taxonomy$Species_in_merwin_et_al) 

#picking one tree to work with
root <- rootTre
root$tip.label <- lastDropper(root$tip.label, seper = "_")

#plotting it with a color for each genus
generaList <- str_split(root$tip.label, pattern  = "_", simplify = T)[,1]
genusColors <- mapvalues(generaList, from = unique(generaList), to =  distinctColorPalette(100))
plot(root, show.tip.label = T, cex = .2, tip.color = genusColors)
}

#import and clean color data
{
##load in the data
#setwd(paste(dropboxPointer, "/Dropbox (AMNH)/AllParrotHomeAnalysis/Scripts/AllParrotAnalysis/Data/RawColorData", sep = ""))
temp = list.files(path = "2_processing_color_and_morphology_data/data/RawColorData/", pattern="*.csv", full.names = T)
print(temp)

myfiles = lapply(temp, read.delim, sep = ",", header = T)

#combine the days data into one data frame
l <- myfiles
colnames(l[[23]]) <- colnames(l[[20]])

#values below 0 can result from digital noise and cause later software errors 
count = 1
for(bb in l)
{
  if(min(bb$lwMean) < 0){l[[count]]$lwMean <- l[[count]]$lwMean - min(l[[count]]$lwMean); print(min(l[[count]]$lwMean))}
  if(min(bb$mwMean) < 0){l[[count]]$mwMean <- l[[count]]$mwMean - min(l[[count]]$mwMean)}
  if(min(bb$swMean) < 0){l[[count]]$swMean <- l[[count]]$swMean - min(l[[count]]$swMean)}
  if(min(bb$uvMean) < 0){l[[count]]$uvMean <- l[[count]]$uvMean - min(l[[count]]$uvMean)}
  print(count)
  print(colnames(bb))
  print(temp[count])
  count <- count+1
}

#total_data <- do.call(rbind, lapply(l, function(x) x[match(names(l[[1]]), names(x))]))

total_data <- bind_rows(l, .id = "source")

summary(total_data$lwMean)
windows()
plot(log((total_data$mwMean)+(total_data$uvMean)+(total_data$lwMean)+(total_data$swMean)), pch = 16, cex = .02, col = total_data$source)

#ggplot(total_data, aes(y = ((mwMean-uvMean)), group = source)) + geom_boxplot()

#extract genus, spp, specimen number, other metadata from data labels
library(stringr)
headerSplit <- str_split_fixed(total_data$Label, pattern = "_", n = Inf)
colnames(headerSplit) <- c("Genus", "Species", "SpecimenNumber", "Resolution", "ImgNum", "View", "Patch", "Subsample", "end", "end2")

total_data_with_metadata <- cbind(headerSplit, total_data)
total_data_metadata_cleaned <- (total_data_with_metadata)

#this section removes data from template mspecs and connects data from 
#larger birds which had to be separated into top and bottom photos 
#remove test rows
testRows <- grep(pattern = "Test", x = total_data_metadata_cleaned$Genus, fixed = F, ignore.case = T)

#rowsForRemove
#find wrong tops
topPatches <- grep(pattern = "top", x = total_data_metadata_cleaned$Label, fixed = F, ignore.case = T)
bottoms <- c("UpperTailCoverts", "Rump", "Crissum", "CruralFeathers", "LowerAbdomen", "UpperAbdomen", "FoldedPrimaries", "FoldedSecondaries")
whichBottoms <- which(apply(total_data_metadata_cleaned, 1, function(r) any(r %in% bottoms)))
topRemoves <- intersect(topPatches, whichBottoms)

#find wrong bottoms
botPatches <- grep(pattern = "bot", x = total_data_metadata_cleaned$Label, fixed = F, ignore.case = T)
botKeeps <- intersect(botPatches, whichBottoms)
botDrops <- setdiff(botPatches, botKeeps)

#remove both

total_data_metadata_cleaned_dropped <- total_data_metadata_cleaned[-c(testRows, topRemoves, botDrops),]
}

#swap bad numbers, typos from the image titles. 
{
badNumber_hotSwap <- read.csv("data/tip_hotswap_7.15.25.csv")

parrotData <- total_data_metadata_cleaned_dropped
totColsAllPoints <- rgb(red = (parrotData$lwMean), green = (parrotData$mwMean), blue = (parrotData$swMean), maxColorValue = 1.5)
totColsFixed <- lighten(totColsAllPoints, .2, fixup = T, method = "relative", space = "HLS")

UVtoMW <- (parrotData$uvMean-parrotData$mwMean)
SWtoLW <- (parrotData$swMean-parrotData$lwMean)
BrightnessAllDat <- (parrotData$uvMean+parrotData$mwMean+parrotData$swMean+parrotData$lwMean)
SDAllDat <- (parrotData$uvSD+parrotData$mwSD+parrotData$swSD+parrotData$lwSD)

#these diagnostic color opponent plots take a long time to plot. run at own risk
#plot(UVtoMW~SWtoLW, col= totColsFixed, pch = 16, cex = .5, xlim = c(-.5,.2), ylim = c(-.5,.2))
#plot(UVtoMW~SWtoLW, col= totColsAllPoints, pch = 16, cex = .5)
#plot(SDAllDat~BrightnessAllDat, col= totColsFixed, pch = 16, cex = 1)

parrotData$species_identity <- mapvalues(x=parrotData$SpecimenNumber, from = metadata_taxonomy$Database_num, to = metadata_taxonomy$Species_in_merwin_et_al)

#summarize the multiples into the final full color dataset
as_tibble(parrotData) %>% group_by(species_identity, SpecimenNumber, Patch) %>% summarise_at(.vars = c("uvMean", "swMean", "mwMean", "lwMean"), .funs = mean) -> ParrotTib

#export color dataset for each bird. This creates a supplementary dataset. 
#write.csv(ParrotTib, file = paste(dropboxPointer, "Dataset_Sx_color_data_by_patch", sep = "")) 

#calculate alpha and non-alpha color volume percentages
{
full_color_space <- colspace(vismodeldata = data.frame(u = c(1,0,0,0), s = c(0,1,0,0), m = c(0,0,1,0), l = c(0,0,0,1)), space = "tcs")
spaceAll <- pavo::colspace(data.frame(ParrotTib[,c("uvMean", "swMean", "mwMean", "lwMean")]), space = "tcs", qcatch = "Qi")

#color space calculations here with alphashape3d package which must be manually installed
#percentageAlpha <- alphashape3d::volume_ashape3d(alpha_allParrots) / total_volume_colspace
#percentageNo_Alpha <- alphashape3d::volume_ashape3d(no_alpha_allParrots) / total_volume_colspace
#print(paste("Alpha percentage (a = .079): ", percentageAlpha))
#print(paste("Non-Alpha (convex hull) percentage: ", percentageNo_Alpha))
}

#create hex colors for each patch for each bird
totColsAllPoints_tib <- rgb(red = (ParrotTib$lwMean), green = (ParrotTib$mwMean), blue = (ParrotTib$swMean), maxColorValue = 1.5)
totColsFixed_tib <- lighten(totColsAllPoints_tib, .2, fixup = T, method = "relative", space = "HLS")


#below make exploratory plots of color volumes
#plot(spaceAll, col = totColsFixed_tib, cex = 1)
#plot(spaceAll, col = "transparent", cex = 1)
#vol(spaceAll, grid = T, fill = T, type = "alpha", avalue = .07, alpha = .2, col = "gold4")
#plot(spaceAll, col = "transparent", cex = 1)
#vol(spaceAll, grid = T, fill = F, alpha = .2, col = "red4")


#png(filename = "c:/Users/jtm378/Dropbox (AMNH)/AllParrotHomeAnalysis/Scripts/AllParrotAnalysis/FigureOutput/allParrotColorHex.png", width = 5000, height = 5000)
#pavo::hexplot(colspace(abs(parrotData[,c("uvMean", "swMean", "mwMean", "lwMean")]), qcatch = "Qi", space = "tcs"), col = totColsFixed, pch = 15, cex = 1)
#dev.off()

#separate Male and Female Color Data
#import taxonomic key
morphData <- read.delim("data/raw_morphology_data.csv")
sampList <- read.delim("data/PGLS_comparative_dataset.csv",sep = ",")
maleNums <- unique(metadata_taxonomy$Database_num[which(metadata_taxonomy$Sex == "M")])
femaleNums <- unique(metadata_taxonomy$Database_num[which(metadata_taxonomy$Sex == "F")])

maleCols <- total_data_metadata_cleaned_dropped[total_data_metadata_cleaned_dropped$SpecimenNumber %in% maleNums, ]
femaleCols <- total_data_metadata_cleaned_dropped[total_data_metadata_cleaned_dropped$SpecimenNumber %in% femaleNums, ]

total_data_metadata_cleaned_dropped$taxName <- paste(total_data_metadata_cleaned_dropped$Genus, total_data_metadata_cleaned_dropped$Species, sep = "_")

data <- total_data_metadata_cleaned_dropped  
}

#continue with color_make_tile_matrix.R to output figures and tables for manuscript.
###STOP FOR Merwin et al., 2026###