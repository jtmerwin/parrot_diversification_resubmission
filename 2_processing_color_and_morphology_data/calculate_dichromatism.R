#getting dimorphism by patch

genDimorphPatchScores <- function(sppTest = NULL, colDatForAnalysis = parrotData)
{
  sppIndx  <- which(paste(colDatForAnalysis$Genus,colDatForAnalysis$Species, sep = "_")==sppTest) 
  sppDat <- colDatForAnalysis[sppIndx,]
  nums <- unique(sppDat$SpecimenNumber)
  
  Dat1 <- sppDat[grep(pattern = nums[1], x = sppDat$SpecimenNumber),]
  Dat2 <- sppDat[grep(pattern = nums[2], x = sppDat$SpecimenNumber),]
  
  patchList <- c("Lores", "Back", "Scapulars", "Crissum")
  patchList <- unique(parrotData$Patch)[1:35]
  outputScore <- NULL
  outputScore <- data.frame(patch = patchList, score = NA)
  
  for(thisPatch in patchList)
  {
    if(length(patchDatSp1)==16 && length(patchDatSp2)==16){
    
    patchDatSp1 <- pavo::colspace(Dat1[grep(pattern = thisPatch, x = Dat1$Label),c("uvMean", "swMean", "mwMean", "lwMean")], space = "tcs", "Qi")
    patchDatSp2 <-  pavo::colspace(Dat2[grep(pattern = thisPatch, x = Dat2$Label),c("uvMean", "swMean", "mwMean", "lwMean")], space = "tcs", "Qi")
    
    
    distMx <- dist(rbind(patchDatSp1[,1:4], patchDatSp2[,1:4]))
    outputScore[which(outputScore$patch==thisPatch),"score"] <- mean(distMx)
    }
  }
  
  return(as.numeric(t(outputScore$score)))
}

allDimorphismPatchDist <- data.frame(spp = unique(paste(parrotData$Genus, parrotData$Species,sep = "_")))

genDimorphPatchScores(sppTest = "Eclectus_roratus")
genDimorphPatchScores("Ara_macao")
genDimorphPatchScores("Pyrrhura_peruviana")

testOutDim <- sapply(allDimorphismPatchDist$spp, FUN = function(X) genDimorphPatchScores(sppTest = X))
rownames(testOutDim) <- unique(parrotData$Patch)[1:35]

write.csv(testOutDim, file = "c:/Users/jtm378/Dropbox (AMNH)/AllParrotHomeAnalysis/Scripts/AllParrotAnalysis/Data/outputDichromatismScoresByPatch.7.21.2022.csv")


testOutDim <- read.csv("c:/Users/jtm378/Dropbox (AMNH)/AllParrotHomeAnalysis/Scripts/AllParrotAnalysis/Data/outputDichromatismScoresByPatch.7.21.2022.csv")

plot(testOutDim[,1], col = "transparent", pch = 16, ylim = c(0,1))
for(i in 1:389)
{
  points(testOutDim[,i], pch = 16, cex = 2, col = alpha("red", .01))

}


par(mar = c(5,5,5,5))
#boxplot(t(testOutDim), names = unique(parrotData$Patch)[1:35])

library(tidyverse)

dfDich <- data.frame((testOutDim))
rownames(dfDich) <- unique(parrotData$Patch)[1:35]
dfDich <- apply(dfDich, MARGIN = 2, function(X) as.numeric(X))
dichromPC <- prcomp(na.omit(t(dfDich)))
rownames(dichromPC$rotation) <- unique(parrotData$Patch)[1:35]
summary(dichromPC)


pcWeightChart <- data.frame(dichromPC$rotation[,1:5]) %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
head(pcWeightChart)


#Now we are ready to use the ggplot2 to plot the heatmap using geom_tile.
pdf(file = "c:/Users/jonme/Dropbox (AMNH)/AllParrotHomeAnalysis/Scripts/AllParrotAnalysis/FigureOutput/Dichrom.PCWeights.pdf", width = 30, height = 30)
ggplot(pcWeightChart, aes(x = rowname, y = colname, fill = (value)), pal = "pal12") +
  geom_tile()
dev.off()

dotplot(dichromPC$rotation[,1], col = c("red"), pch = 16, cex = 2)
dotplot(dichromPC$rotation[,2], col = "blue", pch = 16)
rownames(dichromPC$x)[order(dichromPC$x[,1])]

dt2 <- rbind(data.frame(dfDich),colMeans((dfDich))) %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
head(dt2)

thresh_dichrom <- mean(log(na.omit(dt2$value)))+ 2*sd(log(na.omit(dt2$value)))
dt2$significant_dichrom <- (log(dt2$value) > thresh_dichrom)

hist(dt2$value^2) 
levels(dt2$rowname) <- as.character(c(1:36))
#Now we are ready to use the ggplot2 to plot the heatmap using geom_tile.
pdf(file = "c:/Users/jtm378/Dropbox (AMNH)/AllParrotHomeAnalysis/Scripts/AllParrotAnalysis/FigureOutput/Dichrom.Heatmap.pdf", width = 10, height = 60)
ggplot(dt2, aes(x = rowname, y = colname, fill = (value)^2)) +
  geom_tile() + theme(panel.background = element_rect(fill = "white", colour = "transparent"))
dev.off()

png(file = "c:/Users/jtm378/Dropbox (AMNH)/AllParrotHomeAnalysis/Scripts/AllParrotAnalysis/FigureOutput/Dichrom.Heatmap.png", width = 10, height = 60, units = "in", res = 300)
ggplot(dt2, aes(x = rowname, y = colname, fill = (value)^2)) +
  geom_tile() + theme_minimal()
dev.off()

#Now we are ready to use the ggplot2 to plot the heatmap using geom_tile.
pdf(file = "c:/Users/jtm378/Dropbox (AMNH)/AllParrotHomeAnalysis/Scripts/AllParrotAnalysis/FigureOutput/sigDichrom.pdf", width = 10, height = 60)
ggplot(dt2, aes(x = rowname, y = colname, fill = significant_dichrom)) +
  geom_tile() + theme(panel.background = element_rect(fill = "white", colour = "transparent"))
dev.off()

hist((dt2$value))





