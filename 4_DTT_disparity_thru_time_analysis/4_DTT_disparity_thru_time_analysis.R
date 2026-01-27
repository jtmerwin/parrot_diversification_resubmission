#make DTT plots for color and morphology
#make single column matrices for UV reflectance, 

#manuscript PCA plots
#library and custom function loading
{
  library(ape)
  library(geiger)
  library(phytools)
  library(caper)
  library(picante)
  library(plyr)
  library(dplyr)
  library(dispRity)
  
  #lastDropper
  lastDropper <- function(string = NULL, seper = " ")
  {
    require("stringr")
    lab <- word(string, 1, -2, sep = seper)
    return(lab)
  }
  
  genRandColors <- function(n = 10, plot = F)
  {
    library(RColorBrewer)
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'div',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    randCols <- sample(col_vector, n)
    if(plot == T){pie(rep(1,n), col=randCols)}
    return(randCols)
  }
}

#data loading and diagnostic plot
{
  
  #set working directory to parent directory
  #data_pointer <- <YOUR PATH>
  data_pointer <- "c:/Users/jonme/AMNH Dropbox/Jon Merwin/AllParrotHomeAnalysis/Scripts/GITHUB_UPLOAD/"
  setwd(data_pointer)
  
  # treeRead <- read.nexus(file = "data/parrot_timetree.tre")
  # treeRead$tip.label <- lastDropper(treeRead$tip.label, "_")
  # 
  # datRead <- read.csv(file = "data/PGLS_comparative_dataset.csv")
  # datRead$colsPlot <- rgb(red = datRead$Average.of.centroid.l, green = datRead$Average.of.centroid.m, blue = datRead$Average.of.centroid.s)
  # datRead$contCols <- mapvalues(x = datRead$Continent, from = c("Oceania ", "Oceania", "Eurasia", "North America", "South America", "Africa"), to = c("blue", "blue", "green", "red", "orange", "pink"))

  # rownames(datRead) <- datRead$TreeNames
  # 
  # allMatched <- match.phylo.data(phy = treeRead, data = datRead)
  # allMatched$data$Continent <- gsub(allMatched$data$Continent, pattern = "Oceania ", replacement = "Oceania")
  # plot(as.numeric(allMatched$data$PC2)~log(110+(as.numeric(allMatched$data$PC1)*-1)), col = allMatched$data$colsPlot, pch = 16, cex = 2)
  
  allMatched <- readRDS("data/musscrat_matched.rdj")
}
#source envelope method
{
  #This is a modification of the MDI function in Geiger
  #Here the p-value is computed as in Slater GJ, Price, SA, Santini, F, Alfaro, MA. 2010. Proceedings of the Royal Society B. 277: 3097 -3104.
  #However, the original code is for a one-tailed test, and the modification below is for a two-tailed test.
  
  #The user needs to realise that since the test is two sided the reported p-value needs to be less than 0.025 for significance at the 5% level of significance
  
  #This function is called by dtt1 (which is turn a modification from Geiger).
  #
  #
  
  getMDIp<-function(dttRes) {
    foo<-function(x) {
      return(geiger:::.area.between.curves(x= dttRes$times, f1=x, f2=dttRes$dtt))
    }
    mdis<-apply(dttRes$sim,2,foo)
    
    #Two sided test
    p1<-length(which(mdis>=0))/length(mdis)
    p2<-length(which(mdis<=0))/length(mdis)
    print(min(p1))
    pVal<-min(p1,p2)
    return(pVal)
  }
}

#method to run DTT, and run time-jenk p-value MDI analysis to test for divergence from BM
makeDTT_parrotData <- function(traitName = NULL, nSim = 10, rescale = T, jenk_count = 1, calculate_branching_jenks = T)
{
  if(!rescale){hueDisp_Mean <- (allMatched$data[,traitName])}
  if(rescale){hueDisp_Mean <- as.numeric(na.omit(as.numeric(allMatched$data[,traitName])))}
  names(hueDisp_Mean) <- rownames(allMatched$data)[c(which(!is.na(as.numeric(allMatched$data[,traitName]))))]
  if(!rescale){names(hueDisp_Mean) <- rownames(allMatched$data)[c(which(!is.na((allMatched$data[,traitName]))))]}
  
  if(rescale == T){hueDisp_mean <- log(scales::rescale(hueDisp_Mean, to = c(1,2)))}
  matchedTree_DTT <- match.phylo.data(allMatched$phy, (hueDisp_Mean))
  
  if(rescale){dtt_calculated <- geiger::dtt(phy = matchedTree_DTT$phy, data = na.omit(matchedTree_DTT$data), nsim = nSim, calculateMDIp = T)}
  if(!rescale){dtt_calculated <- geiger::dtt(phy = matchedTree_DTT$phy, data = na.omit(matchedTree_DTT$data), nsim = 0, calculateMDIp = F, index = "num.states")}
  
  time_jenks <- c(0:(jenk_count))/jenk_count
  if(calculate_branching_jenks)
  {
    jenk_count <- jenk_count+1
    rescaledBranchTimes <- scales::rescale(branching.times(matchedTree_DTT$phy), from = c(max(branching.times(matchedTree_DTT$phy)), 0), to = c(0,1))
    #cluster branching time to 
    clusters <- kmeans(rescaledBranchTimes, centers = jenk_count)
    #boxplot(rescaledBranchTimes~clusters$cluster, horizontal = T)
    time_jenks <- c(0,clusters$centers[order(clusters$centers)],1)
    jenk_count <- length(time_jenks-1)
    print(jenk_count)
  }
  
  print(time_jenks[1:(jenk_count-1)])
  print(time_jenks[2:(jenk_count)])
  jenk_df <- data.frame(start = time_jenks[1:(jenk_count-1)], end = time_jenks[2:(jenk_count)], MDI = rep(NA, jenk_count-1), pVal = rep(NA, jenk_count-1))
  print(jenk_df)
  for(i in 1:(length(time_jenks)-1))
  {
    this_jenk <- geiger::dtt(phy = matchedTree_DTT$phy, data = na.omit(matchedTree_DTT$data), nsim = nSim, CI = .95, calculateMDIp = T, mdi.range = c(time_jenks[i], time_jenks[i+1]), plot = F)
    jenk_df[i,"MDI"] <- this_jenk$MDI
    #2 tailed correction
    firstP <- getMDIp(this_jenk)*2
    #dunn-sidak correction, add this to plotting method
    #jenk_df[i,"pVal"] <- 1-((1-(firstP))^(1/jenk_count))
    jenk_df[i,"pVal"] <- firstP
    print(getMDIp(this_jenk))
  }
  
  dtt_calculated$jenkStats <- jenk_df
  title(main = paste(traitName," p = ", dtt_calculated$MDIpVal, sep = " "))
  return(dtt_calculated)
}

makeDTT_parrotData_Guillerme <- function(traitName = NULL, nSim = 10, jenk_count = 1)
{
  matchedTree_DTT <- clean.data(tree = allMatched$phy, data = (allMatched$data))
  matchedTree_DTT$tree$root.time <- 37.455
  guiTestDTT <-  dispRity::dispRity.through.time(tree = matchedTree_DTT$tree, data = matchedTree_DTT$data["PC1"], time = 10)

  dtt_calculated$jenkStats <- jenk_df
  title(main = paste(traitName," p = ", dtt_calculated$MDIpVal, sep = " "))
  return(dtt_calculated)
}

#function to take the table of sig. values and plot them as colored rectangles
#takes a dtt out object, a color for the rectangle, width, where on the yaxis to place the rect. and newplot
plotSigJenks <- function(dttOut = NULL, colSig = "red", width = 1, placement = 0, newplot = F)
{
  if(newplot){plot.new()}
  jenkTab <- dttOut$jenkStats
  jenkNum <- length(jenkTab[,1])
  #sigThresh_DunnSidak <- 1-((1-(.05))^(1/jenkNum))
  sigThresh_DunnSidak <- .05
  print(paste("Significant Threshold = ", sigThresh_DunnSidak))
  
  for(i in 1:jenkNum)
  {
    if(!is.na(jenkTab[i, "pVal"]))
    {
    #dunn-sidak correction
    firstP <- p.adjust(jenkTab[i, "pVal"], method = "fdr") < sigThresh_DunnSidak
    if(firstP){
    try(rect(xleft = jenkTab[i,"start"], xright = jenkTab[i,"end"], ybottom = placement-(width), ytop = placement, col = colSig, border = NA))
    }
    }
  }
}

#conversions of significance windows 
hue_SigWindow <- scales::rescale(x = c(0.7783615, 0.8781289), from = c(0,1), to = c(37.45, 0))
saturation_SigWindow <- scales::rescale(x = c(0.9479673, 1.0000000), from = c(0,1), to = c(37.45, 0))

names(allMatched$data)

#pdf("FigureOutput/morphoPC_DTT_plots_7.17.23.pdf", width = 15, height = 15, bg = "transparent")
#par(mfrow = c(2,1))
pickSims <- 1000
#4 conssitently picked as best jenk number
pickJenks <- 4

#run DTT simulation tests for morphology
PC1_dtt <- makeDTT_parrotData(traitName = "PC1", nSim = pickSims, jenk_count = pickJenks)
PC2_dtt <-makeDTT_parrotData(traitName = "PC2", nSim = pickSims, jenk_count = pickJenks)
#graphics.off()

#pdf("FigureOutput/plumage_DTT_plots_7.17.23.pdf",  width = 5, height = 10, bg = "transparent")
#par(mfrow = c(2,1))
HueDisp_dtt <- makeDTT_parrotData(traitName = "huedisp.m", nSim = pickSims, jenk_count = pickJenks)
MaxRA_dtt <- makeDTT_parrotData(traitName = "mean.ra", nSim = pickSims, jenk_count = pickJenks)

#graphics.off()
#uncomment ^^ depending upon your preferred plotting options.
climPC1_dtt <- makeDTT_parrotData(traitName = "XCentroid0.5", nSim = pickSims, jenk_count = pickJenks)
climPC2_dtt <- makeDTT_parrotData(traitName = "yCentroid0.5", nSim = pickSims, jenk_count = pickJenks)

#make combined plot of DTT for different traits.
pdf("compositeColored_DTT_plots_8.20.25.pdf", width = 5, height = 5, bg = "transparent")
  dttLine_cols <- c("gold2", "gold4", "purple4","purple", "aquamarine2", "aquamarine4")
  setWidth = .01
  
  makeDTT_parrotData(traitName = "mean.ra", nSim = 5, rescale = T, jenk_count = 1)

  lines(PC1_dtt$dtt~PC1_dtt$times, col = "gold2", lwd = 3)
  lines(PC2_dtt$dtt~PC2_dtt$times, col = "gold4", lwd = 3)
  lines(HueDisp_dtt$dtt~HueDisp_dtt$times, col = "purple4", lwd = 3)
  lines(MaxRA_dtt$dtt~MaxRA_dtt$times, col = "purple", lwd = 3)
  lines(climPC1_dtt$dtt~climPC1_dtt$times, col = "aquamarine2", lwd = 3)
  lines(climPC2_dtt$dtt~climPC2_dtt$times, col = "aquamarine4", lwd = 3)
  
  setWidth = -.05
  placeStart = 0
  plotSigJenks(PC1_dtt, col = "gold2", placement = placeStart, width = setWidth)
  plotSigJenks(PC2_dtt, col = "gold4", placement = placeStart-(setWidth), width = setWidth)
  plotSigJenks(HueDisp_dtt, col = "purple4", placement = placeStart-(setWidth*2), width = setWidth)
  plotSigJenks(MaxRA_dtt, col = "purple", placement = placeStart-(setWidth*3), width = setWidth)
  plotSigJenks(climPC1_dtt, col = "aquamarine2", placement = placeStart-(setWidth*4), width = setWidth)
  plotSigJenks(climPC2_dtt, col = "aquamarine4", placement = placeStart-(setWidth*5), width = setWidth)

  legend("topright", legend = c("Body Size", "Body Shape", "Hue Disparity", "Saturation", "Climate PC1", "Climate PC2"),
       fill = c(dttLine_cols), border = c(dttLine_cols), bty = 'n')

graphics.off()

system("open compositeColored_DTT_plots_8.20.25.pdf")

combined_jenk_stats <- cbind.data.frame(PC1_dtt$jenkStats, PC2_dtt$jenkStats[,c(3,4)], HueDisp_dtt$jenkStats[,c(3,4)], MaxRA_dtt$jenkStats[,c(3,4)], climPC1_dtt$jenkStats[,c(3,4)], climPC2_dtt$jenkStats[,c(3,4)])
colnames(combined_jenk_stats) <- c("Time Start", "Time End", "Body Size MDI", "Body Size p", "Body Shape MDI", "Body Shape p", "Hue Disparity MDI", "Hue Disparity p", "Saturation MDI", "Saturation p", "Climate PC1 MDI", "Climate PC1 p", "Climate PC2 MDI", "Climate PC2 p")

#output table for supplement
PVALUE_threshold <- 0.0085
dataTab <- (combined_jenk_stats)
dataTab |>
  gt() |>
  tab_header(
    title = paste("Disparity through time window statistics"),
    subtitle = paste("Dunn-sidak corrected p threshold =", PVALUE_threshold)
  ) -> combined_jenk_stats_plot





