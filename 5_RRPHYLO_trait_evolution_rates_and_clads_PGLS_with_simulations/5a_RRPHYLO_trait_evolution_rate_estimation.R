##script to estimate rates of trait evolution with RRphylo

#load dependencies, datasets
{
#utilities
library(devtools)

#phylo tools
library(phytools)
library(motmot)
library(ape)
library(geiger)
library(RRphylo)
library(ggpmisc)
library(caper)
library(stringdist)
  
#plotting tools
library(ggplot2)
library(gridExtra)

#data management
library(tidyr)
library(dplyr)

#for making tables
library(gt)

library(stringr)
  
#useful method for matching later 
lastDropper <- function(string = NULL, seper = " ")
{
  require("stringr")
  lab <- word(string, 1, -2, sep = seper)
  return(lab)
}

#set working directory to parent directory
#data_pointer <- <YOUR PATH>
data_pointer <- "c:/Users/jtm378/AMNH Dropbox/Jon Merwin/AllParrotHomeAnalysis/Scripts/GITHUB_UPLOAD/"
#data_pointer <- "c:/Users/jonme/AMNH Dropbox/Jon Merwin/AllParrotHomeAnalysis/Scripts/GITHUB_UPLOAD/"
setwd(data_pointer)

#match tree and data
matchForCont <- readRDS("data/tipmatched_parrotdata_obj_with_phyPCA.rds")
View(matchForCont$data)

#load safe_merge function from before 
# Function to normalize species names (convert underscores to spaces, clean up)
normalize_species_name <- function(name) {
  name <- gsub("_", " ", name)  # Convert underscores to spaces
  name <- gsub("\\s+", " ", name)  # Remove multiple spaces
  name <- trimws(name)  # Remove leading/trailing whitespace
  return(name)
}

# Function for fuzzy matching species names
fuzzy_match_species <- function(target_names, reference_names, max_dist = 1) {
  # Normalize all names first
  target_norm <- normalize_species_name(target_names)
  reference_norm <- normalize_species_name(reference_names)
  
  matches <- character(length(target_names))
  
  for (i in seq_along(target_norm)) {
    # First try exact match on normalized names
    exact_match <- which(reference_norm == target_norm[i])
    
    if (length(exact_match) > 0) {
      matches[i] <- reference_names[exact_match[1]]
    } else {
      # Try fuzzy matching with string distance
      distances <- stringdist(target_norm[i], reference_norm, method = "lv")
      min_dist <- min(distances)
      
      if (min_dist <= max_dist) {
        best_match <- which.min(distances)[1]
        matches[i] <- reference_names[best_match]
        cat("Fuzzy match:", target_names[i], "->", reference_names[best_match], 
            "(distance:", min_dist, ")\n")
      } else {
        matches[i] <- NA
        cat("No match found for:", target_names[i], "\n")
      }
    }
  }
  
  return(matches)
}

# Function to safely merge datasets with fuzzy matching
safe_merge <- function(base_df, new_df, by_rownames = TRUE, max_dist = 2, append = "_new") {
  if (by_rownames) {
    # Get species names
    base_species <- rownames(base_df)
    new_species <- rownames(new_df)
    
    # Perform fuzzy matching
    matches <- fuzzy_match_species(base_species, new_species, max_dist = max_dist)
    
    # Create mapping of matched species
    valid_matches <- !is.na(matches)
    if (sum(valid_matches) > 0) {
      base_matched <- base_species[valid_matches]
      new_matched <- matches[valid_matches]
      
      # Add columns from new_df, excluding species identifier columns
      exclude_cols <- c("spp", "SPECIES", "Sp", "DataNames")
      new_cols <- new_df[new_matched, !names(new_df) %in% exclude_cols, drop = FALSE]
      
      # Make sure column names don't conflict
      existing_cols <- intersect(names(new_cols), names(base_df))
      if (length(existing_cols) > 0) {
        cat("Warning: Overlapping column names found:", paste(existing_cols, collapse = ", "), "\n")
        # Rename overlapping columns
        names(new_cols)[names(new_cols) %in% existing_cols] <- paste0(names(new_cols)[names(new_cols) %in% existing_cols], append)
      }
      
      base_df[base_matched, names(new_cols)] <- new_cols
    }
    
    cat("Matched", sum(valid_matches), "out of", length(base_species), "species\n\n")
  }
  return(base_df)
}
}

#Trait evolution rates with MOTMOT, generate summary stats
{
  
  #this transforms data into motmot format, runs the tm1 model, and saves the resulting plot PDF
  motmot_run <- function(tree = matchForCont$phy, data = matchForCont$data, trait_column = "PC1", plot_medusa = T)
  {
    trait_vector <- (as.numeric(data[,trait_column]))
    names(trait_vector) <- (data$species)
    motmot_matched <- picante::match.phylo.data(tree, unlist(as.list(na.omit(trait_vector))))
    motmot_matched$data_matrix <- as.matrix(motmot_matched$data)
    rownames(motmot_matched$data_matrix) <- names(motmot_matched$data)
    # This tests every internal node for potential rate shifts, do for color and morph
    #morph rates
    
    medusa_shifts <- transformPhylo.ML(phy = motmot_matched$phy, y = (motmot_matched$data_matrix),
                                       model = "tm1", minCladeSize = length(motmot_matched$phy$tip.label)/20)
    if(plot_medusa)
    {
      plot_name <- paste("medusa_shifts_", trait_column, ".pdf", sep = "")
      #plot the medusa shifts for body size
      try(pdf(plot_name, width = 20, height = 60))
      motmot::plot.traitMedusa.model(summary(medusa_shifts), reconType = "rates", edge.width=3)
      graphics.off()
      system(paste("open", plot_name))
    }
    return(medusa_shifts)
  }
  
  motmot_run_results <- list()
  
  motmot_run_results$PC1 <- motmot_run(trait_column = "PC1")  
  motmot_run_results$PC2 <- motmot_run(trait_column = "PC2")
  
  motmot_run_results$huedisp.m <- motmot_run(trait_column = "huedisp.m")
  motmot_run_results$mean.ra <- motmot_run(trait_column = "mean.ra")
  
  motmot_run_results$XCentroid0.5 <- motmot_run(trait_column = "XCentroid0.5")
  motmot_run_results$yCentroid0.5 <- motmot_run(trait_column = "yCentroid0.5")
  
  motmot_run_results$PC1_climate <- motmot_run(trait_column = "PC1_climate")
  motmot_run_results$PC2_climate <- motmot_run(trait_column = "PC2_climate")
  
  #now print the resulting statistics
  summary(motmot_run_results$PC1)$Rates
  summary(motmot_run_results$PC2)$Rates
  
  summary(motmot_run_results$huedisp.m)$Rates
  summary(motmot_run_results$mean.ra)$Rates
  
  summary(motmot_run_results$XCentroid0.5)$Rates
  summary(motmot_run_results$yCentroid0.5)$Rates
  
  summary(motmot_run_results$PC1_climate)$Rates
  summary(motmot_run_results$PC2_climate)$Rates
  
}

#prep data for RRphylo (convert to numeric, check branch lengths)
{
#check for edge lengths 0 after matching (0 branch lengths can cause issues)
which(matchForCont$phy$edge.length <= 0)

thisDatForRR <- apply(matchForCont$data, MARGIN = 2, function(X) as.numeric(X))
#thisDatForRR <- data.frame(scale(matchForCont$data, scale = T, center = 5))
thisDfForRR <- data.frame(thisDatForRR, row.names = rownames(thisDatForRR))
thisDfForRR$Continent <- matchForCont$data$Continent
thisDfForRR$Biogeographic.Region <- matchForCont$data$Biogeographic.Region
rownames(thisDfForRR) <- rownames(matchForCont$data)
sampFract <- length(thisDfForRR[,1]) / 385
print(sampFract)
}

##now get trait evolution rates using RRphylo

#esimtate RRphylo rates with morphometric data 
{
  cladPhy <- matchForCont$phy
  which(cladPhy$edge.length <= 0)
  vars <- c("PC1", "PC2","PC3")
  vars_raw <- c("BILL.LENGTH", "BILL.WIDTH","BILL.HEIGHT","TARSUS.LENGTH","WING.LENGTH.TO.PRIMARIES","WING.LENGTH.TO.SECONDARIES","TAIL")
  #vars <- c("BILL.LENGTH", "BILL.WIDTH","BILL.HEIGHT","TARSUS.LENGTH","WING.LENGTH.TO.PRIMARIES","TAIL")
  
  #vars <- c("PC1")
  
  firstcov <- picante::match.phylo.data(phy = (cladPhy), data = na.omit(thisDfForRR[vars]))
  
  which(firstcov$phy$edge.length <= 0)
  
  rrInputDat <- data.frame(scale(firstcov$data,center = F,scale = T))
  
  hist(rrInputDat[,1])
  hist(rrInputDat[,2])
  hist(rrInputDat[,3])
  
  # corrplot::corrplot(cor(t(rrInputDat_bodyscaled)))
  
  
  # rrInputDat <- log(2.1-firstcov$data)
  rrInputDat$PC1_trans <- (log(abs(min(rrInputDat$PC1))-rrInputDat$PC1))
  hist(rrInputDat$PC1_trans)
  # par(mfrow = c(2,2))
  # hist((rrInputDat[,"PC1"]))
  # hist(rrInputDat[,"PC2"])
  # hist(rrInputDat[,"PC3"]) 
  # plot(rrInputDat[,2:1])
  
  firstCovMorph <- RRphylo::RRphylo(tree = (firstcov$phy), y = setNames(rrInputDat[,"PC1_trans"], nm = rownames(rrInputDat)))
  firstCovMorph <- RRphylo::RRphylo(tree = (firstcov$phy), y = rrInputDat[,c("PC1_trans","PC2","PC3")])
  search.shift(RR=firstCovMorph,status.type="clade", f = length(firstCovMorph$tree$tip.label)/20)->shiftSearch
  # Visualize search.shift results 
  shiftSearch$single.clades
  plotShift(RR=(firstCovMorph),SS=shiftSearch)->plotSS
  plotSS$plotClades()
  
  plot(log(as.numeric(firstCovMorph$rates))~as.numeric(c(firstCovMorph$aces[,1],rrInputDat[,"PC1_trans"])))
  test_mod <- lm(log(as.numeric(firstCovMorph$rates))~as.numeric(c(firstCovMorph$aces[,3],rrInputDat[,"PC3"])))
  summary(test_mod)
  
  #run morphology rates using body size as a covariate
  covDat <- c(firstCovMorph$aces[,1], setNames(rrInputDat[,"PC1_trans"], nm = rownames(rrInputDat)))
  #covDat <- rowMeans(rbind.data.frame(firstCovMorph$aces[,1], rrInputDat[,"PC1_trans"]))
  #covDat <- apply(X= covDat, MARGIN = 1, FUN = function(X) sd(X), simplify = T)
  
  View(firstcov$data)
  
  #withCov <- RRphylo(tree = firstcov$phy, y = firstcov$data, cov = (covDat))
  withCov <- RRphylo(tree = firstcov$phy, y = rrInputDat[,c("PC1_trans","PC2","PC3")])
  
  plot(log(as.numeric(withCov$rates[,1]))~(as.numeric(c(withCov$aces[,1], rrInputDat[,1]))))
  test_mod <- lm(log(as.numeric(withCov$rates[,1]))~(as.numeric(c(withCov$aces[,1],rrInputDat[,"PC1_trans"])))^2)
  summary(test_mod)
  
  logRatesWithCov <- log(withCov$rates) 
  limsMorphCov <- (c(mean(logRatesWithCov)-2*sd(logRatesWithCov), mean(logRatesWithCov)+2*sd(logRatesWithCov)))
  
  withCov$rates
  
  png("pc_morph_rrphylo_nocovariance.png", height = 4000, width = 1000, units = "px")
  plot(firstCovMorph$tree)
  addShift(search.shift(RR=firstCovMorph,status.type="clade", f = length(firstCovMorph$tree$tip.label)/20))
  #addShift(search.shift(RR=withCov,status.type="clade", f = 1))
  graphics.off()
  system("open pc_morph_rrphylo_nocovariance.png")
  
  png("pc_morph_rrphylo_withcovariance.png", height = 4000, width = 1000, units = "px")
  plot(withCov$tree)
  addShift(search.shift(RR=withCov,status.type="clade", f = length(withCov$tree$tip.label)/20))
  #addShift(search.shift(RR=withCov,status.type="clade", f = 1))
  graphics.off()
  system("open pc_morph_rrphylo_withcovariance.png")
  
  
  
  #make plots of trait evolution rates
  #png(file = here::here("FigureOutput/withCovMorphPCsRRphylo.png"), width = 5, height = 5, units = "in", res = 300, bg = "transparent")
  #contMap(method = "user", tree = withCov$tree, x = (withCov$rates[-c(1:365),1]), anc.states  = (withCov$rates[c(1:365),1]),lwd = 3, outline = F, fsize = .2, type = "fan", ftype = "off")
  #dev.off()
  
  #png(file = here::here("FigureOutput/withCovMorphRates2SDLimsNoLimsNOLab.png"), width = 5, height = 5, units = "in", res = 300, bg = "transparent")
  #contMap(method = "user", tree = withCov$tree, x = log(withCov$rates[-c(1:length(withCov$tree$tip.label)-1),1]), anc.states  = log(withCov$rates[c(1:363),1]),lwd = 3, outline = F, fsize = .2, type = "fan", ftype = "off", lims = limsMorphCov, legend = F)
  #dev.off() 
  
  #not scaled for extreme values
  png(file = ("firstCovRateShifts.png"), width = 20, height = 40, units = "in", res = 300, bg = "white")
  search.shift(RR=firstCovMorph,status.type="clade", f = length(withCov$tree$tip.label)/20)->shiftSearch
  # Visualize search.shift results 
  plot(firstCovMorph$tree, cex = .5)
  addShift(SS=shiftSearch)
  #plotShift(RR=firstCovMorph,SS=shiftSearch)->plotSS
  #plotSS$plotClades()
  graphics.off()
  system("open firstCovRateShifts.png")
  
  #scaled for extreme values
  #png(file = ("WithCovRateShifts.png"), width = 20, height = 40, units = "in", res = 300, bg = "white")
  pdf(file = ("WithCovRateShifts.pdf"), width = 10, height = 10, bg = "white")
  search.shift(RR=withCov,status.type="clade", f = length(withCov$tree$tip.label)/20)->shiftSearch_cov
  # Visualize search.shift results 
  plot(firstCovMorph$tree, cex = .1, edge.width = 3, edge.color = "gray20")
  addShift(SS=shiftSearch_cov)
  #plotSS$plotClades()
  graphics.off()
  system("open WithCovRateShifts.pdf")
}

##esimtate RRphylo rates with color data 
{
  cladPhy <- matchForCont$phy
  which(cladPhy$edge.length <= 0)
  vars <- c("huedisp.m", "mean.ra","colspan.m")
  #vars <- c("centroid.m","centroid.s","centroid.l","centroid.u","huedisp.m", "mean.ra", "a.vol","colspan.m")
  #vars <- c("centroid.l","centroid.u")
  firstcov <- picante::match.phylo.data(phy = drop.tip(cladPhy, tip = "Neophema_splendida"), data = data.frame(scale(na.omit(thisDfForRR[vars]),center = F,scale = T)))
  which(firstcov$phy$edge.length <= 0)
  rrInputDat <- firstcov$data[vars]
  View(rrInputDat)
  hist((rrInputDat[,1]))
  hist(rrInputDat[,2])
  rrInputDat$Average.of.huedisp.m <- sqrt(as.numeric(rrInputDat$huedisp.m))
  par(mfrow = c(2,2))
  hist((firstcov$data[,"huedisp.m"]))
  hist(firstcov$data[,"mean.ra"])
  #rrInputDat$a.vol <- log(firstcov$data$a.vol)
  #hist((firstcov$data[,"a.vol"]))
  plot(rrInputDat[,2:1])
  
  firstCovCol <- RRphylo::RRphylo(tree = (firstcov$phy), y = rrInputDat)
  search.shift(RR=firstCovCol,status.type="clade", f = length(firstCovCol$tree$tip.label)/20)->shiftSearchCol
  # Visualize search.shift results 
  plotShift(RR=(firstCovCol),SS=shiftSearchCol)->plotSS
  #png(filename = here::here("noCovColorRateShiftsRR.png"), width = 20, height = 40, units = "in", res = 300)
  plotSS$plotClades(symbols.args=list(fg=c(pos="blue",neg="red")))
  #graphics.off()
  
  covColDat <- (rbind.data.frame(firstCovCol$aces, rrInputDat))
  
  #removing cov because it doesn't make biological sense...
  withCovCol <- RRphylo::RRphylo(tree = (firstcov$phy), y = rrInputDat)
  search.shift(RR=withCovCol,status.type="clade", f = length(withCovCol$tree$tip.label)/20)->shiftSearchCol
  #search.shift(RR=withCovCol,status.type="clade", f = 1)->shiftSearchCol
  
  # Visualize search.shift results 
  #plotShift(RR=(withCovCol),SS=shiftSearchCol)->plotSS
  
  #png(filename = here::here("withCovColorRateShiftsRR.png"), width = 20, height = 40, units = "in", res = 300)
  pdf(file = ("withCovColorRateShiftsRR.pdf"), width = 10, height = 10, bg = "white")
  # Visualize search.shift results 
  plot(withCovCol$tree, cex = .1, edge.width = 3, edge.color = "gray20")
  addShift(SS=shiftSearchCol)
  #plotSS$plotClades()
  graphics.off()
  colRRLen <- length(withCovCol$tree$tip.label)
  system("open withCovColorRateShiftsRR.pdf")
}

##now the raw color analysis...this does not seem to make a huge difference!
{
  cladPhy <- matchForCont$phy
  which(cladPhy$edge.length <= 0)
  raw_col_vars <- c(
    grep(x = colnames(matchForCont$data), pattern = c("huedisp.m_")),
    grep(x = colnames(matchForCont$data), pattern = c("colspan.m_")),
    grep(x = colnames(matchForCont$data), pattern = c("mean.ra_")))
  
  
  
  
  
  raw_col_vars <- c(
    grep(x = colnames(matchForCont$data), pattern = c("BRIGHT_")))
  
  
  firstcov <- picante::match.phylo.data(phy = drop.tip(cladPhy, tip = "Neophema_splendida"), data = data.frame(scale(na.omit(thisDfForRR[vars]),center = F,scale = T)))
  
  
  raw_dat_numeric <- apply(data.frame((matchForCont$data[raw_col_vars])),2,as.numeric)
  rownames(raw_dat_numeric) <- rownames(matchForCont$data)
  raw_dat_numeric <- data.frame(na.omit(raw_dat_numeric))
  
  firstcov <- picante::match.phylo.data(phy = drop.tip(matchForCont$phy, tip = "Neophema_splendida"), data =raw_dat_numeric)
  
  firstCovCol_raw <- RRphylo::RRphylo(tree = (firstcov$phy), y = firstcov$data)
  search.shift(RR=firstCovCol_raw,status.type="clade", f = length(firstCovCol_raw$tree$tip.label)/20)->shiftSearchCol_raw
  # Visualize search.shift results 
  plotShift(RR=(firstCovCol_raw),SS=shiftSearchCol_raw)->plotSS
  #png(filename = here::here("noCovColorRateShiftsRR.png"), width = 20, height = 40, units = "in", res = 300)
  plotSS$plotClades(symbols.args=list(fg=c(pos="blue",neg="red")))
  #graphics.off()
  
  covColDat <- (rbind.data.frame(firstCovCol$aces, rrInputDat))
  
  #removing cov because it doesn't make biological sense...
  withCovCol <- RRphylo::RRphylo(tree = (firstcov$phy), y = rrInputDat)
  search.shift(RR=withCovCol,status.type="clade", f = length(withCovCol$tree$tip.label)/20)->shiftSearchCol
  #search.shift(RR=withCovCol,status.type="clade", f = 1)->shiftSearchCol
  
  # Visualize search.shift results 
  #plotShift(RR=(withCovCol),SS=shiftSearchCol)->plotSS
  
  #png(filename = here::here("withCovColorRateShiftsRR.png"), width = 20, height = 40, units = "in", res = 300)
  pdf(file = ("withCovColorRateShiftsRR.pdf"), width = 10, height = 10, bg = "white")
  # Visualize search.shift results 
  plot(withCovCol$tree, cex = .1, edge.width = 3, edge.color = "gray20")
  addShift(SS=shiftSearchCol)
  #plotSS$plotClades()
  graphics.off()
  colRRLen <- length(withCovCol$tree$tip.label)
  system("open withCovColorRateShiftsRR.pdf")
}

##now the climate analysis
{
  cladPhy <- matchForCont$phy
  which(cladPhy$edge.length <= 0)
  vars <- c("Total.niche","Niche.over.50.","XCentroid0.0","yCentroid0.0","XCentroid0.5", "yCentroid0.5")
  #vars <- c("PC1_climate", "PC2_climate")
  firstcov <- picante::match.phylo.data(phy = (cladPhy), data = data.frame(scale(center = T, scale = T, na.omit(thisDfForRR[vars]))))
  which(firstcov$phy$edge.length <= 0)
  rrInputDat <- firstcov$data[,vars]

  firstCovClim <- RRphylo::RRphylo(tree = (firstcov$phy), y = rrInputDat)
  search.shift(RR=firstCovClim,status.type="clade",f = 1)->shiftSearchClim
  # Visualize search.shift results 
  plotShift(RR=(firstCovClim),SS=shiftSearchClim)->plotSS
  #png(filename = here::here("noCovColorRateShiftsRR.png"), width = 20, height = 40, units = "in", res = 300)
  plotSS$plotClades()
  #graphics.off()
  
  covclimDat <- rowMeans(rbind.data.frame(firstCovClim$aces, rrInputDat))
  withCovClim <- RRphylo::RRphylo(tree = (firstcov$phy), y = rrInputDat)
  
  all(is.na(rrInputDat))
  all(is.na(covclimDat))
  
  search.shift(RR=withCovClim,status.type="clade", f = length(withCovClim$tree$tip.label)/20)->shiftSearchClim
  #search.shift(RR=withCovClim,status.type="clade", f = 1)->shiftSearchClim
  
  # Visualize search.shift results 
  #plotShift(RR=(withCovClim),SS=shiftSearchClim)->plotSS
  #png(filename = here::here("withCovClimRateShiftsRR.png"), width = 20, height = 40, units = "in", res = 300)
  pdf(file = ("withCovClimRateShiftsRR.pdf"), width = 10, height = 10, bg = "white")
  #plotSS$plotClades()
  plot(withCovClim$tree, cex = .1, edge.width = 3, edge.color = "gray20")
  addShift(SS=shiftSearchClim)
  graphics.off()
  colRRLen <- length(withCovClim$tree$tip.label)
  system("open withCovClimRateShiftsRR.pdf")
}

#additional RRphylo hypothesis test (auto-test is not full hyp. test)
#this generates the statistics in the manuscript table
{
shiftSearch_morph_fulltree <- search.shift(RR=withCov,status.type="clade", nrep = 10000, f = length(withCov$tree$tip.label)/20)
shiftSearch_morph_hypTest <- search.shift(RR=withCov,status.type="clade", node = c(as.numeric(rownames(shiftSearch_morph_fulltree$single.clades))), nrep = 10000, f = length(withCov$tree$tip.label)/20)
#shiftSearch_morph_hypTest <- search.shift(RR=withCov,status.type="clade", node = 470, nrep = 100000)

add_descendent_column <- function(tree = NA, shift_nodes = NA)
{
  print(paste("Getting descendents for:", length(shift_nodes), "nodes"))
  res_list <- c(rep("",length(shift_nodes)))
  for(i in 1:length(shift_nodes))
  {
  print(i)
  try(these_taxa <- na.omit(tree$tip.label[getDescendants(tree, node = shift_nodes[i])])) #decrease,lorikeets
  (res_list[i]<-paste(unique(str_split_fixed(c(these_taxa), pattern = "_", n = 2)[,1]), collapse = ";"))
  }
  print(res_list)
  if(length(res_list)==length(shift_nodes)){return(res_list)}else{
  print("Error! Input and output length do not match.")}
}

#print descendents
shiftSearch_morph_hypTest$single.clades$singles
shift_nodes_morph <-as.numeric(rownames(shiftSearch_morph_hypTest$single.clades$singles))
# na.omit(withCov$tree$tip.label[getDescendants(withCov$tree, node = shift_nodes[1])]) #decrease,lorikeets
# na.omit(withCov$tree$tip.label[getDescendants(withCov$tree, node = shift_nodes[2])]) #decrease, amazona
# na.omit(withCov$tree$tip.label[getDescendants(withCov$tree, node = shift_nodes[3])]) #increase, Macaws
shiftSearch_morph_hypTest_table <- data.frame(shiftSearch_morph_hypTest$single.clades$singles)
shiftSearch_morph_hypTest_table$trait <- "Body Shape / Size"
shiftSearch_morph_hypTest_table$taxa <- add_descendent_column(tree = withCov$tree, shift_nodes = shift_nodes_morph)

search_clade_size = length(withCovCol$tree$tip.label)/20
shiftSearch_col_fulltree <- search.shift(RR=withCovCol,status.type="clade", nrep = 10000,f=search_clade_size)
col_shift_nodes <- as.numeric(rownames(shiftSearch_col_fulltree$single.clades))
shiftSearch_col_hypTest <- search.shift(RR=withCovCol,status.type="clade", f=search_clade_size, node = col_shift_nodes, nrep = 10000)
shiftSearch_col_hypTest_table <- data.frame(shiftSearch_col_hypTest$single.clades$singles)
shiftSearch_col_hypTest_table$trait <- "Plumage Color"
col_shift_nodes <- as.numeric(rownames(shiftSearch_col_hypTest_table))
shiftSearch_col_hypTest_table$taxa <- add_descendent_column(withCovCol$tree,shift_nodes = col_shift_nodes)


# withCovCol$tree$tip.label[getDescendants(withCovCol$tree, node = col_shift_nodes[1])] #cockatoos, decrease
# withCovCol$tree$tip.label[getDescendants(withCovCol$tree, node = col_shift_nodes[2])] #prioniturus, psittacula, decrease
# withCovCol$tree$tip.label[getDescendants(withCovCol$tree, node = col_shift_nodes[3])] #lorikeets trichoglossus, eos, increase
#withCovCol$tree$tip.label[getDescendants(withCovCol$tree, node = col_shift_nodes[4])] #arinae, macaws, pyrrhura, increase

#now climate 
search_clade_size = length(withCovClim$tree$tip.label)/20
shiftSearch_clim_fulltree <- search.shift(RR=withCovClim,status.type="clade", nrep = 10000,f=search_clade_size)
shift_nodes_clim <- c(as.numeric(rownames(shiftSearch_clim_fulltree$single.clades)))
shiftSearch_clim_hypTest <- search.shift(RR=withCovClim,status.type="clade", f=search_clade_size, node = shift_nodes_clim, nrep = 10000)
shiftSearch_clim_hypTest_table <- data.frame(shiftSearch_clim_hypTest$single.clades)
shiftSearch_clim_hypTest_table$trait <- "Climatic Niche"
colnames(shiftSearch_clim_hypTest_table)[1:2] <- c("rate.difference","p.value")
shiftSearch_clim_hypTest_table$taxa <- add_descendent_column(tree = withCovClim$tree, shift_nodes = shift_nodes_clim)

RRphylo_results_table <- rbind.data.frame(shiftSearch_col_hypTest_table, shiftSearch_morph_hypTest_table, shiftSearch_clim_hypTest_table)
RRphylo_results_table$node <- rownames(RRphylo_results_table)
min_dist <- pmin(abs(RRphylo_results_table$p.value - 0), abs(RRphylo_results_table$p.value - 1))
RRphylo_results_table_reordered <- RRphylo_results_table[c("taxa", "trait", "rate.difference", "p.value")]
RRphylo_results_table_reordered$corrected_p <- round(p.adjust(min_dist, "fdr"),digits = 4)

#save these results as table
gt::gt(RRphylo_results_table_reordered)

#dotchart(shiftSearch$all.clades$rate.difference[order(shiftSearch$all.clades$rate.difference)], pch = 16, col = as.factor(shiftSearch$all.clades$p.value[order(shiftSearch$all.clades$rate.difference)] > .975))
#dotchart(shiftSearchCol$all.clades$rate.difference, pch = 16, labels = rownames(shiftSearch$all.clades))
}

#get EVO rates with RRphylo for ventral and dorsal coloration
{
  windows()
  pairs(data.frame(lapply(na.omit(matchForCont$data[,c(3:16,19:20)]), as.numeric)))
  cor_test <- cor(data.frame(lapply(na.omit(matchForCont$data[,c(3:16,19:20)]), as.numeric)))
  corrplot::corrplot(cor_test)
  
  #first ventral data
  {
  cladPhy <- matchForCont$phy
  which(cladPhy$edge.length <= 0)
  vars <- c("huedisp.m_ventral", "mean.ra_ventral", "a.vol_ventral")
  firstcov_ventral <- picante::match.phylo.data(phy = (cladPhy), data = na.omit(thisDfForRR[vars]))
  which(firstcov_ventral$phy$edge.length <= 0)
  rrInputDat <- firstcov_ventral$data[vars]
  hist((rrInputDat[,1]))
  hist(rrInputDat[,2])
  rrInputDat$huedisp.m <- sqrt(as.numeric(rrInputDat$huedisp.m))
  par(mfrow = c(2,2))
  hist((firstcov_ventral$data[,"huedisp.m_ventral"]))
  hist(firstcov_ventral$data[,"mean.ra_ventral"])
  firstcov_ventral$data$a.vol <- log(1+firstcov_ventral$data$a.vol)
  hist((firstcov_ventral$data[,"a.vol"]))
  plot(rrInputDat[,2:1])
  
  firstCovCol_ventral <- RRphylo::RRphylo(tree = (firstcov_ventral$phy), y = rrInputDat)
  search.shift(RR=firstCovCol_ventral,status.type="clade")->shiftSearchCol_ventral
  # Visualize search.shift results 
  plotShift(RR=(firstCovCol_ventral),SS=shiftSearchCol_ventral)->plotSS
  #png(filename = here::here("noCovColorRateShiftsRR.png"), width = 20, height = 40, units = "in", res = 300)
  plotSS$plotClades()
  #graphics.off()
  
  covColDat_ventral <- c(firstCovCol_ventral$aces[,1], firstCovCol_ventral$predicted.phenotype[,1])
  withCovCol_ventral <- RRphylo::RRphylo(tree = (firstCovCol_ventral$tree), y = rrInputDat, x1 = covColDat_ventral)
  }
  
  #now dorsal data
  {
    cladPhy <- matchForCont$phy
    which(cladPhy$edge.length <= 0)
    vars <- c("huedisp.m_dorsal", "mean.ra_dorsal", "a.vol_dorsal")
    firstcov_dorsal <- picante::match.phylo.data(phy = (cladPhy), data = na.omit(thisDfForRR[vars]))
    which(firstcov_dorsal$phy$edge.length <= 0)
    rrInputDat <- firstcov_dorsal$data[vars]
    hist((rrInputDat[,1]))
    hist(rrInputDat[,2])
    rrInputDat$huedisp.m <- sqrt(as.numeric(rrInputDat$huedisp.m))
    par(mfrow = c(2,2))
    hist((firstcov_dorsal$data[,"huedisp.m_dorsal"]))
    hist(firstcov_dorsal$data[,"mean.ra_dorsal"])
    firstcov_dorsal$data$a.vol <- log(1+firstcov_dorsal$data$a.vol)
    hist((firstcov_dorsal$data[,"a.vol"]))
    plot(rrInputDat[,2:1])
    
    firstCovCol_dorsal <- RRphylo::RRphylo(tree = (firstcov_dorsal$phy), y = rrInputDat)
    search.shift(RR=firstCovCol_dorsal,status.type="clade")->shiftSearchCol_dorsal
    # Visualize search.shift results 
    plotShift(RR=(firstCovCol_dorsal),SS=shiftSearchCol_dorsal)->plotSS
    #png(filename = here::here("noCovColorRateShiftsRR.png"), width = 20, height = 40, units = "in", res = 300)
    plotSS$plotClades()
    #graphics.off()
    
    covColDat_dorsal <- c(firstCovCol_dorsal$aces[,1], firstCovCol_dorsal$predicted.phenotype[,1])
    withCovCol_dorsal <- RRphylo::RRphylo(tree = (firstCovCol_dorsal$tree), y = rrInputDat, x1 = covColDat_dorsal)
  
    plot(log(withCovCol_dorsal$rates)~log(withCovCol_ventral$rates)) 
    abline(a = 0, b = 1, col = "red")
    
    plot(log(firstCovCol_dorsal$rates)~log(firstCovCol_ventral$rates)) 
    abline(a = 0, b = 1, col = "red")
    
  }
  
}

#resave the RDS, adding in trait diversification rates from medusa, rrphylo
#this requires sourcing the safe_merge function from "collate_datasets_for_comparative.R"

#the reads the data output from the collate_datasets_for_comparative to attatch RRphylo rates 
master_dataset_tipmatched <- readRDS("data/tipmatched_parrotdata_obj_with_phyPCA.rds")

#convert RRphylo rates to data frames for matchign
morph_rates_df <- data.frame(withCov$rates)
color_rates_df <- data.frame(withCovCol$rates)
clim_rates_df <-  data.frame(withCovClim$rates)

#setd column names
colnames(morph_rates_df) <- "morph_rates"
colnames(color_rates_df) <- "color_rates"
colnames(clim_rates_df) <- "clim_rates"

#merge rates to the comparative dataset
test_merge_RR <- safe_merge(master_dataset_tipmatched$data, new_df = morph_rates_df, by_rownames = T, max_dist = 3, append = "_unscaled")
test_merge_RR <- safe_merge(test_merge_RR, new_df = color_rates_df, by_rownames = T, max_dist = 3, append = "_unscaled")
test_merge_RR <- safe_merge(test_merge_RR, new_df = clim_rates_df, by_rownames = T, max_dist = 3, append = "_unscaled")

#match comparative dataset to tree
rds_with_rates <- match.phylo.data(data = test_merge_RR,master_dataset_tipmatched$phy)

#save tree-matched object
saveRDS(rds_with_rates, file = "5_RRPHYLO_trait_evolution_rates_and_clads_PGLS_with_simulations/PGLS_RDS_with_rates.rds")

#fix a continent naming issue
rds_with_rates$data$Continent_clean <- gsub(x = rds_with_rates$data$Continent, pattern = "Oceania ", replacement = "Oceania")

#basic plots of trait rates vs. continents
ggplot(rds_with_rates$data[which(!is.na(rds_with_rates$data$Continent_clean)),],aes(x=log(as.numeric(color_rates)))) + geom_histogram(aes(fill=Continent_clean), color = "transparent", alpha = .5) + theme_classic()
ggplot(rds_with_rates$data[which(!is.na(rds_with_rates$data$Continent_clean)),],aes(x=log(as.numeric(morph_rates)))) + geom_histogram(aes(fill=Continent_clean), color = "transparent", alpha = .5) + theme_classic()
ggplot(rds_with_rates$data[which(!is.na(rds_with_rates$data$Continent_clean)),],aes(x=log(as.numeric(clim_rates)))) + geom_histogram(aes(fill=Continent_clean), color = "transparent", alpha = .5) + theme_classic()





















