#testing ratematrix

#install.packages("ratematrix")
library(ratematrix)
library(phytools)
library(picante)
library(geiger)
library(corrplot)

setwd("c:/Users/jtm378/AMNH Dropbox/Jon Merwin/AllParrotHomeAnalysis/Scripts/GITHUB_UPLOAD/")

rds_ratematrix <- readRDS("5_RRPHYLO_trait_evolution_rates_and_clads_PGLS_with_simulations/PGLS_RDS_with_rates.rds")
rds_ratematrix$phy <- root(rds_ratematrix$phy,outgroup = c("Calyptomena_viridis","Caracara_cheriway","Icterus_cucullatus"),resolve.root = T)
#ratematrix_names <- names(rds_ratematrix$data)[c(9,10,11,12,13,14,21,22,23,24,35,36,37,38,50:56,63:68)]
#ratematrix_names <- names(rds_ratematrix$data)[c(79:(79+34))]
#ratematrix_names <- names(rds_ratematrix$data)[c(114:(114+34))]
#ratematrix_names <- names(rds_ratematrix$data)[c(149:(149+34))]


dir.create("mcmcburn")

#ratematrix_names <- c("huedisp.m","mean.ra", "PC1","PC2","XCentroid0.0","yCentroid0.0")
ratematrix_names <- c("PC1", "PC2","huedisp.m","mean.ra","Total.niche","XCentroid0.0","yCentroid0.0","genus")

rds_ratematrix_matched <- match.phylo.data(rds_ratematrix$phy, na.omit(rds_ratematrix$data[,ratematrix_names]))

rds_ratematrix_matched$genus <- rds_ratematrix_matched$data$genus
rds_ratematrix_matched$num_data <- as.matrix(apply(rds_ratematrix_matched$data, 2,function(X) scales::rescale(as.numeric(X),to=c(1,5)),simplify = T))
rds_ratematrix_matched$num_data[,1] <- (rds_ratematrix_matched$num_data[,1])

rownames(rds_ratematrix_matched$num_data) <- rownames(rds_ratematrix_matched$data)
colnames(rds_ratematrix_matched$num_data) <- colnames(rds_ratematrix_matched$data)

calc_matrix <- geiger::ratematrix(rds_ratematrix_matched$phy, rds_ratematrix_matched$num_data)
corrplot((calc_matrix),is.corr = F,order = "hclust")

ratematrixMCMC(phy = anoles$phy,data = anoles$data,dir="mcmcburn",gen = 1000)

mcmc_out <- ratematrixJointMCMC(phy = rds_ratematrix_matched$phy, data_BM = (rds_ratematrix_matched$num_data),data_Mk = setNames(rds_ratematrix_matched$data[,"genus"],rownames(rds_ratematrix_matched$data)),dir="mcmcburn/",gen = 1000)

conv_check <- checkConvergence(readMCMC(mcmc_out))
conv_check$heidel

plotRatematrix(readMCMC(mcmc_out),colors = c("purple","gold","cyan3"),show.zero = T, n.lines = 100,ell.wd = 2,p = c(1))

