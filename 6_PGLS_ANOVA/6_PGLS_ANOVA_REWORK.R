##revised rrphylo

#load packages and data
{
library(devtools)
library(phytools)
library(motmot)
library(ape)
library(geiger)
library(RRphylo)
library(ggplot2)
library(ggpmisc)
library(caper)
library(gridExtra)
#for making tables
library(gt)
library(ggpubr)
library(picante)
library(mvMORPH)


#useful method for matching later 
lastDropper <- function(string = NULL, seper = " ")
{
  require("stringr")
  lab <- word(string, 1, -2, sep = seper)
  return(lab)
}

#set working directory to parent directory
#data_pointer <- <YOUR PATH>
data_pointer <- "c:/Users/jonme/AMNH Dropbox/Jon Merwin/AllParrotHomeAnalysis/Scripts/GITHUB_UPLOAD/"
#data_pointer <- "c:/Users/jtm378/AMNH Dropbox/Jon Merwin/AllParrotHomeAnalysis/Scripts/GITHUB_UPLOAD/"
setwd(data_pointer)

#match tree and data
allMatched <- readRDS("5_RRPHYLO_trait_evolution_rates_and_clads_PGLS_with_simulations/PGLS_RDS_with_rates.rds")

#fix data type (extra space after "Oceania")
allMatched$data$Continent_clean <- gsub(allMatched$data$Continent, pattern = "Oceania ", replacement = "Oceania")
}

#boxplots of continent-scale data, phylogenetic ANOVAs
{
clads_vs_continent <- ggplot(allMatched$data, aes(x = Continent_clean, y = log(as.numeric(ClaDS)), fill = Continent_clean)) + geom_boxplot() + theme_bw()
morphrate_vs_continent <- ggplot(allMatched$data, aes(x = Continent_clean, y = log(as.numeric(morph_rates)), fill = Continent_clean)) + geom_boxplot() + theme_bw()
colrate_vs_continent <- ggplot(allMatched$data, aes(x = Continent_clean, y = log(as.numeric(color_rates)), fill = Continent_clean)) + geom_boxplot() + theme_bw()
clim_vs_continent <- ggplot(allMatched$data, aes(x = Continent_clean, y = log(as.numeric(clim_rates)), fill = Continent_clean)) + geom_boxplot() + theme_bw()
ggplot(allMatched$data, aes(x = Continent_clean, y = log(as.numeric(color.rate)), fill = Continent_clean)) + geom_boxplot() + theme_bw() + ylim(c(-5,-4.8))
ggplot(allMatched$data, aes(x = Continent_clean, y = log(as.numeric(raw.morpho.rate)), fill = Continent_clean)) + geom_boxplot() + theme_bw() + ylim(c(2,5))

ggarrange(ncol =2, nrow = 2, clads_vs_continent, morphrate_vs_continent, colrate_vs_continent, clim_vs_continent)

#phylogenetic anova, nothing is near significance
phylANOVA_morph_match <- match.phylo.data(phy = allMatched$phy, data = na.omit(allMatched$data[,c("Continent_clean", "ClaDS", "morph_rates", "color_rates","color.rate","raw.morpho.rate")]))
phylANOVA_clim_match  <-  match.phylo.data(phy = allMatched$phy, data = na.omit(allMatched$data[,c("Continent_clean", "ClaDS", "clim_rates")]))
phylANOVA(phylANOVA_morph_match$phy, x = as.factor(phylANOVA_morph_match$data$Continent_clean), y = as.numeric(phylANOVA_morph_match$data$ClaDS), nsim = 1000, posthoc = T, p.adj = "fdr")
phylANOVA(phylANOVA_morph_match$phy, x = as.factor(phylANOVA_morph_match$data$Continent_clean), y = as.numeric(phylANOVA_morph_match$data$morph_rates), nsim = 1000, posthoc = T, p.adj = "fdr")
phylANOVA(phylANOVA_morph_match$phy, x = as.factor(phylANOVA_morph_match$data$Continent_clean), y = as.numeric(phylANOVA_morph_match$data$color_rates), nsim = 1000, posthoc = T, p.adj = "fdr")
phylANOVA(phylANOVA_clim_match$phy, x = as.factor(phylANOVA_clim_match$data$Continent_clean), y = as.numeric(phylANOVA_clim_match$data$clim_rates), nsim = 1000, posthoc = T, p.adj = "fdr")

#no pattern w. the bayesian rates either
phylANOVA(phylANOVA_morph_match$phy, x = as.factor(phylANOVA_morph_match$data$Continent_clean), y = as.numeric(phylANOVA_morph_match$data$color.rate), nsim = 1000, posthoc = T, p.adj = "fdr")
phylANOVA(phylANOVA_morph_match$phy, x = as.factor(phylANOVA_morph_match$data$Continent_clean), y = as.numeric(phylANOVA_morph_match$data$raw.morpho.rate), nsim = 1000, posthoc = T, p.adj = "fdr")
}

#test whether specific morphologies are associated with specific climates



#root the tree
rooted_phy <- root(allMatched$phy, outgroup = c("Caracara_cheriway","Icterus_cucullatus", "Calyptomena_viridis"), resolve.root = T)

#Question 1: Is parrot body size associated with climate (larger parrots in colder areas, etc...)
#the answer seems to be no... R2:-0.001689, p: 0.4744 
#however, whne not accounting for phylogeny, we did find a significant relationship...
#felsenstein strikes again!
{
  dataset1 <- comparative.data(rooted_phy,data = na.omit(allMatched$data[,c("species","PC1", "PC2","KIPPS_RATIO","WING.LENGTH.TO.PRIMARIES","WING.LENGTH.TO.SECONDARIES","TARSUS.LENGTH", "BILL.WIDTH","BILL.HEIGHT", "XCentroid0.0", "yCentroid0.0","Total.niche")]),names.col = "species")
  save_rownames <- rownames(dataset1$data)
  save_colnames <- colnames(dataset1$data)
  dataset1$data <- data.frame(apply(dataset1$data, MARGIN = 2, FUN = function(X) as.numeric(X)))
  dataset1$data$HWI <- 100*((dataset1$data$WING.LENGTH.TO.PRIMARIES-dataset1$data$WING.LENGTH.TO.SECONDARIES)/dataset1$data$WING.LENGTH.TO.PRIMARIES)
  
  rownames(dataset1$data) <- save_rownames
  colnames(dataset1$data) <- c(save_colnames,"HWI")
  
  fitc_list <- (dataset1$data[,"BILL.WIDTH"])
  names(fitc_list) <- rownames(dataset1$data)
  bm_fit <- fitContinuous(phy = dataset1$phy, dat = fitc_list, model = c("BM"))
  ou_fit <- fitContinuous(phy = dataset1$phy, dat = fitc_list, model = c("OU"))
  wt_fit <- fitContinuous(phy = dataset1$phy, dat = fitc_list, model = c("white"))
  bm_fit$opt$aicc;bm_fit$opt$sigsq
  ou_fit$opt$aicc;ou_fit$opt$alpha
  wt_fit$opt$aicc;wt_fit$opt$sigsq
  
  #check data for skew
  hist(sqrt(max(dataset1$data$PC1)-dataset1$data$PC1))
  hist(dataset1$data$XCentroid0.0)
  hist(dataset1$data$yCentroid0.0)
  hist(sqrt(dataset1$data$Total.niche))
  hist((dataset1$data[,"HWI"]))
  formula1 <- formula("sqrt(max(PC1)-PC1)~XCentroid0.0:yCentroid0.0")
  formula1 <- formula("(HWI)~XCentroid0.0")
  #formula1 <- formula("sqrt(max(PC1)-PC1)~(XCentroid0.0)")
  model1 <- caper::pgls(data =dataset1, formula = formula1, lambda = "ML")
  par(mfrow = c(2,2));plot(model1)
  model1_summary <- summary(model1)
  print(model1_summary)
  ggplot(dataset1$data, aes(y=sqrt(max(PC1)-PC1), x=((XCentroid0.0)))) + geom_point() + geom_smooth(method = "lm", lwd=2, col = "black") + xlab("Climate (PC1)") + ylab("Body Size") + theme_classic()
  ggplot(dataset1$data, aes(y=HWI, x=((XCentroid0.0)))) + geom_point() + geom_smooth(method = "lm", lwd=2, col = "black") + xlab("Climate (PC1)") + ylab("Kipps Index") + theme_classic() + geom_text(aes(label=rownames(dataset1$data)))
}

#Question 2: Is parrot body shape associated with climate
#no...R2: 0.002974, p: 0.2724, lambda = 0.996
#for kipps vs temp seasonality though... Adjusted R-squared: 0.04685 
#F-statistic:  8.57 on 2 and 306 DF,  p-value: 0.0002392  
{
  dataset2 <- comparative.data(rooted_phy,data = na.omit(allMatched$data[,c("species","PC2", "PC3", "KIPPS_RATIO", "XCentroid0.0", "yCentroid0.0","Total.niche")]),names.col = "species")
  save_rownames <- rownames(dataset2$data)
  save_colnames <- colnames(dataset2$data)
  dataset2$data <- data.frame(apply(dataset2$data, MARGIN = 2, FUN = function(X) as.numeric(X)))
  rownames(dataset2$data) <- save_rownames
  colnames(dataset2$data) <- save_colnames
  #check data for skew
  par(mfrow = c(2,2))
  hist((100+dataset2$data$PC2))
  hist(dataset2$data$XCentroid0.0)
  hist(dataset2$data$yCentroid0.0)
  hist(sqrt(dataset2$data$Total.niche))
  #dataset2$data$PC2 <- dataset2$data$PC2^.6
  #formula2 <- formula("PC2^.6~XCentroid0.0+yCentroid0.0+Total.niche")
  #cutting out insignificant variables
  #formula2 <- formula("PC2~XCentroid0.0+sqrt(Total.niche)")
  formula2 <- formula("PC2~XCentroid0.0")
  
  #test kipps ratio...comment to just model PC2, stronger signal w. kipps
  formula2 <- formula("KIPPS_RATIO~XCentroid0.0+yCentroid0.0")
  
  model2   <- caper::pgls(data =dataset2, formula = formula2, lambda = "ML")
  par(mfrow = c(2,2));plot(model2)
  model2_summary <- summary(model2)
  print(model2_summary)
  #rownames(dataset2$data)[order(dataset2$data$PC2, decreasing = T)]
}

#Question 3: Are parrots in certain climates more colorful? In certain climates less colorful? (color span)
#seems like a slight effect where birds in dryer envs are less colorful, less saturated,  but super weak
# Multiple R-squared: 0.01446,	Adjusted R-squared: 0.01124 
# F-statistic:   4.6 on 1 and 307 DF,  p-value: 0.03469 
{
  dataset3 <- comparative.data(rooted_phy,data = na.omit(allMatched$data[,c("species","huedisp.m","centroid.l","centroid.s","colspan.m","mean.ra", "rel.c.vol","a.vol","XCentroid0.0", "yCentroid0.0","Total.niche","PC1_climate","PC2_climate","color_PC1","color_PC1","color_PC2")]),names.col = "species")
  save_rownames <- rownames(dataset3$data)
  save_colnames <- colnames(dataset3$data)
  dataset3$data <- data.frame(apply(dataset3$data, MARGIN = 2, FUN = function(X) as.numeric(X)))
  #dataset3$data <- data.frame(scale(dataset3$data, center = F, scale = T))
  rownames(dataset3$data) <- save_rownames
  colnames(dataset3$data) <- save_colnames
  #check data for skew
  hist((dataset3$data$mean.ra))
  hist(dataset3$data$XCentroid0.0)
  hist(dataset3$data$yCentroid0.0)
  hist(sqrt(dataset3$data$Total.niche))
  dataset3$data$Total.niche_sqrt <- sqrt(dataset3$data$Total.niche)
  dataset3$data <- data.frame(scale(dataset3$data,scale = T))
  #formula3 <- formula("(colspan.m)~XCentroid0.0*yCentroid0.0+(Total.niche_sqrt)")
  #formula3 <- formula("(huedisp.m+mean.ra)~XCentroid0.0:yCentroid0.0")
  formula3 <- formula("(huedisp.m)~XCentroid0.0:yCentroid0.0")
  
  ggplot(dataset3$data) + geom_text(aes(x=huedisp.m, y = mean.ra, label=rownames(dataset3$data), ))
  
  model3   <-caper::pgls(data =dataset3, formula = formula3, lambda = "ML")
  par(mfrow = c(2,2));plot(model3)
  profile_lambda=pgls.profile(model3, which="lambda") # vary lambda
  plot(profile_lambda)
  model3_summary <- summary(model3)
  print(model3_summary)
  
  
}

#Question 4: Does larger body size lead to reduced selection for crypsis? (more green)
#yes! Larger parrots are less likely to be green! (R2: .11, p val: 6.474e-11, coefficient= -0.02)
{
  dataset4 <- comparative.data(rooted_phy,data = na.omit(allMatched$data[,c("species","centroid.u","centroid.m","centroid.l","centroid.s","PC1", "centroid.u_ventral","centroid.s_ventral","centroid.m_ventral","centroid.l_ventral","centroid.u_dorsal","centroid.s_dorsal","centroid.m_dorsal","centroid.l_dorsal", "mean.ra_dorsal", "mean.ra_ventral")]),names.col = "species")
  
  save_rownames <- rownames(dataset4$data)
  save_colnames <- colnames(dataset4$data)
  dataset4$data <- data.frame(apply(dataset4$data, MARGIN = 2, FUN = function(X) as.numeric(X)))
  
  rownames(dataset4$data) <- save_rownames
  colnames(dataset4$data) <- save_colnames
  
  #Ventral Centroid Green Max Proportion
  majority_table <- tabulate(apply(X = dataset4$data[,1:4], MARGIN = 1, function(X) which.max(X)))
  majority_table[2] / length(rooted_phy$tip.label)
  
  print("Ventral Centroid Green Max Proportion")
  majority_table_dorsal <- tabulate(apply(X = dataset4$data[,c("centroid.u_dorsal", "centroid.s_dorsal", "centroid.m_dorsal", "centroid.l_dorsal")], MARGIN = 1, function(X) which.max(X)))
  majority_table_dorsal[3] / length(rooted_phy$tip.label)
  
  print("Ventral Centroid Green Max Proportion")
  majority_table_ventral <- tabulate(apply(X = dataset4$data[,c("centroid.u_ventral", "centroid.s_ventral", "centroid.m_ventral", "centroid.l_ventral")], MARGIN = 1, function(X) which.max(X)))
  majority_table_ventral[3] / length(rooted_phy$tip.label)
  
  #check data for skew
  hist(sqrt(max(dataset4$data$PC1)-dataset4$data$PC1))
  hist(dataset4$data$centroid.u)
  hist(dataset4$data$centroid.s)
  hist(dataset4$data$centroid.m)
  hist((dataset4$data$centroid.l))
  formula4 <- formula("centroid.m/(centroid.u+centroid.s+centroid.l)~sqrt(max(PC1)-PC1)")
  model4   <-caper::pgls(data =dataset4, formula = formula4, lambda = "ML")
  par(mfrow = c(2,2));plot(model4)
  model4_summary <- summary (model4)
  print(model4_summary)
  #plot the pattern
  greenColor_solo_scatterplot <- ggplot(dataset4$data, aes(x = log(260-as.numeric(PC1)), y = (as.numeric(centroid.m)/(as.numeric(centroid.u)+as.numeric(centroid.s)+as.numeric(centroid.l))))) + geom_smooth(method = "lm",col="black",fill="black") + geom_point(size = 6,aes(col=I(rgb(red=centroid.l,green=centroid.m,blue=centroid.s)))) + theme_bw() + xlab("Body Size (Small->Large)") + ylab("Green / Red+Blue+UV")
}

#plots of larger parrots being less green...
{
ventral_ventral_green <- na.omit(allMatched$data[,c("species", "PC1", "max.ra", "max.ra_dorsal","color_rates","centroid.m","centroid.u_ventral","centroid.s_ventral","centroid.m_ventral","centroid.l_ventral","centroid.u_dorsal","centroid.s_dorsal","centroid.m_dorsal","centroid.l_dorsal", "mean.ra_dorsal", "mean.ra_ventral")])
save_spp <- ventral_ventral_green$species
ventral_ventral_green <- data.frame(apply(ventral_ventral_green, MARGIN = 2, FUN = function(X) as.numeric(X)))
ventral_ventral_green$species <- save_spp

ventral_ventral_green$dorsal_rgb <- rgb(red=ventral_ventral_green$centroid.l_dorsal,green=ventral_ventral_green$centroid.m_dorsal, blue=ventral_ventral_green$centroid.s_dorsal*1.1, maxColorValue = 1)
ventral_ventral_green$ventral_rgb <- rgb(red=ventral_ventral_green$centroid.l_ventral,green=ventral_ventral_green$centroid.m_ventral, blue=ventral_ventral_green$centroid.s_ventral*1.1, 1)

ggplot(ventral_ventral_green) + geom_point(size = 3, shape = 16,aes(x=log(260-PC1),y= centroid.m_ventral, color=I(dorsal_rgb))) + geom_abline(slope = 1,intercept = 0)
ggplot(ventral_ventral_green) + geom_point(size = 3, shape = 16,aes(x=log(260-PC1),y= centroid.m_ventral, color=I(ventral_rgb))) + geom_abline(slope = 1,intercept = 0)


ggplot(ventral_ventral_green) + geom_point(size = 6, shape = 15,aes(x=log(260-PC1),y= centroid.m, color=I(ventral_rgb))) + 
  geom_point(size = 6, shape = 15,aes(x=log(260-PC1),y= centroid.m+.008, color=I(dorsal_rgb))) + geom_abline(slope = 1,intercept = 0)

ggplot(ventral_ventral_green) + geom_point(size = 6, shape = 16,aes(x=centroid.m_dorsal/centroid.u_dorsal,y= centroid.m_ventral/centroid.u_ventral, color=I(ventral_rgb))) + 
  geom_point(size = 6, shape = 16,aes(x=centroid.m_dorsal/centroid.u_dorsal,y= centroid.m_ventral/centroid.u_ventral+.2, color=I(dorsal_rgb))) + geom_abline(slope = 1,intercept = 0)

ggplot(ventral_ventral_green) + geom_point(size = 6, shape = 15,aes(x=centroid.s_dorsal/centroid.l_dorsal,y= centroid.s_ventral/centroid.l_ventral, color=I(ventral_rgb))) + 
  geom_point(size = 6, shape = 15,aes(x=centroid.s_dorsal/centroid.l_dorsal,y=.06+(centroid.s_ventral/centroid.l_ventral), color=I(dorsal_rgb))) + geom_abline(slope = 1,intercept = 0)
  #scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10")
 # geom_label(size = 6, shape = 16,aes(x=log(260-PC1),y= centroid.m+.6, label=species), size = .6)

ggplot(ventral_ventral_green,aes(x=log(260-PC1),y= centroid.s_ventral/centroid.l_ventral)) + geom_point(size = 6, shape = 15,aes(color=I(ventral_rgb))) + 
  geom_point(size = 6, shape = 15,aes(x=log(260-PC1),y=.06+(centroid.s_ventral/centroid.l_ventral), color=I(dorsal_rgb))) + geom_abline(slope = 1,intercept = 0) + scale_x_continuous(trans="log10") +
  geom_text(aes(label=species), size = 2)
}

#Question 5: Are rates of color evolution more tied to diversification rates or rates of niche evolution
#Answer: Color is tied to both, but more closely with diversification rates.
# Coefficients:
#   Estimate Std. Error t value  Pr(>|t|)    
# (Intercept)     -2.293636   0.497776 -4.6078 6.974e-06 ***
#   log(clim_rates)  0.211800   0.086114  2.4696   0.01446 *  
#   log(ClaDS)       0.320036   0.108273  2.9668   0.00336 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.06 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2026 on 307 degrees of freedom
# Multiple R-squared: 0.0632,	Adjusted R-squared: 0.04703 
# F-statistic: 8.624 on 2 and 307 DF,  p-value: 0.0002269 
{
  dataset5 <- comparative.data(rooted_phy,data = na.omit(allMatched$data[,c("species", "color_rates", "clim_rates","morph_rates","DR","ClaDS")]),names.col = "species")
  save_rownames <- rownames(dataset5$data)
  save_colnames <- colnames(dataset5$data)
  dataset5$data <- data.frame(apply(dataset5$data, MARGIN = 2, FUN = function(X) as.numeric(X)))
  rownames(dataset5$data) <- save_rownames
  colnames(dataset5$data) <- save_colnames
  #check data for skew
  hist(log(dataset5$data$color_rates))
  hist(log(dataset5$data$morph_rates))
  hist(log(dataset5$data$clim_rates))
  hist(log(dataset5$data$color_rates))
  hist(log(dataset5$data$DR))
  hist(log(dataset5$data$ClaDS))
  
  pairs(log(dataset5$data))
  
  plot(log(dataset5$data$color_rates)~log(dataset5$data$clim_rates)/log(dataset5$data$DR))
  
  #formula5 <- formula("log(color_rates)~log(clim_rates)+log(ClaDS)+log(morph_rates)")
  #formula5 <- formula("log(color_rates)~log(ClaDS)")
  formula5 <- formula("log(color_rates)~log(clim_rates)")
  model5   <- caper::pgls(data =dataset5, formula = formula5, lambda = "ML")
  par(mfrow = c(2,2));plot(model5)
  model5_summary <- summary (model5)
  print(model5_summary)
}

#Question 6: Are rates of morphological evolution more tied to diversification rates, rates of niche evolution, rates of color evolution
# Branch length transformations:
#   
#   kappa  [Fix]  : 1.000
# lambda [ ML]  : 0.136
# lower bound : 0.000, p = 0.19138
# upper bound : 1.000, p = < 2.22e-16
# 95.0% CI   : (NA, 0.465)
# delta  [Fix]  : 1.000
# 
# Coefficients:
#   Estimate Std. Error t value  Pr(>|t|)    
# (Intercept)       0.918550   0.163752  5.6094 4.544e-08 ***
#   log(color_rates)  0.092769   0.031988  2.9001  0.004001 ** 
#   log(clim_rates)  -0.049939   0.048811 -1.0231  0.307071    
# log(ClaDS)        0.057150   0.060064  0.9515  0.342114    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.08427 on 306 degrees of freedom
# Multiple R-squared: 0.03466,	Adjusted R-squared: 0.0252 
# F-statistic: 3.663 on 3 and 306 DF,  p-value: 0.01278 
{
  dataset6 <- comparative.data(rooted_phy,data = na.omit(allMatched$data[,c("species", "PC1","color_rates", "clim_rates","morph_rates","DR","ClaDS","PC1_climate","PC2_climate")]),names.col = "species")
  save_rownames <- rownames(dataset6$data)
  save_colnames <- colnames(dataset6$data)
  dataset6$data <- data.frame(apply(dataset6$data, MARGIN = 2, FUN = function(X) as.numeric(X)))
  rownames(dataset6$data) <- save_rownames
  colnames(dataset6$data) <- save_colnames
  #check data for skew
  hist(log(dataset6$data$color_rates))
  hist(log(dataset6$data$morph_rates))
  hist(log(dataset6$data$clim_rates))
  hist(log(dataset6$data$color_rates))
  hist(log(dataset6$data$DR))
  hist(log(dataset6$data$ClaDS))
  #switch between models to test morph relationship with clads and clim rates
  #formula6 <- formula("log(morph_rates)~log(clim_rates)+log(ClaDS)+log(color_rates)")
  #formula6 <- formula("log(morph_rates)~log(ClaDS)+log(color_rates)")
  formula6 <- formula("log(morph_rates)~log(ClaDS)")
  formula6 <- formula("log(morph_rates)~log(clim_rates)")
  formula6 <- formula("log(clim_rates)~log(ClaDS)")
  model6   <-caper::pgls(data =dataset6, formula = formula6, lambda = "ML")
  par(mfrow = c(2,2));plot(model6)
  model6_summary <- summary(model6)
  print(model6_summary)
}

#adjust p values for multiple testing
p.adjust(c(8e-4, 5.58e-11, .005, 7.1e-13, 4.0e-3, .0003), method = "bonferroni")<.05

#Body size ~ Climate
model1_summary

#Body shape ~ Climate
#or Kipps ~ Temp Seasonality
model2_summary

#Color span ~ Climate
model3_summary

# G/RBU ~ Body Size
model4_summary

#Color Rates ~ Niche Rates + CLADS
model5_summary

#Morph Rates ~ Niche Rates + CLADS
model6_summary












