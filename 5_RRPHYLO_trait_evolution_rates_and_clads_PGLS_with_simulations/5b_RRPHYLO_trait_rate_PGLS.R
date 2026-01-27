##revised rrphylo
##source functions and load data
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
library(caper)
library(phytools)
library(ggpubr)

#set working directory to parent directory
#data_pointer <- <YOUR PATH>
data_pointer <- "c:/Users/jtm378/AMNH Dropbox/Jon Merwin/AllParrotHomeAnalysis/Scripts/GITHUB_UPLOAD/"
setwd(data_pointer)

PGLS_data_with_traitrates <- readRDS(file = "5_RRPHYLO_trait_evolution_rates_and_clads_PGLS_with_simulations/PGLS_RDS_with_rates.rds")
PGLS_data_with_traitrates$phy <- root(PGLS_data_with_traitrates$phy,outgroup = c("Calyptomena_viridis","Caracara_cheriway","Icterus_cucullatus"),resolve.root = T)

plot(log(as.numeric(PGLS_data_with_traitrates$data$morph_rates)),(as.numeric(PGLS_data_with_traitrates$data$PC1)))
}

#reformat and transform data for PGLS analysis
{
morph_color_clads_PGLS_data <- cbind.data.frame(names = PGLS_data_with_traitrates$data$species, 
                                          body_size =  scales::rescale(as.numeric(PGLS_data_with_traitrates$data$PC1),to=c(.01,1)), 
                                          body_shape =  scales::rescale(as.numeric(PGLS_data_with_traitrates$data$PC2),to=c(.01,1)),
                                          hue_disp =  scales::rescale(as.numeric(PGLS_data_with_traitrates$data$huedisp.m),to=c(.01,1)),
                                          mean.ra =  scales::rescale(as.numeric(PGLS_data_with_traitrates$data$huedisp.m),to=c(.01,1)),
                                          rr_morph_cov =  scales::rescale(as.numeric(PGLS_data_with_traitrates$data$morph_rates),to=c(.01,1)), 
                                          rr_color_cov =  scales::rescale(as.numeric(PGLS_data_with_traitrates$data$color_rates),to=c(.01,1)),
                                          clads = as.numeric(PGLS_data_with_traitrates$data$ClaDS),
                                          dr_stat =  as.numeric(PGLS_data_with_traitrates$data$DR),
                                          bayes_color = as.numeric(PGLS_data_with_traitrates$data$color.rate),
                                          bayes_morph_raw= as.numeric(PGLS_data_with_traitrates$data$raw.morpho.rate),
                                          bayes_morph_pc= as.numeric(PGLS_data_with_traitrates$data$pcs.morpho.rate))
                                          
morph_color_clads_clim_PGLS_data <- cbind.data.frame(names = PGLS_data_with_traitrates$data$species, 
                                                rr_morph_cov =  scales::rescale(as.numeric(PGLS_data_with_traitrates$data$morph_rates),to=c(.01,1)), 
                                                rr_color_cov =  scales::rescale(as.numeric(PGLS_data_with_traitrates$data$color_rates),to=c(.01,1)),
                                                rr_clim_cov =  scales::rescale(as.numeric(PGLS_data_with_traitrates$data$clim_rates),to=c(.01,1)),
                                                clads = as.numeric(PGLS_data_with_traitrates$data$ClaDS),
                                                dr_stat =  as.numeric(PGLS_data_with_traitrates$data$DR),
                                                bayes_color = as.numeric(PGLS_data_with_traitrates$data$color.rate),
                                                bayes_morph_raw= as.numeric(PGLS_data_with_traitrates$data$raw.morpho.rate),
                                                bayes_morph_pc= as.numeric(PGLS_data_with_traitrates$data$pcs.morpho.rate))

morph_color_clads_comparative <- caper::comparative.data(data =na.omit(morph_color_clads_PGLS_data),phy = PGLS_data_with_traitrates$phy, names.col = "names", na.omit = T, vcv = T)
morph_color_clim_clads_comparative <- caper::comparative.data(data =na.omit(morph_color_clads_clim_PGLS_data),phy = PGLS_data_with_traitrates$phy, names.col = "names", na.omit = T, vcv = T)
}


#exploratory plots
# {
# pairs(morph_color_clads_comparative$data, log = c(1,2,3,4,5,6,7,8,9,10,11))
# pairs(morph_color_clads_comparative$data[,c("rr_morph_cov","rr_color_cov","bayes_color","bayes_morph_raw","bayes_morph_pc")], log = c(1,2,3,4,5,6,7,8,9,10,11))
# pairs(morph_color_clim_clads_comparative$data, log = c(1,2,3,4,5,6))
# 
# ggplot(morph_color_clads_comparative$data) + geom_point(aes(x=sqrt(3+log(5+log(bayes_color))),y=log(rr_color_cov)),col="purple2") + xlim(c(.95,1.1))
# ggplot(morph_color_clads_comparative$data) + geom_point(aes(x=log(bayes_morph_raw),y=log(rr_morph_cov)),col="gold4")
# ggplot(morph_color_clads_comparative$data) + geom_point(aes(x=log(bayes_morph_pc),y=log(rr_morph_cov)),col="red")
# }

clads_regressions <- ggplot(morph_color_clads_comparative$data) + theme_bw() +
  geom_point(alpha = 1, size = 1, color = "gold2",aes(y = log(as.numeric(rr_morph_cov)), x = log(as.numeric(clads)))) +
  geom_point(alpha = 1, size = 1, color = "purple2",  aes(y = log(as.numeric(rr_color_cov)),x = log(as.numeric(clads)))) +
  geom_point(alpha = 1, data = morph_color_clads_clim_PGLS_data, size = 1, color = "cyan3",  aes(y = log(as.numeric(rr_clim_cov)),x = log(as.numeric(clads)))) +
  geom_smooth(color = "gold2", aes(y = log(as.numeric(rr_morph_cov)), x = log(as.numeric(clads))), method = "lm", fill = "gold2") + 
  geom_smooth(color = "purple2", aes(y = log(as.numeric(rr_color_cov)), x = log(as.numeric(clads))), method = "lm", fill = "purple2") +
  geom_smooth(data = morph_color_clads_clim_PGLS_data, color = "cyan3", aes(y = log(as.numeric(rr_clim_cov)), x = log(as.numeric(clads))), method = "lm", fill = "cyan3") +
  xlab("CLADS Rate") + ylab("Trait Evolution Rate")

drstat_regressions <- ggplot(morph_color_clads_comparative$data) + theme_bw() +
  geom_point(alpha = 1, size = 1, color = "gold2",aes(y = log(as.numeric(rr_morph_cov)), x = log(as.numeric(dr_stat)))) +
  geom_point(alpha = 1, size = 1, color = "purple2",  aes(y = log(as.numeric(rr_color_cov)),x = log(as.numeric(dr_stat)))) +
  geom_point(alpha = 1, data = morph_color_clads_clim_PGLS_data, size = 1, color = "cyan3",  aes(y = log(as.numeric(rr_clim_cov)),x = log(as.numeric(clads)))) +
  geom_smooth(color = "gold2", aes(y = log(as.numeric(rr_morph_cov)), x = log(as.numeric(dr_stat))), method = "lm", fill = "gold3") + 
  geom_smooth(color = "purple2", aes(y = log(as.numeric(rr_color_cov)), x = log(as.numeric(dr_stat))), method = "lm", fill = "purple2") +
  geom_smooth(data = morph_color_clads_clim_PGLS_data, color = "cyan3", aes(y = log(as.numeric(rr_clim_cov)), x = log(as.numeric(dr_stat))), method = "lm", fill = "cyan3") +
  xlab("DR Statistic") + ylab("Trait Evolution Rate")


clads_regressions_bayes <- ggplot(morph_color_clads_comparative$data) + theme_bw() +
  geom_point(alpha = 1, size = 1, color = "gold2",aes(y = log(as.numeric(bayes_morph_raw)), x = log(as.numeric(clads)))) +
  geom_point(alpha = 1, size = 1, color = "purple2",  aes(y = log(as.numeric(bayes_color)),x = log(as.numeric(clads)))) +
  geom_smooth(color = "gold2", aes(y = log(as.numeric(bayes_morph_raw)), x = log(as.numeric(clads))), method = "lm", fill = "gold3") + 
  geom_smooth(color = "purple2", aes(y = log(as.numeric(bayes_color)), x = log(as.numeric(clads))), method = "lm", fill = "purple2") +
  xlab("CLADS") + ylab("Bayes Trait Evolution Rate")

clads_regressions_bayes <- ggplot(morph_color_clads_comparative$data) + theme_bw() +
  geom_point(alpha = 1, size = 3, color = "purple2",  aes(y = log(5+log(as.numeric(bayes_color))),x = log(as.numeric(clads)))) +
  geom_smooth(color = "purple2", aes(y = log(5+log(as.numeric(bayes_color))), x = log(as.numeric(clads))), method = "lm", fill = "purple2") +
  xlab("CLADS") + ylab("Bayes Trait Evolution Rate") #+ ylim(c(-2,-1.85))

clads_regressions_bayes <- ggplot(morph_color_clads_comparative$data) + theme_bw() +
  geom_point(alpha = 1, size = 3, color = "gold2",aes(y = log(as.numeric(bayes_morph_raw)), x = log(as.numeric(clads)))) +
 geom_smooth(color = "gold2", aes(y = log(as.numeric(bayes_morph_raw)), x = log(as.numeric(clads))), method = "lm", fill = "gold3") + 
  xlab("CLADS") + ylab("Bayes Trait Evolution Rate")

clads_regressions_bayes

#generates supplemental figure SX?
ggarrange(clads_regressions, drstat_regressions, nrow = 2, ncol =1)

#diagnostic plots of PC space
# pdf("test_morph_pc_space.pdf")
# ggplot(PGLS_data_with_traitrates$data) + theme_bw() +
#   geom_label(aes(label=species, ,y = (as.numeric(PC2)), x = log(250-(as.numeric(PC1)))), size = .3)
# graphics.off()
# system("open test_morph_pc_space.pdf")
# 
# pdf("test_color_coordinate_space.pdf")
# ggplot(PGLS_data_with_traitrates$data) + theme_bw() +
#   geom_label(aes(label=species,y = (as.numeric(huedisp.m)), x = (as.numeric(mean.ra))), size = .3)
# graphics.off()
# system("open test_color_coordinate_space.pdf")

#run PGLS analysis of clads vs. color, plot the quality control plots
{
par(mfrow=c(3,1))
#run PGLS on morphology and color rates. 
#are rates of color and morphology evolution correlated with diversification rates?
summary_morphRates_pgls_clads <- summary(caper::pgls(data = (morph_color_clads_comparative), formula = sqrt(rr_morph_cov)~log(clads), lambda = "ML"))
summary_morphRates_pgls_clads
#check normalcy of residuals
hist(summary_morphRates_pgls_clads$residuals)

summary_colRates_pgls_clads <- summary(caper::pgls(data = morph_color_clads_comparative, formula = log(rr_color_cov)~log(clads), lambda = "ML"))
summary_colRates_pgls_clads
#check normalcy of residuals
hist(summary_colRates_pgls_clads$residuals)

summary_climRates_pgls_clads <- summary(caper::pgls(data = morph_color_clim_clads_comparative, formula = log(rr_clim_cov)~log(clads), lambda = "ML"))
summary_climRates_pgls_clads
#check normalcy of residuals
hist(summary_climRates_pgls_clads$residuals)

par(mfrow=c(2,2))
plot(caper::pgls(data = (morph_color_clads_comparative), formula = sqrt(rr_morph_cov)~log(clads), lambda = "ML"))
plot(caper::pgls(data = morph_color_clads_comparative, formula = log(rr_color_cov)~log(clads), lambda = "ML"))
plot(caper::pgls(data = morph_color_clim_clads_comparative, formula = log(rr_clim_cov)~log(clads), lambda = "ML"))
}

#run PGLS analysis of drstat vs. color, plot the quality control plots
{
  par(mfrow=c(3,1))
  #run PGLS on morphology and color rates. 
  #are rates of color and morphology evolution correlated with diversification rates?
  summary_morphRates_pgls_drstat <- summary(caper::pgls(data = (morph_color_clads_comparative), formula = sqrt(rr_morph_cov)~log(dr_stat), lambda = "ML"))
  summary_morphRates_pgls_drstat
  #check normalcy of residuals
  hist(summary_morphRates_pgls_clads$residuals)
  
  summary_colRates_pgls_drstat <- summary(caper::pgls(data = morph_color_clads_comparative, formula = log(rr_color_cov)~log(dr_stat), lambda = "ML"))
  summary_colRates_pgls_drstat
  #check normalcy of residuals
  hist(summary_colRates_pgls_clads$residuals)
  
  summary_climRates_pgls_drstat <- summary(caper::pgls(data = morph_color_clim_clads_comparative, formula = log(rr_clim_cov)~log(dr_stat), lambda = "ML"))
  summary_climRates_pgls_drstat
  #check normalcy of residuals
  hist(summary_climRates_pgls_drstat$residuals)
  
  par(mfrow=c(2,2))
  plot(caper::pgls(data = (morph_color_clads_comparative), formula = sqrt(rr_morph_cov)~log(dr_stat), lambda = "ML"))
  plot(caper::pgls(data = morph_color_clads_comparative, formula = log(rr_color_cov)~log(dr_stat), lambda = "ML"))
  plot(caper::pgls(data = morph_color_clim_clads_comparative, formula = log(rr_clim_cov)~log(dr_stat), lambda = "ML"))
}

#summarize model output using clads and DR as predictors
summary_morphRates_pgls_clads
summary_morphRates_pgls_drstat

summary_colRates_pgls_clads
summary_colRates_pgls_drstat

summary_climRates_pgls_clads
summary_climRates_pgls_drstat

summary(caper::pgls(data = morph_color_clim_clads_comparative, formula = log(rr_color_cov)~log(rr_clim_cov)+log(clads), lambda = "ML"))
summary(caper::pgls(data = morph_color_clim_clads_comparative, formula = log(rr_morph_cov)~log(rr_clim_cov)+log(clads), lambda = "ML"))

summary(caper::pgls(data = morph_color_clim_clads_comparative, formula = log(rr_color_cov)~log(rr_clim_cov)+log(dr_stat), lambda = "ML"))
summary(caper::pgls(data = morph_color_clim_clads_comparative, formula = log(rr_morph_cov)~log(rr_clim_cov)+log(dr_stat), lambda = "ML"))
#yes! for all

#randomized PGLS trials
#make empty list for R-squared values
#need column named "names", rr_color or rr_morph, dr_stat, or clads
#run this now
run_bootstrapped_PGLS <- function(tree, data, simNum = 100, percentage = 10, divStat = "clads", trait = "morph"){
  R2_list <- list()
  R2_banana_list <- list()
  R2_banana_rate_list <- list()
  trialNum <- simNum
  #run simulations dropping 10 taxa for both real data and randomized data
  for(i in 1:trialNum){
    #scaledDat <- apply(X = (matched_plus_geoClim), MARGIN = 2, FUN = function(X) scaleOnlyNumeric(X))
    bootNum <- round(length(data[,1])*(percentage/100))
    randDrop <- sample(x = 1:length(data[,1]), size = bootNum)
    #comparativeCulled <- comparative.data(phy = withCov$tree, data = (morph_clads_PGLS_data[-c(randDrop),]), names.col = "names")
    
    comparativeCulled <- comparative.data(phy = tree, data = (data[-c(randDrop),]), names.col = "names")
    
    datRandomized <-(data[-c(randDrop),])
    #datRandomized$rr_morph <- comparativeCulled$data$rr_morph[sample(1:length(comparativeCulled$data$rr_morph),)]
    #brownian simulation with bounds set to real empirical bounds
    
    colnameGrab <- paste("rr_", trait, "_cov", sep = "")
    if(verbatim_column){colnameGrab <- trait}
    
    datRandomized$rr_morph <- fastBM(tree, bounds = c(min(comparativeCulled$data[,colnameGrab]),max(comparativeCulled$data[,colnameGrab])))[datRandomized$names]
    datRandomized$rr_color <- datRandomized$rr_morph
    rownames(datRandomized) <- datRandomized$names
    
    comparativeCulledRandomized <- comparative.data(phy = tree, data = datRandomized, names.col = "names")
    
    thisY <- comparativeCulledRandomized$data[comparativeCulledRandomized$phy$tip.label,colnameGrab]
    names(thisY) <- comparativeCulledRandomized$phy$tip.label
    
    rr_rate_randomized <- RRphylo::RRphylo(comparativeCulledRandomized$phy, y = (thisY))
    comparativeCulledRandomized$data$random_rates <- log(scales::rescale(rr_rate_randomized$rates[comparativeCulledRandomized$phy$tip.label,1], to = c(1,2)))
    
    if(divStat == "clads"){
      if(trait == "morph"){thisFormula <- formula("(rr_morph_cov)~clads")}
      if(trait == "color"){thisFormula <- formula("(rr_color_cov)~clads")}
      if(trait == "clim"){thisFormula <- formula("(rr_clim_cov)~clads")}
      thisFormula_rates <- formula("random_rates~clads")
    }
    if(divStat == "dr_stat"){
      if(trait == "color"){thisFormula <- formula("(rr_col_cov)~dr_stat")}
      if(trait == "morph"){thisFormula <- formula("(rr_morph_cov)~dr_stat")}
      if(trait == "climate"){thisFormula <- formula("(rr_morph_cov)~dr_stat")}
      thisFormula_rates <- formula("random_rates~dr_stat")
    }
    
    #thisFormula <- formula("sqrt(1+numLouseGen)~sqrt(1+Mass)+log(1800+X_precipitationPC1_mean)+log(350+X_temperaturePC1_mean)+log(380+X_elevation_max)+log(1+sympatryStats_max)")
    resultPGLS <- summary(caper::pgls(data = (comparativeCulled), formula = thisFormula), lambda = .9)
    resultPGLS_bananas <- summary(caper::pgls(data = comparativeCulledRandomized, formula = thisFormula), lambda = .9)
    resultPGLS_bananas_rates <- summary(caper::pgls(data = comparativeCulledRandomized, formula = thisFormula_rates), lambda = .9)
    #print(resultPGLS)
    #plot(caper::pgls(data = (comparativeCulled), formula = thisFormula))
    
    R2_list[i] <- paste(resultPGLS$adj.r.squared)
    R2_banana_list[i] <- paste(resultPGLS_bananas$adj.r.squared)
    R2_banana_rate_list[i] <- paste(resultPGLS_bananas_rates$adj.r.squared)
    cat("\014")
    print(paste(i, "/", trialNum))
  }
  
  #convert r-squared lists to data frames and bind them together
  bb <- as.data.frame(as.numeric(unlist(R2_list)))
  bb_bananas <- as.data.frame(as.numeric(unlist(R2_banana_list)))
  bb_bananas_rates <- as.data.frame(as.numeric(unlist(R2_banana_rate_list)))
  names(bb) <- "r2"
  names(bb_bananas) <- "r2"
  names(bb_bananas_rates) <- "r2"
  
  bb$trial <- "real_bootstrap"
  bb_bananas$trial <- "bananas_bootstrap"
  bb_bananas_rates$trial <- "bananas_bootstrap_rates"
  
  bb_plus_bananas <- rbind.data.frame(bb, bb_bananas, bb_bananas_rates)
  
  return(bb_plus_bananas)
}

simple_permutation_PGLS <- function(tree, data, simNum = 100, percentage = 10, divStat = "clads", trait = "rr_morph_cov", set_lambda = 1)
{
  require(picante)
  require(caper)
  if(!trait %in% c("rr_morph_cov", "rr_clim_cov", "rr_color_cov")){print("Pick one: rr_morph_cov, rr_clim_cov, or rr_color_cov");return()}
  if(!divStat %in% c("clads", "dr_stat")){print("Pick one: clads, or dr_stat");return()}
  
  #data <- morph_color_clads_clim_PGLS_data
  #tree <- morph_color_clads_comparative$phy
  
  R2_permutation_list <- list()
  comparative_for_bootstrap <- comparative.data(tree, data,names.col = "names",vcv = T)

  if(divStat == "clads"){
    if(trait == "rr_morph_cov"){thisFormula <- formula("log(rr_morph_cov)~log(clads)")}
    if(trait == "rr_color_cov"){thisFormula <- formula("log(rr_color_cov)~log(clads)")}
    if(trait == "rr_clim_cov"){thisFormula <- formula("log(rr_clim_cov)~log(clads)")}
    rand_formula <- formula("log(trait_randomized)~log(clads)")
  }
  if(divStat == "dr_stat"){
    if(trait == "rr_color_cov"){thisFormula <- formula("log(rr_color_cov)~log(dr_stat)")}
    if(trait == "rr_morph_cov"){thisFormula <- formula("log(rr_morph_cov)~log(dr_stat)")}
    if(trait == "rr_clim_cov"){thisFormula <- formula("log(rr_clim_cov)~log(dr_stat)")}
    rand_formula <- formula("log(trait_randomized)~log(dr_stat)")
  }
  
  fullmodel_result_summary <- summary(caper::pgls(data = (comparative_for_bootstrap), formula = thisFormula), lambda = "ML")
  
  for(i in 1:simNum)
  {
    #thisFormula <- formula("sqrt(1+numLouseGen)~sqrt(1+Mass)+log(1800+X_precipitationPC1_mean)+log(350+X_temperaturePC1_mean)+log(380+X_elevation_max)+log(1+sympatryStats_max)")
    data <- na.omit(data)
    random_df <- data.frame(names= data$names, trait_randomized = sample((data[,trait]), size = length(data[,1])), clads = data$clads, dr_stat = data$dr_stat)
    randomized_comparative_data <- comparative.data(tree, random_df,names.col = "names",vcv = T, na.omit = T)
    randomized_comparative_PGLS_summary <- summary(caper::pgls(data = randomized_comparative_data, formula = rand_formula, lambda = set_lambda))
    #print(resultPGLS)
    #plot(caper::pgls(data = (comparativeCulled), formula = thisFormula))
    
    R2_permutation_list[i] <- paste(randomized_comparative_PGLS_summary$adj.r.squared)
    print(paste(i, "/", simNum))
  }
  return_perm_results <- data.frame(R2 = as.numeric(unlist(R2_permutation_list)), trait = trait, stat = divStat, model="permutation")
  return(return_perm_results)
}

#run permutation for clads and plot (this takes a LONG time, try with low simNum)
{
simple_bootstrapped_clads_color_rates <- simple_permutation_PGLS(tree = morph_color_clads_comparative$phy, data = morph_color_clads_PGLS_data, simNum = 1000, divStat = "clads", trait = "rr_color_cov")
simple_bootstrapped_clads_morph_rates <- simple_permutation_PGLS(tree = morph_color_clads_comparative$phy, data = morph_color_clads_PGLS_data, simNum = 1000, divStat = "clads", trait = "rr_morph_cov")
simple_bootstrapped_clads_climate_rates <- simple_permutation_PGLS(tree = morph_color_clads_comparative$phy, data = morph_color_clads_clim_PGLS_data, simNum = 1000, divStat = "clads", trait = "rr_clim_cov")
simple_bootstrapped_clads_color_rates$plot_color   <- "purple2"
simple_bootstrapped_clads_morph_rates$plot_color   <-  "gold2"
simple_bootstrapped_clads_climate_rates$plot_color <-  "cyan3"

all_permutation_data_clads <- rbind.data.frame(simple_bootstrapped_clads_color_rates, simple_bootstrapped_clads_morph_rates, simple_bootstrapped_clads_climate_rates)

permutation_histogram_clads <- ggplot() +
  geom_histogram(data = all_permutation_data_clads, color = "transparent", alpha = .5, aes(fill = I(plot_color), x = (R2))) +
  geom_vline(xintercept = (summary_colRates_pgls_clads$adj.r.squared), size = 2, color = "purple2") +
  geom_vline(xintercept = (summary_morphRates_pgls_clads$adj.r.squared), size = 2, color = "gold2") +
  geom_vline(xintercept = (summary_climRates_pgls_clads$adj.r.squared), size = 2, color = "cyan3") +
  ggplot2::ggtitle("Permutation Test Clads Rates")+ theme(legend.position="none") + theme_bw() + xlab("Adjusted R-squared") + ylab("Count")

# Calculate p-value (one-tailed)
p_value_col_clads <- sum(((simple_bootstrapped_clads_color_rates$R2)) >= (c((summary_colRates_pgls_clads$adj.r.squared)))) / length((simple_bootstrapped_clads_color_rates$R2))
p_value_clim_clads <- sum((simple_bootstrapped_clads_climate_rates$R2) >= (c(summary_climRates_pgls_clads$adj.r.squared))) / length(simple_bootstrapped_clads_climate_rates$R2)
p_value_morph_clads <- sum((simple_bootstrapped_clads_morph_rates$R2) >= abs(c(summary_morphRates_pgls_clads$adj.r.squared))) / length(simple_bootstrapped_clads_morph_rates$R2)

print(paste("Color", p_value_col_clads))
print(paste("Morph", p_value_morph_clads))
print(paste("Climate", p_value_col_clads))


print(paste("Upper tail p-value:", p_value_col_clads))
print(paste("Lower tail p-value:", p_value_col_clads))
}

#run permutation for dr_stat and plot (this takes a LONG time, try with low simNum)
{
  simple_bootstrapped_dr_stat_color_rates <- simple_permutation_PGLS(tree = morph_color_clads_comparative$phy, data = morph_color_clads_PGLS_data, simNum = 1000, divStat = "dr_stat", trait = "rr_color_cov")
  simple_bootstrapped_dr_stat_morph_rates <- simple_permutation_PGLS(tree = morph_color_clads_comparative$phy, data = morph_color_clads_PGLS_data, simNum = 1000, divStat = "dr_stat", trait = "rr_morph_cov")
  simple_bootstrapped_dr_stat_climate_rates <- simple_permutation_PGLS(tree = morph_color_clads_comparative$phy, data = morph_color_clads_clim_PGLS_data, simNum = 1000, divStat = "dr_stat", trait = "rr_clim_cov")
  simple_bootstrapped_dr_stat_color_rates$plot_color   <- "purple2"
  simple_bootstrapped_dr_stat_morph_rates$plot_color   <-  "gold2"
  simple_bootstrapped_dr_stat_climate_rates$plot_color <-  "cyan3"
  
  all_permutation_data_dr_stat <- rbind.data.frame(simple_bootstrapped_dr_stat_color_rates, simple_bootstrapped_dr_stat_morph_rates, simple_bootstrapped_dr_stat_climate_rates)
  
  permutation_histogram_dr_stat <- ggplot() +
    geom_histogram(data = all_permutation_data_dr_stat, color = "transparent", alpha = .5, aes(fill = I(plot_color), x = (R2))) +
    geom_vline(xintercept = (summary_colRates_pgls_drstat$adj.r.squared), size = 2, color = "purple2") +
    geom_vline(xintercept = (summary_morphRates_pgls_drstat$adj.r.squared), size = 2, color = "gold2") +
    geom_vline(xintercept = (summary_climRates_pgls_drstat$adj.r.squared), size = 2, color = "cyan3") +
    ggplot2::ggtitle("Permutation Test DR statistic rates") + theme(legend.position="none") + theme_bw() +
    xlab("Adjusted R-squared") + ylab("Count")
  
  # Calculate p-value (one-tailed)
  p_value_col_dr_stat <- 100*sum((simple_bootstrapped_dr_stat_color_rates$R2) >= (c(summary_colRates_pgls_clads$adj.r.squared))) / length(simple_bootstrapped_clads_color_rates$R2)
  p_value_clim_dr_stat <- 100*sum((simple_bootstrapped_dr_stat_climate_rates$R2) >= (c(summary_climRates_pgls_clads$adj.r.squared))) / length(simple_bootstrapped_dr_stat_climate_rates$R2)
  p_value_morph_dr_stat <- 100*sum((simple_bootstrapped_dr_stat_morph_rates$R2) >= (c(summary_morphRates_pgls_clads$adj.r.squared))) / length(simple_bootstrapped_dr_stat_morph_rates$R2)

  print(paste("Color", p_value_col_dr_stat))
  print(paste("Morph", p_value_morph_dr_stat))
  print(paste("Climate", p_value_clim_dr_stat))
}

ggarrange(clads_regressions, drstat_regressions, permutation_histogram_clads, permutation_histogram_dr_stat,ncol = 2, nrow = 2)

color_rates_histogram_clads

#essim returned similar results
#source essim method (downloaded from github)
#Harvey, MG and Rabosky, DL. 2018. Continuous traits and speciation rates: Alternatives to state-dependent diversification models. Methods in Ecology and Evolution 9: 984-993. 
#https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12949
#download available at https://github.com/mgharvey/ES-sim
source("5_RRPHYLO_trait_evolution_rates_and_clads_PGLS_with_simulations/essim.R")
{
nsims_param = 1000
  
es_dat <- c(morph_color_clads_comparative$data$rr_morph_cov)
named_es_clads <- c(morph_color_clads_comparative$data$clads)
names(es_dat) <- rownames(morph_color_clads_comparative$data)
names(named_es_clads) <- rownames(morph_color_clads_comparative$data)
essim_morph_es <- essim(phy = morph_color_clads_comparative$phy, trait = es_dat,nsim = nsims_param)
essim_morph_clads <- essim(phy = morph_color_clads_comparative$phy, trait = es_dat, es = named_es_clads,nsim = nsims_param)

es_dat <- c(morph_color_clads_comparative$data$rr_color_cov)
named_es_clads <- c(morph_color_clads_comparative$data$clads)
names(es_dat) <- rownames(morph_color_clads_comparative$data)
names(named_es_clads) <- rownames(morph_color_clads_comparative$data)
essim_color_es <- essim(phy = morph_color_clads_comparative$phy, trait = es_dat,nsim = nsims_param)
essim_color_clads <- essim(phy = morph_color_clads_comparative$phy, trait = es_dat, es = named_es_clads, nsim = nsims_param)

es_dat <- c(morph_color_clim_clads_comparative$data$rr_clim_cov)
named_es_clads <- c(morph_color_clim_clads_comparative$data$clads)
names(es_dat) <- rownames(morph_color_clim_clads_comparative$data)
names(named_es_clads) <- rownames(morph_color_clim_clads_comparative$data)
essim_clim_es <- essim(phy = morph_color_clim_clads_comparative$phy, trait = es_dat,nsim = nsims_param)
essim_clim_clads <- essim(phy = morph_color_clim_clads_comparative$phy, trait = es_dat,es = named_es_clads,nsim = nsims_param)

bound_essim_results <- rbind.data.frame(essim_morph_es,essim_morph_clads,essim_color_es,essim_color_clads,essim_clim_es,essim_clim_clads)
colnames(bound_essim_results) <- c("rho", "P_val")
bound_essim_results$div_metric <- c("ES", "CLADS")
bound_essim_results$trait_rate <- c("Morph", "Morph", "Color", "Color", "Climate", "Climate")

ggplot(bound_essim_results) +
  geom_point(data = bound_essim_results[which(p.adjust(bound_essim_results$P_val, method = "fdr")<.05),],aes(x=div_metric,y=rho, color = "FDR-corrected P < .05"), size = 10, shape = 5) +
  geom_point(aes(x=div_metric,y=rho, color = trait_rate), size = 5) +
  theme_bw() + scale_color_viridis_d(begin = .1, end = .8)
}

morph_color_clads_PGLS_data
morph_color_clads_clim_PGLS_data


