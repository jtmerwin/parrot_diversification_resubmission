#get colors for all points
#note, pavo and alphashape3d can cause issues with R 4.4, make sure you are on 4.5
#you may need to manually delete and re-install xfun and rgl.

#initialize functions and library, load plumage data
{
library(plot.matrix)
library(circlize)
library(plyr)
library(dplyr)
library(tibble)
library(colorspace)
library(tibble)
library(stringr)
library(pavo)
library(picante)
library(alphashape3d)
library(purrr)

#useful method for matching later 
lastDropper <- function(string = NULL, seper = " ")
{
  require("stringr")
  lab <- word(string, 1, -2, sep = seper)
  return(lab)
}

#set data pointer
data_pointer <- "c:/Users/jtm378/AMNH Dropbox/Jon Merwin/AllParrotHomeAnalysis/Scripts/GITHUB_UPLOAD/"
# data_pointer <- <YOUR DIRECTORY HERE>
setwd(data_pointer)

#if you have not run color_load_collate.R, uncomment line below 
#total_data_metadata_cleaned <- read.csv("data/new_specimens_total_color_data.Jul.2025")

total_data_metadata_cleaned <- total_data_metadata_cleaned_dropped

#requires objects from color_load_collate.R in memory, or uncomment and load
head(total_data_metadata_cleaned)

}

#clean data and add relevant useful columns (i.e. sex, hex values, color opponent coords.)
#also start fixing taxonomy / naming issues
#load tree and reshape
{
#generate all hex colors
totColsAllPoints <- rgb(red = scales::rescale(x = total_data_metadata_cleaned$lwMean, to = c(0,256)), green = scales::rescale(x = total_data_metadata_cleaned$mwMean, to = c(0, 256)), blue = scales::rescale(x = total_data_metadata_cleaned$swMean, to = c(0, 256)), maxColorValue = 256)

#quality control histograms
# windows()
# par(mfrow = c(2, 2))
# hist(total_data_metadata_cleaned$lwMean)
# hist(total_data_metadata_cleaned$mwMean)
# hist(total_data_metadata_cleaned$swMean)
# hist(total_data_metadata_cleaned$uvMean)

#rotate data to get all colors for each taxon 
#badNumber_hotSwap <- read.csv("data/number_hotfixes.csv")
#try this sheet instead

#this sheet helps correct data mislabeling, typos, etc.
badNumber_hotSwap <- read.csv("data/tip_hotswap_7.16.25.csv")
metadata_taxonomy <- read.csv("data/tip_hotswap_7.16.25.csv")
#specimen_numbers_cleaned <- as.numeric(stringr::str_trim(stringr::str_extract_all(totDatTib$SpecimenNumber, "\\d+", simplify = T)))

#total_data_metadata_cleaned$SpecimenNumber <- mapvalues(total_data_metadata_cleaned$SpecimenNumber, from = badNumber_hotSwap$PhotoNum, to = badNumber_hotSwap$SpecNum)

totDatTib <- tibble(total_data_metadata_cleaned)
specimen_numbers_cleaned <- totDatTib$SpecimenNumber

# print(total_data_metadata_cleaned[grep(total_data_metadata_cleaned$Label, pattern= "281404"),])
# print(total_data_metadata_cleaned[grep(total_data_metadata_cleaned$Label, pattern= "Psittacella"),])

results <- data.frame()
for (i in seq_along(metadata_taxonomy$Number_in_color_data)) {
  num <- metadata_taxonomy$Number_in_color_data[i]
  matching_rows <- totDatTib %>%
    filter(str_detect(Label, as.character(num))) %>%
    mutate(number_index = i,
           matched_number = num)
  results <- bind_rows(results, matching_rows)
}

results <- results %>%
  mutate(taxName = mapvalues(matched_number, 
                             from = metadata_taxonomy$Number_in_color_data,
                             to = metadata_taxonomy$Species_in_merwin_et_al,
                             warn_missing = FALSE))


results <- results %>%
  mutate(sex = mapvalues(matched_number, 
                             from = metadata_taxonomy$Number_in_color_data,
                             to = metadata_taxonomy$Sex,
                             warn_missing = FALSE))

results$sex[which(is.na(results$sex))] <- "M_infer"

root <- read.nexus(file = "data/parrot_timetree.tre")
root$tip.label <- lastDropper(root$tip.label, "_")
root <- ape::root(root, outgroup = c("Caracara_cheriway", "Icterus_cucullatus", "Calyptomena_viridis"))

#okay now i have birds that probably need to be matched up...try tomorrow
intersect(gsub(root$tip.label, pattern="_", replacement=" "), unique(results$taxName))
setdiff(gsub(root$tip.label, pattern="_", replacement=" "), unique(results$taxName))
setdiff(y=gsub(root$tip.label, pattern="_", replacement=" "), x=unique(results$taxName))

ab <- pavo::colspace(abs(results[, c("uvMean", "swMean", "mwMean", "lwMean")]/results$dblMean), space = "tcs", "Qi")
plot(ab)
abCols <- rgb(red = ab$l, green = ab$m, blue=ab$s)
abColsReal <- rgb(red = scales::rescale(results$lwMean, to = c(.3,1)), blue = scales::rescale(results$swMean, to = c(.3,1)), green = scales::rescale(results$mwMean, to = c(.3,1)), maxColorValue = 1)

root <- read.nexus(file = "data/parrot_timetree.tre")
root$tip.label <- lastDropper(root$tip.label, "_")
root <- ape::root(root, outgroup = c("Caracara_cheriway", "Icterus_cucullatus", "Calyptomena_viridis"))

intersect(unique(results$taxName), gsub(root$tip.label, pattern ="_",replacement = " "))

in_data_not_tree <- setdiff(x=results$taxName, y=gsub(root$tip.label, pattern ="_",replacement = " "))
in_tree_not_data <- setdiff(y=results$taxName, x=gsub(root$tip.label, pattern ="_",replacement = " "))

#should be missing Pionus sordidus, Pyrrhura caeruleiceps, and the 3 outgroups
print(in_tree_not_data)


#attach to data
totDatTib <- as_tibble(results)
totDatTib <- add_column(totDatTib, abColsReal)
sppNames <- gsub(totDatTib$taxName, pattern = " ", replacement = "_")
totDatTibSpp <- add_column(totDatTib, sppNames)
totDatTibGrouped <- group_by(totDatTibSpp, taxName)

#nameSwapSheet <- read.csv("data/colDatNames_7.14.23.csv")
#nameSwapSheet1 <- nameSwapSheet[which(str_length(nameSwapSheet$ColorNames) > 1), ]

#testing
totDatTib[grep(x=totDatTib$taxName, pattern = "Pyrrhura amazonum"),]

unique(totDatTibGrouped$taxName)

setdiff(gsub(root$tip.label,pattern="_",replacement=" "),totDatTib$taxName)

#root2 <- root
#root2$tip.label <- gsub(pattern = "_", replacement = " ", x = root2$tip.label)
#newLabels <- plyr::mapvalues(x = root2$tip.label, from = nameSwapSheet1$TreeNames, to = nameSwapSheet1$ColorNames)

#newSpp_forColorData <- plyr::mapvalues(totDatTibGrouped$sppNames, from = nameSwapSheet1$ColorNames, to = nameSwapSheet1$TreeNames)
#totDatTibGrouped_TipSpp <- add_column(totDatTibGrouped, newSpp_forColorData)

totDatTibGrouped_TipSpp <- totDatTibGrouped

#root$tip.label <- newLabels 

#get tree into right shape
graphics.off()
plot(root, cex = .2)
root <- ladderize(root, right = F)
root <- root(root, c("Caracara_cheriway", "Icterus_cucullatus", "Calyptomena_viridis"))
root$tip.label

is_tip <- root$edge[,2] <= length(root$tip.label)
ordered_tips <- root$edge[is_tip, 2]
root$tip.label[ordered_tips]

plot(root, cex = .2)
nodelabels()
rootRotated <- ape::rotate(phy = root, node = 387)
plot(ladderize(rootRotated), cex = .2)
rootRotated <- ladderize(rootRotated)

is_tip <- rootRotated$edge[,2] <= length(rootRotated$tip.label)
ordered_tips <- rootRotated$edge[is_tip, 2]
rootRotated$tip.label[ordered_tips]

#list giving the plotting order for measured patches
orderedPatchList <- c("Lores","SuperciliaryLine","Eyeline","AuricularFeathers","MalarRegion",
                      "LateralNeck","LateralBreast","Wrist","Shoulder","MarginalWingCoverts","LesserCoverts",
                      "MedianCoverts","GreaterSecondaryCoverts","FoldedSecondaries","FoldedPrimaries","UpperFlank","LowerFlank",
                      'ThroatSides',"UpperBreast","LowerBreast","SideBreast","UpperAbdomen","LowerAbdomen","CruralFeathers","Crissum",
                      "Collar","Forehead","Crown","Occiput","Nape","NucalCollar","Scapulars",'Back',"Rump","UpperTailCoverts")

# Function to find patch in string (expects only one patch per string)
find_patch <- function(text, patch_list) {
  if (is.na(text) || text == "") return(NA)
  
  # Create regex pattern that matches any of the patches
  pattern <- paste(patch_list, collapse = "|")
  match <- str_extract(text, pattern)
  return(match)
}

# Add the real_patch column, by querying the "label" column of the dataset for each patch.
totDatTibGrouped_TipSpp$real_patch <- sapply(totDatTibGrouped_TipSpp$Label, find_patch, patch_list = orderedPatchList)

#this function pulls the 35 measurements for each tip.
grabColorForTip_rework <- function(tipLab = NULL, UV_return = FALSE, real_dat_return = FALSE, sex_filter = NA)
{
  #pick the color dataset that we will pull color from
  colorDataset_toQuery <- totDatTibGrouped_TipSpp
  theseCols <- which(adist(colorDataset_toQuery$taxName, tipLab) <= 3)
  
  allow_partials=T
  print(length(theseCols))
  if((length(theseCols) < 1) && allow_partials){return(t(rep("transparent", 35)))}
  if((length(theseCols) < 35) && !allow_partials){return(t(rep("transparent", 35)))}
  
  #trying to clip out the super-bright white standard bad patches...
  tipColors <- colorDataset_toQuery[theseCols, ]
  
  bright_threshold <- .4
  uv_thresh <- .6
  
  tipColors_uv_filt <- tipColors[-c(which(tipColors$uvMean > uv_thresh)),]
  if(length(tipColors_uv_filt[,1])>5){tipColors <- tipColors_uv_filt}
  
  #edge case 
  cockatoo_genera <- c("Cacatua", "Callocephalon", "Eolophus", 
              "Lophochroa", "Nymphicus")
  
  if(str_detect(tipLab, paste(cockatoo_genera, collapse = "|"))){bright_threshold <- .9}

  #edge case for big birds
  
  bottom_rows <- grep(tipColors$Label, pattern = "_bot", fixed = F, ignore.case = T)
  top_rows <- grep(tipColors$Label, pattern = "_top", fixed = F, ignore.case = T)
  
  bottom_patches_acceptable <- c("UpperTailCoverts", "Rump", "Crissum", "CruralFeathers", "UpperAbdomen", "LowerAbdomen", "FoldedPrimaries", "FoldedSecondaries")

  which_bottom_patches <- which(apply(tipColors, 1, function(r) any(r %in% bottoms)))
  
  topDrops <- intersect(top_rows, which_bottom_patches)
  bottomDrops <- setdiff(bottom_rows, which_bottom_patches)
  
  unique(tipColors[bottomDrops,"real_patch"])
  unique(tipColors[topDrops,"real_patch"])

  print(paste("Dropping rows because of big birds:", length(c(bottomDrops, topDrops))))
  
  if(length(c(bottomDrops,topDrops))>0){
  tipColors <- tipColors[-c(bottomDrops, topDrops),]}
  
  print(paste("Remaining data rows: ", length(tipColors[,1])))
  
  hist((tipColors$dblMean))
  
  if(!is.na(sex_filter))
  {
    if((sex_filter == "M") && real_dat_return){
    tipColors <- tipColors %>% group_by(sex) %>% filter((sex != "F"))}
    if((sex_filter == "F") && real_dat_return){
      tipColors <- tipColors %>% group_by(sex) %>% filter((sex == "F"))}
    print(unique(tipColors$sex))
  }
  
  tipColors <- tipColors[which(tipColors$dblMean<bright_threshold),]
  mean_lum <- mean(tipColors$lumMean, trim = .25, na.rm = T)
  
  if(length(grep(tipColors$Label, pattern = "whole_image"))){print("Missing ROIS for: ");print(unique(tipColors$taxName));print(unique(tipColors$SpecimenNumber))}
  
  if(length(unique(tipColors$real_patch)) > 6)
  {
    summarize_at(group_by(tipColors, real_patch), .vars = c("lwMean", "swMean", "mwMean", "uvMean"), .funs = function(X) mean(X, trim=.25)) -> x
  } else {
    return(t(rep("transparent", 35)))
  }
  
  print(paste(tipLab))
  print(x)
  thisColMax = 1
  if((max(x[,2:5])) > 1){thisColMax <- (2*max(x[,2:5]))}
  x <- na.omit(x)
  
  #if lum adjust?
  lum_adjust = F
  if(lum_adjust){x[,2:5] <- x[,2:5]*mean_lum}
  
  rownames(x) <- x$real_patch
  available_patches <- intersect(orderedPatchList, x$real_patch)
  x <- x[available_patches,]
  hist(log(unlist(x[,2:5])), main = tipLab)
  
  returnTipColors <- data.frame(color = rep("transparent", 35))
  rownames(returnTipColors) <- orderedPatchList
  tipColorsSummarized <- rgb(red = log(1+abs(x$lwMean)), green =log(1+abs(x$mwMean)), blue = log(1+abs(x$swMean)), maxColorValue = thisColMax)
  
  if(real_dat_return){
    tipColorsSummarized_realdat <- data.frame(tipLab = tipLab, patch = paste(x$real_patch),l = log(1+abs(x$lwMean)), m =log(1+abs(x$mwMean)), s = log(1+abs(x$swMean)), u = log(1+abs(x$uvMean)), hex = tipColorsSummarized)
    return(tipColorsSummarized_realdat)
    }
  
  
  #turn only based on UV reflectance
  if(UV_return == T){
    print(x)
    xMean <- mean(rowMeans(x[,-c(1)]))
    xUVnegatives <- which(x$uvMean < 0);
    xUVmaxes <- which(x$uvMean > .3);
    x$uvMean[xUVnegatives] <- 0;
    x$uvMean[xUVmaxes] <- .3;
    tipColorsSummarized <- rgb(red = x$uvMean^2/xMean, green = 0, blue = x$uvMean^2/xMean, maxColorValue = .3^2/xMean)}
  
  
  tipColorsSummarized_df <- data.frame(tipColorsSummarized)
  rownames(tipColorsSummarized_df) <- available_patches
  
  returnTipColors[available_patches,] <- tipColorsSummarized_df[available_patches,1]
  
  #colSpaceCols
  
  return(t(as.data.frame(returnTipColors)[1:35,"color"]))
  
}

unique(totDatTibGrouped_TipSpp$sppNames)
unique(totDatTibGrouped_TipSpp$sppNames)

intersect(totDatTibGrouped_TipSpp$taxName, root$tip.label)

colorName_matchTest <- c(unique(totDatTibGrouped_TipSpp$taxName))
names(colorName_matchTest) <- gsub(x = colorName_matchTest, pattern = " ", replacement = "_")
testinMismatches <- match.phylo.data(phy = rootRotated, colorName_matchTest)
}

#make the tile matrices, grab patchwise summarized data...
{
#get correct order of tips for the clads rate figure to plot the tile matrix
orderedLabels_clads <- read.csv("data/tipOrder_ElkinCladsTree_7.14.2023.csv")

#grab the data for the tile matrix
colList <- data.frame(V2 = orderedLabels_clads$Spp_Name_noNum_spaces)
length(unique(colList$V2))

#this is where the color values for each species are pulled from the data
#should take about 5-10 minutes...
cc <- sapply(X = unique(colList$V2), simplify = T, FUN = function(X) grabColorForTip_rework(X))
cc<- t(cc)
ccDF <- data.frame(cc)
rownames(ccDF) <- gsub(colList$V2, pattern=" ", replacement = "_")

#fixing known naming issues
rootRotated$tip.label[grep(rootRotated$tip.label, pattern = "luchsi")] <- "Myiopsitta_fuchsi"
rootRotated$tip.label[grep(rootRotated$tip.label, pattern = "melanonotus")] <- "Touit_melanotus"
rootRotated$tip.label[grep(rootRotated$tip.label, pattern = "cornutus")] <- "Eunymphicus_cornotus"

b <- match.phylo.data(data =ccDF, phy = rootRotated)
b$data <-  ccDF

bb<- b$data[which(b$data[,5] != "transparent"),]

#plot the tile matrix as an image, with each pixel representing a patch for a single taxon.
##fullMatrix
bMatrix <- as.matrix(b$data)
allDatWithBlanks <- bMatrix
which(allDatWithBlanks[,1] == "transparent")
mxSize <- length(allDatWithBlanks[,1])*35
dummyMatrixBlanks <- matrix(1:mxSize, nrow = length(allDatWithBlanks[,1]), ncol = 35)
png(filename = "allColsRGB_Blanks.png", height = 35*100, width = length(allDatWithBlanks[,1])*100, units = "px", bg = "transparent")
image(dummyMatrixBlanks, col = lighten(allDatWithBlanks[,], amount = .15, method = "absolute", space = "HLS"), axes = F)
dev.off()
#open tile matrix
system("open allColsRGB_Blanks.png")

#now grab the raw data for males, females, and a combined dataset.
#this will take about 15 minutes to run the below
cc_realdat_nosex <- sapply(X = unique(colList$V2), simplify = T, FUN = function(X) grabColorForTip_rework(X, real_dat_return = T))
cc_realdat_males <- sapply(X = unique(colList$V2), simplify = T, FUN = function(X) grabColorForTip_rework(X, real_dat_return = T,sex_filter = "M"))
cc_realdat_females <- sapply(X = unique(colList$V2), simplify = T, FUN = function(X) grabColorForTip_rework(X, real_dat_return = T,sex_filter = "F"))

# Remove NULL, empty, and problematic data frames
clean_list_males <- cc_realdat_males %>%
  purrr::discard(is.null) %>%
  purrr::discard(~is.null(.x) || nrow(.x) == 0) %>%
  purrr::keep(is.data.frame)  # Ensure they're actually data frames

# If you want to bind them
cc_realdat_males_df <- clean_list_males %>% bind_rows(.id = "source")
cc_realdat_males_df$sex <- "M"

# Remove NULL, empty, and problematic data frames
clean_list_females <- cc_realdat_females %>%
  purrr::discard(is.null) %>%
  purrr::discard(~is.null(.x) || nrow(.x) == 0) %>%
  purrr::keep(is.data.frame)  # Ensure they're actually data frames

# If you want to bind them
cc_realdat_females_df <- clean_list_females %>% bind_rows(.id = "source")
cc_realdat_females_df$sex <- "F"

cc_all_dat_combined <- rbind.data.frame(cc_realdat_males_df, cc_realdat_females_df)
PC_result <-  prcomp(cc_all_dat_combined[,c("u","s","m","l")], retx = T, center = T, scale. = T)
cc_all_dat_combined$PC1 <- PC_result$x[,1]
cc_all_dat_combined$PC2 <- PC_result$x[,2]
cc_all_dat_combined$PC3 <- PC_result$x[,3]
cc_all_dat_combined$MW_UV_oppo <- (cc_all_dat_combined$m-cc_all_dat_combined$u)/(cc_all_dat_combined$m+cc_all_dat_combined$u+cc_all_dat_combined$l+cc_all_dat_combined$s)
cc_all_dat_combined$LW_SV_oppo <- (cc_all_dat_combined$l-cc_all_dat_combined$s)/(cc_all_dat_combined$m+cc_all_dat_combined$u+cc_all_dat_combined$l+cc_all_dat_combined$s)
cc_all_dat_combined$BRIGHTNESS <- cc_all_dat_combined$m+cc_all_dat_combined$u+cc_all_dat_combined$l+cc_all_dat_combined$s
}

#now map to body regions, views, make summary datasets for these regions
{
#group patches into views and body parts for export and plots
patch_groups <- data.frame(patch = c("Lores","SuperciliaryLine","Eyeline","AuricularFeathers","MalarRegion",
                                                 "LateralNeck","LateralBreast","Wrist","Shoulder","MarginalWingCoverts","LesserCoverts",
                                                 "MedianCoverts","GreaterSecondaryCoverts","FoldedSecondaries","FoldedPrimaries","UpperFlank","LowerFlank",
                                                 'ThroatSides',"UpperBreast","LowerBreast","SideBreast","UpperAbdomen","LowerAbdomen","CruralFeathers","Crissum",
                                                 "Collar","Forehead","Crown","Occiput","Nape","NucalCollar","Scapulars",'Back',"Rump","UpperTailCoverts"))

patch_groups$view <- c("ventral","ventral","ventral","ventral","ventral",
                       "ventral","ventral","dorsal","dorsal","dorsal","dorsal",
                       "dorsal","dorsal","dorsal","dorsal","ventral","ventral",
                       'ventral',"ventral","ventral","ventral","ventral","ventral","ventral","ventral",
                       "ventral","dorsal","dorsal","dorsal","dorsal","dorsal","dorsal",'dorsal',"dorsal","dorsal")

patch_groups$grouped_bodypart <- c("head","head","head","head","head",
                                   "head","head","wing","back","wing","wing",
                                   "wing","wing","wing","wing","body","body",
                                   'head',"body","body","body","body","body","body","body",
                                   "head","head","head","head","head","back","back",'back',"back","back")

library(ggpubr)

cc_all_dat_combined$view <- mapvalues(x=cc_all_dat_combined$patch, from = patch_groups$patch, to = patch_groups$view)
cc_all_dat_combined$body_part <- mapvalues(x=cc_all_dat_combined$patch, from = patch_groups$patch, to = patch_groups$grouped_bodypart)

cc_all_dat_body_part <- cc_all_dat_combined %>%
  filter(body_part == "head")

ggplot(cc_all_dat_body_part, aes(x = (MW_UV_oppo), y = (LW_SV_oppo))) + #xlim(c(-3,2.5)) + ylim(c(-2,2)) +
  geom_density_2d_filled(show.legend = F) +
  scale_fill_brewer(palette = "Blues") + theme_pubclean()# + xlab(label = "Precipitation Distance") + ylab("Phenotype Distance")  #+ geom_smooth(method = "lm", color = "lightblue", fill = "black", alpha = .1)

cc_all_dat_body_part <- cc_all_dat_combined %>%
  filter(body_part == "body")
ggplot(cc_all_dat_body_part, aes(x = MW_UV_oppo, y = LW_SV_oppo)) + #xlim(c(-3,2.5)) + ylim(c(-2,2)) +
  geom_density_2d_filled(show.legend = F) +
  scale_fill_brewer(palette = "Reds") + theme_pubclean()# + xlab(label = "Precipitation Distance") + ylab("Phenotype Distance")  #+ geom_smooth(method = "lm", color = "lightblue", fill = "black", alpha = .1)

cc_all_dat_body_part <- cc_all_dat_combined %>%
  filter(body_part == "wing")
ggplot(cc_all_dat_body_part, aes(x = MW_UV_oppo, y = LW_SV_oppo)) + #+ xlim(c(-3,2.5)) + ylim(c(-2,2)) +
  geom_density_2d_filled(show.legend = F) +
  scale_fill_brewer(palette = "Greens") + theme_pubclean()# + xlab(label = "Precipitation Distance") + ylab("Phenotype Distance")  #+ geom_smooth(method = "lm", color = "lightblue", fill = "black", alpha = .1)

ggplot(cc_all_dat_combined) + geom_boxplot(aes(group=body_part,x=(MW_UV_oppo), fill=body_part))
ggplot(cc_all_dat_combined) + geom_boxplot(aes(group=body_part,x=LW_SV_oppo, fill=body_part))
ggplot(cc_all_dat_combined) + geom_boxplot(aes(group=body_part,x=BRIGHTNESS, fill=body_part))

cc_all_dat_combined$hex_brightscaled <- lighten(cc_all_dat_combined$hex, amount = ((cc_all_dat_combined$BRIGHTNESS))/2, method = "relative", space = "HLS")

head_cols <- ggplot(cc_all_dat_combined%>%filter(body_part=="head")) + geom_point(aes(x=MW_UV_oppo, y = LW_SV_oppo, col = I(hex_brightscaled)),alpha = .01,size = 20, shape = 16) + theme_classic() + ggtitle("Head") + xlim(c(-.8, .8))+ ylim(c(-.8, .8))
body_cols <- ggplot(cc_all_dat_combined%>%filter(body_part=="body")) + geom_point(aes(x=MW_UV_oppo, y = LW_SV_oppo, col = I(hex_brightscaled)),alpha = .01,size = 20, shape = 16) + theme_classic() + ggtitle("Body") + xlim(c(-.8, .8))+ ylim(c(-.8, .8))
wing_cols <- ggplot(cc_all_dat_combined%>%filter(body_part=="wing")) + geom_point(aes(x=MW_UV_oppo, y = LW_SV_oppo, col = I(hex_brightscaled)),alpha = .01,size = 20, shape = 16) + theme_classic() + ggtitle("Wing") + xlim(c(-.8, .8))+ ylim(c(-.8, .8))
back_cols <- ggplot(cc_all_dat_combined%>%filter(body_part=="back")) + geom_point(aes(x=MW_UV_oppo, y = LW_SV_oppo, col = I(hex_brightscaled)),alpha = .01,size = 20, shape = 16) + theme_classic() + ggtitle("Back") + xlim(c(-.8, .8))+ ylim(c(-.8, .8))

ggpubr::ggarrange(head_cols,body_cols,wing_cols,back_cols)

cc_all_dat_combined$genus <- str_split_i(cc_all_dat_combined$source, pattern = " ",i = 1)

rand_genus <- sample(cc_all_dat_combined$genus)
rand_genus <- "Polytelis"
this_plot_dat <- cc_all_dat_combined %>% filter(genus == rand_genus)
ggplot(this_plot_dat) + geom_point(aes(y=MW_UV_oppo, x = LW_SV_oppo, col = I(hex_brightscaled)),alpha = 1,size = 10, shape = 16) + theme_classic() + ggtitle(rand_genus)
}

#make color plots for each genus...also runs a pairwise T-test to compare 
#brightness of dorsal and ventral (see Fig 1)
{
genus_plot_list <- list()
for(i in 1:length(unique(cc_all_dat_combined$genus)))
{
rand_genus <- unique(cc_all_dat_combined$genus)[i]
this_plot_dat <- cc_all_dat_combined %>% filter(genus == rand_genus)
this_genus_plot <- ggplot(this_plot_dat) + geom_point(aes(y=MW_UV_oppo, x = LW_SV_oppo, col = I(hex_brightscaled)),alpha = 1,size = 2, shape = 16) + theme_void() + ggtitle(rand_genus) +xlab(label = "") +ylab(label="") + xlim(c(-.6,.6)) + ylim(c(-.6,.6))
genus_plot_list[[i]] <- this_genus_plot
}
names(genus_plot_list) <- unique(cc_all_dat_combined$genus)

png("genus_spaces.png", width = 1000, height = 1000, units = "px")
ggarrange(plotlist =genus_plot_list, nrow=10,ncol=10)
graphics.off()
system("open genus_spaces.png")

patch_plot_list <- list()
for(i in 1:length(unique(cc_all_dat_combined$patch)))
{
  rand_patch <- unique(cc_all_dat_combined$patch)[i]
  this_patch_dat <- cc_all_dat_combined %>% filter(patch == rand_patch)
  this_patch_plot <- ggplot(this_patch_dat) + geom_point(aes(y=MW_UV_oppo, x = LW_SV_oppo, col = I(hex_brightscaled)),alpha = 1,size = 2, shape = 16) + theme_void() + ggtitle(rand_patch) +xlab(label = "") +ylab(label="")
  patch_plot_list[[i]] <- this_patch_plot
}
names(patch_plot_list) <- unique(cc_all_dat_combined$patch)

png("patch_spaces.png", width = 1400, height = 1000, units = "px")
ggarrange(plotlist =patch_plot_list, nrow=5,ncol=7)
graphics.off()
system("open patch_spaces.png")

ggplot(cc_all_dat_combined) + geom_point(aes(x=MW_UV_oppo, y = LW_SV_oppo, col = body_part),alpha = .01,size = 30) + theme_bw()

write.csv(file ="2_processing_color_and_morphology_data/data/datasets_to_collate/sexes_combined_raw_color_longways.csv", cc_all_dat_combined)

ggplot(cc_all_dat_combined, aes(x = factor(patch, level = orderedPatchList), y = (l-s)/(l+s), col = I(lighten(hex, .2, method = "relative", space = "HLS")))) + geom_quasirandom(width =.5, size = 1, alpha = 1) + theme_bw()

pdf("brightness_countershading.pdf", width = 8, height = 10)
ggplot(cc_all_dat_combined, aes(x = factor(view), y = log(BRIGHTNESS), col = I(hex_brightscaled))) + geom_quasirandom(width =.5, size = 1, alpha = 1) + theme_bw()
graphics.off()
system("open brightness_countershading.pdf")

ggplot(cc_all_dat_combined, aes(x = factor(body_part), y = (log(u+s+m+l)), col = I(lighten(hex, .3, method = "relative", space = "HLS")))) + geom_quasirandom(width =.5, size = 1, alpha = 1) + theme_bw()
ggplot(cc_all_dat_combined, aes(x = factor(view), y = (log(u/(s+m+l))), col = I(lighten(hex, .3, method = "relative", space = "HLS")))) + geom_quasirandom(width =.5, size = 1, alpha = 1) + theme_bw()
ggplot(cc_all_dat_combined, aes(x = factor(body_part), y = (log(u/(s+m+l))), col = I(lighten(hex, .3, method = "relative", space = "HLS")))) + geom_quasirandom(width =.5, size = 1, alpha = 1) + theme_bw()

pairwise_results <- pairwise.t.test(
  x = log(cc_all_dat_combined$BRIGHTNESS),
  g = factor(cc_all_dat_combined$view),
  p.adjust.method = "bonferroni"  # or "holm", "BH", etc.
)
}

#now we use PAVO to generate summary statistics for each species...
#data for use in rest of the analysis...
{
cc_all_colspace <- colspace(cc_all_dat_combined[,c("u","s","m","l")], space = "tcs", qcatch = "Qi")
plot(cc_all_colspace$x~cc_all_colspace$y, col = "gray", pch = 16, cex = .6, alpha = .2)
by_spp_summary <- summary(cc_all_colspace, by = cc_all_dat_combined$source)
by_spp_summary$spp <- rownames(by_spp_summary)

#calculate PAVO summary stats for dorsal
cc_vent_colspace <- colspace(cc_all_dat_combined[which(cc_all_dat_combined$view=="ventral"),c("u","s","m","l")], space = "tcs", qcatch = "Qi")
by_spp_summary_vent <- summary(cc_vent_colspace, by = cc_all_dat_combined[which(cc_all_dat_combined$view=="ventral"),"source"])
by_spp_summary_vent$spp <- rownames(by_spp_summary_vent)

#calculate PAVO summary stats for ventral
cc_dorsal_colspace <- colspace(cc_all_dat_combined[which(cc_all_dat_combined$view=="dorsal"),c("u","s","m","l")], space = "tcs", qcatch = "Qi")
by_spp_summary_dorsal <- summary(cc_vent_colspace, by = cc_all_dat_combined[which(cc_all_dat_combined$view=="dorsal"),"source"])
by_spp_summary_dorsal$spp <- rownames(by_spp_summary_dorsal)

#calculate PAVO summary stats for body regions
cc_head_colspace <- colspace(cc_all_dat_combined[which(cc_all_dat_combined$body_part=="head"),c("u","s","m","l")], space = "tcs", qcatch = "Qi")
by_spp_summary_head <- summary(cc_head_colspace, by = cc_all_dat_combined[which(cc_all_dat_combined$body_part=="head"),"source"])
by_spp_summary_head$spp <- rownames(by_spp_summary_head)

cc_body_colspace <- colspace(cc_all_dat_combined[which(cc_all_dat_combined$body_part=="body"),c("u","s","m","l")], space = "tcs", qcatch = "Qi")
by_spp_summary_body <- summary(cc_body_colspace, by = cc_all_dat_combined[which(cc_all_dat_combined$body_part=="body"),"source"])
by_spp_summary_body$spp <- rownames(by_spp_summary_body)

cc_wing_colspace <- colspace(cc_all_dat_combined[which(cc_all_dat_combined$body_part=="wing"),c("u","s","m","l")], space = "tcs", qcatch = "Qi")
by_spp_summary_wing <- summary(cc_wing_colspace, by = cc_all_dat_combined[which(cc_all_dat_combined$body_part=="wing"),"source"])
by_spp_summary_wing$spp <- rownames(by_spp_summary_wing)

cc_back_colspace <- colspace(cc_all_dat_combined[which(cc_all_dat_combined$body_part=="back"),c("u","s","m","l")], space = "tcs", qcatch = "Qi")
by_spp_summary_back <- summary(cc_back_colspace, by = cc_all_dat_combined[which(cc_all_dat_combined$body_part=="back"),"source"])
by_spp_summary_back$spp <- rownames(by_spp_summary_back)

# Add body_part identifier to each summary
by_spp_summary$body_part <- "whole_body"
by_spp_summary_head$body_part <- "head"
by_spp_summary_body$body_part <- "body"
by_spp_summary_wing$body_part <- "wing"
by_spp_summary_back$body_part <- "back"

# Combine all summaries
all_summaries <- bind_rows(
  by_spp_summary,
  by_spp_summary_head,
  by_spp_summary_body,
  by_spp_summary_wing,
  by_spp_summary_back
)

# Pivot to wide format: one row per species, body parts as column prefixes
combined_wide_bodypart_PAVO <- all_summaries %>%
  pivot_wider(
    id_cols = spp,
    names_from = body_part,
    values_from = -c(spp, body_part),  # All columns except spp and body_part
    names_sep = "_"
  )

ggplot(cc_all_colspace, aes(x = log(u/m), y = log(l/s), col = I(rgb(l,m,s)))) + geom_point(alpha = 1, size = 1)

#save those datasets as CSV
 write.csv(file ="2_processing_color_and_morphology_data/data/datasets_to_collate/sexes_combined_spp_color_spaces.csv", by_spp_summary)
 write.csv(file ="2_processing_color_and_morphology_data/data/datasets_to_collate/sexes_combined_spp_color_spaces_ventral.csv", by_spp_summary_vent)
 write.csv(file ="2_processing_color_and_morphology_data/data/datasets_to_collate/sexes_combined_spp_color_spaces_dorsal.csv", by_spp_summary_dorsal)
 write.csv(file ="2_processing_color_and_morphology_data/data/datasets_to_collate/sexes_combined_spp_color_spaces_back.csv", by_spp_summary_back)
 write.csv(file ="2_processing_color_and_morphology_data/data/datasets_to_collate/sexes_combined_spp_color_spaces_body.csv", by_spp_summary_body)
 write.csv(file ="2_processing_color_and_morphology_data/data/datasets_to_collate/sexes_combined_spp_color_spaces_wing.csv", by_spp_summary_wing)
 write.csv(file ="2_processing_color_and_morphology_data/data/datasets_to_collate/sexes_combined_spp_color_spaces_head.csv", by_spp_summary_head)
 write.csv(file ="2_processing_color_and_morphology_data/data/datasets_to_collate/sexes_combined_spp_color_spaces_chunked.csv", combined_wide_bodypart_PAVO)
 
col_summary_pca <- prcomp(by_spp_summary[,1:13], center = T, scale. = T)
biplot(col_summary_pca)
dotchart(col_summary_pca$rotation[,1:3], pch = 15, cex = 1) 
by_spp_summary[,c("PAVO_SUMM_PC1", "PAVO_SUMM_PC2", "PAVO_SUMM_PC3")] <- col_summary_pca$x[,1:3]

pdf("pavo_summary_pca_plot.pdf", width = 30, height = 30)
ggplot(by_spp_summary)+ geom_label(aes(x=(mean.ra),y=(a.vol),label = spp)) + theme_bw()
graphics.off()

savenames <- combined_wide_bodypart_PAVO$spp
rownames(combined_wide_bodypart_PAVO) <- savenames
combined_wide_bodypart_PAVO<-data.frame(combined_wide_bodypart_PAVO[,-c(1)])
rownames(combined_wide_bodypart_PAVO) <- savenames
col_summary_pca <- prcomp(na.omit(combined_wide_bodypart_PAVO), center = T, scale. = T)

ggplot(all_summaries) + geom_point()

View(col_summary_pca$x)
dotchart(col_summary_pca$rotation[,1:3], pch = 15, cex = 1) 
by_spp_summary[,c("PAVO_SUMM_PC1", "PAVO_SUMM_PC2", "PAVO_SUMM_PC3")] <- col_summary_pca$x[,1:3]

pdf("pavo_summary_pca_plot.pdf", width = 30, height = 30)
ggplot(by_spp_summary)+ geom_label(aes(label = spp)) + theme_bw()
graphics.off()

system("open pavo_summary_pca_plot.pdf")
diagnostic_data <-cc_realdat$`Eclectus roratus`
pavo_hex <- rgb(pavo::colspace(diagnostic_data[,c("u","s","m","l")], space = "tcs", qcatch = "Qi")[,c("l","m","s")])
pavo::tcsplot(pavo::colspace(diagnostic_data[,c("u","s","m","l")], space = "tcs", qcatch = "Qi"), col = lighten(diagnostic_data$hex, .2, method = "relative", space = "HLS"))
pavo::tcsplot(pavo::colspace(diagnostic_data[,c("u","s","m","l")], space = "tcs", qcatch = "Qi"), col = pavo_hex)
}















