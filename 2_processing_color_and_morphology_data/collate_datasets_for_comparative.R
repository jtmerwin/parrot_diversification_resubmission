#make comparative dataset with color,morphology,niche,habitat,clads rate, DR stat, named with tree names

#load dependencies
{
library(ape)
library(geiger)
library(phytools)
library(caper)
library(picante)
library(plyr)
library(dplyr)
library(stringdist)
library(ggplot2)
library(tidyr)
library(corrplot)

#set working directory
data_pointer <- "c:/Users/jonme/AMNH Dropbox/Jon Merwin/AllParrotHomeAnalysis/Scripts/GITHUB_UPLOAD/"
#data_pointer <- "c:/Users/jtm378/AMNH Dropbox/Jon Merwin/AllParrotHomeAnalysis/Scripts/GITHUB_UPLOAD/"

setwd(data_pointer)

#lastDropper
lastDropper <- function(string = NULL, seper = " ")
{
  require("stringr")
  lab <- word(string, 1, -2, sep = seper)
  return(lab)
}
}

#load datasets
{
#load taxonomy key
tip_hotswap <- read.csv("2_processing_color_and_morphology_data/data/datasets_to_collate/tip_hotswap_7.16.25.csv")

#load color
color_collate <- read.csv("2_processing_color_and_morphology_data/data/datasets_to_collate/sexes_combined_spp_color_spaces.csv")
rownames(color_collate) <- color_collate$spp

color_collate_dorsal <- read.csv("2_processing_color_and_morphology_data/data/datasets_to_collate/sexes_combined_spp_color_spaces_dorsal.csv")
rownames(color_collate_dorsal) <- color_collate_dorsal$spp

color_collate_ventral <- read.csv("2_processing_color_and_morphology_data/data/datasets_to_collate/sexes_combined_spp_color_spaces_ventral.csv")
rownames(color_collate_ventral) <- color_collate_ventral$spp

color_collate_chunkedPAVO <- read.csv("2_processing_color_and_morphology_data/data/datasets_to_collate/sexes_combined_spp_color_spaces_chunked.csv")
rownames(color_collate_chunkedPAVO) <- color_collate_chunkedPAVO$X

#load raw color measures 
color_raw <- read.csv("2_processing_color_and_morphology_data/data/datasets_to_collate/sexes_combined_raw_color_longways.csv")

summary_by_patch <- color_raw %>%
  group_by(tipLab, patch) %>%
  summarise(
    n = n(),
    L_S = mean(LW_SV_oppo, na.rm = TRUE),
    M_U = mean(MW_UV_oppo, na.rm = TRUE),
    BRIGHT = mean(BRIGHTNESS, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(summary_by_patch,aes(y=(L_S), x=(M_U),fill = patch,color=patch)) + geom_point(size = 1,alpha=1)

# Pivot wider: each species gets one row, patches become column prefixes
raw_patchwise_wide <- summary_by_patch %>%
  pivot_wider(
    id_cols = tipLab,
    names_from = patch,
    values_from = c(L_S, M_U, BRIGHT),
    names_sep = "_"
  )

rownames(raw_patchwise_wide) <- raw_patchwise_wide$tipLab
names(raw_patchwise_wide)

c(2:36,72:(72+34),142:(142+34))

corrplot(cor(na.omit(raw_patchwise_wide[,c(2:36)])), method = "square", order = "hclust", diag = F, addrect = 5)

#load morphology
morph_collate <- read.csv("2_processing_color_and_morphology_data/data/datasets_to_collate/Morph_collate.csv")
morph_collate_sppsummed <- data.frame(morph_collate %>%
  group_by(Species) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)))
rownames(morph_collate_sppsummed) <- morph_collate_sppsummed$Species   

morph_collate_sppsummed$TAIL_SUM <- as.numeric(paste(morph_collate_sppsummed$TAIL))+as.numeric(paste(morph_collate_sppsummed$TailOrnament.))
morph_collate_sppsummed$KIPPS_RATIO <- NA
morph_collate_sppsummed$KIPPS_RATIO <- (morph_collate_sppsummed$WING.LENGTH.TO.PRIMARIES-morph_collate_sppsummed$WING.LENGTH.TO.SECONDARIES)/morph_collate_sppsummed$WING.LENGTH.TO.PRIMARIES

#load niche and clads
niche_collate <- read.delim("2_processing_color_and_morphology_data/data/datasets_to_collate/Climate_collate.txt")
rownames(niche_collate) <- lastDropper(niche_collate$Sp, seper = "_")  

#load burgio habitat
burgio_collate <- read.csv("2_processing_color_and_morphology_data/data/datasets_to_collate/burgioCanopyScoring_collate.csv")
rownames(burgio_collate) <- make.names(gsub(burgio_collate$TreeNames, pattern = " ", replacement = "_"), unique = T)
burgio_collate <- burgio_collate[,-c(2)]

#load biogeography
load_biogeo <- read.csv("2_processing_color_and_morphology_data/data/datasets_to_collate/Geo_area_collate.csv")
rownames(load_biogeo) <- make.names(load_biogeo$TaxonNames.9.13.22Update, unique = T)
merge_biogeo <- load_biogeo[,c("TaxonNames.9.13.22Update", "Continent", "Biogeographic.Region", "graze")]

# #load updated biogeo for musscrat
# load_musscrat_biogeo <- read.nexus.data("2_processing_color_and_morphology_data/data/datasets_to_collate/parrot_continents_musscrat.nex")
# df_musscrat_biogeo <- data.frame(t(as.data.frame(load_musscrat_biogeo)))
# names(df_musscrat_biogeo) <- "musscrat_biogeo"

 #load updated biogeo for musscrat 9.11.2025
 load_musscrat_biogeo <- read.csv("2_processing_color_and_morphology_data/data/datasets_to_collate/10-Sept-2025-parrot-pheno-data.csv")
 df_musscrat_biogeo <- data.frame((as.data.frame(load_musscrat_biogeo["New.Biogeographic.Region"])))
 rownames(df_musscrat_biogeo) <- (load_musscrat_biogeo$species)
 names(df_musscrat_biogeo) <- "musscrat_biogeo"

#load tree
load_timetree <- read.nexus("2_processing_color_and_morphology_data/data/datasets_to_collate/parrot_timetree.tre")
load_timetree$tip.label <- lastDropper(load_timetree$tip.label, "_")
rooted_timetree <- ape::root(load_timetree, outgroup = c("Caracara_cheriway", "Icterus_cucullatus", "Calyptomena_viridis"))
#rooted_timetree$tip.label[which(rooted_timetree$tip.label=="Amazona_tresmariae")] <- "Amazona_oratrix"

#load bayestraits rates
load_bayesrates <- read.csv("2_processing_color_and_morphology_data/data/datasets_to_collate/species_tip_rates_REVBAYES_BRIAN.csv")
rownames(load_bayesrates) <- load_bayesrates$species
}

#set up the merge, custom function sourcing
{
# Create master dataset by merging all data sources
# Start with tree tips as the base
master_dataset <- data.frame(
  species = rooted_timetree$tip.label,
  stringsAsFactors = FALSE
)
rownames(master_dataset) <- make.names(rooted_timetree$tip.label, unique = T)

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

#merge datasets, save as an RDS (tipmatched_parrotdata_obj.rds)
{
# Merge color data
master_dataset <- safe_merge(master_dataset, color_collate, append = "_wholebody", max_dist = 2)
master_dataset <- safe_merge(master_dataset, color_collate_dorsal, append = "_dorsal")
master_dataset <- safe_merge(master_dataset, color_collate_ventral, append = "_ventral")
#379

# Merge morphology data
master_dataset <- safe_merge(master_dataset, morph_collate_sppsummed)
#379 matched

# Merge niche/climate data
master_dataset <- safe_merge(master_dataset, niche_collate)
#Matched 379 out of 384 species

# Merge habitat data
master_dataset <- safe_merge(master_dataset, burgio_collate)
#matched 363

master_dataset <- safe_merge(master_dataset, merge_biogeo)
#matched 382

setdiff(rownames(df_musscrat_biogeo), master_dataset$species)
setdiff(master_dataset$species, rownames(df_musscrat_biogeo))

master_dataset <- safe_merge(master_dataset, df_musscrat_biogeo)
#matched 382

master_dataset <- safe_merge(master_dataset,raw_patchwise_wide)

master_dataset <- safe_merge(master_dataset,color_collate_chunkedPAVO)

master_dataset <- safe_merge(master_dataset,load_bayesrates,append = "_bayes")



# Optional: Use tip_hotswap for name conversions if needed
# This depends on the structure of the tip_hotswap file
# You might need to adjust species names before merging
rownames(master_dataset) <- master_dataset$species

# Summary statistics
tipmatched_superset <- match.phylo.data(phy = rooted_timetree, data = master_dataset)

# Save any R object (preserves exact structure)
saveRDS(tipmatched_superset, "data/tipmatched_parrotdata_obj.rds")
}

#now below, we are running a phyloPCA to attatch to the dataset...
{
#If starting from here, load back
tipmatched_superset <- readRDS("data/tipmatched_parrotdata_obj.rds")

#phypca for morphology
colnames(tipmatched_superset$data)

morph_variables <- c("BILL.LENGTH","BILL.WIDTH","BILL.HEIGHT","TARSUS.LENGTH","WING.LENGTH.TO.PRIMARIES",
                     "WING.LENGTH.TO.SECONDARIES","TAIL","KIPPS_RATIO")

#option to exclude kipps...
#morph_variables <- c("BILL.LENGTH","BILL.WIDTH","BILL.HEIGHT","TARSUS.LENGTH","WING.LENGTH.TO.PRIMARIES",
#                     "WING.LENGTH.TO.SECONDARIES","TAIL")

comp_morph <- (tipmatched_superset$data[,c(morph_variables)])
bill_morph <- (tipmatched_superset$data[,c(morph_variables[1:3])])
bill_morph$BILL.VOL <- bill_morph$BILL.LENGTH

matched_pca_morph <- match.phylo.data(tipmatched_superset$phy, data = na.omit(comp_morph))
matched_pca_bill <- match.phylo.data(tipmatched_superset$phy, data = na.omit(bill_morph))

# Save row names if they exist
matched_pca_morph_rownames <- rownames(matched_pca_morph$data)
matched_pca_bill_rownames <- rownames(matched_pca_bill$data)
# Convert to numeric
matched_pca_morph$data <- data.frame(lapply(matched_pca_morph$data, as.numeric))
matched_pca_bill$data <- data.frame(lapply(matched_pca_bill$data, as.numeric))

# Restore row names
rownames(matched_pca_morph$data) <- matched_pca_morph_rownames
rownames(matched_pca_bill$data) <- matched_pca_bill_rownames

#run phylo pca for morphology, and for bill morphology specifically
phylo_pca_morph <- phytools::phyl.pca(tree = matched_pca_morph$phy, Y = matched_pca_morph$data, method = "lambda")
phylo_pca_bill <- phytools::phyl.pca(tree = matched_pca_bill$phy, Y = matched_pca_bill$data[,1:3], method = "lambda")

summary(phylo_pca_morph)
summary(phylo_pca_bill)

dotchart(phylo_pca_morph$L[,1])
dotchart(phylo_pca_morph$L[,2])
dotchart(phylo_pca_morph$L[,3])
dotchart(phylo_pca_morph$L[,4])
dotchart(phylo_pca_morph$L[,5])

PCA_for_merge <- data.frame(phylo_pca_morph$S[,1:3])
colnames(PCA_for_merge) <- c("PC1", "PC2", "PC3")

master_dataset <- safe_merge(master_dataset, PCA_for_merge)

#phypca for color
colnames(tipmatched_superset$data)[3:15]

color_variables <- c("centroid.u","centroid.l","rel.c.vol","colspan.m", "huedisp.m",  "mean.ra",  "a.vol")  

comp_color <- (tipmatched_superset$data[,c(color_variables)])


matched_pca_color <- match.phylo.data(tipmatched_superset$phy, data = data.frame(na.omit(comp_color)))

# Save row names if they exist
matched_pca_color_rownames <- rownames(matched_pca_color$data)
# Convert to numeric
matched_pca_color$data <- data.frame(lapply(matched_pca_color$data, as.numeric))
# Restore row names
rownames(matched_pca_color$data) <- matched_pca_color_rownames

phylo_pca_color <- phytools::phyl.pca(tree = matched_pca_color$phy, Y = matched_pca_color$data, method = "lambda")

summary(phylo_pca_color)

dotchart(phylo_pca_color$L[,1])
dotchart(phylo_pca_color$L[,2])
dotchart(phylo_pca_color$L[,3])

color_PCA_for_merge <- data.frame(phylo_pca_color$S[,1:3])
colnames(color_PCA_for_merge) <- c("color_PC1", "color_PC2", "color_PC3")

master_dataset <- safe_merge(master_dataset, color_PCA_for_merge)


#now do the PCA for climate
clim_pca_data <- na.omit(master_dataset[,c("XCentroid0.0","yCentroid0.0")])
save_names <- rownames(clim_pca_data)
clim_pca_data <- apply(clim_pca_data, MARGIN = 2, FUN = function(X) as.numeric(X))
rownames(clim_pca_data) <- save_names
clim_pca <- prcomp(clim_pca_data, center = T, scale. = T)
biplot(clim_pca, pc.biplot = T)
clim_pca_df <- data.frame(clim_pca$x[,1:2])
clim_pca_df$species <- rownames(clim_pca_df)
clim_pca_df$genus <- str_split_fixed(string = clim_pca_df$species, pattern = "_", n = 2)[,1]
ggplot(clim_pca_df, aes(x=PC1, y=PC2))  + geom_label(aes(label=species), size = 3)  +theme_bw()
hist(clim_pca_df$PC1)
hist(clim_pca_df$PC2) 
clim_pca$rotation

parrot_rds_climdat <- safe_merge(master_dataset, clim_pca_df, append = "_climate")

master_dataset_tipmatched <- match.phylo.data(data = parrot_rds_climdat, phy=tipmatched_superset$phy)

#raw_color_pc

# Save any R object (preserves exact structure)
saveRDS(master_dataset_tipmatched, "data/tipmatched_parrotdata_obj_with_phyPCA.rds")

# Read it back
tipmatched_superset_with_PC_matched <- readRDS("data/tipmatched_parrotdata_obj_with_phyPCA.rds")
tipmatched_superset_with_PC <- tipmatched_superset_with_PC_matched$data

pdf("test_morph_pc_space.pdf")
ggplot(tipmatched_superset_with_PC) + theme_bw() +
  geom_label(aes(label=species, y = sqrt(log(abs(400)-as.numeric(PC1))), x = (as.numeric(PC2))), size = 1)
graphics.off()
system("open test_morph_pc_space.pdf")

pdf("test_color_pc_space.pdf", width = 30, height = 30)
ggplot(tipmatched_superset_with_PC, aes(y = (as.numeric(color_PC1)), x = (as.numeric(color_PC2)))) + theme_bw() +
  #geom_image(aes(image = jpg_location), size = 0.01, alpha = 1) 
  geom_label(aes(label=species,y = (as.numeric(color_PC2)), x = (as.numeric(color_PC3))), size = 1)
graphics.off()
system("open test_color_pc_space.pdf")

pdf("test_clim_pc_space.pdf", width = 30, height = 30)
ggplot(tipmatched_superset_with_PC, aes(y = (as.numeric(PC2_climate)), x = (as.numeric(PC1_climate)))) + theme_bw() +
  #geom_image(aes(image = jpg_location), size = 0.01, alpha = 1) 
  geom_label(aes(label=species,y = (as.numeric(PC2_climate)), x = (as.numeric((PC1_climate)))), size = 1)
graphics.off()
system("open test_clim_pc_space.pdf")

ggplot(tipmatched_superset_with_PC, aes(y = (as.numeric(PC2_climate)), x = (as.numeric(PC1_climate)))) + theme_bw() + geom_point(size = 3, aes(col=Continent))
}

#below here writes data in a format for musscrat state dependent analysis
#variables included:
#  "TreeNames"            "graze"                "DR"                  
# "PC1"                  "PC2"                  "PC3"                 
# "Average.of.huedisp.m" "Average.of.mean.ra"   "mean.Dichrom"        
# "sd.Dichrom"           "Continent"            "Biogeographic.Region"
# "XCentroid0.5.y"       "yCentroid0.5.y"    
#data for musccrat
{
musscrat_data <- tipmatched_superset_with_PC[,c("species", "graze", "DR", "ClaDS", "PC1", "PC2", "PC3", "huedisp.m", "mean.ra", "a.vol", "huedisp.m_dorsal", "mean.ra_dorsal", "a.vol_dorsal","huedisp.m_ventral", "mean.ra_ventral", "a.vol_ventral", "Continent", "musscrat_biogeo", "XCentroid0.0", "yCentroid0.0","PC1_climate","PC2_climate","nectar","arboreal","ground.dwelling")]
musscrat_tree <- tipmatched_superset$phy

musscrat_matched <- match.phylo.data(data = musscrat_data, phy = musscrat_tree)

# Save any R object (preserves exact structure)
saveRDS(musscrat_matched, "data/musscrat_matched.rdj")
musscrat_matched <- readRDS("data/musscrat_matched.rdj")

hist(log(as.numeric(musscrat_matched$data$a.vol)))
musscrat_matched$data$genus <- str_split(musscrat_matched$data$species, pattern = "_", n = 2, simplify = T)[,1]
}

#script to generate scatterplots colored by habitat
#also does scatterplots colored by continents
#with marginal histograms showing distributions
{
#we need the new biogeobears areas
dec.25_biogeo_areas <- read.csv("transition_matrix_untrimmed.csv")
merged_biogeo_areas <- merge.data.frame(x = musscrat_matched$data, y = dec.25_biogeo_areas[,c("node","to_map")], all.x = T, all.y=F, by.x="species",by.y="node")

library(readODS)
#we need the new habitat designation
dec.25_habitat_score <- read_ods("8_ancestral_areas_biogeobears/data/3-nov-2025-parrot-continents.ods")
dec.25_habitat_score$PC1_new_trans <- sqrt(max(na.omit(as.numeric(dec.25_habitat_score$PC1_new)))-as.numeric(dec.25_habitat_score$PC1_new))

ggplot(dec.25_habitat_score, aes(x=sqrt(max(na.omit(as.numeric(PC1_new)))-as.numeric(PC1_new)),y=as.numeric(PC2_new))) + geom_point(size = 2, alpha =1,aes(col=factor(habitat, levels = c(0,1,2)))) + theme_classic() #+ geom_text(aes(label=species)) + geom_density()

ggplot(dec.25_habitat_score, aes(x=sqrt(max(na.omit(as.numeric(PC1_new)))-as.numeric(PC1_new)))) + 
  geom_histogram(alpha=1,
                 aes(fill=factor(habitat, levels = c(0,1,2))),  # Reversed
                 position = "dodge") + 
  theme_classic()

ggplot(dec.25_habitat_score, aes(x=as.numeric(PC2_new))) + 
  geom_histogram(alpha=1,
                 aes(fill=factor(habitat, levels = c(0,1,2))),  # Reversed
                 position = "dodge") + 
  theme_classic()

library(ggplot2)
library(ggExtra)

# Create the morph scatterplot
p <- ggplot(dec.25_habitat_score, 
            aes(x=sqrt(max(na.omit(as.numeric(PC1_new)))-as.numeric(PC1_new)),
                y=as.numeric(PC2_new))) + 
  # geom_point(stroke = .8, size = 1, alpha = 1, pch = 1, 
  #            aes(col=factor(habitat, levels = c(0,1,2)))) +
  geom_point(size = 1, alpha = .5, pch = 16, 
             aes(col=factor(habitat, levels = c(0,1,2)))) +
  scale_color_manual(values = c("purple", "skyblue2", "gold2")) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Body Size (Phylogenetic PC1)") +
  ylab("Body Shape (Phylogenetic PC2)")

# Add marginal histograms
morph_pc <- ggMarginal(p, type = "histogram", 
           color = "transparent", 
           groupFill = T,
           position = "dodge", alpha = 1, lwd = 1, size = 3, bins = 15)

ggMarginal(p, type = "density", 
           fill = "transparent", 
           groupColour = T,
           position = "dodge", alpha = .4, lwd = 1, size = 4)

# Now do niche plots?
p <- ggplot(dec.25_habitat_score, 
            aes(x=as.numeric(XCentroid0.5),
                y=as.numeric(yCentroid0.5))) + 
  # geom_point(stroke = .8, size = 1, alpha = 1, pch = 1, 
  #            aes(col=factor(habitat, levels = c(0,1,2)))) +
  geom_point(size = 1, alpha = .5, pch = 16, 
             aes(col=factor(habitat, levels = c(0,1,2)))) +
  scale_color_manual(values = c("purple", "skyblue2", "gold2")) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Temperature (Niche PC1)") +
  ylab("Precipitation (Niche PC2)")

# Add marginal histograms
niche_pc <- ggMarginal(p, type = "histogram", 
           color = "transparent", 
           groupFill = T,
           position = "dodge", alpha = 1, lwd = 1, size = 3, bins = 15)

ggMarginal(p, type = "density", 
           fill = "transparent", 
           groupColour = T,
           position = "dodge", alpha = .4, lwd = 1, size = 3)

# Now do color plots
p <- ggplot(dec.25_habitat_score, 
            aes(x=as.numeric(huedisp.m),
                y=as.numeric(mean.ra))) + 
  #geom_point(stroke = .8, size = 1, alpha = 1, pch = 1, 
  #           aes(col=factor(habitat, levels = c(0,1,2)))) +
  geom_point(size = 1, alpha = .5, pch = 16, 
             aes(col=factor(habitat, levels = c(0,1,2)))) +
  scale_color_manual(values = c("purple", "skyblue2", "gold2")) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Hue Disparity") +
  ylab("Mean Saturation") #+ 
  #geom_text(aes(label = factor(habitat, levels = c(0,1,2))))

# Add marginal histograms
color_pc <- ggMarginal(p, type = "histogram", 
           color = "transparent", 
           groupFill = T,
           position = "dodge", alpha = 1, lwd = 1, size = 3, bins = 15)

ggMarginal(p, type = "density", 
           fill = "transparent", 
           groupColour = T,
           position = "dodge", alpha = .4, lwd = 1, size = 8)

library(gridExtra)

pdf("pc_plots_histograms_axes.pdf", width = 9, height = 3)
grid.arrange(morph_pc, color_pc, niche_pc, ncol = 3, nrow = 1)
graphics.off()
system("open pc_plots_histograms_axes.pdf")
system("open ../")

ggplot(dec.25_habitat_score, aes(x=as.numeric(huedisp.m),y=as.numeric(mean.ra))) + geom_point(size = 6, alpha =1,aes(col=factor(habitat, levels = c(2,1,0)))) + theme_bw() #+ geom_text(aes(label=species))

# Now do color plots
p <- ggplot(dec.25_habitat_score[which(!is.na(dec.25_habitat_score$New.Continent)),], 
            aes(x=(as.numeric(huedisp.m)),
                y=as.numeric(mean.ra))) + 
  geom_point(stroke = .8, size = 3, alpha = 1, pch = 1, 
             aes(col=factor(New.Continent))) +
  geom_point(size = 3, alpha = .8, pch = 16, 
             aes(col=factor(New.Continent))) +
  #scale_color_manual(values=c("#5aba58",  "#64cdda",, "#3557a0", "#2e9ccc",  "#64cdda", "#2a8270", "#5aba58", "#f58643", "#bc529b", "#ea2636") +
  scale_color_manual(values = c("#5aba58", "#bc529b","#3557a0","#2e9ccc", "#f58643","#64cdda")) +
  theme_bw() +
  theme(legend.position = "off") +
  xlab("Hue Disparity") +
  ylab("Mean Saturation")

# Add marginal histograms
continent_color <- ggMarginal(p, type = "histogram", 
                              color = "transparent", 
                              groupFill = T, alpha = 1, size = 4,position = "dodge",bins = 10)

ggMarginal(p, type = "density", 
           fill = "transparent", 
           groupColour = T,
          alpha = 1, lwd = 1, size = 2)

# Now do morph plots
p <- ggplot(dec.25_habitat_score[which(!is.na(dec.25_habitat_score$New.Continent)),], 
            aes(x=(as.numeric(PC1_new_trans)),
                y=as.numeric(PC2_new))) + 
  geom_point(stroke = .8, size = 3, alpha = 1, pch = 1, 
             aes(col=factor(New.Continent))) +
  geom_point(size = 3, alpha = .8, pch = 16, 
             aes(col=factor(New.Continent))) +
  scale_color_manual(values = c("#5aba58", "#bc529b","#3557a0","#2e9ccc", "#f58643","#64cdda")) +
  #scale_color_manual(values = c("green4", "skyblue", "orange2", "darkblue", "skyblue4", "maroon")) +
  theme_bw() +
  theme(legend.position = "off") +
  xlab("Body Size (Phylogenetic PC1)") +
  ylab("Body Shape (Phylogenetic PC2)")

# Add marginal histograms
continent_morph <- ggMarginal(p, type = "histogram", 
           color = "transparent", 
           groupFill = T, alpha = 1, size = 4,position = "dodge",bins = 10)

ggMarginal(p, type = "density", 
           fill = "transparent", 
           groupColour = T,
           
           alpha = 1, lwd = 2, size = 4)

grid.arrange(ncol = 2, nrow = 1, continent_morph, continent_color)
}

#below generates more simple plots of hulls for genera and continents for morph and color
#genus-level morph PC plot convex hulls
{
library(ggplot2)
library(dplyr)

# First, clean and transform the data
clean_data <- musscrat_matched$data %>%
  mutate(
    PC1_num = as.numeric(PC1),
    PC2_num = as.numeric(PC2),
    genus = as.factor(genus)
  ) %>%
  # Remove rows with NA values in PC1 or PC2
  filter(!is.na(PC1_num), !is.na(PC2_num)) %>%
  mutate(
    x = log(250 - PC1_num),
    y = sqrt(150 - PC2_num)
  ) %>%
  # Remove rows where transformations create invalid values
  filter(is.finite(x), is.finite(y)) %>%
  # Also check that we don't have negative values under sqrt or log
  filter(250 - PC1_num > 0, 150 - PC2_num >= 0)

# Check if we have data left and enough points per group
clean_data %>% 
  group_by(genus) %>% 
  summarise(n = n()) %>%
  print()

# Compute convex hulls (need at least 3 points per group)
hull_data <- clean_data %>%
  group_by(genus) %>%
  filter(n() >= 3) %>%  # Need at least 3 points for a polygon
  slice(chull(x, y)) %>%
  ungroup()

# Option 5: Using distinct colors optimized for many categories
library(randomcoloR)
colors_100 <- distinctColorPalette(100)

ggplot() +
  geom_polygon(data = hull_data, 
               aes(x = x, y = y, color = genus, fill = genus), 
               size = 1, alpha = .1) +
  scale_fill_manual(values = colors_100) +
  scale_color_manual(values = colors_100) +
  theme_bw()
}

#continent-level morph PC plot convex hulls
{
  library(ggplot2)
  library(dplyr)
  
  # First, clean and transform the data
  clean_data <- musscrat_matched$data %>%
    mutate(
      PC1_num = as.numeric(PC1),
      PC2_num = as.numeric(PC2),
      genus = as.factor(musscrat_biogeo)
    ) %>%
    # Remove rows with NA values in PC1 or PC2
    filter(!is.na(PC1_num), !is.na(PC2_num)) %>%
    mutate(
      x = log(250 - PC1_num),
      y = sqrt(150 - PC2_num)
    ) %>%
    # Remove rows where transformations create invalid values
    filter(is.finite(x), is.finite(y)) %>%
    # Also check that we don't have negative values under sqrt or log
    filter(250 - PC1_num > 0, 150 - PC2_num >= 0)
  
  # Check if we have data left and enough points per group
  clean_data %>% 
    group_by(musscrat_biogeo) %>% 
    summarise(n = n()) %>%
    print()
  
  # Compute convex hulls (need at least 3 points per group)
  hull_data <- clean_data %>%
    group_by(musscrat_biogeo) %>%
    filter(n() >= 3) %>%  # Need at least 3 points for a polygon
    slice(chull(x, y)) %>%
    ungroup()
  
  # Option 5: Using distinct colors optimized for many categories
  library(randomcoloR)
  colors_100 <- distinctColorPalette(100)
  
  ggplot() +
    geom_polygon(data = hull_data, 
                 aes(x = x, y = y, color = musscrat_biogeo, fill = musscrat_biogeo), 
                 size = 1, alpha = .1) +
    scale_fill_manual(values = colors_100) +
    scale_color_manual(values = colors_100) +
    theme_bw()
}

#genus-level color PC plot convex hulls
{
  library(ggplot2)
  library(dplyr)
  
  # First, clean and transform the data
  clean_data <- musscrat_matched$data %>%
    mutate(
      var1 = as.numeric(huedisp.m),
      var2 = as.numeric(mean.ra),
      genus = as.factor(genus)
    ) %>%
    # Remove rows with NA values in PC1 or PC2
    filter(!is.na(var1), !is.na(var2)) %>%
    mutate(
      x = (var1),
      y = (var2)
    ) %>%
    # Remove rows where transformations create invalid values
    filter(is.finite(x), is.finite(y))
    # Also check that we don't have negative values under sqrt or log

  # Check if we have data left and enough points per group
  clean_data %>% 
    group_by(genus) %>% 
    summarise(n = n()) %>%
    print()
  
  # Compute convex hulls (need at least 3 points per group)
  hull_data <- clean_data %>%
    group_by(genus) %>%
    filter(n() >= 3) %>%  # Need at least 3 points for a polygon
    slice(chull(x, y)) %>%
    ungroup()
  
  # Option 5: Using distinct colors optimized for many categories
  library(randomcoloR)
  colors_100 <- distinctColorPalette(100)
  
  ggplot() +
    geom_polygon(data = hull_data, 
                 aes(x = x, y = y, color = genus, fill = genus), 
                 size = 1, alpha = .4) +
    scale_fill_manual(values = colors_100) +
    scale_color_manual(values = colors_100) +
    theme_bw()
}

#continent-level color PC plot convex hulls
{
  library(ggplot2)
  library(dplyr)
  
  # First, clean and transform the data
  clean_data <- musscrat_matched$data %>%
    mutate(
      PC1_num = as.numeric(huedisp.m),
      PC2_num = as.numeric(mean.ra),
      genus = as.factor(musscrat_biogeo)
    ) %>%
    # Remove rows with NA values in PC1 or PC2
    filter(!is.na(PC1_num), !is.na(PC2_num)) %>%
    mutate(
      x = (PC1_num),
      y = (PC2_num)
    ) %>%
    # Remove rows where transformations create invalid values
    filter(is.finite(x), is.finite(y))
  
  # Check if we have data left and enough points per group
  clean_data %>% 
    group_by(musscrat_biogeo) %>% 
    summarise(n = n()) %>%
    print()
  
  # Compute convex hulls (need at least 3 points per group)
  hull_data <- clean_data %>%
    group_by(musscrat_biogeo) %>%
    filter(n() >= 3) %>%  # Need at least 3 points for a polygon
    slice(chull(x, y)) %>%
    ungroup()
  
  # Option 5: Using distinct colors optimized for many categories
  library(randomcoloR)
  colors_100 <- distinctColorPalette(100)
  
  ggplot() +
    geom_polygon(data = hull_data, 
                 aes(x = x, y = y, color = musscrat_biogeo, fill = musscrat_biogeo), 
                 size = 1, alpha = .1) +
    scale_fill_manual(values = colors_100) +
    scale_color_manual(values = colors_100) +
    theme_bw()
}


