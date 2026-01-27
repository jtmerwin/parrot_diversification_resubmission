#trying to understand the climate variables

setwd("c:/Users/jtm378/AMNH Dropbox/Jon Merwin/AllParrotHomeAnalysis/Scripts/GITHUB_UPLOAD/")

parrot_rds <- readRDS("5_RRPHYLO_trait_evolution_rates_and_clads_PGLS_with_simulations/PGLS_RDS_with_rates.rds")

names(parrot_rds$data)

ggplot(parrot_rds$data, aes(x=as.numeric(XCentroid0.0), y=as.numeric(yCentroid0.0))) + geom_label(aes(label=species))
ggplot(parrot_rds$data, aes(x=as.numeric(Total.niche), y=as.numeric(Niche.over.50.))) + geom_label(aes(label=species),size=3) +
  scale_x_continuous("log10") + scale_y_continuous("log10")

ggplot(parrot_rds$data, aes(x=as.numeric(yCentroid0.0)*as.numeric(XCentroid0.0), y=as.numeric(yCentroid0.0)/as.numeric(XCentroid0.0))) + geom_label(aes(label=species))

clim_pca_data <- na.omit(parrot_rds$data[,c("XCentroid0.0","yCentroid0.0")])
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

parrot_rds_climdat <- safe_merge(parrot_rds$data, clim_pca_df, append = "_climate")

body_size_model <- lm(data = parrot_rds_climdat, formula = this_PGLS_Formula)
summary(body_size_model)

ggplot(parrot_rds_climdat, aes(x=(as.numeric(PC2_climate)), y=log(250-as.numeric(PC1)))) + geom_point(size=3) + theme_bw() + geom_smooth(method = "lm")

body_size_model <- lm(data = parrot_rds_climdat, formula = this_PGLS_Formula)
summary(body_size_model)
plot(body_size_model)

clim_pca_comparative <- comparative.data(phy = root(parrot_rds$phy,outgroup=c("Caracara_cheriway"), resolve.root = T),parrot_rds_climdat[,c("species","PC1_climate","PC2_climate","PC1")], names.col = "species")
this_PGLS_Formula <- formula(log(250-as.numeric(PC1))~(as.numeric(PC1_climate)))
this_PGLS_Formula <- formula((as.numeric(PC1))~(as.numeric(PC2_climate)))

pca_climate_mod<- caper::pgls(data = clim_pca_comparative, formula = this_PGLS_Formula)
summary(pca_climate_mod)
plot(pca_climate_mod)

