#plot a bird?
library(grid)
library(scales)
library(colorspace)
library(stringi)
library(stringr)

data_poistringidata_pointer <- "c:/Users/jtm378/AMNH Dropbox/Jon Merwin/AllParrotHomeAnalysis/Scripts/GITHUB_UPLOAD/"
setwd(data_pointer)

PGLS_data_with_traitrates <- readRDS(file = "5_RRPHYLO_trait_evolution_rates_and_clads_PGLS_with_simulations/PGLS_RDS_with_rates.rds")
PGLS_data_with_traitrates$phy <- root(PGLS_data_with_traitrates$phy,outgroup = c("Calyptomena_viridis","Caracara_cheriway","Icterus_cucullatus"),resolve.root = T)

plot_bird <- function(dat = PGLS_data_with_traitrates$data[5,])
{
  back_col <- lighten(method = "relative",space="HLS",rgb(red = dat$centroid.l_back, blue = dat$centroid.s_back, green = dat$centroid.m_back),amount = as.numeric(dat$BRIGHT_Back)/2)
  body_col <- lighten(method = "relative",space="HLS",rgb(red = dat$centroid.l_body, blue = dat$centroid.s_body, green = dat$centroid.m_body),amount = as.numeric(dat$BRIGHT_UpperAbdomen)/2)
  head_col <- lighten(method = "relative",space="HLS",rgb(red = dat$centroid.l_head, blue = dat$centroid.s_head, green = dat$centroid.m_head),amount = as.numeric(dat$BRIGHT_Crown)/2)
  wing_col <- lighten(method = "relative",space="HLS",rgb(red = dat$centroid.l_wing, blue = dat$centroid.s_wing, green = dat$centroid.m_wing),amount = as.numeric(dat$BRIGHT_GreaterSecondaryCoverts)/2)
  
  # Create a 3x3 matrix:
  # [transparent] [back] [back]
  # [head]        [wing] [wing]
  # [transparent] [body] [transparent]
  
  bird_mat <- matrix(c(
    head_col,head_col,back_col,body_col,wing_col,body_col,wing_col,
    "transparent"
  ), nrow = 2, ncol = 4)
  
  grid.raster((bird_mat), interpolate = F)
}

plot_bird(PGLS_data_with_traitrates$data[sample(1:379,1),])
graphics.off()

output_dir <- "bird_images"
if(!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Modified plot_bird function to work with jpg output
plot_bird <- function(dat, species_name, output_dir) {
  back_col <- lighten(method = "relative", space = "HLS", 
                      rgb(red = dat$centroid.l_back, blue = dat$centroid.s_back, green = dat$centroid.m_back),
                      amount = as.numeric(dat$BRIGHT_Back)/2)
  body_col <- lighten(method = "relative", space = "HLS",
                      rgb(red = dat$centroid.l_body, blue = dat$centroid.s_body, green = dat$centroid.m_body),
                      amount = as.numeric(dat$BRIGHT_UpperAbdomen)/2)
  head_col <- lighten(method = "relative", space = "HLS",
                      rgb(red = dat$centroid.l_head, blue = dat$centroid.s_head, green = dat$centroid.m_head),
                      amount = as.numeric(dat$BRIGHT_Crown)/2)
  wing_col <- lighten(method = "relative", space = "HLS",
                      rgb(red = dat$centroid.l_wing, blue = dat$centroid.s_wing, green = dat$centroid.m_wing),
                      amount = as.numeric(dat$BRIGHT_GreaterSecondaryCoverts)/2)
  
  bird_mat <- matrix(c(
    head_col, head_col, back_col, body_col, 
    wing_col, body_col, wing_col, "transparent"
  ), nrow = 2, ncol = 4)
  
  # Clean species name for filename (remove special characters)
  clean_name <- gsub("[^A-Za-z0-9_-]", "_", species_name)
  filename <- file.path(output_dir, paste0(clean_name, ".jpg"))
  
  # Save to file
  jpeg(filename, width = 400, height = 200, quality = 95, bg = "white")
  grid.raster(bird_mat, interpolate = FALSE)
  dev.off()
  
  return(filename)
}

# Loop through all birds

# library(grid)  # Make sure grid is loaded
# for(i in 1:nrow(PGLS_data_with_traitrates$data)) {
#   dat <- PGLS_data_with_traitrates$data[i, ]
#   
#   # Get species name (adjust column name as needed)
#   species_name <- dat$species  # Change this to match your actual column name
#   
#   # Plot and save
#   tryCatch({
#     plot_bird(dat, species_name, output_dir)
#     cat("Saved:", species_name, "(", i, "of", nrow(PGLS_data_with_traitrates$data), ")\n")
#   }, error = function(e) {
#     cat("Error with species", species_name, ":", e$message, "\n")
#   })
# }

cat("Done! All images saved to:", output_dir, "\n")


system("open bird_images")


library(ggplot2)
library(ggimage)  # You may need to install this: install.packages("ggimage")

# Prepare data with image paths
plot_data <- PGLS_data_with_traitrates$data
plot_data$species_name <- plot_data$Species  # Adjust column name as needed
plot_data$clean_name <- gsub("[^A-Za-z0-9_-]", "_", plot_data$species_name)
plot_data$image_path <- file.path(output_dir, paste0(plot_data$clean_name, ".jpg"))
plot_data$genus <- str_split_i(plot_data$species, pattern = "_", i=1)


plot_data$norm_tars_len <- scales::rescale(as.numeric(plot_data$TARSUS.LENGTH), to = c(0.02, 0.1))

plot_data_clean <- plot_data[!is.na(plot_data$TARSUS.LENGTH), ]
plot_data_clean <- plot_data_clean[plot_data_clean$TARSUS.LENGTH > 0, ]

# Create the plot with bird images as points
# Replace 'variable1' and 'variable2' with your actual column names
windows()
color_var_explanatory <- ggplot(plot_data_clean[!is.na(plot_data$Species),], aes(x = (as.numeric(mean.ra)), y = (as.numeric(huedisp.m)))) +
  geom_image(aes(image = image_path), size = .03) +  # Adjust size as needed
  labs(x = "SATURATION", 
       y = "HUEDISP",
       title = "Bird Color Patterns") +
  theme_bw() +
  geom_text(aes(label=species),size = 2,nudge_y = -.2)+
  theme(panel.grid.minor = element_blank())

ggplot(plot_data_clean[!is.na(plot_data$species),], aes(x = log(as.numeric(centroid.l)/as.numeric(centroid.s)), y = log(as.numeric(centroid.u)/as.numeric(centroid.m)))) +
  geom_image(aes(image = image_path), size = .03) +  # Adjust size as needed
  labs(x = "<-BLUE RED->", 
       y = "<-GREEN UV->",
       title = "Bird Color Patterns") +
  theme_bw() +
  #geom_text(aes(label=species),size = 2,nudge_y = -.2)+
  theme(panel.grid.minor = element_blank())


ggplot(plot_data_clean[(!is.na(plot_data$species)),], aes(size=(400-as.numeric(PC1)),y = (as.numeric(yCentroid0.0)), x = (as.numeric(XCentroid0.0)))) + geom_point()
  
ggplot(plot_data_clean[(!is.na(plot_data$species)),], aes(y = (as.numeric(yCentroid0.0)), x = (as.numeric(XCentroid0.0)))) +
  geom_image(aes(image = image_path),size=.03) +  # Adjust size as needed
  labs(x = "Temperature", 
       y = "Precipitation",
       title = "Bird Color Patterns") +
  theme_bw() +
  geom_text(aes(label=species),size = 1.5,nudge_y = -.2)+
  theme(panel.grid.minor = element_blank())

ggplot(plot_data_clean[!is.na(plot_data$species),], aes(x = (as.numeric(yCentroid0.0)/as.numeric(XCentroid0.0)), y = sqrt(as.numeric(huedisp.m)))) +
  geom_image(aes(image = image_path), size = .03) +  # Adjust size as needed
  labs(x = "Precipitation/Temperature", 
       y = "Hue Disparity",
       title = "Bird Color Patterns") +
  theme_bw() +
  geom_text(aes(label=species),size = 1.5,nudge_y = -.02)+
  theme(panel.grid.minor = element_blank())

ggplot(plot_data_clean[!is.na(plot_data$species),], aes(x = (as.numeric(yCentroid0.0)/as.numeric(XCentroid0.0)), y = sqrt(as.numeric(mean.ra)))) +
  geom_image(aes(image = image_path), size = .03) +  # Adjust size as needed
  labs(x = "Precipitation/Temperature", 
       y = "Saturation",
       title = "Bird Color Patterns") +
  theme_bw() +
  geom_text(aes(label=species),size = 1.5,nudge_y = -.02)+
  theme(panel.grid.minor = element_blank())

ggplot(plot_data_clean[!is.na(plot_data$species),], aes(x = (as.numeric(XCentroid0.0)), y = (as.numeric(color_PC2)))) +
  geom_image(aes(image = image_path), size = .03) +  # Adjust size as needed
  labs(x = "Color PC1", 
       y = "Color PC2",
       title = "Bird Color Patterns") +
  theme_bw() +
  geom_text(aes(label=species),size = 1.5,nudge_y = -.02)+
  theme(panel.grid.minor = element_blank())


clim_var_explanatory <- ggplot(plot_data_clean[!is.na(plot_data$species),], aes(y = (as.numeric(yCentroid0.5)), x = (as.numeric(XCentroid0.5)))) +
  geom_image(aes(image = image_path), size = .03) +  # Adjust size as needed
  labs(x = "Temperature", 
       y = "Precipitation",
       title = "Bird Color Patterns") +
  theme_bw() +
  geom_text(aes(label=species),size = 2,nudge_y = -.2)+
  theme(panel.grid.minor = element_blank())


morph_var_explanatory <- ggplot(plot_data_clean[!is.na(plot_data$species),], aes(x = log(250-as.numeric(PC1)), y = (as.numeric(PC2)))) +
  geom_image(aes(image = image_path), size = .03) +  # Adjust size as needed
  labs(x = "Body Size", 
       y = "Body Shape",
       title = "Bird Color Patterns") +
  theme_bw() +
  geom_text(aes(label=species),size = 2,nudge_y = -.2)+
  theme(panel.grid.minor = element_blank())

png("explanatory_cartoons.png",width = 6000, height = 2000, units = "px")
ggpubr::ggarrange(morph_var_explanatory,color_var_explanatory,clim_var_explanatory,nrow = 1,ncol = 3)
graphics.off()
system("open explanatory_cartoons.png")

ggplot(plot_data_clean[!is.na(plot_data$species),], aes(y = log(as.numeric(centroid.l)/as.numeric(centroid.s)),x = log(as.numeric(Total.niche)))) +
  geom_image(aes(image = image_path), size = .03) +  # Adjust size as needed
  labs(x = "Niche", 
       y = "<-BLUE RED->",
       title = "Bird Color Patterns") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())


ggplot(plot_data_clean[!is.na(plot_data$species),], aes(x = log((as.numeric(clim_rates))), y = log(as.numeric(color_rates)))) +
  geom_image(aes(image = image_path), size = .03) +  # Adjust size as needed
  labs(y = "Color Rates", 
       x = "Clim Rates",
       title = "Bird Color Patterns") +
  geom_text(aes(label=species),size = 2,nudge_y = -.2)+
  theme_bw() +
  theme(panel.grid.minor = element_blank()) #+ geom_smooth(method = "lm")

ggplot(plot_data_clean[!is.na(plot_data$species),], aes(x = log((as.numeric(clim_rates))), y = log(as.numeric(color_rates)))) +
  geom_image(aes(image = image_path), size = .03) +  # Adjust size as needed
  labs(y = "Color Rates", 
       x = "Clim Rates",
       title = "Bird Color Patterns") +
  geom_text(aes(label=species),size = 2,nudge_y = -.2)+
  theme_bw() +
  theme(panel.grid.minor = element_blank()) #+ geom_smooth(method = "lm")


ggplot(plot_data_clean[!is.na(plot_data$species),], aes(x = log((as.numeric(clim_rates))), y = log(as.numeric(morph_rates)))) +
  geom_image(aes(image = image_path), size = .03) +  # Adjust size as needed
  labs(y = "Morph Rates", 
       x = "CLADS / Clim Rates",
       title = "Bird Color Patterns") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) + geom_smooth(method = "lm")





ggplot(plot_data_clean[!is.na(plot_data$species),], aes(y = log(as.numeric(color_rates)), x = log(as.numeric(clim_rates)))) +
  geom_image(aes(image = image_path), size = .03) +  # Adjust size as needed
  labs(y = "Color Rates", 
       x = "Morph Rates*Clim Rates",
       title = "Bird Color Patterns") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) + geom_smooth(method = "lm")


ggplot(plot_data_clean[!is.na(plot_data$species),], aes(y = log(as.numeric(raw.morpho.rate)), x = log(as.numeric(morph_rates)))) +
  geom_image(aes(image = image_path), size = .03) +  # Adjust size as needed
  labs(y = "BAYES Morph Rates", 
       x = "RR Morph Rates",
       title = "Bird Color Patterns") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) + geom_smooth(method = "lm")

ggplot(plot_data_clean[!is.na(plot_data$species),], aes(y = log(2+log(5+log(as.numeric(color.rate)))), x = log(as.numeric(color_rates)))) +
  geom_image(aes(image = image_path), size = .03) +  # Adjust size as needed
  labs(y = "BAYES Color Rates", 
       x = "RR Color Rates",
       title = "Bird Color Patterns") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) + geom_smooth(method = "lm") +
  geom_text(aes(label=species),size = 2,nudge_y = -.2)




brians_tree <- read.nexus(file = "data/color-multivariate_relaxed_BM.tre")
pdf("replot_relax_tree.pdf", width = 5, height = 10)
plot(brians_tree, cex = .1)
graphics.off()
system("open replot_relax_tree.pdf")




















