library(circlize)
library(plyr)
library(dplyr)

#make circlize diagam showing transitions in Rbiogeobears matrix
setwd("c:/Users/jtm378/AMNH Dropbox/Jon Merwin/AllParrotHomeAnalysis/Scripts/GITHUB_UPLOAD/8_ancestral_areas_biogeobears/")

paletteContinents <- setNames(c("#3557a0", "#2e9ccc", "#64cdda", "#2a8270", "#5aba58", "#f58643", "#bc529b", "#ea2636", "gray", "gray", "gray", "gray", "gray"), c("A", "B", "C", "D", "E", "F", "G", "H", "GH","AC", "FG", "FGH"))
#colors_list_for_states  <- c("#3557a0", "#2e9ccc",  "#64cdda", "#2a8270", "#5aba58", "#f58643", "#bc529b", "#ea2636", rep("grey50", times = (92-8)))

geo_trans_matrix <- read.csv("transition_matrix_trim.csv")

#weird issue i need to figure out tomorrow
which(geo_trans_matrix$node<4)

#geo_trans_matrix$from_map[1]<-"GH"
#geo_trans_matrix$from_map[c(32,33,56)]<-"GH"

geo_trans_matrix_plot <- geo_trans_matrix[,6:7]

colors <- plyr::mapvalues(x = geo_trans_matrix_plot$from_map, from = names(paletteContinents), to = paletteContinents)

circlize::chordDiagramFromDataFrame(geo_trans_matrix_plot, col = colors, grid.col = paletteContinents, annotationTrack = c("name", "grid"), order = rev(c("A", "AC", "B", "C", "D", "E", "F", "FG", "FGH", "G", "GH", "H")), link.border = "transparent")
system("open ./")
