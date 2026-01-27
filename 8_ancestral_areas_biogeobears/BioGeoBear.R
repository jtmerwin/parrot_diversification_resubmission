neededPackages <- c("snow", "devtools", "FD", "GenSA", "parallel", "rexpokit", 
                    "cladoRcpp", "BioGeoBEARS", "ape", "phytools", "here", "dispRity", "phylobase", "ggpmisc")
lapply(neededPackages, require, character.only = TRUE)

#useful method for matching later 
lastDropper <- function(string = NULL, seper = " ")
{
  require("stringr")
  lab <- word(string, 1, -2, sep = seper)
  return(lab)
}

setwd("c:/Users/jtm378/AMNH Dropbox/Jon Merwin/AllParrotHomeAnalysis/Scripts/GITHUB_UPLOAD/8_ancestral_areas_biogeobears/data/")
getwd()

extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
extdata_dir
list.files(extdata_dir)
system("open extdata_dir")


colors_list_for_states  <- c("#EB8645", "#3557a0", "#2e9ccc", "#bc529b", "#64cdda", "#2a8270", "#5aba58", "#f58643", "#ebe654", "#ea2636", "#327ab6", rep("grey50", times = 45))


colors_list_for_states  <- c("#3557a0", "#2e9ccc",  "#64cdda", "#2a8270", "#5aba58", "#f58643", "#bc529b", "#ea2636", rep("grey50", times = (93-8)))


length("possible_ranges_list_txt")

#### MAKE PARROT TREE PDF ####
#### NOT NEEDED FOR CODE ####

trfn <- "parrot_timetree.tre"
moref(trfn)

tr = read.nexus(trfn)
tr$tip.label <- lastDropper(tr$tip.label,seper = "_")
bio_geo_trim_data <- read.delim("3-nov-2025-biogeographic-areas-reassigned.txt",sep = "\t", colClasses = "character")
View(bio_geo_trim_data)
#bio_geo_trim_data$X.A.B.C.E.F.G.H.I. <- as.character(bio_geo_trim_data$X.A.B.C.E.F.G.H.I.)
rownames(bio_geo_trim_data) <- bio_geo_trim_data$X381
trim_tr <- picante::match.phylo.data(tr, data = bio_geo_trim_data)
write.tree(phy = trim_tr$phy, file = "tree.newick")

pdffn = "tree.pdf"
pdf(file=pdffn, width=30, height= 100)

plot(trim_tr$phy)
title("Parrot Tree")
axisPhylo() # plots timescale

graphics.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)

#### Set up BioGeoBears ####
tr <- read.tree("tree.newick")
trfn <- "tree.newick"
#geogfn <- "RangeMetadataParrots2020 - BioGeo and TaxonBetter - Copy.txt"
geogfn <- "3-nov-2025-biogeographic-areas-reassigned.txt"
typeof(geogfn)
system("open 3-nov-2025-biogeographic-areas-reassigned.txt")

system("open ./")


BioGeoBEARS_run_object = define_BioGeoBEARS_run(trfn = trfn, geogfn = geogfn)
BioGeoBEARS_run_object$geogfn = geogfn

tipranges = (getranges_from_LagrangePHYLIP(lgdata_fn=geogfn))

# tipranges

max(rowSums(dfnums_to_numeric(tipranges@df)))
numstates_from_numareas(numareas=4, maxareas=4, include_null_range=FALSE)
numstates_from_numareas(numareas=4, maxareas=4, include_null_range=FALSE)
numstates_from_numareas(numareas=4, maxareas=3, include_null_range=FALSE)
numstates_from_numareas(numareas=4, maxareas=2, include_null_range=FALSE)
numstates_from_numareas(numareas=10, maxareas=10, include_null_range=FALSE)
numstates_from_numareas(numareas=8, maxareas=3, include_null_range=FALSE)
BioGeoBEARS_run_object$max_range_size = 3
BioGeoBEARS_run_object$min_branchlength = -0.00000    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = FALSE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 16
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

BioGeoBEARS_run_object
BioGeoBEARS_run_object$BioGeoBEARS_model_object
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table
BioGeoBEARS::check_BioGeoBEARS_run(BioGeoBEARS_run_object)

max_range_size = 3

areas = getareas_from_tipranges_object(tipranges)
states_list_0based = rcpp_areas_list_to_states_list(areas = areas, maxareas = max_range_size, include_null_range = FALSE)
length(states_list_0based)

## Manual Modifications of States ####

max_range_size = 3
numstates_from_numareas(numareas=8, maxareas=3, include_null_range=FALSE)

areas = getareas_from_tipranges_object(tipranges)
states_list_0based = rcpp_areas_list_to_states_list(areas = areas,
                                                    maxareas = max_range_size, 
                                                    include_null_range = FALSE)
length(states_list_0based)

## Make the list of ranges

ranges_list = NULL
for (i in 1:length(states_list_0based))
{
  if ( (length(states_list_0based[[i]]) == 1) && (is.na(states_list_0based[[i]])) )
  {
    tmprange = "_"
  } else {
    tmprange = paste(areas[states_list_0based[[i]]+1], collapse = "")
  }
  ranges_list = c(ranges_list, tmprange)
}
ranges_list

#### Time Strata for Time-Stratified Analysis ####

#in new analysis
#A = North America
#B = Caribbean
#C = South America
#D = Africa
#E = Eurasia
#F = WallaceaPhillippines
#G = AustraliaPapuaNewguinea
#H = Oceania

# STRATUM 1 (3 mya) #
#current day
nonadjacent = c("A", "B", "C", "D", "E", "F",
                "G", "H", "AC", "DE", "EF", "FG","FGH", "GH")

#trying only physical connection
nonadjacent = c("A", "B", "C", "D", "E", "F",
                "G", "H", "AC", "DE")
keepTF1 = ranges_list %in% nonadjacent == TRUE|FALSE
ranges_list_NEW1 = ranges_list#[keepTF1]
ranges_list_NEW1


# STRATUM 2 (10 mya) #
#basically the same as present day
nonadjacent = c("A", "B", "C", "D", "E", "F",
                "G", "H", "AC", "DE", "EF", "FG","FGH","GH")

#citation for sundaland: https://doi.org/10.1016/j.pld.2022.05.006
#trying only physical connection, isthmus and sundaland exist
nonadjacent = c("A", "B", "C", "D", "E", "F",
                "G", "H", "AC", "DE", "EF", "DEF")

keepTF2 = ranges_list %in% nonadjacent == TRUE
ranges_list_NEW2 = ranges_list#[keepTF2]
ranges_list_NEW2

# STRATUM 3 (20 mya) #
#basically the same as 10mya except no isthmus of panama
#citation for wallacea: https://doi.org/10.1016/j.pld.2022.05.006
#trying only physical connection
nonadjacent = c("A", "B", "C", "D", "E", "F",
                "G", "H", "DE", "EF", "DEF")

keepTF3 = ranges_list %in% nonadjacent == TRUE
ranges_list_NEW3 = ranges_list#[keepTF3]
ranges_list_NEW3

# STRATUM 4 (30 mya) 
#NZ quite distant from distant from asia
#sundaland inundated, no isthmus
nonadjacent = c("A", "B", "C", "D", "E", "F",
                "G", "H", "DE", "EF")

#trying only physical connection
nonadjacent = c("A", "B", "C", "D", "E", "F",
                "G", "H", "DE", "EF")

keepTF4 = ranges_list %in% nonadjacent == TRUE
ranges_list_NEW4 = ranges_list#[keepTF4]
ranges_list_NEW4

# STRATUM 5 (before 30 mya)
#australia + ng + nz were quite distant from asia, no adjacency
# nonadjacent = c("A", "B", "C", "E", "F", "G",
#                 "H", "I", "EF", "FG","HI")

nonadjacent = c("A", "B", "C", "D", "E", "F",
                "G", "H", "DE", "EF","GH")

#trying only physical connection
nonadjacent = c("A", "B", "C", "D", "E", "F",
                "G", "H", "DE", "EF")

keepTF5 = ranges_list %in% nonadjacent == TRUE
ranges_list_NEW5 = ranges_list#[keepTF5]
ranges_list_NEW5

states_list_0based_NEW1 = states_list_0based[keepTF1]
states_list_0based_NEW2 = states_list_0based[keepTF2]
states_list_0based_NEW3 = states_list_0based[keepTF3]
states_list_0based_NEW4 = states_list_0based[keepTF4]
states_list_0based_NEW5 = states_list_0based[keepTF5]

BioGeoBEARS_run_object$lists_of_states_0based[[1]] = states_list_0based_NEW1
BioGeoBEARS_run_object$lists_of_states_0based[[2]] = states_list_0based_NEW2
BioGeoBEARS_run_object$lists_of_states_0based[[3]] = states_list_0based_NEW3
BioGeoBEARS_run_object$lists_of_states_0based[[4]] = states_list_0based_NEW4
BioGeoBEARS_run_object$lists_of_states_0based[[5]] = states_list_0based_NEW5

# ######################################################
# Run DEC
# ######################################################

# Initialize a default model (DEC model)
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = FALSE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
#  Jeremy M.; Matzke, Nicholas J.; O'Meara, Brian C. (2015). Non-null Effects of 
#  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
#  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change

# Set up a time-stratified analysis:
# 1. Here, un-comment ONLY the files you want to use.
# 2. Also un-comment "BioGeoBEARS_run_object = section_the_tree(...", below.
# 3. For example files see (a) extdata_dir, 
#  or (b) http://phylo.wikidot.com/biogeobears#files
#  and BioGeoBEARS Google Group posts for further hints)
#
# Uncomment files you wish to use in time-stratified analyses:
BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
BioGeoBEARS_run_object$areas_allowed_fn = NA
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

max_range_size

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 16
# (use more cores to speed it up; this requires
# library(parallel) and/or library(snow). The package "parallel" 
# is now default on Macs in R 3.0+, but apparently still 
# has to be typed on some Windows machines. Note: apparently 
# parallel works on Mac command-line R, but not R.app.
# BioGeoBEARS checks for this and resets to 1
# core with R.app)

# Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
# I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
# but the results are imprecise and so I haven't explored it further.
# In a Bayesian analysis, it might work OK, but the ML point estimates are
# not identical.
# Also, I have not implemented all functions to work with force_sparse=TRUE.
# Volunteers are welcome to work on it!!
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DEC model
# (nothing to do; defaults)

# Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
BioGeoBEARS_run_object

# This contains the model object
BioGeoBEARS_run_object$BioGeoBEARS_model_object

# This table contains the parameters of the model 
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

# Run this to check inputs. Read the error messages if you get them!
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = TRUE
resfn = "restable_v1.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDEC = res
} else {
  # Loads to "res"
  load(resfn)
  resDEC = res
}

#######################################################
# Run DEC+J
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = FALSE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
#  Jeremy M.; Matzke, Nicholas J.; O'Meara, Brian C. (2015). Non-null Effects of 
#  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
#  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change

# Set up a time-stratified analysis:
BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
BioGeoBEARS_run_object$areas_allowed_fn = NA
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 16
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DEC+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resDEC$outputs@params_table["d","est"]
estart = resDEC$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "Psittaciform_DEC+J_M0_unconstrained_v1.Rdata"
runslow = TRUE
if (runslow)
{
  #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
  
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDECj = res
} else {
  # Loads to "res"
  load(resfn)
  resDECj = res
}

#######################################################
# PDF plots
#######################################################
pdffn = "Psittaciformes_DEC_vs_DEC+J_M0_unconstrained_v1.pdf"
pdf(pdffn, width=30, height=100)

#######################################################
# Plot ancestral states - DEC
#######################################################
analysis_titletxt ="BioGeoBEARS DEC on Psittaciform M0_unconstrained"

# Setup
results_object = resDEC
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

colors_list_for_states  <- c("#3557a0", "#2e9ccc", "#64cdda","#2a8270", "#5aba58", "#f58643", "#bc529b", "#ea2636", randomColor(count = 92-8, hue = "monochrom"))
names(colors_list_for_states)
#colors_list_for_states  <- randomColor(count = 92)

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=1, statecex=1, splitcex=1, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, colors_list_for_states = colors_list_for_states,include_null_range=FALSE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=.3, statecex=.3, splitcex=.3, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, colors_list_for_states = colors_list_for_states,include_null_range=FALSE, tr=tr, tipranges=tipranges, plotlegend = T)

#######################################################
# Plot ancestral states - DECJ
#######################################################
analysis_titletxt ="BioGeoBEARS DEC+J on Psittaciform M0_unconstrained"

# Setup
results_object = resDECj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# length(colors_list_for_states)
# length(unique(geogfn[2]))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=1, statecex=1, splitcex=1, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, colors_list_for_states = colors_list_for_states, include_null_range=FALSE, tr=trim_tr$phyim_tr$phy, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=.3, statecex=.3, splitcex=.3, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, colors_list_for_states = colors_list_for_states,include_null_range=FALSE, tr=trim_tr$phyim_tr$phy, tipranges=tipranges, plotlegend = T)

graphics.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

#now make the painted simmap tree!

tip_length <- 381

#start here tomorrow to start making the biogeobears plots and pulling the data out of the biogeobears object...
colors_list_for_states  <- c("#3557a0", "#2e9ccc",  "#64cdda", "#2a8270", "#5aba58", "#f58643", "#bc529b", "#ea2636", rep("grey50", times = (92-8)))

modelSelectPlot <- resDEC
#modelSelectPlot <- resDECj

#get plots out of biogeobears ecosystem
summarizedBestNodes <- apply(modelSelectPlot$ML_marginal_prob_each_state_at_branch_top_AT_node, MARGIN = 1, FUN = function(X) which.max(X))
rootState <- which.max(modelSelectPlot$ML_marginal_prob_each_state_at_branch_top_AT_node[tip_length+1,])
summarizedTipStates <- apply(tipranges@df, MARGIN = 1, FUN = function(X) paste(names(tipranges@df)[which(as.numeric(X)==1)],collapse = ""))
resWithStates <- cbind.data.frame(prt(trim_tr$phy), summarizedBestNodes = c(unlist(summarizedBestNodes)[1:(tip_length)],rootState, unlist(summarizedBestNodes)[(tip_length+2):(tip_length*2-1)]))

resWithStates[length(resWithStates$node),]

#get the plotting order of the tips
plot(trim_tr$phy, cex = .1)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tipsInYorder <- rev(paintedTreeDec$tip.label[order(lastPP$yy[1:Ntip(paintedTreeDec)])])

#need to rerun this
for(i in 1:length(resWithStates$node))
{
try(if(i == 1){paintedTree <- paintSubTree(trim_tr$phy, node = resWithStates$node[i], state = resWithStates$summarizedBestNodes[i], stem =T)})
try(if((i != 382)&&(i!= 760)){paintedTree <- paintSubTree(paintedTree, node = resWithStates$node[i], state = resWithStates$summarizedBestNodes[i], stem = T)})
try(if((i == 382)||(i == 760)){paintedTree <- paintSubTree(paintedTree, node = resWithStates$node[i], state = resWithStates$summarizedBestNodes[i], stem = F)})
print(resWithStates$summarizedBestNodes[i])
print(i)
}

#branchCols <- (rep("black", count = length(unique(resWithStates$summarizedBestNodes))))
#names(branchCols) <- unique(resWithStates$summarizedBestNodes) 

#turn these on and off to switch which model plot you are running
paintedTreeDec <- paintedTree
#paintedTreeDecJ <- paintedTree

colorAreaMatch <- data.frame(colors = colors_list_for_states, ranges = ranges_list)

branchCols <- colorAreaMatch$colors
tipCols <- colorAreaMatch$colors
names(branchCols) <- rownames(colorAreaMatch)
names(tipCols) <- colorAreaMatch$ranges

#plotting script to generate a comparison of the two models, using plotSimmap
png(filename = "biogeobearsModelComparison.png", width = 5000, height = 5000, units = "px")
par(mfrow = c(1,2))
tipRangeColors <- plyr::mapvalues(x = summarizedTipStates, to = colorAreaMatch$colors, from = colorAreaMatch$ranges)
plotSimmap(paintedTreeDec,colors = branchCols,lwd=5,pts=F, fsize = F, ftype = "off", mar = c(2,2,2,0))
title(main = paste("DEC: lnl = ", resDEC$grand_total_likelihood), cex.main = 10)
#legend("bottomleft", fill = branchCols[order(as.numeric(names(branchCols)))], legend = ranges_list[order(as.numeric(names(branchCols)))], ncol = 3)
points(1:381~rep(37.5-h, 381), pch = 15, col = tipRangeColors[rev(tipsInYorder)], cex = 2)
legend("bottomleft", fill = branchCols[1:8], legend = colorAreaMatch$ranges[1:8], ncol = 2, border = "transparent", bg = "transparent", box.col = "transparent", cex = 4)
plotSimmap(paintedTreeDecJ,colors = branchCols,lwd=5,pts=F, fsize = 0, direction = "leftwards", mar = c(2,0,2,2))
h <- min(nodeHeights(paintedTreeDecJ))
points(1:381~rep(h, 381), pch = 15, col = tipRangeColors[rev(tipsInYorder)], cex = 2)
obj<-markChanges(tree = paintedTreeDecJ, colors = branchCols, plot = F)
from <- stringr::str_split(rownames(obj), pattern = "->", simplify = T)[,1]
to <- stringr::str_split(rownames(obj), pattern = "->", simplify = T)[,2]
points(obj[which(from != to),c(1,2)], cex = 10, col = "red", pch = 0)
title(main = paste("DEC+J: lnl = ", resDECj$grand_total_likelihood), cex.main = 10)
#nodelabels(frame = "none", adj = -.1, cex = 2)
graphics.off()
system("open biogeobearsModelComparison.png") # Plot it

#plot DecJ for manuscript
pdf("painted_decJ_manuscript.pdf", width = 10, height = 20)
plotSimmap(paintedTreeDecJ,colors = branchCols,lwd=5,pts=F,ftype = "i",fsize=.5, direction = "rightwards", mar = c(2,0,2,2))
h <- min(nodeHeights(paintedTreeDecJ))
points(1:381~rep(37.8-h, 381), pch = 15, col = tipRangeColors[rev(tipsInYorder)], cex = .5)
obj<-markChanges(tree = paintedTreeDecJ, colors = branchCols, plot = F)
from <- stringr::str_split(rownames(obj), pattern = "->", simplify = T)[,1]
to <- stringr::str_split(rownames(obj), pattern = "->", simplify = T)[,2]
#points(obj[which(from != to),c(1,2)], cex = 10, col = "red", pch = 0)
#title(main = paste("DEC+J: lnl = ", resDECj$grand_total_likelihood), cex.main = 10)
graphics.off()
system("open ./")

#get the internal changes
internal_transitions <- data.frame(str_split(rownames(obj),pattern = "->",simplify = T))
colnames(internal_transitions) <- c("from","to")
internal_transitions$type <- "internal"
internal_transitions$node <- rownames(internal_transitions)

#get tip changes
tipStates <- plyr::mapvalues(summarizedTipStates, from = colorAreaMatch$ranges, to = rownames(colorAreaMatch))
tipParents <- castor::get_ancestral_nodes(tree = trim_tr$phy, descendants = names(summarizedTipStates), Nsplits = 1)
parentStates <- getStates(paintedTreeDecJ, "nodes")[tipParents]

tip_transitions <- data.frame(
  node = names(tipStates),
  from = parentStates,
  to = tipStates,
  type = "tip"
)

transition_matrix <- rbind(internal_transitions, tip_transitions)
transition_matrix_trim <- transition_matrix[which(transition_matrix$from!=transition_matrix$to),]
transition_matrix_trim$from_map <- plyr::mapvalues(x = transition_matrix_trim$from, to = colorAreaMatch$ranges, from = rownames(colorAreaMatch))
transition_matrix_trim$to_map <-  plyr::mapvalues(x = transition_matrix_trim$to, to = colorAreaMatch$ranges, from = rownames(colorAreaMatch))
transition_matrix$from_map <- plyr::mapvalues(x = transition_matrix$from, to = colorAreaMatch$ranges, from = rownames(colorAreaMatch))
transition_matrix$to_map <-  plyr::mapvalues(x = transition_matrix$to, to = colorAreaMatch$ranges, from = rownames(colorAreaMatch))

write.csv(x=transition_matrix_trim, file = "transition_matrix_trim.csv")
write.csv(x=transition_matrix, file = "transition_matrix_untrimmed.csv")

#finish here to make circos plot...


#continue below to summarize values below each biogeographic transition (not in manuscript)

rateData <- read.csv(file = "C:/Users/jtm378/Dropbox (AMNH)/AllParrotHomeAnalysis/Scripts/AllParrotAnalysis/Data/fullBoundDataListGeoColorMorphWithRR.csv")
rateData$names <- rateData$TreeNames

nodeAges <- dispRity::tree.age(tree = tr, order = "past")
nodeSheet <- read.csv("NewAreasNodes.csv")
climRates <- data.frame(names = rownames(rateData), climRate = rateData$climRate)

#climRates <- data.frame(names = rownames(firstRRClim$rates), climRate = firstRRClim$rates)

View(nodeSheet)
#plotting the relationship between jumps to new areas and changes in trait and div. rates

#lastDropper
lastDropper <- function(string = NULL, seper = " ")
{
  require("stringr")
  lab <- word(string, 1, -2, sep = seper)
  return(lab)
}

#rm(phy); rm(node); rm(toRegion)
#phy = tr; node = 384; toRegion = "H"
#phy = tr; node = NULL; tip = "Conuropsis_carolinensis"; toRegion = "A"

newAreaDataReturn <- function(phy = tr, node = NA, toRegion = NA, tip = "")
{
  sisterTree <- phy
  sisterTree$tip.label <- lastDropper(phy$tip.label, seper = "_")
  
  print(node)
  print(tip)
  
  if(tip != ""){
    sisterTax <- getSisters(tree = sisterTree, node = tip, mode = "label")
    if(is.numeric(sisterTax$nodes)){sisterTax <- getDescendants(sisterTree, node = sisterTax$nodes)}
    node <- as.numeric(phylobase::ancestors(phy = phylo4(sisterTree), node = tip, type = "parent"))
  }

  print(node); print(tip); print(toRegion)
  
  output <- list()
  
  if(is.null(node) || is.na(toRegion)){
    print("No node and/or region specified! Returning NA"); 
    output$nodeNum <- NA
    output$morphRateMeanClade <- NA
    output$morphRateSEClade <- NA
    output$tipNum <- NA
    output$genus <- NA
    return(output)}
  
  #get the age of the node of interest
  output$nodeAge <- nodeAges$ages[node]
  print(output$nodeAge)
  
  #get tips of tree downstream of the shift node   
  descentTips <- (sisterTree$tip.label[getDescendants(tree = sisterTree, node = node)[getDescendants(tree = sisterTree, node = node) < 383]])
  print(descentTips)
  fullDescentList <- paste(descentTips, sep = ";")
  
  #pull out just the tips which are still in the ancestral area
  inRegion <- which(rateData$Biogeographic.Region == c(toRegion))
  inData <- which(rateData$TreeNames %in% descentTips)
  descentTips <- rateData$TreeNames[base::intersect(inRegion, inData)]
    
  print(descentTips)
  
  if(tip != ""){descentTips[1] <- tip}

  #pull out clads div. rate data from the rate sheet from those tips
  cladsNode <- na.omit(rateData[which(rateData$TreeNames %in% descentTips), "ClaDS"])
  output$cladsMeanClade <- mean(cladsNode)
  output$cladsRateSEClade <- sd(cladsNode) / sqrt(length(cladsNode))
  
  #pull out climate rate change data from the rate sheet from those tips
  climNode <- na.omit(rateData[which(rateData$TreeNames %in% descentTips), "Total.niche.y"])
  output$nichesizeMeanClade <- mean(climNode)
  output$nicheSizeSEClade <- sd(climNode) / sqrt(length(climNode))
  
  #pull out climate rate change data from the rate sheet from those tips
  (climRRNode <- na.omit(rateData[which(rateData$names %in% descentTips), "climRate"]))
  output$nicheRRMeanClade <- mean(climRRNode)
  output$nicheRRSEClade <- sd(climRRNode) / sqrt(length(climRRNode))
  
  #pull out color rate data from the rate sheet from those tips
  colRatesNode <- as.numeric(na.omit(rateData[which(rateData$TreeNames %in% descentTips), c("colRate")]))
  output$colorRateMeanClade <- mean(colRatesNode)
  output$colorRateSEClade <- sd(colRatesNode) / sqrt(length(colRatesNode))
  
  #pull out greenness data from the rate sheet from those tips
  greennessNode <- as.numeric(na.omit(rateData[which(rateData$TreeNames %in% descentTips), c("Average.of.centroid.m")]))
  rednessNode <- as.numeric(na.omit(rateData[which(rateData$TreeNames %in% descentTips), c("Average.of.centroid.l")]))
  everythingNode <- as.numeric(unlist(na.omit(rateData[which(rateData$TreeNames %in% descentTips), c("Average.of.centroid.u", "Average.of.centroid.s", "Average.of.centroid.l")])))
  output$greens <- mean(greennessNode)-mean(everythingNode)

  #pull out saturation data from the rate sheet from those tips
  saturationNode <- as.numeric(na.omit(rateData[which(rateData$TreeNames %in% descentTips), c("Average.of.mean.ra")]))
  output$Saturation <- mean(saturationNode)
  
  #pull out climate measures from the rate sheet from those tips
  climateNodeX <- as.numeric(na.omit(rateData[which(rateData$TreeNames %in% descentTips), c("XCentroid0.0.y")]))
  output$climateMeanX <- mean(climateNodeX)
  
  #pull out climate measures from the rate sheet from those tips
  climateNodeY <- as.numeric(na.omit(rateData[which(rateData$TreeNames %in% descentTips), c("yCentroid0.0.y")]))
  output$climateMeanY <- mean(climateNodeY)
  
  #pull out color rate data from the rate sheet from those tips
  bodysizeNode <- as.numeric(na.omit(rateData[which(rateData$TreeNames %in% descentTips), c("PC1")]))
  output$BodySize <- mean(bodysizeNode)
  
  
  
  #pull out morph rate data from the rate sheet from those tips
  morphRatesNode <- as.numeric(na.omit(rateData[which(rateData$TreeNames %in% descentTips), "morphRate"]))
  output$nodeNum <- node
  output$morphRateMeanClade <- mean(morphRatesNode)
  output$morphRateSEClade <- sd(morphRatesNode) / sqrt(length(morphRatesNode))
  output$tipNum <- length(descentTips)
  print(fullDescentList)
  #output$genus <- fullDescentList
  output$genus <- paste(descentTips[1], sep = ";")
  head(output)
  return(output)
}

#nodeSheet <- read.csv("Data/NewAreasNodes.csv")
nodeSheet[33,]

for(i in 1:60)
{
  print(i);print(nodeSheet$Node.Num[i]);nodeSheet$To.Area[i]
  thisResult <- newAreaDataReturn(phy = tr, tip = nodeSheet$Tip[i], node = nodeSheet$Node.Num[i], toRegion = nodeSheet$To.Area[i])
  print("Thisresult: ")
  print(thisResult)
  if(!is.na(thisResult$tipNum) && thisResult$nodeNum > 0)
  {
  nodeSheet$Node.Num[i] <- thisResult$nodeNum
  nodeSheet$Node.Age[i] <- thisResult$nodeAge
  nodeSheet$Mean.CLADS[i] <- thisResult$cladsMeanClade
  #nodeSheet$SE.CLADS[i] <- thisResult$cladsRateSEClade
  nodeSheet$Mean.Plum.RR[i] <- thisResult$colorRateMeanClade
  nodeSheet$SE.Plum.RR[i] <- thisResult$colorRateSEClade
  nodeSheet$Mean.Morph.RR[i] <- thisResult$morphRateMeanClade
  nodeSheet$SE.Morph.RR[i] <- thisResult$morphRateSEClade
  nodeSheet$Mean.Clim.RR[i] <- thisResult$nicheRRMeanClade
  nodeSheet$SE.Morph.RR[i] <- thisResult$nicheRRSEClade
  nodeSheet$Clade.Description[i] <- thisResult$genus
  nodeSheet$Greenness[i] <- thisResult$greens
  nodeSheet$BodySize[i] <- thisResult$BodySize
  nodeSheet$ClimX[i] <- thisResult$climateMeanX
  nodeSheet$ClimY[i] <- thisResult$climateMeanY
  nodeSheet$Saturation[i] <- thisResult$Saturation
  }
}

write.csv(nodeSheet, file = "Data/nodeSheet_transition_data_12.17.25.csv")

nodeSheet <- read.csv("Data/nodeSheet_transition_data_7.25.2023.csv")


node_hist <- ggplot(nodeSheet, aes(x = log(.5+1/Node.Age), fill = To.Area, colour = To.Area)) +
 geom_density(position = "stack", alpha = 1) + theme_minimal()+
 scale_color_manual(values=c("#3557a0", "#2e9ccc", "#2a8270", "#5aba58", "#f58643", "#ebe654", "#ea2636", "gray")) +
  scale_fill_manual(values=c("#3557a0", "#2e9ccc", "#2a8270", "#5aba58", "#f58643", "#ebe654", "#ea2636", "gray")) 
  
plumRR_hist <- ggplot(nodeSheet, aes(x = log(Mean.Plum.RR), fill = To.Area, colour = To.Area)) +
  geom_density(position = "stack", alpha = 1) + theme_minimal()+
  scale_color_manual(values=c("#3557a0", "#2e9ccc", "#2a8270", "#5aba58", "#f58643", "#ebe654", "#ea2636", "gray")) +
  scale_fill_manual(values=c("#3557a0", "#2e9ccc", "#2a8270", "#5aba58", "#f58643", "#ebe654", "#ea2636", "gray")) 

morphRR_hist <- ggplot(nodeSheet, aes(x = log(Mean.Morph.RR), fill = To.Area, colour = To.Area)) +
  geom_density(position = "stack", alpha = 1) + theme_minimal()+
  scale_color_manual(values=c("#3557a0", "#2e9ccc", "#2a8270", "#5aba58", "#f58643", "#ebe654", "#ea2636", "gray")) +
  scale_fill_manual(values=c("#3557a0", "#2e9ccc", "#2a8270", "#5aba58", "#f58643", "#ebe654", "#ea2636", "gray")) 


clads_hist <- ggplot(nodeSheet, aes(x = (log(Mean.CLADS)), fill = To.Area, colour = To.Area)) +
  geom_density(position = "stack", alpha = 1) + theme_minimal()+
  scale_color_manual(values=c("#3557a0", "#2e9ccc", "#2a8270", "#5aba58", "#f58643", "#ebe654", "#ea2636", "gray")) +
  scale_fill_manual(values=c("#3557a0", "#2e9ccc", "#2a8270", "#5aba58", "#f58643", "#ebe654", "#ea2636", "gray")) 


climRR_hist <- ggplot(nodeSheet, aes(x = log(Mean.Clim.RR), fill = To.Area, colour = To.Area)) +
  geom_density(position = "stack", alpha = 1) + theme_minimal()+
  scale_color_manual(values=c("#3557a0", "#2e9ccc",  "#2a8270", "#5aba58", "#f58643", "#ebe654", "#ea2636", "gray")) +
  scale_fill_manual(values=c("#3557a0", "#2e9ccc", "#2a8270", "#5aba58", "#f58643", "#ebe654", "#ea2636", "gray")) 


pdf("FigureOutput/biogeo_histograms.pdf", width = 20, height = 10)
gridExtra::grid.arrange(node_hist, plumRR_hist, morphRR_hist, clads_hist, climRR_hist, ncol = 3)
graphics.off()


#violinplots start here
ggplot(nodeSheet, aes(x = factor(To.Area), y = log(Mean.Clim.RR), fill = To.Area)) +
  geom_violin(color = "transparent") + theme_minimal()+  geom_jitter(height = 0, width = 0.1, colour = "red") +
  scale_color_manual(values=c("#3557a0", "#2e9ccc", "#2a8270", "#5aba58", "#f58643", "#ebe654", "#ea2636", "gray")) +
  scale_fill_manual(values=c("#3557a0", "#2e9ccc", "#2a8270", "#5aba58", "#f58643", "#ebe654", "#ea2636", "gray")) 

ggplot(nodeSheet, aes(x = factor(To.Area), y = log(Mean.Plum.RR), fill = To.Area)) +
  geom_violin(color = "transparent") + theme_minimal()+  geom_jitter(height = 0, width = 0.1, colour = "red") +
  scale_color_manual(values=c("#3557a0", "#2e9ccc","#2a8270", "#5aba58", "#f58643", "#ebe654", "#ea2636", "gray")) +
  scale_fill_manual(values=c("#3557a0", "#2e9ccc", "#2a8270", "#5aba58", "#f58643", "#ebe654", "#ea2636", "gray")) 

ggplot(nodeSheet, aes(x = factor(To.Area), y = log(Mean.Morph.RR), fill = To.Area)) +
  geom_violin(color = "transparent") + theme_minimal()+  geom_jitter(height = 0, width = 0.1, colour = "red") +
  scale_color_manual(values=c("#3557a0", "#2e9ccc", "#2a8270", "#5aba58", "#f58643", "#ebe654", "#ea2636", "gray")) +
  scale_fill_manual(values=c("#3557a0", "#2e9ccc", "#2a8270", "#5aba58", "#f58643", "#ebe654", "#ea2636", "gray")) 


earth_plotmodsel(cladsStat, col.grsq = "red")
earth_plotmodsel(colorStat, add = T, col.grsq = "goldenrod")
earth_plotmodsel(morphStat, add = T, col.grsq = "purple")
earth_plotmodsel(climStat, add = T, col.grsq = "cyan")


cladsStat <- earth(log(Mean.CLADS) ~ log(.5+1/Node.Age), data = na.omit(nodeSheet[,c("Mean.CLADS", "Node.Age")]), nk = 50, degree = 1, )
colorStat <- earth(log(Mean.Plum.RR) ~  log(.5+1/Node.Age), data = na.omit(nodeSheet[,c("Mean.Plum.RR", "Node.Age")]), nk = 50, degree = 1)
morphStat <- earth(log(Mean.Morph.RR) ~ log(.5+1/Node.Age), data = na.omit(nodeSheet[,c("Mean.Morph.RR", "Node.Age")]), nk = 50, degree = 1)
climStat <- earth(log(Mean.Clim.RR) ~  log(.5+1/Node.Age), data = na.omit(nodeSheet[,c("Mean.Clim.RR", "Node.Age")]), nk = 50, degree = 1)

nodeSheet$uniqueNames <- make.names(unique = T, names = nodeSheet$Clade.Description)

cladsFormula <- formula(log(Mean.CLADS) ~ log(.5+1/Node.Age))
colorFormula <- formula(log(Mean.Plum.RR) ~ log(.5+1/Node.Age))
morphFormula <- formula(log(Mean.Morph.RR) ~ log(.5+1/Node.Age))
climFormula <- formula(log(Mean.Clim.RR) ~ log(.5+1/Node.Age))
combinedFormula <-  formula(log(Mean.Clim.RR)+log(Mean.Plum.RR) ~ log((.5+1/Node.Age))+log(Mean.Clim.RR)+log(Mean.CLADS))

tree_file <- treeRead
tree_file$tip.label <- lastDropper(treeRead$tip.label, "_")
rownames(nodeSheet) <- nodeSheet$uniqueNames
pglsDat_earthAlternative <- comparative.data(phy = tree_file, data = nodeSheet[,-c(14)], names.col = "uniqueNames") 
clads_bio_PGLS <- summary(caper::pgls(data = pglsDat_earthAlternative, formula = cladsFormula, lambda = "ML", kappa = "ML", delta = "ML"))
color_bio_PGLS<- summary(caper::pgls(data = pglsDat_earthAlternative, formula = colorFormula, lambda = "ML", kappa = "ML", delta = "ML"))
morph_bio_PGLS<- summary(caper::pgls(data = pglsDat_earthAlternative, formula = morphFormula, lambda = "ML", kappa = "ML", delta = "ML"))
clim_bio_PGLS<- summary(caper::pgls(data = pglsDat_earthAlternative, formula = climFormula, lambda = "ML", kappa = "ML", delta = "ML"))

summary(caper::pgls(data = pglsDat_earthAlternative, formula = combinedFormula, lambda = "ML", kappa = "ML", delta = "ML"))

par(mfrow = c(2,2))
plot(caper::pgls(data = pglsDat_earthAlternative, formula = combinedFormula, lambda = "ML", kappa = "ML", delta = "ML"))
dev.off()

cladsPlot <- ggplot(nodeSheet, aes(y=(Mean.CLADS), x=log(.5+1/nodeSheet$Node.Age))) + 
    geom_point(aes(col = To.Area), size = 4)+
    geom_text(aes(x = -.3, y = .4, label = paste("Earth R-sq =", cladsStat$rsq,"\n PGLS R-sq =", clads_bio_PGLS$adj.r.squared))) + 
    #geom_text(aes(label = Clade.Description), hjust = 0, nudge_x = 0.05,  size = 2) +
    geom_smooth(stat= "smooth", method = "gam", color="blue", size = 1, formula = y ~ s(x)) +
    scale_color_manual(values=c("#3557a0", "#2e9ccc", "#64cdda", "#2a8270", "#5aba58", "#f58643", "#ebe654", "#ea2636", "gray")) +
    theme_minimal() + xlab("log Recency of Colonization") + ylab("CLADS Div. Rate")

colPlot <- ggplot(nodeSheet, aes(y=log(Mean.Plum.RR), x=log(.5+1/nodeSheet$Node.Age))) + 
  geom_point(aes(col = To.Area), size = 4)+
  geom_text(aes(x = -.3, y = .4, label = paste("Earth R-sq =", colorStat$rsq, "\n PGLS R-sq =", color_bio_PGLS$adj.r.squared))) + 
  #geom_text(aes(label = nodeSheet$Clade.Description), hjust = 0, nudge_x = 0.05,  size = 2) +
  geom_smooth(stat= "smooth", method = "gam", color="blue", size = 1, formula = y ~ s(x)) +
  scale_color_manual(values=c("#3557a0", "#2e9ccc", "#64cdda", "#2a8270", "#5aba58", "#f58643", "#ebe654", "#ea2636", "gray")) +
  theme_minimal() + xlab("log Recency of Colonization") + ylab("log Plumage Div. Rate")

morphPlot <- ggplot(nodeSheet, aes(y=log(Mean.Morph.RR), x=log(.5+1/nodeSheet$Node.Age))) + 
  geom_point(aes(col = To.Area), size = 4)+
  #geom_text(aes(label = nodeSheet$Clade.Description), hjust = 0, nudge_x = 0.05,  size = 2) +
  geom_text(aes(x = 0, y = .4, label = paste("R-sq =", morphStat$rsq, "\n PGLS R-sq =", morph_bio_PGLS$adj.r.squared))) + 
  geom_smooth(stat= "smooth", method = "gam", color="blue", size = 1, formula = y ~ s(x)) +
  scale_color_manual(values=c("#3557a0", "#2e9ccc", "#64cdda", "#2a8270", "#5aba58", "#f58643", "#ebe654", "#ea2636", "gray")) +
  theme_minimal() + xlab("log Recency of Colonization") + ylab("log Morphometric Div. Rate")

climPlot <- ggplot(nodeSheet, aes(y=log(Mean.Clim.RR), x=log(.5+1/nodeSheet$Node.Age))) + 
  geom_point(aes(col = To.Area), size = 4)+
  geom_text(aes(x = -.3, y = 2, label = paste("R-sq =", climStat$rsq, "\n PGLS R-sq =", clim_bio_PGLS$adj.r.squared))) + 
  #geom_text(aes(label = nodeSheet$Clade.Description), hjust = 0, nudge_x = 0.05, size = 2) +
  geom_smooth(stat= "smooth", method = "gam", color="blue", size = 1, formula = y ~ s(x)) +
  scale_color_manual(values=c("#3557a0", "#2e9ccc", "#64cdda", "#2a8270", "#5aba58", "#f58643", "#ebe654", "#ea2636", "gray")) +
  theme_minimal() + xlab("log Recency of Colonization") + ylab("log Niche Div. Rate")



pdf("FigureOutput/node_age_plots.pdf", width = 10, height = 10)
gridExtra::grid.arrange(colPlot,cladsPlot, morphPlot, climPlot, ncol=2)
graphics.off()



colClimPlot <- ggplot(nodeSheet, aes(x=log(Mean.Clim.RR), y=log(Mean.Plum.RR))) + 
  geom_point(aes(col = To.Area), size = 4)+
  geom_text(aes(x = -.5, y = -2, label = paste("R-sq =", "X"))) + 
  geom_text(aes(label = nodeSheet$Clade.Description), hjust = 0, nudge_x = 0.05, size = 2) +
  geom_smooth(stat= "smooth", method = "gam", color="blue", size = 1, formula = y ~ s(x)) +
  scale_color_manual(values=c("#3557a0", "#2e9ccc", "#64cdda", "#bc529b", "#2a8270", "#5aba58", "#f58643", "#ebe654", "#ea2636", "gray")) +
  theme_minimal() + xlab("log Niche Div. Rate") + ylab("log Plumage Div. Rate")

colCladsPlot <- ggplot(nodeSheet, aes(x=log(Mean.CLADS), y=log(Mean.Plum.RR))) + 
  geom_point(aes(col = To.Area), size = 4)+
  geom_text(aes(x = -.5, y = -2, label = paste("R-sq =", "X"))) + 
  geom_text(aes(label = nodeSheet$Clade.Description), hjust = 0, nudge_x = 0.05, size = 2) +
  geom_smooth(stat= "smooth", method = "gam", color="blue", size = 1, formula = y ~ s(x)) +
  scale_color_manual(values=c("#3557a0", "#2e9ccc", "#64cdda", "#bc529b", "#2a8270", "#5aba58", "#f58643", "#ebe654", "#ea2636", "gray")) +
  theme_minimal() + xlab("log CLADS Rate") + ylab("log Plumage Div. Rate")

climCladsPlot <- ggplot(nodeSheet, aes(x=log(Mean.CLADS), y=log(Mean.Clim.RR))) + 
  geom_point(aes(col = To.Area), size = 4)+
  geom_text(aes(x = -.5, y = -2, label = paste("R-sq =", "X"))) + 
  geom_text(aes(label = nodeSheet$Clade.Description), hjust = 0, nudge_x = 0.05, size = 2) +
  geom_smooth(stat= "smooth", method = "gam", color="blue", size = 1, formula = y ~ s(x)) +
  scale_color_manual(values=c("#3557a0", "#2e9ccc", "#64cdda", "#bc529b", "#2a8270", "#5aba58", "#f58643", "#ebe654", "#ea2636", "gray")) +
  theme_minimal() + xlab("log CLADS Rate") + ylab("log Mean Clim RR")

gridExtra::grid.arrange(colClimPlot,colCladsPlot, climCladsPlot, ncol=3)















