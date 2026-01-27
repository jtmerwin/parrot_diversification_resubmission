List of analyses in methods
Updated by JTM on 20.January.2026

1) Estimation of div rates (ELKIN TENORIO)
	This script runs in JULIA and estimates
	diversification rates using CLADS.

2) Processing color and morph data
	color_load_collate
		this take raw color data output from imagej and
		transforms it into usable data.
	color_make_tile_matrix
		this takes the organized data from color_load_collate and outputs an "image"
		representing our color measurements across the phylogeny, as well as tables
		and statistics.
	morphology_load_collate.R
		this script takes raw
		morphology measurements
		and summarized them
		by species / specimen
	calculate_dichromatism.R
		this takes the organized data from color_load_collate and outputs Euclidean
		distances between male and female individuals for each patch and each taxon.
		Also plots this data as a tile matrix.

3) Niche estimation (ELKIN TENORIO)
	3_NicheEstimations.R
		This script estimates environmental Niches from eBird data.

4) DTT plotting
	4_DTT_disparity_thru_time_analysis
		This script runs a windowed DTT analysis to show trait disparity 
		in time windows across the tree
	
5) RRphylo for trait evo rates
	5_RRPHYLO_trait_evolution_rates.R
		This script estimates trait evolution rates using RR_phylo
		Then it runs PGLS analysis of the relationships between
		estimated trait evolution rates and CLADS diversification rates

6) PGLS, ANOVA, simulations
	6_PGLS_ANOVA.R
		This script runs a phylogenetic anova to test whether different continents have
		different diversification or trait evolution rates.
		Then it runs multivariate PGLS models to test the relationship of climate and div. rates
		rates on morphological and plumage color evolution. 
	
7) MUSSCRAT (BRIAN)
	musscrat_commands-24-may-2024.txt
		Shell scripts to run MUSSCRAT on linux with associated nexus files

8) BIOGEOBEARS
	BioGeoBear.R
		Modified template script to run ancestral area reconstruction
		Additional code here paints a phylogeny with most likely ranges.
	ChordgramForBiogeoBears.R
		This script builds a chordgram figure from a matrix of area transitions estimated
		using biogeobears.

		