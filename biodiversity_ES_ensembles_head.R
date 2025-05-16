## ---------------------------
##
## Script name: biodiversity_ES_ensembles_head.R
##
## Purpose of script: This is the first script for the ensembles uncertainty project.
##                    This script will set up all the remaining scripts
##                    with a consistent date and structure
##                    
## ------------ Notes --------------  ##
## Below, the user should choose which aspect of the model they want to run
## ------------ ----- --------------  ##

rm(list = ls())

## ---------------------------
options(scipen = 6, digits = 4) # for non-scientific notation
## ---------------------------

#### User input ####
##### Which elements should be run? #####
## ------------ Notes --------------  ##
## The below lines should be assigned 'TRUE' (or 'T') if you want them to run
## ------------ ----- --------------  ##
# Step: create continents into the form for the analysis
createContinents <- F
# Step: create countries into the form for the analysis
createCountries <- F
# Step: take the plant data from Cai et al. (2023) and convert it to 1 km
preparePlantDataCai <- F
# Step: create raster templates for the biodiversity data
rasterTemplates <- F
# Step: convert all data to 100 km2, which is a compromise between 7,774 km2 and 1 km2
data_to_continent_extent <- F
## for 'data_to_100km2' do you need global and / or continental results
globalAnalysis <- F
continentalAnalysis <- F
# Step: correlate data, to find least correlated
correlateData <- F
# Step: run hotspot analysis using percentiles
hotspotPercentiles <- F
# Step: create the final maps and / or final graphs
graph_creation_final <- T
map_creation_final <- F
# Step: create final figures for main and / or SI
finalFigures <- F
# choose which figures to run
mainFigures <- F
SI <- F
singlePercGraphs <- F
## which percentiles to include?
spg <- c(
  spg70 <- F
  , spg80 <- T
  , spg90 <- F
)

#### 0 - paths ####
pathBioD <- "D:/biodiversity"
pathSaveBio <- "data/intermediate_outputs/biodiversity"
path.esData <- "D:/Carbon Model Layers/Normalised and clipped (43200 x 18600)"
pathDataInter <- "data/intermediate_outputs"

#### load settings ####
## automatic install of packages if they are not installed already
list.of.packages <- c(
  "terra", "raster", "dplyr", "tidyr"
  , "ggplot2", "tidyterra", "gridExtra"
  , "patchwork"
  # , "scater" # for multiplot displays
  , "sf", "stars", "fasterize", "exactextractr"
  , "tictoc", "data.table", "parallel", "pbapply"
  , "matrixStats" # for 'weightedMedian' in dplyr
  , "DescTools" # for use of the 'Winsorize' function [https://search.r-project.org/CRAN/refmans/DescTools/html/Winsorize.html]
  , "beepr" # for sounds when things are finished
  , "rnaturalearth" # for shapefiles
  , "foreach", "doParallel" # for parallel processing
)
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0){
  install.packages(new.packages, dep=TRUE)
}
# loading packages
for(package.i in list.of.packages){
  suppressPackageStartupMessages(library(package.i, character.only = TRUE)
  )
}
# tidy
rm(list.of.packages, new.packages, package.i)
terraOptions(progress=10)
options(dplyr.summarise.inform = F)
rasterOptions(progress = 'text', timer=TRUE)
terraOptions(memfrac=0.5, tempdir = "N:/temp")
tempFilePat <- dirname(rasterTmpFile())
removeTmpFiles(h=.1)
gc()

#### set the structure for the whole project ####
files.toInclude <- c("data/raw", "data/intermediate_outputs"
                     # non-data architecture
                     , "code/r_scripts", "code/python_scripts"
                     , "results/tables", "results/reports"
                     , "images"
                     , "docs/references", "docs/spreadsheets", "docs/presentations"
                     
                     # continental saves
                     , "data/intermediate_outputs/ecosystem_services/normalised/res100/continents"
                     , "data/intermediate_outputs/ecosystem_services/normalised/res100/countries"
                     , "data/intermediate_outputs/biodiversity/normalised/res100/continents"
                     , "data/intermediate_outputs/biodiversity/normalised/res100/countries"
                     , "data/checks"
                     
                     # global saves
                     , "data/intermediate_outputs/ecosystem_services/normalised/res100/global"
                     , "data/intermediate_outputs/biodiversity/normalised/res100/global"
                     
                     # continent shapefiles
                     , "data/intermediate_outputs/continents/rasters"
                     # country
                     , "data/intermediate_outputs/countries/rasters"
                     
                     # saving figures
                     , "results/figures", "results/figures/single_item", "results/figures/combined"
                     , "results/figures/maps"
)
pblapply(files.toInclude, function(x) {dir.create(x, showWarnings = F, recursive = T)})
rm(files.toInclude)

#### run elements ####
if(createContinents){
  source(file.path("code", "r_scripts", "create_continents.R"))
}
if(createCountries){
  source(file.path("code", "r_scripts", "create_countries.R"))
}
if(preparePlantDataCai){
  source(file.path("code", "r_scripts", "Cai_data_download.R"))
}
if(rasterTemplates){
  source(file.path("code", "r_scripts", "template_and weightings.R"))
}
if(data_to_continent_extent){
  source(file.path("code", "r_scripts", "data_to_continent_extent.R"))
}
if(correlateData){
  source(file.path("code", "r_scripts", "correlation_test.R"))
}
if(hotspotPercentiles){
  source(file.path("code", "r_scripts", "hotspot_analysis.R"))
}
if(graph_creation_final){
  source(file.path("code", "r_scripts", "uncert_graphs.R"))
}
if(map_creation_final){
  source(file.path("code", "r_scripts", "map_creation_final.R"))
}
if(finalFigures){
  source(file.path("code", "r_scripts", "paper_graphs.R"))
}
