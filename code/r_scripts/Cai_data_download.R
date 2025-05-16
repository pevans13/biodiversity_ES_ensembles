## ---------------------------
##
## Script name: Cai_data_download.R
##
## Purpose of script: download spatial data from Cai et al (2023) paper
##                    
## ------------ Notes --------------  ##
## The original data were actually downloaded from
## https://gift.uni-goettingen.de/shiny/predictions/
## They were imported as SF data files, saved in the form of RData,
## and in the LongLat projection
## They were saved: .../data/raw/vascular_plants/Cai_2023
## ------------ ----- --------------  ##
##
## Run after: continent and country creation scripts 
##
## Run before: any analysis scripts
##
## Specific numbered tasks:
## 1 - convert downloaded Cai data from RData to raster
##
## list of final outputs:
##    Geopackages of data downloaded from Cai et al (2023)
##
## Author: Dr. Paul M. Evans
##
## Date Created: 2024-04-18
##
## Copyright (c) Paul M. Evans, 2024
## Email: paueva@ceh.ac.uk
## Git: https://github.com/pevans13

cat("\n\n####################################"
    , "\nstarting Cai_data_download.R"
    , "\n####################################\n\n"
    , sep = "\n")
# stop("Cai_data_download.R")

#### 0 - paths ####
path.data <- "data/raw/vascular_plants/Cai_2023"

#### 1 - biodiversity files ####
# list biodiversity files downloaded from Cai
biodList <- list.files(path.data
                       , pattern = ".RData$"
                       , full.names = T)
biodList

# open and save all as sf files (specifically '.gpkg')
for(i in 1:length(biodList)){
  
  cat("working with", basename(biodList[i]), "| number", i, "\n")
  
  # open file and save as gpkg
  load(biodList[i]) # this should produce an object called 'predictions_grid'
  head(predictions_grid)
  class(predictions_grid) # should be sf
  plot(predictions_grid[11])
  
  # convert to WGS84 to match ES models
  predictions_grid.wgs <- predictions_grid %>%
    st_transform(4326)
  
  # get new name (with projection changed)
  nn <- gsub("LongLat", ""
             , paste0(gsub(".RData", "", basename(biodList[i])), "wgs"))
  # save the sf file
  cat("\twriting", paste0(nn, ".gpkg"), "\n")
  st_write(predictions_grid.wgs
           , file.path(path.data
                       , paste0(nn, ".gpkg"))
           , append = F, quiet = T)
}

# write readme
#### x - write readme ####
readmePath <- file.path(path.data, "readme.md")
# Create and open the file connection
fileConn <- file(readmePath)
writeLines(c(
  "Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13"
  , "Date created: 2024-06-05"
  , paste0("Last update:",  format(Sys.Date()))
  , "Produced by 'Cai_data_download.R'

The files in this folder were all downloaded from https://gift.uni-goettingen.de/shiny/predictions
The file name is in the form '[unit]_[method]_prediction_[resolution]_[projection]_wgs.RData' where:
	[unit] = 'sr' (species richness), or 'pr' (phylogenetic richness)
	[method] = 	'GAM' (Generalised Additive Model)
				'GAM-geographic-coordinates' (GAM with spline of geographic coordinates)
				'GLMsimplified-interactions' (Generalised Linear Model that included interaction terms)
				'Neural-networks'
				'Neural-networks-trend-surface' (NN with cubic polynomial trend surface)
				'Random-Forest'
				'Random-Forest-Floristic-Kingdom' (RF with floristic kingdom [regions of plants])
				'Random-Forest-trend-surface' (RF with cubic polynomial trend surface)
				'XGBoost' (eXtreme Gradient Boosting)
				'XGBoost-Floristic-Kingdom' (XGBoost with floristic kingdom)
				'XGBoost-trend-surface' (XGBoost with cubic polynomial trend surface)
				'Ensemble' (ensemble predictions weighted by cross-validation values)
	[resolution] = 7774 km2
	[projection] = wgs (World Geodetic System 1984; EPSG: 4326)")
  , fileConn)
close(fileConn)