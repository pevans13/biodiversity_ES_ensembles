## ---------------------------
##
## Script name: correlation_test.R
##
## Purpose of script: To determine which combination of continent and ES has the
##                    lowest correlation, for the purpose of using that information
##                    to determine the likely number of runs required for all the 
##                    research.
##                    
## ------------ Notes --------------  ##
## This is based on ES at the continent level
## ------------ ----- --------------  ##
##
## Run after: data_to_continent_extent.R
##
## list of final outputs:
##    correlation_summary_BD.csv <- summary of the correlations for species richness
##    correlation_summary_ES.csv <- summary of the correlations for ecosystem services
##
## Author: Dr. Paul M. Evans
##
## Date Created: 2024-07-25
##
## Copyright (c) Paul M. Evans, 2024
## Email: paueva@ceh.ac.uk
## Git: https://github.com/pevans13

cat('\n\n####################################'
    , '\n\n      correlation_test.R'
    , '\n\n####################################\n\n')
# stop("correlation_test.R")

#### 0 - paths ####
bioNormPath <- "data/intermediate_outputs/biodiversity/normalised"
esNormPath <- "data/intermediate_outputs/ecosystem_services/normalised"

#### 0 - load info ####
resolutions <- c(100) # resolution of 100 km

#### 1 - list all the files ####
# biodiversity
bioNormFiles <- list.files(bioNormPath
                           , pattern = ".tif$", full.names = T
                           , recursive = T)
bioNormFiles

# ecosystem services
esNormFiles <- list.files(esNormPath
                          , recursive = T
                          , pattern = ".tif$", full.names = T)
esNormFiles

#### 2 - determine correlations between all of the layers - species ####
## ------------ Notes --------------  ##
## Go through all the continents and different resolutions
## ------------ ----- --------------  ##
continents <- unique(gsub(".*norm_(.+).tif$", "\\1", bioNormFiles))
stopifnot(length(continents) == 6)

for(rr in resolutions){
  for(co in continents){
    
    cat(co, rr, sep = " | ")
    cat("\n")
    
    # create names
    saveName <- file.path("images", paste0("hist_corr_Spp_", rr, "_", co,".jpg"))
    graphTitle <- paste0("Pearson correlations histogram between all\nthe winsorised, normalised Cai et al.\n(2023) plant species richness data, at ", rr, " km2")
    
    corrList <- bioNormFiles[grepl(co, bioNormFiles) & grepl(rr, bioNormFiles) & grepl("^sr_", basename(bioNormFiles))]
    # stack the resolution and continent matches
    bioStack <- rast(corrList) %>% c()
    str(bioStack)
    
    # Compute the correlation matrix
    corMatrix <- layerCor(bioStack, fun=cor, use = "complete.obs")
    corMatrix[lower.tri(corMatrix, diag = TRUE)] <- NA
    # print(corMatrix)
    
    jpeg(saveName, res = 200
         , height = 20, width = 30, units = "cm")
    hist(corMatrix
         , main = graphTitle
         , xlab = "Pearson correlation coefficient"
         , breaks = seq(min(corMatrix, na.rm = T), max(corMatrix, na.rm = T)
                        , (max(corMatrix, na.rm = T) - min(corMatrix, na.rm = T))/50))
    dev.off()
    
    if(rr == resolutions[1] & co == continents[1]){
      # save the table
      outTable <- bind_cols(continent = co
                            , resoution = rr
                            , metric = "species richness"
                            , lowest_pearson = min(corMatrix, na.rm = T)
                            , mean_pearson = mean(corMatrix, na.rm = T)
                            , median_pearson = median(corMatrix, na.rm = T)
                            , highest_pearson = max(corMatrix, na.rm = T)
      )
    } else {
      outTable <- bind_rows(outTable, 
                            bind_cols(continent = co
                                      , resoution = rr
                                      , metric = "species richness"
                                      , lowest_pearson = min(corMatrix, na.rm = T)
                                      , mean_pearson = mean(corMatrix, na.rm = T)
                                      , median_pearson = median(corMatrix, na.rm = T)
                                      , highest_pearson = max(corMatrix, na.rm = T)
                            ))
    }
  }
}
## save the correlation summary table
fwrite(outTable
       , file.path("results", "tables"
                   , paste0("correlation_summary_BD.csv"))
       , row.names = F)


#### 3 - determine correlations between all of the layers - ES ####
## ------------ Notes --------------  ##
## Go through all the continents and different resolutions, and ES
## ------------ ----- --------------  ##
ESs <- esNormFiles %>% dirname() %>% basename() %>% unique()

for(rr in resolutions){
  for(co in continents){
    for(esx in ESs){
      
      cat(co, rr, esx, sep = " | ")
      cat("\n")
      
      # create names
      saveName <- file.path("images", paste0("hist_corr_", esx, "_", rr, "_", co,".jpg"))
      graphTitle <- paste0("Pearson correlations histogram between all\nthe winsorised, normalised ", esx, " data, at ", rr, " km2")
      
      corrList <- esNormFiles[grepl(co, esNormFiles) & grepl(rr, esNormFiles) & grepl(esx, esNormFiles)]
      corrList
      # stack the resolution and continent matches
      esStack <- rast(corrList) %>% c()
      str(esStack)
      
      # Compute the correlation matrix
      corMatrix <- layerCor(esStack, fun=cor, use = "pairwise.complete.obs")
      corMatrix[lower.tri(corMatrix, diag = TRUE)] <- NA
      # print(corMatrix)
      
      jpeg(saveName, res = 200
           , height = 20, width = 30, units = "cm")
      hist(corMatrix
           , main = graphTitle
           , xlab = "Pearson correlation coefficient"
           , breaks = seq(min(corMatrix, na.rm = T), max(corMatrix, na.rm = T)
                          , (max(corMatrix, na.rm = T) - min(corMatrix, na.rm = T))/50))
      dev.off()
      
      if(rr == resolutions[1] & co == continents[1]){
        # save the table
        outTable <- bind_cols(continent = co
                              , resoution = rr
                              , metric = esx
                              , lowest_pearson = min(corMatrix, na.rm = T)
                              , mean_pearson = mean(corMatrix, na.rm = T)
                              , median_pearson = median(corMatrix, na.rm = T)
                              , highest_pearson = max(corMatrix, na.rm = T)
        )
      } else {
        outTable <- bind_rows(outTable, 
                              bind_cols(continent = co
                                        , resoution = rr
                                        , metric = esx
                                        , lowest_pearson = min(corMatrix, na.rm = T)
                                        , mean_pearson = mean(corMatrix, na.rm = T)
                                        , median_pearson = median(corMatrix, na.rm = T)
                                        , highest_pearson = max(corMatrix, na.rm = T)
                              ))
      }
    }
  }
}
## save the correlation summary table
fwrite(outTable
       , file.path("results", "tables"
                   , paste0("correlation_summary_ES.csv"))
       , row.names = F)

#### 4 - write readme ####
readmePath <- file.path("results", "tables", "readme_correlation_summary.md")
# Create and open the file connection
fileConn <- file(readmePath)
writeLines(c(
  "Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13",
  "Date created: 2024-07-26",
  paste0("-----------------------------"),
  paste0("All files in the below section were created using the 'correlation_test.R' script"
         , "\nSection last updated: ", format(Sys.Date())),
  paste0("-----------------------------\n"),
  
  paste0(
    "Files:
  'correlation_summary_ES.csv'
      Contains the results of pairwise correlations when all of the models for a single ecsoystem service (ES) where correlated with each other. The file gives summary information of min, max, mean, and median of the correlation per ES, continent, and resolution
  'correlation_summary_BD.csv'
      Contains the results of pairwise correlations when all of the models for species richness where correlated with each other. The file gives summary information of min, max, mean, and median of the correlation per continent, and resolution.

Columns of both files:
  'continent' = the continent that was the focus of the analysis
  'resolution' = the resolution of the data that were correlated with each other
  'metric' = the result that was correlated
  '[mm]_pearson', where 'mm' was either:
      'lowest' = minimum Pearson correlation that was obtained after a summary of all parwise correlations
      'highest' = maximum Pearson correlation that was obtained after a summary of all parwise correlations
      'mean' = mean Pearson correlation that was obtained after a summary of all parwise correlations
      'median' = median Pearson correlation that was obtained after a summary of all parwise correlations
  
Data: Original data for biodiversity: Cai et al. (2023)
      Original data for ES: Willcock et al. (2023)
Original citation: Cai, L., Kreft, H., Taylor, A., Denelle, P., Schrader, J., Essl, F., van Kleunen, M., Pergl, J., Pyšek, P., Stein, A., Winter, M., Barcelona, J. F., Fuentes, N., Inderjit, Karger, D. N., Kartesz, J., Kuprijanov, A., Nishino, M., Nickrent, D., … Weigelt, P. (2023). Global models and predictions of plant diversity based on advanced machine learning techniques. New Phytologist, 237(4), 1432–1445. https://doi.org/10.1111/nph.18533
                   Willcock, S., Hooftman, D. A. P., Neugarten, R. A., Chaplin-Kramer, R., Barredo, J. I., Hickler, T., Kindermann, G., Lewis, A. R., Lindeskog, M., Martínez-López, J., & Bullock, J. M. (2023). Model ensembles of ecosystem services fill global certainty and capacity gaps. Science Advances, 9(14), eadf5492. https://doi.org/10.1126/sciadv.adf5492")
)
, fileConn)
close(fileConn)