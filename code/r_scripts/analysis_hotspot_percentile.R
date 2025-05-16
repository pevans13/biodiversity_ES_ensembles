## ---------------------------
##
## Script name: analysis_hotspot_percentile.R
##
## Purpose of script: to conduct the hotspot analysis 
##                    
##                    
## ------------ Notes --------------  ##
## As of 23/07/24: ensemble the models normalised models first, get the median
## from each pixel, renormalise, then get the top percentiles
## ------------ ----- --------------  ##
##
## Run after: ??
##
## Run before: ??
##
## Specific numbered tasks:
## 1 - list files
## 2 - get 80, 85, 90, 95 percentiles for all data
## 3 - hotspot analysis
## 4 
##
## list of final outputs:
##    percentiles_Bio_ES.csv - contains 80 - 95 percentiles of all continents for BD and ES
##    
##
## Author: Dr. Paul M. Evans
##
## Date Created: 2024-06-14
##
## Copyright (c) Paul M. Evans, 2024
## Email: paueva@ceh.ac.uk
## Git: https://github.com/pevans13

cat('\n\n####################################'
    , '\n\nstarting ... analysis_hotspot_percentile.R'
    , '\n\n####################################\n\n')

#### 0 - paths ####
# paths for the results so far, which are stored in the 'Normised' directories
bioNormPath <- "data/intermediate_outputs/biodiversity/normalised"
esNormPath <- "data/intermediate_outputs/ecosystem_services/normalised"

#### 0 - load functions ####
# Function to create new column names based on their position
renameSeqFunc <- function(names) {
  paste0("map", seq_along(names), "_1s")
}

# Function to create new column names based on their position
renameSeqNames <- function(names) {
  paste0("name_map", seq_along(names))
}

# Function to create new column names based on their position, specifically for final map
renameSeqFinal <- function(names) {
  paste0("map_X_", seq_along(names))
}

# Function to create new column names based on their position, specifically for final map (ES)
renameSeqES <- function(names) {
  paste0("map_Xes_", seq_along(names))
}
# Function to create new column names based on their position, specifically for final map (BD)
renameSeqBD <- function(names) {
  paste0("map_Xbd_", seq_along(names))
}

#### 1 - list all the files ####
## ------------ Notes --------------  ##
## Here, only continental data will be assessed
## ------------ ----- --------------  ##
# biodiversity
bioNormFiles <- list.files(bioNormPath
                           , recursive = T
                           , pattern = ".tif$", full.names = T)
bioNormFiles

# ecosystem services
esNormFiles <- list.files(esNormPath
                          , recursive = T
                          , pattern = ".tif$", full.names = T)
esNormFiles

## combine all
allNormFiles <- c(bioNormFiles, esNormFiles)

#### 2 - top quantiles ####
## ------------ Notes --------------  ##
## in this section, the top percentiles (80%, 85%, 90%, amd 95%)
## need to be extracted, and that needs to used to compared Es and biodiversity
## models outputs via hotspot analysis
## ------------ ----- --------------  ##

## first, see if final table already produced
if(file.exists(file.path("data", "intermediate_outputs", "percentiles_Bio_ES.csv"))){
  esBdOut <- fread(file.path("data", "intermediate_outputs", "percentiles_Bio_ES.csv"))
} else {
  
  ##### get percentiles for biodviersity #####
  bioOut <- lapply(1:length(bioNormFiles), function(nm){
    
    tic("One rotation")
    # get name
    nmIn <- basename(bioNormFiles[[nm]])
    cat(nm, "|", nmIn, "\n")
    
    # get resolution
    resNm <- basename(dirname(dirname(bioNormFiles[[nm]])))
    # get continent
    cont <- gsub(".*norm_(.+).tif.*", "\\1", bioNormFiles[[nm]])
    # get richness measure
    richness <- substring(nmIn, 1, 2)
    
    cat("reading in", nmIn, "...\n")
    r <- rast(bioNormFiles[[nm]])
    max(values(r), na.rm = T)
    min(values(r), na.rm = T)
    
    # get the quartiles, every 5%
    cat("currently quantiling", nmIn, "...\n")
    tic("quantile1")
    probabilities <- seq(0, 1, by = 0.05)
    quartiles <- quantile(values(r), probabilities, na.rm = T)
    print(quartiles)
    ## get just the top ones
    perc95 <- quartiles[[20]]
    perc90 <- quartiles[[19]]
    perc85 <- quartiles[[18]]
    perc80 <- quartiles[[17]]
    
    xdf <- bind_cols(name = gsub(".tif", "", nmIn)
                     , resolution = resNm
                     , continent = cont
                     , richness = richness
                     , perc80 = perc80, perc85 = perc85
                     , perc90 = perc90, perc95 = perc95)
    
    return(xdf)
  }) %>% bind_rows()
  head(bioOut)
  
  ##### get percentiles for ecosystem services #####
  esOut <- lapply(1:length(esNormFiles), function(nm){
    
    tic("One rotation")
    # get name
    nmIn <- basename(esNormFiles[[nm]])
    cat(nm, "|", nmIn, "\n")
    
    # get resolution
    resNm <- basename(dirname(dirname(dirname(esNormFiles[[nm]]))))
    # get continent
    cont <- gsub(".*norm_(.+).tif.*", "\\1", esNormFiles[[nm]])
    if(nchar(cont) > 15){
      print("cont being shortened")
      Sys.sleep(5)
      cont <- gsub(".*norm_(.+).tif.*", "\\1", cont)
      stop("continent too vast")
    }
    # get richness measure
    es <- basename(dirname(esNormFiles[[nm]]))
    
    cat("reading in", nmIn, "...\n")
    r <- rast(esNormFiles[[nm]])
    max(values(r), na.rm = T)
    min(values(r), na.rm = T)
    
    # get the quartiles, every 5%
    cat("currently quantiling", nmIn, "...\n")
    tic("quantile1")
    probabilities <- seq(0, 1, by = 0.05)
    quartiles <- quantile(values(r), probabilities, na.rm = T)
    print(quartiles)
    ## get just the top ones
    perc95 <- quartiles[[20]]
    perc90 <- quartiles[[19]]
    perc85 <- quartiles[[18]]
    perc80 <- quartiles[[17]]
    
    xdf <- bind_cols(name = gsub(".tif", "", nmIn)
                     , resolution = resNm
                     , continent = cont
                     , ES = es
                     , perc80 = perc80, perc85 = perc85
                     , perc90 = perc90, perc95 = perc95)
    
    return(xdf)
  }) %>% bind_rows()
  head(esOut)
  
  ## combine both ES and BD outputs, to get an easy reference df for both
  esBdOut <- bioOut %>%
    bind_rows(esOut) %>%
    # and get the metric used into one column
    mutate(metric = ifelse(is.na(ES), richness, ES)) %>%
    dplyr::select(-c(ES, richness)) %>%
    relocate(name, resolution, continent, metric)
  head(esBdOut)
  
  ### save the result
  fwrite(esBdOut
         , file.path("data", "intermediate_outputs", "percentiles_Bio_ES.csv")
         , row.names = F)
} # end of file.exists

#### 3 - task3 ####
## ------------ Notes --------------  ##
## This section focuses on assessing biodiversity and ecosystem services against
## each other. This is currently done just using the crossovers between
## the top percentiles from the normalised versions of all data, for each
## continent. 
## ------------ ----- --------------  ##

# see the dataset
head(esBdOut)

## get all combinations of resolutions, and continents
combos <- esBdOut %>%
  tidyr::nest(.by = c(resolution, continent))
str(combos)

### get all unique metrics, and split them into ES and BD
ug <- unique(esBdOut$metric)
metricsBD <- ug[grepl("sr|pr", ug)] # biodiversity
metricsES <- ug[!grepl("sr|pr", ug)] # non-biodiversity (i.e., ecosystem services)

#### loop through the combos and use them to extract the correct map
for(i in 1:nrow(combos)){
  # if(i == 1){stop("starting loop")}
  
  # get res
  xRes <- combos[[i, "resolution"]]
  # get continent
  xCont <- combos[[i, "continent"]]  
  
  ## get table
  xdf <- combos[[i, "data"]] %>% as.data.frame()
  
  ### further split the data by metric
  xdf2 <- xdf %>% tidyr::nest(.by = c(metric))
  head(xdf2)
  
  for(bd in metricsBD){ # loop through BD metrics
    for(esx in metricsES){  # loop through ES metrics
      
      cat("\nThis location is", xCont)
      cat(" at a resolution of", xRes)  
      cat("\n   with the ES of", esx) 
      cat(" and the BD of", bd, "\n") 
      
      # create a name as the save file path
      saveName <- file.path("results", "tables"
                            , paste0("hs_"
                                     , xCont, gsub("res", "", xRes)
                                     , "_", esx, "_", bd
                                     , ".csv"))
      
      # reduce to ES and BD that will be compared
      xdf3 <- xdf2 %>%
        filter(metric %in% c(bd, esx))
      stopifnot(nrow(xdf3) == 2) # check
      
      #### use all the models of SR, and grazing (as an example)
      xdf2.bd <- xdf3 %>% 
        filter(metric == bd) %>% purrr::pluck("data") %>% as.data.frame()
      head(xdf2.bd)
      
      xdf2.es <- xdf2 %>% 
        filter(metric == esx) %>%
        purrr::pluck("data") %>%
        as.data.frame()
      head(xdf2.es)
      # get length of combined names
      lenlen <- length(c(xdf2.es$name, xdf2.bd$name))
      
      ##### get all combinations of these two model outputs
      l <- rep(list(0:1), lenlen) # get matrix with all possible combinations
      xxx <- expand.grid(l) %>% as.data.frame()
      names(xxx) <- gsub("(.+)_norm.*", "\\1", c(xdf2.es$name, xdf2.bd$name))
      ## get number of models in each df
      nMods1 <- 1:nrow(xdf2.es)
      nMods1
      nMods2 <- (max(nMods1) + 1) : (max(nMods1) + nrow(xdf2.bd))
      nMods2
      
      # get two sums, one of each metric
      xxx.sums <- xxx %>%
        mutate(metEsSum = rowSums(dplyr::select(., nMods1), na.rm = TRUE)
               , metBdSum = rowSums(dplyr::select(., nMods2), na.rm = TRUE)
               , met_ES_BD_Sum = rowSums(dplyr::select(., 1:ncol(.)), na.rm = TRUE),
               .keep = "none")
      head(xxx.sums)
      # keep all rows that sum to above 2, have at least one entry in both metrics, and match
      keepRows <- which(xxx.sums$met_ES_BD_Sum >= 2 & 
                          xxx.sums$metEsSum > 0 & xxx.sums$metBdSum > 0 &
                          xxx.sums$metEsSum == xxx.sums$metBdSum)
      keepRows[1:5]
      xxx.sums[keepRows[1:5], ]
      xxx.keeps <- xxx %>%
        slice(c(keepRows))
      rm(xxx, xxx.sums, keepRows
         , xdf2.es, xdf2.bd) # to save space
      gc()
      
      ## ------------ Notes --------------  ##
      ## go through each row, extracting names of the models to use.
      ## Extract those results for 10 km2 and 7,774 km2 res
      ## Determine the crossover in terms of hotspots for BD and ES
      ## ------------ ----- --------------  ##
      
      stopCluster(cl)
      
      ## use a cluster to do the next bit
      # Set the number of cores
      nCores <- detectCores() - 1
      cl <- makeCluster(nCores)
      registerDoParallel(cl)
      
      # set the change amount
      nr <- 50
      nrSeq <- c(seq(0, nrow(xxx.keeps), nr), nrow(xxx.keeps))
      nrSeq
      
      # for(pseq in 1:length(nrSeq)){
      #   cat(pp1 <- nrSeq[pseq], "to", pp2 <- nrSeq[pseq+1], "\n")
      
      # foreach(pseq = 1:(length(nrSeq) - 1)
      foreach(pseq = 11:20
              , .combine = 'c', .packages = c("data.table", "dplyr", "janitor", "pbapply", "terra")) %dopar% {
                pp1 <- nrSeq[pseq]
                pp2 <- nrSeq[pseq + 1]
                
                # define the next save point
                savePoint <- gsub(".csv", paste0(pp2, ".csv"), saveName)
                
                if(pseq != 1){
                  # define the previous save point
                  prevSavePoint <- gsub(paste0(pp2, ".csv"), paste0(pp2 - nr, ".csv"), savePoint)
                } else {
                  prevSavePoint <- NULL
                }
                
                rwsOut <- pblapply((pp1 + 1):(pp2), function(rw){
                  # rwsOut <- lapply(1:4, function(rw){
                  # if(rw == 1){stop("starting loop [rw]")}
                  
                  cat("\n----------\nThis row is", rw)
                  
                  # extract the names that have '1' in
                  nmExtract <- names(xxx.keeps)[which(xxx.keeps[rw, ] == 1)]
                  nmExtract
                  cat(" which features", length(nmExtract), "maps", "\n")
                  
                  # get those maps - based on model
                  mapsExtract <- allNormFiles[grepl(paste(nmExtract, collapse = "|"), allNormFiles)]
                  ## based on continent / resolution
                  mapsExtract <- mapsExtract[grepl(xCont, mapsExtract)]
                  mapsExtract <- mapsExtract[grepl(xRes, mapsExtract)]    
                  mapsExtract
                  
                  # get the percentiles from df
                  dfExtract <- xdf %>%
                    filter(grepl(paste(nmExtract, collapse = "|"), name))
                  
                  ## load maps in, with percentiles
                  mapsAll <- lapply(nmExtract, function (x){
                    mapX <- rast(mapsExtract[grepl(x, mapsExtract)])
                    return(mapX)
                  })
                  
                  percAll <- lapply(nmExtract, function (x){
                    percX <- dfExtract %>% filter(grepl(x, name)) %>% dplyr::select(contains("perc"))
                    return(percX)
                  })
                  
                  mmOut <- lapply(1:length(nmExtract), function(mm){
                    
                    cat("\nThis nmExtract is", mm, "[", nmExtract[[mm]], "]",  "\n")
                    
                    percApply <- lapply(names(percAll[[1]]), function(perc){
                      
                      # cat(" | This perc is", perc)
                      
                      # adjust the maps based on whether it was > first percentile - i.e., get hotspots
                      table(as.vector(mapsAll[[mm]]))
                      r <- mapsAll[[mm]]
                      r[is.nan(r)] <- 0
                      mapLogical <- as.numeric(r > percAll[[mm]][, perc])
                      # plot(mapLogical)
                      map1hot1s <- table(as.vector(mapLogical))[[2]]
                      
                      return(list(mapLogical, map1hot1s, perc, nmExtract[[mm]]))
                    })
                    
                    return(percApply)
                  }) # end 'mm'
                  
                  result <- lapply(1:4, function(ip) {
                    
                    # check percentiles
                    pipi <- mmOut[[1]][[ip]][[3]]
                    cat("This perc is", pipi, "\n")
                    
                    n <- length(mmOut) # get length of all the maps
                    
                    # see whether BD or ES
                    allFrom <- sapply(1:n, function(i) mmOut[[i]][[ip]][[4]])
                    allFrom
                    
                    nES <- which(!grepl("^sr_|^pr_", allFrom)) # get length of ES maps
                    nBD <- which(grepl("^sr_|^pr_", allFrom)) # get length of BD maps
                    
                    # get the summed result when including all current rasters - BD
                    allRastersBD <- sapply(nBD, function(i) mmOut[[i]][[ip]][[1]])
                    
                    ## Summing all SpatRaster objects in the list
                    combinedMapBD <- Reduce(`+`, allRastersBD)
                    # plot(combinedMapBD)
                    valOutXbd <- table(as.vector(combinedMapBD)) %>% as.data.frame() %>%
                      slice(2:nrow(.)) %>% t() %>%
                      as.data.frame() %>%
                      janitor::row_to_names(row_number = 1) %>%
                      rename_with(renameSeqBD)
                    
                    # get the summed result when including all current rasters - ES
                    allRastersES <- sapply(nES, function(i) mmOut[[i]][[ip]][[1]])
                    ## Summing all SpatRaster objects in the list
                    combinedMapES <- Reduce(`+`, allRastersES)
                    # plot(combinedMapES)
                    valOutXes <- table(as.vector(combinedMapES)) %>% as.data.frame() %>%
                      slice(2:nrow(.)) %>% t() %>%
                      as.data.frame() %>%
                      janitor::row_to_names(row_number = 1) %>%
                      rename_with(renameSeqES)
                    
                    # get the summed result when including all current rasters
                    allRasters <- sapply(1:n, function(i) mmOut[[i]][[ip]][[1]])
                    ## Summing all SpatRaster objects in the list
                    combinedMap <- Reduce(`+`, allRasters)
                    # plot(combinedMap)
                    valOutX <- table(as.vector(combinedMap)) %>% as.data.frame() %>%
                      slice(2:nrow(.)) %>% t() %>%
                      as.data.frame() %>%
                      janitor::row_to_names(row_number = 1) %>%
                      rename_with(renameSeqFinal)
                    
                    # get all the original cell numbers that > percentile cut-off
                    allCutoffs <- sapply(1:n, function(i) mmOut[[i]][[ip]][[2]])
                    valOut1 <- allCutoffs %>% t() %>% as.data.frame() %>%
                      rename_with(renameSeqFunc)
                    
                    # get names
                    valOutNm <- nmExtract %>% t() %>% as.data.frame() %>%
                      rename_with(renameSeqNames)
                    
                    hotSpotsDf <- bind_cols(continent = xCont
                                            , res = xRes
                                            , perc = pipi
                                            , valOutNm
                                            , valOut1
                                            
                                            # final values for different categories
                                            , valOutXbd
                                            , valOutXes
                                            , valOutX)
                    
                    return(hotSpotsDf)
                  }) %>% bind_rows()
                  
                  return(result)
                }) %>% bind_rows()
                
                # convert all columns that begin with 'map' to numeric
                rwsOut <- rwsOut %>%
                  mutate(across(starts_with("map_"), as.numeric))
                # str(rwsOut)
                
                # save the semi-file - first checking the previous one and adding to it
                if (!is.null(prevSavePoint) && file.exists(prevSavePoint)) {
                  prevandNew <- fread(prevSavePoint, stringsAsFactors = FALSE) %>%
                    bind_rows(., rwsOut, .id = "id")
                  fwrite(prevandNew, savePoint, row.names = FALSE)
                  file.remove(prevSavePoint)
                } else {
                  fwrite(rwsOut, savePoint, row.names = FALSE)
                }
                
              } # end 'pseq'
      
      stopCluster(cl)
      
      # save the file
      fwrite(rwsOut, saveName, row.names = F)
      
    }} # end of 'bd' and 'esx'
  stop()
} # 'i' end

