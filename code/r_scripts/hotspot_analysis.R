## ---------------------------
##
## Script name: hotspot_analysis.R
##
## Purpose of script: to run the final hotspot analysis, creating the final CSVs.
##                    
## ------------ Notes --------------  ##
## This script loops over all continents of interest and ecosystem services (ES)
## of interest. The resolution can also be changed, but currently is just 100 km2. 
## ------------ ----- --------------  ##
##  
## Run after: correlation_test.R
##
## Run before:  uncert_graphs.R
##              map_creation_final.R
##
## list of final outputs:
## 	  results/tables/[ecosystem service]100[continent]_sr_rw[row number].csv
##    results/tables/[ecosystem service]100[continent]_sr_END.csv    
##    
## Author: Dr. Paul M. Evans
##
## Date Created: 2024-07-26
##
## Copyright (c) Paul M. Evans, 2024
## Email: paueva@ceh.ac.uk
## Git: https://github.com/pevans13

cat('\n\n####################################'
    , '\n\n   hotspot_analysis.R'
    , '\n\n####################################\n\n')
# stop("hotspot_analysis.R")

#### 0 - load functions ####
# Define a function to normalise raster values
normalise <- function(x) {
  (x - min(values(x), na.rm = TRUE)) / (max(values(x), na.rm = TRUE) - min(values(x), na.rm = TRUE))
}

# Function to create new column names based on their position
renameSeqNames <- function(names) {
  paste0("name_map", seq_along(names))
}

# Function to calculate quartiles for a single layer
calculateQuartiles <- function(raster_layer, probsx) {
  return(quantile(values(raster_layer), probs = c(probsx), na.rm = TRUE))
}

#### 0 - paths ####
## ------------ Notes --------------  ##
## if running on JASMIN the paths are slightly different
## ------------ ----- --------------  ##
# paths for the results so far, which are stored in the 'Normised' directories
bioNormPath <- file.path("data/intermediate_outputs/biodiversity/normalised")
esNormPath <- file.path("data/intermediate_outputs/ecosystem_services/normalised")
# the path to save the results
resTabPath <- file.path("results", "tables")

#### 0 - cores for foreach loop ####
## Set the number of cores
nCores <- detectCores() - 4

#### 0 - running on JASMIN? ####
jasmin <- F
if(jasmin){
  bioNormPath <- file.path("ensembles_uncertainty", bioNormPath)
  esNormPath <- file.path("ensembles_uncertainty", esNormPath)
  resTabPath <- file.path("ensembles_uncertainty", resTabPath)
  nCores <- 12
}

cat("paths for this script:"
    , paste0("bioNormPath = ", bioNormPath)
    , paste0("esNormPath = ", esNormPath)
    , paste0("resTabPath = ", resTabPath)
    , sep = "\n")

#### 0 - choice of ES and continents ####
## ------------ Notes --------------  ##
## in this section you choose which continents and ES to analysis.
## Comment out those you are not interested in, together with their preceding comma
## ------------ ----- --------------  ##

chosenConts <- c(
  "Africa"
  ,
  "Asia"
  ,
  "Europe"
  ,
  "North America"
  ,
  "Oceania"
  ,
  "South America"
)

chosenES <- c(
  "Carbon"
  ,
  "Water"
  ,
  "Recreation"
)

#### 1 - list all the files ####
## ------------ Notes --------------  ##
## Here, only continental data will be assessed
## ------------ ----- --------------  ##
# biodiversity
bioNormFiles <- list.files(bioNormPath
                           , recursive = T
                           , pattern = ".*sr_(.+).tif$", full.names = T)
## remove any global
bioNormFiles <- bioNormFiles[!grepl("global", bioNormFiles)]
# ecosystem services
esNormFiles <- list.files(esNormPath
                          , recursive = T
                          , pattern = ".tif$", full.names = T)
## remove any global
esNormFiles <- esNormFiles[!grepl("global", esNormFiles)]

#### 0 - set defining properties ####
xRes <- "100"
# derive all continents from files
conts <- chosenConts
# derive all individual ES
allES <- chosenES

# set properties for quantiles
probabilities <- seq(0, 1, by = 0.05)

## loop for all continents
for(iCont in conts){
  for(iES in allES){
    
    cat('\n\n####################################---####################################'
        , '\n      ##############################---##############################'
        , paste0('\n              Let ', iES, ' in ', iCont, " [at ", xRes, " km2] begin!")
        , '\n      ##############################---##############################'
        , '\n####################################---####################################\n\n')
    Sys.sleep(5)
    
    # shorten the strings
    iESshort <- substring(iES, 1, 4)
    iContshort <- substring(iCont, 1, 4)
    # create repeating pattern that will be used to save all final data
    finalNameRep <- paste0(iESshort, xRes, iContshort, "_sr")
    
    ## see if the final output has already been created
    if(file.exists(file.path(resTabPath
                             , paste0(finalNameRep, "_END.csv")))){
      cat("The combination of", iES, "and", iCont, "at the resolution of", xRes
          , "already exists...\n")
      Sys.sleep(5)
      
    } else {
      
      # biodiversity
      bioNormFiles <- list.files(bioNormPath
                                 , recursive = T
                                 , pattern = ".*sr_(.+).tif$", full.names = T)
      ## remove any global
      bioNormFiles <- bioNormFiles[!grepl("global", bioNormFiles)]
      # ecosystem services
      esNormFiles <- list.files(esNormPath
                                , recursive = T
                                , pattern = ".tif$", full.names = T)
      ## remove any global
      esNormFiles <- esNormFiles[!grepl("global", esNormFiles)]
      
      ## select only continent
      bioNormFiles <- bioNormFiles[grepl(iCont, bioNormFiles)]
      ## remove the already ensembled one
      bioNormFiles <- bioNormFiles[!grepl("sr_Ensemble", bioNormFiles)]
      ### only 100 resolution
      bioNormFiles <- bioNormFiles[grepl("res100", bioNormFiles)]
      
      ## select only continent
      esNormFiles <- esNormFiles[grepl(iCont, esNormFiles)]
      ### select only recreation
      esNormFiles <- esNormFiles[grepl(paste0("/", iES, "/"), esNormFiles)]
      #### only 100 resolution
      esNormFiles <- esNormFiles[grepl("res100", esNormFiles)]
      ##### see which ES models have values - remove those that are not 
      essx <- pblapply(seq_along(esNormFiles), function(i){
        xIn <- rast(esNormFiles[[i]])
        # print(xIn)
        if(is.infinite(max(values(xIn), na.rm = T))){
          cat("No good")
          keepIn <- "No"
        } else {
          keepIn <- i
        }
        return(keepIn)
      }) %>% unlist() %>% as.vector() %>% as.numeric() %>% na.omit()
      esNormFiles <- esNormFiles[essx]
      
      # combine ES and BD for the for loop later
      allNormFiles <- c(esNormFiles, bioNormFiles)
      
      #### 2 - hotspot analysis ####
      ## ------------ Notes --------------  ##
      ## This section focuses on assessing biodiversity and ecosystem services against
      ## each other. This is currently done just using the crossovers between
      ## the top percentiles from the normalised versions of all data, for each
      ## continent. 
      ## ------------ ----- --------------  ##
      
      ##### 2a - determine all possible combinations of ES and BD models #####
      ## ------------ Notes --------------  ##
      ## This assumes 1 BD x 1 ES, then 2 x 2, then 3 x 3, etc.
      ## ------------ ----- --------------  ##
      
      # get unique models for both
      uBD <- gsub("^sr_(.+)_norm.*", "\\1", basename(bioNormFiles))
      uES <- gsub("^(.+)_norm.*", "\\1", basename(esNormFiles))
      cat("uBD:", uBD, sep = "\n")
      cat("uES:", uES, sep = "\n")
      
      # get all combinations of these two model outputs
      l <- rep(list(0:1), (length(uBD) + length(uES))) # get matrix with all possible combinations
      combos <- expand.grid(l) %>% as.data.frame()
      ## assign names
      names(combos) <- c(uBD, uES)
      ### get number of models in each df
      nMods1 <- 1:length(uBD); nMods1
      nMods2 <- (length(nMods1)+1):(length(nMods1) + length(uES)); nMods2
      
      # get two sums, one for ES and one for BD
      combos.sums <- combos %>%
        mutate(metBdSum = rowSums(dplyr::select(., all_of(nMods1)), na.rm = TRUE)
               , metEsSum = rowSums(dplyr::select(., all_of(nMods2)), na.rm = TRUE)
               , met_ES_BD_Sum = rowSums(dplyr::select(., 1:ncol(.)), na.rm = TRUE),
               .keep = "none")
      head(combos.sums)
      ## keep all rows that sum to above 2, have at least one entry in both metrics, and match
      keepRows <- which(combos.sums$met_ES_BD_Sum >= 2 & 
                          combos.sums$metEsSum > 0 & combos.sums$metBdSum > 0 & # higher than 0
                          combos.sums$metEsSum == combos.sums$metBdSum)
      # keepRows[1:5]
      print(combos.sums[keepRows[1:5], ])
      cdv <- combos.sums[which(combos.sums$metBdSum >= 7 & 
                                 combos.sums$metEsSum == combos.sums$metBdSum), ]
      ### determine the amount of combinations
      combos.sums2 <- combos.sums %>%
        slice(c(keepRows))
      table(as.vector(combos.sums2$metEsSum))
      #### keep the final table that matches the criteria
      combos.keeps <- combos %>%
        slice(c(keepRows))
      ckCheck <- combos.keeps %>%
        mutate(metBdSum = rowSums(dplyr::select(., all_of(nMods1)), na.rm = TRUE)
               , metEsSum = rowSums(dplyr::select(., all_of(nMods2)), na.rm = TRUE)
               , met_ES_BD_Sum = rowSums(dplyr::select(., 1:ncol(.)), na.rm = TRUE),
               .keep = "none")
      ckk <- which(ckCheck$metBdSum >= 7)
      
      # tidy
      rm(combos, combos.sums, keepRows, combos.sums2
         , nMods1, nMods2, l) # to save space
      gc()
      
      ##### 2b - row-by-row analysis #####
      ## ------------ Notes --------------  ##
      ## go through each row, extracting names of the models to use.
      ## Determine the crossover in terms of hotspots for BD and ES.
      ## ------------ ----- --------------  ##
      
      maxRows <- nrow(combos.keeps)
      cat("maxRows:", maxRows, "\n")
      
      # split the analysis, making it save after every user-defined break
      if(iES == "Carbon"){
        userDefineBreak <- ceiling(maxRows/30)
      } else {
        userDefineBreak <- ceiling(maxRows/20)
      }
      cat("userDefineBreak:", userDefineBreak, "\n")
      
      ## create a sequence from it
      userSeq <- c(seq(userDefineBreak, maxRows, userDefineBreak), maxRows)
      
      # stop("before loop")
      
      # determine the highest number that has previously been saved
      prevList <- list.files(file.path(resTabPath)
                             , pattern = finalNameRep)
      ## extract just the numbers
      prevList.no <- gsub(".*rw(.+).csv", "\\1", prevList) %>% as.numeric()
      ### get max
      prevList.max <- unique(prevList.no[which(prevList.no == max(prevList.no, na.rm = T))])
      #### find where that is along the sequence
      if(length(prevList.max) > 0){
        seqStart <- which(userSeq == prevList.max)
        seqStart <- seqStart + 1
      } else {
        seqStart <- 1
      }
      
      for(len in seqStart:length(userSeq)){
        
        # create the start and final numbers
        numStart <- (userSeq[len] - userDefineBreak) + 1
        numEnd <- userSeq[len]
        
        ### check if the final save of this rotation already exists
        saveFile <- file.path(resTabPath
                              , paste0(finalNameRep, "_rw", numEnd, ".csv"))
        saveMap <- file.path(resTabPath
                             , paste0(finalNameRep, "_Maprw", numEnd, ".csv"))
        ### and determine previous
        if(numEnd == maxRows){
          # change start, as to not repeat
          numStart <- userSeq[(length(userSeq)-1)] + 1
          
          # previous table
          saveFilePrev <- file.path(resTabPath
                                    , paste0(finalNameRep, "_rw", userSeq[(length(userSeq)-1)], ".csv"))
          stopifnot(file.exists(saveFilePrev))
          # previous map
          saveMapPrev <- file.path(resTabPath
                                   , paste0(finalNameRep, "_Maprw", userSeq[(length(userSeq)-1)], ".csv"))
          stopifnot(file.exists(saveMapPrev))
          
          # if final, save without row IDs
          saveFile <- file.path(resTabPath
                                , paste0(finalNameRep, "_END.csv"))
          saveMap <- file.path(resTabPath
                               , paste0(finalNameRep, "_Map_END.csv"))
          
        } else if(numStart != 1){
          # previous table
          saveFilePrev <- file.path(resTabPath
                                    , paste0(finalNameRep, "_rw", numEnd - userDefineBreak, ".csv"))
          stopifnot(file.exists(saveFilePrev))
          
          # previous map
          saveMapPrev <- file.path(resTabPath
                                   , paste0(finalNameRep, "_Maprw", numEnd - userDefineBreak, ".csv"))
          stopifnot(file.exists(saveMapPrev))
        }
        
        ## print
        cat("\n")
        cat("Start:", numStart)
        cat(" | End:", numEnd)
        cat("\n")
        # already exist?
        fExists <- file.exists(saveFile)
        #### if it does not already exist, run the loop from this number
        if(!fExists){
          tic("oneRun time for fExist")
          
          # use a cluster to do the next bit
          cl <- makeCluster(nCores)
          cat("Starting foreach loop ...on", nCores, "cores\n")
          cat("Forloop started at", format(Sys.time(), "%Y-%m-%d %H:%M")
              , "[going from", numStart/maxRows * 100, "to", numEnd/maxRows * 100, "%]"
              , "\n")
          
          registerDoParallel(cl)
          on.exit(stopCluster(cl)) # Ensure cluster is stopped if an error occurs
          
          forLoopResults <- foreach(rw = numStart:numEnd
                                    # forLoopResults <- foreach(rw = ckk
                                    , .packages = c("sf", "dplyr", "ggplot2", "terra"
                                                    , "beepr"
                                                    , "DescTools")
                                    # , .combine = rbind
                                    , .errorhandling = "pass") %dopar% {
                                      
                                      # extract the names that have '1' in
                                      nmExtract <- names(combos.keeps)[which(combos.keeps[rw, ] == 1)]
                                      # remove brackets
                                      nmExtract <- gsub("\\(.*", "", nmExtract)
                                      nmExtract
                                      if("GAM" %in% nmExtract){
                                        xw <- which(nmExtract == "GAM")
                                        nmExtract[xw] <- paste0(nmExtract[xw], "_")
                                      }
                                      if("Random-Forest" %in% nmExtract){
                                        xw <- which(nmExtract == "Random-Forest")
                                        nmExtract[xw] <- paste0(nmExtract[xw], "_")
                                      }
                                      if("XGBoost" %in% nmExtract){
                                        xw <- which(nmExtract == "XGBoost")
                                        nmExtract[xw] <- paste0(nmExtract[xw], "_")
                                      }
                                      if("Neural-networks" %in% nmExtract){
                                        xw <- which(nmExtract == "Neural-networks")
                                        nmExtract[xw] <- paste0(nmExtract[xw], "_")
                                      }
                                      if("WaterWorld" %in% nmExtract){
                                        xw <- which(nmExtract == "WaterWorld")
                                        nmExtract[xw] <- paste0(nmExtract[xw], "\\(")
                                      }
                                      nmExtract
                                      ## determine number of pairs
                                      mapsPairs <- length(nmExtract)/2
                                      ### get those maps - based on model
                                      mapsExtract <- allNormFiles[grepl(paste(nmExtract, collapse = "|"), allNormFiles)]
                                      #### separate into BD and ES maps
                                      mapsExtract.bd <- mapsExtract[grepl("^sr_", basename(mapsExtract))]
                                      mapsExtract.es <- mapsExtract[!grepl("^sr_", basename(mapsExtract))]
                                      
                                      stopifnot(length(mapsExtract.bd) == mapsPairs)
                                      stopifnot(length(mapsExtract.es) == mapsPairs)
                                      
                                      ##### deal with raster maps - ensembles #####
                                      # load maps in for biodiversity
                                      mapsBdstack <- rast(mapsExtract.bd)
                                      ## get the ensembled median for each pixel
                                      EnsembleBd <- app(mapsBdstack, fun = median, na.rm = TRUE)
                                      ### normalise to the now-ensembled maps
                                      EnsembleBdNorm <- normalise(EnsembleBd)
                                      #### create another copy with no NANs
                                      EnBDNan <- EnsembleBdNorm
                                      EnBDNan[is.nan(EnBDNan)] <- 0
                                      
                                      # load maps in for biodiversity
                                      mapsESstack <- rast(mapsExtract.es)
                                      ## get the ensembled median for each pixel
                                      EnsembleES <- app(mapsESstack, fun = median, na.rm = TRUE)
                                      ### normalise to the now-ensembled maps
                                      EnsembleESNorm <- normalise(EnsembleES)
                                      #### create another copy with no NANs
                                      EnESNan <- EnsembleESNorm
                                      EnESNan[is.nan(EnESNan)] <- 0
                                      
                                      ##### hotspots #####
                                      # get names
                                      valOutNm <- nmExtract %>% t() %>% as.data.frame() %>%
                                        rename_with(renameSeqNames)
                                      
                                      ## ------------ Notes --------------  ##
                                      ## The next section determines which pixels of both the ensembled BD and
                                      ## ES models overlap.
                                      ## ------------ ----- --------------  ##
                                      qps <- c(70, 80, 90, 95) # list of quartiles to use
                                      qpsOut <- lapply(qps, function(iq) {
                                        qp <- iq
                                        
                                        ## get the current quartile
                                        quartilesES <- quantile(values(EnsembleESNorm), probabilities, na.rm = T)
                                        quartilesES <- quartilesES[which(names(quartilesES) == paste0(qp, "%"))]
                                        ## get the current quartile
                                        quartilesBd <- quantile(values(EnsembleBdNorm), probabilities, na.rm = T)
                                        quartilesBd <- quartilesBd[which(names(quartilesBd) == paste0(qp, "%"))]
                                        
                                        ## adjust the maps based on whether 
                                        ## it was >= current percentile - i.e., get hotspots
                                        mapLogicalES <- as.numeric(EnESNan >= quartilesES)
                                        mapLogicalBD <- as.numeric(EnBDNan >= quartilesBd)
                                        ## get total hotspot pixels in both
                                        tHes <- sum(values(mapLogicalES))
                                        tHbd <- sum(values(mapLogicalBD))
                                        ### sum together
                                        mapLogical <- mapLogicalES + mapLogicalBD
                                        table(values(mapLogical))
                                        
                                        #### determine cross over
                                        mapTab <- table(as.vector(mapLogical))
                                        if(length(mapTab) > 2){
                                          map1hot1s.2 <- mapTab[[which(names(mapTab) == 2)]] # both hotspots overlap
                                        } else {
                                          map1hot1s.2 <- 0
                                        }
                                        map1hot1s.1 <- mapTab[[which(names(mapTab) == 1)]] # one has a hotspot
                                        map1hot1s.0 <- mapTab[[which(names(mapTab) == 0)]] # Neither has a hotspot
                                        map1hot1s <- map1hot1s.1 + map1hot1s.2 # above 0
                                        
                                        # run checks (total of one hotspot, with the overlapping aspect removed from both)
                                        stopifnot(map1hot1s.1 == (tHbd - map1hot1s.2) + (tHes - map1hot1s.2))
                                        
                                        # save as a df
                                        hotSpotsDfqp <- bind_cols(continent = iCont
                                                                  , res = xRes
                                                                  , perc = qp
                                                                  , pairs = mapsPairs
                                                                  , initial_row = rw
                                                                  , valOutNm
                                                                  , Overlap = map1hot1s.2
                                                                  , one_hotpot = map1hot1s.1
                                                                  , bd_hotspots = tHbd
                                                                  , es_hotspots = tHes 
                                                                  , no_hotspots = map1hot1s.0)
                                        
                                        ## ------------ Notes --------------  ##
                                        ## In order to save the final map to obtain spatial overlap info,
                                        ## save each map, one at a time
                                        ## ------------ ----- --------------  ##
                                        # convert the 2s to 1000
                                        mapLogical[mapLogical == 2] <- 1000
                                        ncell(mapLogical)
                                        which(values(mapLogical) == 1000)
                                        # get the coordinates of the overlapping hotspots
                                        coords <- xyFromCell(mapLogical, which(values(mapLogical) == 1000)) %>%
                                          as.data.frame()
                                        coords$x <- as.character(round(coords$x, 2))
                                        coords$y <- as.character(round(coords$y, 2))
                                        coords$xy <- paste(gsub("-", "m", coords$x)
                                                           , gsub("-", "m", coords$y)
                                                           , sep = "_")
                                        percName <- paste0("perc_", qp)
                                        cxqp <- bind_cols(xy = coords$xy
                                                          , continent = iCont
                                                          , res = xRes
                                                          , row = rw
                                                          , pairs = mapsPairs
                                                          , perc = qp) %>%
                                          rename(!!percName := perc)
                                        
                                        return(list(hotSpotsDfqp, cxqp))
                                      })
                                      
                                      ## combine all percentile results
                                      hotSpotsDf <- lapply(qpsOut, function(x) x[[1]]) %>%
                                        bind_rows()
                                      ### spatial overlap
                                      cxx <- lapply(qpsOut, function(x) x[[2]] %>% as.data.frame())
                                      cx2Combine <- Reduce(function(y1, y2) merge(y1, y2, all = TRUE),
                                                           cxx # list of dfs
                                      ) %>%
                                        as.data.frame()
                                      
                                      return(list(hotSpotsDf
                                                  , cx2Combine))
                                      
                                    } # for loop end
          
          toc()
        } # '!fExists' end
        cat("Start:", numStart)
        cat(" | End:", numEnd, "[ENDED]")
        cat("\n")
        stopCluster(cl)
        
        # if(numEnd == 1000) { stop("1000") }
        
        # # extract just the tables
        finalTable <- lapply(forLoopResults, function(x) x[[1]]) %>%
          bind_rows()
        finalMap <- lapply(forLoopResults, function(x) x[[2]]) %>%
          bind_rows() %>%
          # Select columns that contain "perc" in their names
          {
            perc_columns <- dplyr::select(., contains("perc")) %>% colnames()
            group_by_at(., vars(xy, res, continent, pairs, one_of(perc_columns)))
          } %>%
          summarise(count = n())
        head(finalMap)
        
        cat("Forloop finished at", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
        
        if(!fExists){
          # first search for the previous save
          if(numStart != 1){
            
            # load in old tables
            oldf <- fread(saveFilePrev) %>% as.data.frame() %>%
              mutate(res = as.character(res))
            f2 <- oldf %>%
              bind_rows(., finalTable)
            
            # load in old maps
            oldmap <- fread(saveMapPrev) %>% as.data.frame() %>%
              mutate(res = as.character(round(res, 2))
                     , xy = as.character(xy, 2))
            cat("   ... Binding maps\n")
            f3map <- bind_rows(oldmap, finalMap) %>%
              # Select columns that contain "perc" in their names
              {
                perc_columns <- dplyr::select(., contains("perc")) %>% colnames()
                group_by_at(., vars(xy, res, continent, pairs, one_of(perc_columns)))
              } %>%
              summarise(count = n())
            
            cat("... ... Saving df... ... ...\n")
            # save the final file of the df
            fwrite(f2, saveFile, row.names = F)
            file.remove(saveFilePrev)
            
            # save the final map of the df
            fwrite(f3map, saveMap, row.names = F)
            file.remove(saveMapPrev)
            
            cat("Per cent finished at", format(Sys.time(), "%Y-%m-%d %H:%M"), "was"
                , numEnd/maxRows * 100, "%\n")
            
          } else {
            
            cat("... ... Saving df... ... ...\n")
            # save the final file of the df
            fwrite(finalTable, saveFile, row.names = F)
            
            # save the final map 
            fwrite(finalMap, saveMap, row.names = F)
          } 
        } # end of '!fExists'
        
        cat("... ... ... It is finished!! ... ... ... ...\n\n")
        
        # tidy
        # stop()
        rm(forLoopResults)
        gc()
        Sys.sleep(3)
        
        # stopifnot(numStart < 11)
        # stop()
      }
      
    } # end of checking whether the final dataset already existed
  } # end of ES
} # end of continent