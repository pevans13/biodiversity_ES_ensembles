## ---------------------------
##
## Script name: data_to_continent_extent.R
##
## Purpose of script: to convert all data involved in this project to 100 km2.
## The data of interest are vascular plant and ecosystem service data
##                    
## ------------ Notes --------------  ##
## The native resolution of the ES data was 1 km2, whereas it was 7,774 km2
## for the vascular plant species richness data.
## ------------ ----- --------------  ##
##
## Run after: template_and weightings.R
##
## Run before: any analysis (e.g., hotspot_analysis.R)
##             correlation_test.R
##
## Specific numbered tasks:
## 1 - convert all original ES data to 100 km2 [continental]
## 2 - convert all original biodiversity data to 100 km2 [continental]
##
## list of final outputs:
##      data/intermediate_outputs/biodiversity/normalised
##      data/intermediate_outputs/ecosystem_services/normalised/res100/continents/
##    
## Author: Dr. Paul M. Evans
##
## Date Created: 2024-06-13
##
## Copyright (c) Paul M. Evans, 2024
## Email: paueva@ceh.ac.uk
## Git: https://github.com/pevans13

cat("\n\n####################################"
    , "\n\nstarting data_to_continent_extent.R"
    , "\n\n####################################\n\n")
# stop("data_to_continent_extent.R")

#### 0 - paths ####
path.vpData <- "data/raw/vascular_plants/Cai_2023"
path.esData <- "D:/"
# where to save raster biodiversity data
dir.create(path.bdSave <- file.path("data", "intermediate_outputs", "biodiversity", "normalised")
           , showWarnings = F, recursive = T)
# where to save raster ES data
path.esSave <- file.path("data", "intermediate_outputs", "ecosystem_services")

#### 0 - load functions ####
# Function to convert metres to degrees of latitude
metres_to_degrees_latitude <- function(metres) {
  degrees_latitude = metres / 111000
  return(degrees_latitude)
}

# Function to convert metres to degrees of longitude given a latitude
metres_to_degrees_longitude <- function(metres, latitude) {
  # Convert latitude to radians
  latitude_in_radians = latitude * (pi / 180)
  degrees_longitude = metres / (111320 * cos(latitude_in_radians))
  return(degrees_longitude)
}

# Function to check geometries are ok, and then adjust if not
checkValidAdjust <- function(x, conty = co){
  # x2 <- st_make_valid(x)
  x2 <- x
  vy <- st_is_valid(x2)
  vyReason <- st_is_valid(x2, reason = TRUE)
  unique(vyReason)
  
  invalidGeomXY <- ifelse(sum(vy) != nrow(x2), T, F)
  
  # remove ones that are not
  if(invalidGeomXY){
    
    ## determine which rows are invalid
    rwRemove <- which(vy == F)
    cat(paste0("rows removed for ", conty, ":")
        , rwRemove, sep = "\n")
    ## remove invalid rows
    x3 <- x2 %>%
      slice(-c(rwRemove))
    sum(st_is_valid(x3))
    ## if the numbers now match, use those data
    if(sum(st_is_valid(x3)) == nrow(x3)){
      x4 <- x3
      cat("   ... invalid rows fixed for x [input shapefile] ... \n")
    } else {
      stop("Still invalid")
    }
  } else {
    cat("   ... No invalid rows for x [input shapefile] ... \n")
    x4 <- x2
  }
  return(x4)
}

# Define a function to normalise raster values
normalise <- function(x) {
  (x - min(values(x), na.rm = TRUE)) / (max(values(x), na.rm = TRUE) - min(values(x), na.rm = TRUE))
}

# winsorise the data
winsorise <- function(x){
  
  # 'x' = a raster
  xr <- x
  ## ------------ Notes --------------  ##
  ## following Willcock et al. (2023)
  ## the winsor technique should be a 'double-sided Winsorising protocol'
  ## using '2.5 and 97.5% percentiles'
  ## ------------ ----- --------------  ##
  
  ## get the appropriate percentiles
  probabilities <- c(0.025, 0.975)
  quartiles <- quantile(values(xr), probabilities, na.rm = T)
  minq <- quartiles[[1]]; maxq <- quartiles[[2]]
  # cat("......", "2.5% =", minq 
  #     , "| 97.5% =", maxq
  #     , "\n")
  values(xr) <- Winsorize(values(xr)
                          , minval = minq
                          , maxval = maxq
                          , na.rm = T)
  return(xr)
}

# function to load and crop raster
readCropCheckMask <- function(xRast, tempRast, meth){
  
  #' @param xRast A string leading to a raster.
  #' @param tempRast A rast object that provides a template for cropping and masking.
  
  cat("Reading in", xRast, "\n")
  
  # load in rast
  origrast <- rast(xRast) %>%
    ## crop by continent
    crop(., tempRast)
  # plot(origrast)
  
  ## Convert sf object to terra raster - resample to get 7,774 km2 resolution
  origrast.resamp <- terra::resample(x = origrast,
                                     , y = tempRast
                                     , method = meth)
  
  while(max(values(origrast.resamp), na.rm = T) == -Inf){
    
    cat("Still infinite |")
    
    # an error is produced if the values are too big, so they can be reduced
    # and because normalised values are being used, nothing should be affected
    origrast <- origrast / 2
    ## Convert sf object to terra raster - resample to get 7,774 km2 resolution
    origrast.resamp <- terra::resample(x = origrast,
                                       , y = tempRast
                                       , method = "average")
  }
  cat("\n")
  
  ### mask the data
  origrast.masked <- origrast.resamp %>%
    mask(., tempRast)
  # beep(1)
  return(origrast.masked)
}

# Function  to get the weighted average means for each grid cell underlying each hexagon
wmHexFunc <- function(reso = "", conty = co){
  
  # get names to save with
  saveName <- file.path(path.bdSave, paste0("res", reso), "continents"
                        , paste0(bd.name, "_norm", "_", co, ".tif"))
  ## non-normalised name
  saveNameNN <- file.path(dirname(path.bdSave)
                          , paste0(bd.name, "_", reso, "_", co, ".tif"))
  
  if(file.exists(saveName)){
    
    cat("Raster of", basename(saveName), "already saved", "...\n")
    
  } else {
    
    # create hexagons as a table
    hexTable <- predictions_grid %>%
      st_drop_geometry() %>%
      # keep only id and value columns
      dplyr::select(grid_ID, point_x, point_y
                    , value)
    
    head(hexTable)  
    
    if(reso == "7774"){
      ## merge with the proportional hex values
      mergeTableHex <- contTableHex7774 %>%
        merge(., hexTable)
      head(mergeTableHex)
    } else {
      ## merge with the proportional hex values
      mergeTableHex <- contTableHex100 %>%
        merge(., hexTable)
      head(mergeTableHex)
    }
    
    ## ------------ Notes --------------  ##
    ## Now, do the calculation, where propotion is multiplied
    ## by the 'value' (i.e., species richness)
    ## ------------ ----- --------------  ##
    calcTableHex <- mergeTableHex %>%
      mutate(weighted_value = value * proportion) %>%
      # summarise per pixel
      group_by(smlGrd_ID) %>%
      summarise(wmValue = sum(weighted_value)
                , count = n())
    
    if(reso == "7774"){
      ## merge back to spatial, using the small grid IDs
      tableHexSmlGrid <- calcTableHex %>%
        merge(smlGrid7774, .)
    } else {
      ## merge back to spatial, using the small grid IDs
      tableHexSmlGrid <- calcTableHex %>%
        merge(smlGrid100, .)
    }
    
    ### finally, save as a raster - load in template
    rTemp <- rast(file.path(pathDataInter, "continents", "rasters"
                            , paste0(co, reso, "template.tif")))
    ### convert sf object to terra raster based on save raster template
    finalVal.rast <- terra::rasterize(vect(tableHexSmlGrid)
                                      , rTemp
                                      , field = "wmValue")
    
    # save, to check
    # st_write(tableHexSmlGrid, paste0("finalVal_res", reso, ".gpkg"), append = F)
    # save the final values (unnomralised)
    writeRaster(finalVal.rast, saveNameNN, overwrite = T)
    
    ## ------------ Notes --------------  ##
    ## the next section windorises and normalises the data based on
    ## continental minimums and maximums
    ## ------------ ----- --------------  ##
    
    # winsorise based on 2.5 and 97.5 percentiles
    finalVal.w <- winsorise(finalVal.rast)
    finalVal.w
    
    # normalise based on winsorised data
    finalVal.norm <- normalise(finalVal.w)
    finalVal.norm
    plot(finalVal.norm)
    
    stopifnot(res(finalVal.norm)[1] < 0.85) # check res is ok
    if(max(values(finalVal.norm), na.rm = T) != 1){stop("Value NA for norm")}
    
    ## save data
    cat("saving raster of", basename(saveName), "...\n")
    writeRaster(finalVal.norm
                , saveName
                , overwrite = T)
  }
}

# function to get mean for ES data
esDataRastFunc <- function(est = es.type, esn = es.name){
  
  ## save data - create name
  saveName <- file.path(path.esSave, "normalised", "res100", "continents", est
                        , paste0(esn, "_norm", "_", co, ".tif"))
  ### create that path
  dir.create(file.path(path.esSave, "normalised", "res100", "continents", est)
             , showWarnings = F, recursive = T)
  ## non-normalised name
  saveNameNN <- file.path(dirname(path.esSave)
                          , paste0(esn, "_", "100", "_", co, ".tif"))
  
  fexists <- file.exists(saveName)
  
  # does it exist?
  if(fexists){
    mxmx <- max(values(rast(saveName)), na.rm = T)
    if(is.infinite(mxmx)){
      naAdjust <- T
    } else {
      naAdjust <- F
    }
  } else {
    naAdjust <- F
  }
  
  ## add an extra stipulation to only run the next steps for JRC_Resprojected when
  ## it is Europe, as that result only covers Europe
  if(es.name == "JRC_Resprojected" & co != "Europe"){
    fexists <- T
    naAdjust <- F
  }
  
  # if it does not exist, or if it does exist but needs adjusting
  if(fexists & naAdjust || !fexists){
    
    ## ------------ Notes --------------  ##
    ## The next section resamples the rasters for ES. If the ES is largest associated with biomass,
    ## it should use a mean average. These are: Carbon, Grazing, and RawMaterials_Fuelwood.
    ## For the others (i.e., Recreation and Water), sum should be used.
    ## ------------ ----- --------------  ##
    
    if(es.type %in% c("Carbon", "Grazing", "RawMaterials_Fuelwood")){
      # load, crop, check, and then mask, then resample (using average for this group)
      origESrast.resamp <- readCropCheckMask(ies, rTemp, meth = "average")
      plot(origESrast.resamp, main = paste0(basename(ies), " | ", "res7774"))
    } else if(es.type %in% c("Recreation", "Water")){
      # load, crop, check, and then mask, then resample (using sum for this group)
      origESrast.resamp <- readCropCheckMask(ies, rTemp, meth = "sum")
      plot(origESrast.resamp, main = paste0(basename(ies), " | ", "res7774"))
    } else {
      cat("Wrong es.type\n")
      stop("es.type error")
    }
    # save the final values (unnomralised)
    writeRaster(origESrast.resamp, saveNameNN, overwrite = T)
    
    ## ------------ Notes --------------  ##
    ## the next section windorises and normalises the data based on
    ## continental minimums and maximums
    ## ------------ ----- --------------  ##
    
    # winsorise based on 2.5 and 97.5 percentiles
    origESrast.w <- winsorise(origESrast.resamp)
    plot(origESrast.w)
    
    # normalise based on winsorised data
    origESrast.norm <- normalise(origESrast.w)
    plot(origESrast.norm)
    osMax <- max(values(origESrast.norm), na.rm = T)
    stopifnot(round(osMax,2) == 1 | round(osMax,2) == -Inf)
    
    ## save data 
    
    cat("...... Saving here:", saveName, "\n")
    writeRaster(origESrast.norm
                , saveName
                , overwrite = T)
    
    ## ES - see how many currrently exist
    xl <- list.files(file.path(path.esSave, "normalised")
                     , recursive = T)
    xl <- xl[!grepl("global", xl)]
    gxl <- gsub(".*norm_(.+).tif$", "\\1", xl) %>% as.vector() %>%
      table() %>% as.data.frame() %>% rowwise() %>% mutate(count = paste(.[1], Freq)); cat(gxl$count, sep = "\n")
    
  } else {
    
    cat("...... ", basename(saveName), "at"
        , basename(dirname(dirname(dirname(saveName))))
        , "already exists ...\n")
  }
  
}

#### list biodversity and ES files ####
# list vascular plants richness files
vpList <- list.files(path.vpData
                     , pattern = ".*sr_(.+).RData$"
                     , full.names = T) %>%
  # do not want the ensemble
  .[!grepl("ensemble", .)]
# vpList
stopifnot(length(vpList) == 11)

# list ecosystem service files (ones that end in '.tif' and are not normalised)
esList <- list.files(file.path(path.esData)
                     , pattern = ".tif$"
                     , full.names = T
                     , recursive = T)
# esList
## remove any normalised ES data
esList <- esList[!grepl("normalised|biodiversity|recycle.bin", tolower(esList))] %>%
  .[!grepl("stage2", .)]
# removes mask used
esList <- esList[!grepl("masks used", tolower(esList))]
# remove grazing and raw material (not the focus of this study)
esList <- esList[!grepl("RawMaterials_Fuelwood|Grazing", tolower(esList))]
# esList
stopifnot(length(esList) == 24)

# convert metres to degrees for 7,774 km2
resDegreesLat <- metres_to_degrees_latitude(88000)
resDegreesLong <- metres_to_degrees_longitude(88000, resDegreesLat)

#### 1 - conversions to 100 km2 - at the continental level #### 
if(continentalAnalysis){
  
  for(iii in 1){
    # stop("inside iii")
    co.num <- 0; co.list <- list()
    
    # determine the amount that currently exists of all
    ## ES
    xl <- list.files(file.path(path.esSave, "normalised")
                     , recursive = T); xl <- xl[!grepl("global", xl)]
    gxl <- gsub(".*norm_(.+).tif$", "\\1", xl) %>% as.vector() %>%
      table(); gxl
    gxl <- basename(dirname(xl)) %>% as.vector() %>%
      table(); print(gxl)
    ## BD
    xlb <- list.files(file.path(path.bdSave), recursive = T)
    gxlb <- gsub("^(.+?)_.*", "\\1", basename(xlb)) %>% as.vector() %>%
      table()
    print(gxlb)
    
    ## get the continents that were saved in 'create_continents.R'
    conts <- list.files(file.path(pathDataInter, "continents")
                        , pattern = "_final.gpkg")
    cont.unique <- gsub("_final.gpkg", "", conts)
    stopifnot(length(cont.unique) == 6)
    
    #### BD extraction ####
    for(co in cont.unique){
      # if(co == cont.unique[1]){stop("1st co")}
      
      cat("This continent:", co, " \n")
      continents.one <- st_read(file.path("data", "intermediate_outputs", "continents"
                                          , paste0(co, "_final.gpkg")), quiet = T)
      plot(continents.one[1])
      
      ## ------------ Notes --------------  ##
      ## For each continent, load the table that was created in 'template_rasters.R'
      ## This contains each of the hexagons of a region, and its weighting
      ## giving the information on how to calculated a weighted mean per 
      ## overlapping grid cells and hexagons
      ## ------------ ----- --------------  ##
      
      # load in the continent saved table - for the hexagon weightings
      contTableHex100 <- fread(paste0("data/"
                                      , paste0(co, "_biodiversity_100_weighting.csv")))
      head(contTableHex100)
      
      # load in the continent saved table - for the small grid IDs (7774 km)
      smlGrid7774 <- st_read(file.path(pathDataInter, "biodiversity"
                                       , paste0(substring(co, 1, 4)
                                                , "_BD_7774_smallGrid.gpkg")), quiet = T)
      # load in the continent saved table - for the small grid IDs (7774 km)
      smlGrid100 <- st_read(file.path(pathDataInter, "biodiversity"
                                      , paste0(substring(co, 1, 4)
                                               , "_BD_100_smallGrid.gpkg")), quiet = T)
      
      ###### BD extraction ######
      bd.num <- 0; bd.list <- list()
      ## use one of Cai rasters
      for(iv in vpList){
        
        # load in the hexagons, and create as a table
        load(iv)
        
        # get save details
        bd.num <- bd.num + 1 # number
        bd.name <- gsub(".*23/(.+)_Pred.*", "\\1", iv)

        ##### BD extraction (100 x 100 km) #####
        wmHexFunc(reso = "100", conty = co)
        
      } # end of 'vplist'
    }
    
    # does it need to be run as a parallelised process?
    parrellelProcess <- F
    
    #### ES extraction ####
    if(!parrellelProcess){
      for(co in cont.unique){
        # if(co == cont.unique[1]){stop("1st co")}
        
        cat("This continent:", co, " \n")
        continents.one <- st_read(file.path("data", "intermediate_outputs", "continents"
                                            , paste0(co, "_final.gpkg")), quiet = T) 
        
        # load in appropriate raster template
        rTemp <- rast(file.path(pathDataInter, "continents", "rasters"
                                , paste0(co, "100", "template.tif")))
        plot(rTemp)
        
        es.num <- 0; es.list <- list()
        ## use one of Cai rasters
        for(ies in esList){
          # if(ies == esList[1]){ stop("1st es") }
          # if(co == "Europe"){ stop("WD") }
          es.num <- es.num + 1 # number
          es.type <- gsub(".*D:/([^/]+) Model Layers/.*", "\\1", ies)
          es.name <- gsub(".*/(.+).tif*", "\\1", ies)
          cat("...", es.num, "|", es.type, "|", es.name, "| es at 100 km2 for", co, "\n")
          
          # run the function to get the final normalised raster
          esDataRastFunc(est = es.type, esn = es.name)
          
        } # eslist end
      } # continent end
      
    } else {
      
      ## use a cluster to do the next bit
      # Set the number of cores
      nCores <- detectCores() - 4
      
      cl <- makeCluster(nCores)
      registerDoParallel(cl)
      # stopCluster(cl)
      on.exit(stopCluster(cl)) # Ensure cluster is stopped if an error occurs
      
      beep(12)
      # stop("before cluster")
      cat("Starting foreach loop ...on", nCores, "cores\n")
      # results <- foreach(co = cont.unique[2:6]
      results <- foreach(co = cont.unique
                         , .packages = c("sf", "dplyr", "ggplot2", "pbapply", "terra", "beepr", "DescTools")
                         , .errorhandling = "pass") %dopar% {
                           tryCatch({
                             
                             cat("This continent:", co, " \n")
                             continents.one <- st_read(file.path("data", "intermediate_outputs", "continents"
                                                                 , paste0(co, "_final.gpkg")), quiet = T) 
                             
                             # load in appropriate raster template
                             rTemp <- rast(file.path(pathDataInter, "continents", "rasters"
                                                     , paste0(co, "100", "template.tif")))
                             plot(rTemp)
                             
                             es.num <- 0; es.list <- list()
                             ## use one of Cai rasters
                             for(ies in esList){
                               # if(ies == esList[1]){ stop("1st es") }
                               # if(co == "Europe"){ stop("WD") }
                               es.num <- es.num + 1 # number
                               es.type <- gsub(".*D:/([^/]+) Model Layers/.*", "\\1", ies)
                               es.name <- gsub(".*/(.+).tif*", "\\1", ies)
                               cat("...", es.num, "|", es.type, "|", es.name, "| es at 100 km2 for", co, "\n")
                               
                               # run the function to get the final normalised raster
                               esDataRastFunc(est = es.type, esn = es.name)
                               
                             } # eslist end
                             
                             return(list(success = TRUE))
                           }, error = function(e) {
                             return(list(success = FALSE, error = e))
                           })
                         } # foreach end ('co')
      stopCluster(cl)
    }
  } # 'iii' end
} # end of continentalAnalysis
