## ---------------------------
##
## Script name: template_and weightings.R
##
## Purpose of script: create raster templates for the biodiversity 
##                    data, and determine the proportion
##                    overlap of hexagonal polygons and pixels in
##                    the raster templates
##                    
## ------------ Notes --------------  ##
## The biodiversity data initial come in hexagonal polygon form
## ------------ ----- --------------  ##
##
## Run after: Cai_data_download.R
##
## Run before: data_to_continent_extent.R
##
## list of final outputs:
##    data/[continent]_biodiversity_[resolution]_weighting.csv
##    data/intermediate_outputs/continents/[continent][resolution]template.tif
##
## Author: Dr. Paul M. Evans
##
## Date Created: 2024-08-02
##
## Copyright (c) Paul M. Evans, 2024
## Email: paueva@ceh.ac.uk
## Git: https://github.com/pevans13

cat('\n\n####################################'
    , '\n\n starting template_raster.R'
    , '\n\n####################################\n\n')
# stop("template_raster.R")

#### 0 - paths ####
path.vpData <- "data/raw/vascular_plants/Cai_2023"

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
checkValidAdjust <- function(x){
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
    cat(paste0("rows removed for ", ":")
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

# Function to eventually calculate which parts overlap (i.e., intersect)
smallGridFunc <- function(tRast
                          , resLong, resLat # longitude (x), lat (y) resolutions, in degrees
) {
  # create a grid from the raster template  - the point of this is to be able to get a weighted mean
  # using the polygons area and raster pixel overlap.
  rGrid <- st_make_grid(tRast
                        , cellsize = c(resLong # xres / longitude
                                       , resLat) # yres / latitude
  ) %>% st_as_sf() %>%
    mutate(smlGrd_ID = paste0(substring(co,1,4)
                              , "_smGr_", 1:nrow(.)))
  return(rGrid)
}

largeGridFunc <- function(inGrid
                          , nChoice = 10 # choice of how to slice the larger grid
){  
  # split grid for less memory use
  ## choose how to split the grid, both vertical and horizontally
  n <- nChoice
  gridLarger <- st_make_grid(inGrid
                             # how many splits
                             , n = c(n, n)) # this splits the grid
  # Convert the grid to an sf object
  largeGridFunc <- st_as_sf(gridLarger) %>%
    mutate(lrgGrd_ID = paste0(substring(co,1,4)
                              , "_lrgGr_", 1:nrow(.)))
  plot(largeGridFunc)
  return(largeGridFunc)
}

doesIntersectFunc <- function(){
  ## first, find which larger intersect at all
  for(di in 1:nrow(largeGrid)) {
    dfx <- largeGrid[di, ]
    df <- dfx[st_intersects(dfx, hexesBDCo) %>% lengths > 0, ]
    if(nrow(df) > 0){
      xCC <- bind_cols(x = di, intersect = "Yes")
    } else {
      xCC <- bind_cols(x = di, intersect = "No")
    }
    if(di == 1){
      doesIntersect <- xCC
    } else {
      doesIntersect <- doesIntersect %>%
        bind_rows(xCC)
    }
  }
  
  ### get a list of those that do intersect
  doesIntersect <- doesIntersect %>% as.data.frame() %>%
    filter(intersect == "Yes") %>%
    dplyr::select(x) %>% unlist() %>% as.numeric()
  return(doesIntersect)
}

# function for the foreach loop to determine intersections
# Set up parallel backend to use multiple processors
foreachIntersectFunc <- function(inGridSmall = rGrid
                                 , inGridLarge = largeGrid
                                 , hexIn = hexesBDCo){
  numCores <- detectCores() - 3
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl)) # Ensure cluster is stopped if an error occurs
  cat("\nStarting a parallel backend on", numCores, "cores at"
      , format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  
  # Perform intersection on each chunk in parallel
  resultsList <- foreach(iGrs = doesIntersect
                         , .packages = c('sf', 'dplyr', 'ggplot2')
                         , .export=c("checkValidAdjust")) %dopar% {
                           
                           # Print iteration
                           cat("\nStarting ")
                           cat(iGrs, "|", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
                           
                           # Extract just one grid square from the larger grid
                           chunk <- inGridLarge[iGrs, ]
                           
                           # Subset the part of inGridSmall that intersects with the current chunk
                           rLargerGrid_chunk <- st_intersection(inGridSmall, chunk) %>%
                             st_make_valid()
                           rLargerGrid_chunk <- checkValidAdjust(rLargerGrid_chunk)
                           head(rLargerGrid_chunk)
                           
                           # ggplot() +
                           #   geom_sf(data = hexesBDCo, fill = "red") +
                           #   geom_sf(data = chunk, fill = "blue") +
                           #   geom_sf(data = rLargerGrid_chunk, colour = "black", alpha = 0.7)
                           
                           # Perform the intersection
                           intersectionResult <- st_intersection(rLargerGrid_chunk, hexIn)
                           intersectionResult <- checkValidAdjust(intersectionResult)
                           intersectionResult <- intersectionResult %>%
                             # Get area
                             mutate(pixArea = st_area(.))
                           
                           return(intersectionResult)
                         } # 'foreach' end
  
  cat("... Hurra! Parallel backend closed at"
      , format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  
  return(resultsList)
}

#### starting global maps ####
## ------------ Notes --------------  ##
## In this section, you need to load in a global map of both BD and ES
## ------------ ----- --------------  ##

##### Biodiversity #####
# list vascular plants richness files - just the first
vpList <- list.files(path.vpData
                     , pattern = ".RData$"
                     , full.names = T)
vpListEnsem <- vpList[1]
## load ensemble model
load(vpListEnsem)
### convert to WGS84 to match ES models
hexesBD <- predictions_grid %>%
  st_transform(4326)
### check geometries
hexesBD <- checkValidAdjust(hexesBD)
rm(predictions_grid)

# convert metres to degrees for 7,774 km2
resDegreesLat <- metres_to_degrees_latitude(88000)
resDegreesLong <- metres_to_degrees_longitude(88000, resDegreesLat)

##### Ecosystem service #####
# load in a template that is the same extent and resolution as in Willcock et al. (2023)
# this ensures that they can keep the same results when coarsened
resres <- c(0.008333, 0.008333) # set the res which is 1 km (in EPSG 4326)
### get rounded extent
xmin.new <- -180
xmax.new <- 180
ymin.new <- -65
ymax.new <- 90
## create as a new raster
templateRaster.res100 <- rast(ext=c(xmin.new, xmax.new
                                    , ymin.new, ymax.new)
                              , res = resres # copy resolution
                              , vals = 1) %>%
  ## coarsen the data by a factor of 100, making it 100 x 100 km2
  aggregate(., fact=100, fun=mean) 
templateRaster.res100

##### Continent shapefiles #####
# Download continent shapefiles
continents <- list.files(file.path(pathDataInter
                                   , "continents")
                         , pattern = "_final.gpkg"
                         , full.names = T, recursive = T)

# get unique continents
cont.unique <- unique(gsub("_final.gpkg", "", basename(continents)))
cont.unique

## ------------ Notes --------------  ##
## for all continents, crop the extent of a sample BD dataset (they are all the same format output).
## For this, the 'ensemble' BD model can be used.
## Then, convert to 7,774 km2 raster
## Finally, get the overlap between all pixels and polygons, to inform the future analysis.
## ------------ ----- --------------  ##

co.num <- 0; co.list <- list()
for(co in cont.unique){
  co.num <- co.num + 1
  cat("This continent is:", co, " \n")
  # read in continent
  continents.one <- st_read(continents[[co.num]], quiet = T)
  # stop()
  
  # crop the ensemble model to the continent
  # cut the biodiversity output by continent
  hexesBDCo <- hexesBD %>%
    st_intersection(., continents.one)
  hexesBDCo <- checkValidAdjust(hexesBDCo)
  plot(hexesBDCo %>% st_geometry())
  
  ## use the extent to make a new raster
  r <- rast(resolution = c(resDegreesLong # xres / longitude
                           , resDegreesLat) # yres / latitude
            , extent = ext(hexesBDCo))
  # set the values, to see the differentiation between pixels
  set.seed(42); values(r) <- sample(1:10, ncell(r), replace = T)
  plot(r)
  res(r)
  # save, to check
  # writeRaster(r, file.path("data", "checks", "rast.tif"), overwrite = T)
  
  # crop to the shapefile
  rCrop <- crop(r, hexesBDCo) %>%
    mask(., hexesBDCo)
  plot(rCrop)
  ## save for future use
  writeRaster(rCrop
              , file.path(pathDataInter, "continents", "rasters"
                          , paste0(co, "7774template.tif"))
              , overwrite = T)
  
  #### calculate crossover between pixels and BD hexagons ####
  ## ------------ Notes --------------  ##
  ## To get a value of the mean per pixel (i.e., using the 
  ## all polygons it overlaps by area), this has to be done
  ## in several steps, as below.
  ## First, get the intersections between the square grid and hexagons
  ## get a weighted mean value based on the area of crossover
  ## then, create final raster
  ## ------------ ----- --------------  ##
  
  ##### 100 grid #####
  
  cat("\n\n")
  cat("----------- xxx -----------")
  cat("\nStarting the 100 * 100 km section for", co, "\n")
  cat("----------- ooo -----------")
  
  ## ------------ Notes --------------  ##
  ## Repeat the above actions, but for 100 * 100 km pixels of ES
  ## ------------ ----- --------------  ##
  
  resultsSaveName <- paste0("data/"
                            , co
                            , "_biodiversity_100_weighting.csv")
  
  if(file.exists(resultsSaveName)){
    cat("\n", resultsSaveName, "already exists ...\n")
  } else {
    
    # crop the ES 100 km raster to the continent
    rastESCo <- templateRaster.res100 %>%
      crop(., continents.one) %>%
      mask(., continents.one)
    plot(rastESCo)
    # save, to check
    # writeRaster(rastESCo, "data/checks/ES100.tif", overwrite = T)
    
    ## save for future use
    writeRaster(rastESCo, file.path(pathDataInter, "continents", "rasters"
                                    , paste0(co, "100template.tif"))
                , overwrite = T)
    ### get the resolution, in degrees, of the raster
    res100m <- res(rastESCo)
    
    # label the polygons, for ID later
    hexesBDCo <- hexesBDCo %>%
      mutate(hex_id_7774 = paste0(co, "_", 1:n())) %>%
      # keep only certain cols
      dplyr::select(hex_id_7774, grid_ID, point_x, point_y)
    head(hexesBDCo)
    hexesBDCo <- checkValidAdjust(hexesBDCo)
    # determine which parts of the grids overlap
    rGrid <- smallGridFunc(tRast = rastESCo
                           , resLong = res100m, resLat = res100m)
    head(rGrid)
    ## save the small grid, once per continent
    coSaveName <- file.path(pathDataInter, "biodiversity"
                            , paste0(substring(co, 1, 4)
                                     , "_BD_100_smallGrid.gpkg"))
    if(!file.exists(coSaveName)){
      st_write(rGrid, coSaveName, append = F)
    }
    
    largeGrid <- largeGridFunc(rGrid)
    head(largeGrid)
    doesIntersect <- doesIntersectFunc()
    
    ## ------------ Notes --------------  ##
    ## To get a value of the mean per pixel (i.e., using the 
    ## all polygons it overlaps by area), this has to be done
    ## in several steps, as below.
    ## First, get the intersections between the square grid and hexagons
    ## get a weighted mean value based on the area of crossover
    ## then, create final raster
    ## ------------ ----- --------------  ##
    
    # run parallel function that determines intersections, and then saves the output
    resultsListOut <- foreachIntersectFunc(inGridSmall = rGrid)
    head(resultsListOut)
    
    # bind the output
    resultsListBind <- bind_rows(resultsListOut) %>%
      st_drop_geometry()
    head(resultsListBind)
    ## calculate the proportions
    rlProport <- resultsListBind %>%
      group_by(smlGrd_ID) %>%
      summarise(totArea = sum(pixArea)) %>%
      ## merge back
      merge(resultsListBind, .) %>%
      mutate(proportion = pixArea/totArea)
    head(rlProport)
    ### save the information
    fwrite(rlProport
           , resultsSaveName
           , row.names = F)
  }
  # beep(8)
  Sys.sleep(1)
  
} # 'co' (continent) end
