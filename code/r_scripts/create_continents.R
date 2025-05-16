## ---------------------------
##
## Script name: create_continents.R
##
## Purpose of script: create the continent shapefiles that will be used for the 
##                    final analysis.
##                    
## ------------ Notes --------------  ##
## The biggest change here is that Russia will be split into its Asian and
## European constituent parts
## ------------ ----- --------------  ##
##
## Run before: Cai_data_download.R
##             template_and weightings.R
##
## list of final outputs:
##    North America_final.gpkg
##    South America_final.gpkg
##    Africa_final.gpkg
##    Europe_final.gpkg
##    Asia_final.gpkg
##    Oceania_final.gpkg
##
## Author: Dr. Paul M. Evans
##
## Date Created: 2024-08-02
##
## Copyright (c) Paul M. Evans, 2024
## Email: paueva@ceh.ac.uk
## Git: https://github.com/pevans13

cat('\n\n####################################'
    , '\n\n create_continents.R'
    , '\n\n####################################\n\n')
# stop("start of create_continents.R")

#### 1 - task1 ####
# Download continent shapefiles
continents <- ne_countries(scale = "medium", type = "countries"
                           , continent = c("Africa", "Asia", "Europe"
                                           , "North America", "Oceania"
                                           , "South America")
                           , returnclass = "sf")

# get unique continents
cont.unique <- unique(continents$continent)
co.num <- 0; co.list <- list()

# save each continent shapefile
if(file.exists(file.path("data", "intermediate_outputs", "continents"
                         , paste0(cont.unique[length(cont.unique)], "_final.gpkg")))){
  cat("final contients exist ... \n")
} else {
  # get Russia separately
  poo_tin <- geodata::gadm(country="RUS", level=1, path=tempdir()) %>%
    st_as_sf()
  # st_write(poo_tin
  #          , file.path("data", "intermediate_outputs", "continents", "poo_tin.gpkg")
  #          , append = F)
  # plot(poo_tin[1])
  ## select only certain regions
  rusEurRegion <- poo_tin %>%
    filter(!GID_1 %in% c(
      "RUS.2_1", "RUS.3_1", "RUS.6_1", "RUS.9_1", "RUS.11_1"
      , "RUS.12_1", "RUS.16_1", "RUS.18_1", "RUS.24_1", "RUS.27_1"
      , "RUS.28_1", "RUS.29_1", "RUS.30_1", "RUS.35_1", "RUS.36_1"
      , "RUS.40_1", "RUS.50_1", "RUS.51_1", "RUS.53_1", "RUS.55_1"
      , "RUS.56_1", "RUS.60_1", "RUS.61_1", "RUS.66_1", "RUS.69_1"
      , "RUS.71_1", "RUS.73_1", "RUS.80_1", "RUS.82_1", "RUS.83_1"
    )) %>%
    st_union() %>% st_as_sf()
  rename_geometry <- function(g, name){
    current = attr(g, "sf_column")
    names(g)[names(g)==current] = name
    st_geometry(g)=name
    g
  }
  rusEurRegion <- rename_geometry(rusEurRegion, "geometry")
  
  plot(rusEurRegion[1])
  # Asia part
  rusAsiaRegion <- poo_tin %>%
    filter(GID_1 %in% c(
      "RUS.2_1", "RUS.3_1", "RUS.6_1", "RUS.9_1", "RUS.11_1"
      , "RUS.16_1", "RUS.18_1", "RUS.24_1", "RUS.27_1"
      , "RUS.28_1", "RUS.29_1", "RUS.30_1", "RUS.35_1", "RUS.36_1"
      , "RUS.40_1", "RUS.50_1", "RUS.51_1", "RUS.53_1", "RUS.55_1"
      , "RUS.56_1", "RUS.60_1", "RUS.61_1", "RUS.66_1", "RUS.69_1"
      , "RUS.71_1", "RUS.73_1", "RUS.80_1", "RUS.82_1", "RUS.83_1"
    ))  %>%
    st_make_valid() %>%  # Fix invalid geometries
    st_union() %>% st_as_sf()
  rusAsiaRegion <- rename_geometry(rusAsiaRegion, "geometry")
  plot(rusAsiaRegion[1])
  
  for(co in cont.unique){
    cat("This continent:", co, " \n")
    
    ## ------------ Notes --------------  ##
    ## use continent outline to cut the number of pixels to extent
    
    ## adjust Europe and Asia so that the majority of Russia is part of Asia
    ## ------------ ----- --------------  ##
    # create shapefile from continent
    continents.regions <- continents %>%
      filter(continent == co) %>%
      # Split multipolygons into individual polygons
      st_cast(., "POLYGON") %>%
      # get sizes of individual units
      mutate(area = st_area(.)) %>%
      mutate(area = units::set_units(area, "km2")) %>%
      filter(area > units::set_units(7774 * 2, "km2")) %>%
      dplyr::select(continent, name_sort)
    plot(st_geometry(continents.regions))
    
    ## from Europe, remove 9 and 35
    if(co == "Europe"){
      continents.regions <- continents.regions %>%
        # remove extra
        slice(-c(9, 35)) %>%
        # remove all 'Russian Federation'
        filter(name_sort != "Russian Federation") %>%
        # create a union for the european part of russia
        bind_rows(., rusEurRegion)
      plot(continents.regions[1])
    }
    
    ## to Asia, add Russian part
    if(co == "Asia"){
      continents.regions <- continents.regions %>%
        # create a union for the european part of russia
        bind_rows(., rusAsiaRegion)
      plot(continents.regions[1])
    }
    
    # make one outline
    continents.one <- continents.regions %>%
      st_union()
    # check geometry
    valid_continents.one2 <- st_is_valid(continents.one)
    
    ggplot() +
      geom_sf(data = continents.regions, fill = "blue") +
      geom_sf(data = continents.one, fill = "red", alpha = 0.5)
    
    # save the continent as its shapefile
    st_write(continents.one
             , file.path("data", "intermediate_outputs", "continents"
                         , paste0(co, "_final.gpkg"))
             , append = F)
    beep(12)
  }}

#### x - write readme ####
readmePath <- file.path("data", "intermediate_outputs", "continents", "readme.md")
# Create and open the file connection
fileConn <- file(readmePath)
writeLines(c(
"Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13"
, "Date created: 2024-08-02"
, paste0("Last update: ",  format(Sys.Date()))
, "Produced by 'create_continents.R'

Description of files in directory:
The files in this directory contain the final shapefiles for the continents that will be used in the analysis for the Uncertainty project.
The initial continental shapefiles were downloaded from the rnaturalearth in R using the 'ne_countries' function. Russia was considered to be totally in Europe from that dataset. However, for our analysis we split it between its Asian and European constituent parts.

Files:
    North America_final.gpkg
    South America_final.gpkg
    Africa_final.gpkg
    Europe_final.gpkg
    Asia_final.gpkg
    Oceania_final.gpkg")
           , fileConn)
close(fileConn)
