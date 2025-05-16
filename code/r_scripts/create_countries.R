## ---------------------------
##
## Script name: create_countries.R
##
## Purpose of script: create the country shapefiles that will be used for the 
##                    final analysis.
##
## list of final outputs:
##    data/intermediate_outputs/[country]_final.gpkg
##
## Author: Dr. Paul M. Evans
##
## Date Created: 2024-10-02
##
## Copyright (c) Paul M. Evans, 2024
## Email: paueva@ceh.ac.uk
## Git: https://github.com/pevans13
##
## Run before: Cai_data_download.R
##             template_and weightings.R

cat('\n\n####################################'
    , '\n\n create_countries.R'
    , '\n\n####################################\n\n')
# stop("start of create_countries.R")

#### 1 - task1 ####
# Download country shapefiles
countries <- ne_countries(scale = "medium", type = "countries"
                          , continent = c("Africa", "Asia", "Europe"
                                          , "North America", "Oceania"
                                          , "South America")
                          , returnclass = "sf") %>%
  # get the size of each separate shapefile
  mutate(land_area = st_area(.)) %>%
  relocate(formal_en, land_area) %>%
  # Only include countries bigger than 20,000 square kilometres (2 x resolution of this study)
  filter(land_area > units::set_units(10000 * 2, km^2))

# get unique countries
cont.unique <- unique(countries$brk_name)
co.num <- 0; co.list <- list()

for(co in cont.unique){
  cat("This country:", co, " \n")
  
  # save each country shapefile
  if(file.exists(file.path("data", "intermediate_outputs", "countries"
                           , paste0(co, "_final.gpkg")))){
    cat("already exists ... \n")
    
  } else {
    
    # skip Solomon Islands
    if(co != "Solomon Is."){
      
      ## ------------ Notes --------------  ##
      ## use country outline to cut the number of pixels to extent
      ## ------------ ----- --------------  ##
      
      # create shapefile from country
      countrys.regions <- countries %>%
        filter(brk_name == co) %>%
        # Split multipolygons into individual polygons
        st_cast(., "POLYGON") %>%
        # get sizes of individual units
        mutate(area = st_area(.)) %>%
        mutate(area = units::set_units(area, "km2")) %>%
        filter(area > units::set_units(7774 * 2, "km2")) %>%
        dplyr::select(brk_name, name_sort)
      plot(st_geometry(countrys.regions))
      
      # remove most easterly parts of Russia ['1' in the attirbute table]
      if(co == "Russia"){
        countrys.regions <- countrys.regions %>%
          slice(2:nrow(countrys.regions))
        plot(st_geometry(countrys.regions))
      }
      
      # make one outline
      countrys.one <- countrys.regions %>%
        st_union()
      # check geometry
      valid_countrys.one2 <- st_is_valid(countrys.one)
      
      ggplot() +
        geom_sf(data = countrys.regions, fill = "blue") +
        geom_sf(data = countrys.one, fill = "red", alpha = 0.5)
      
      # save the country as its shapefile
      st_write(countrys.one
               , file.path("data", "intermediate_outputs", "countries"
                           , paste0(co, "_final.gpkg"))
               , append = F)
    }
  }
} # co in cont.unique

#### x - write readme ####
readmePath <- file.path("data", "intermediate_outputs", "countries", "readme.md")
# Create and open the file connection
fileConn <- file(readmePath)
writeLines(c(
  "Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13"
  , "Date created: 2024-10-02"
  , paste0("Last update:",  format(Sys.Date()))
  , "Produced by 'create_countries.R'

Description of files in directory:
The files in this directory contain the final shapefiles for the countries that will be used in the analysis for the Uncertainty project.
The initial countryal shapefiles were downloaded from the rnaturalearth in R using the 'ne_countries' function. 

Files:
    data/intermediate_outputs/[country].gpkg")
  , fileConn)
close(fileConn)
