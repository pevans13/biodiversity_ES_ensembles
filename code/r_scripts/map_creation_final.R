## ---------------------------
##
## Script name: map_creation_final.R
##
## Purpose of script: to create the final hotspot maps after running the analysis
##
## Run after: uncert_graphs.R
##
## Run before: paper_graphs.R
##
## list of final outputs:
##    Maps_[ES]_SR_[continent].pdf = maps of spatial agreement
##    spatial_agreement[ES][continent].pdf = graphs of spatial agreement
##
## Author: Dr. Paul M. Evans
##
## Date Created: 2024-07-29
##
## Copyright (c) Paul M. Evans, 2024
## Email: paueva@ceh.ac.uk
## Git: https://github.com/pevans13

cat("\n\n####################################"
    , "\n\n map_creation_final.R"
    , "\n\n####################################\n\n")

#### 0 - paths ####
pathTables <- file.path("results", "tables")
options(dplyr.summarise.inform = F)

#### 0 - load functions ####
# Define a function to normalise raster values
normalise <- function(x, mx) {
  (x - min(values(x), na.rm = TRUE)) / (mx - min(values(x), na.rm = TRUE))
}

# function for separating a graphs legend
get_legend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#### 1 - load ####
# load in results that have 'END' in 
resultsList <- list.files(file.path("results", "tables")
                          , pattern = "END.csv"
                          , full.names = T)
## remove maps
resultsListMaps <- resultsList[grepl("map", tolower(resultsList))]
resultsListMaps <- resultsListMaps[grepl("100", resultsListMaps)]
resultsListMaps

### read in 
results <- lapply(resultsListMaps, function(x) { fread(x) %>%
    mutate(es = substring(basename(x), 1, 4))}) %>%
  bind_rows() %>%
  ### convert back to spatial
  separate(xy, c("x", "y"), "_") %>%
  mutate(x = gsub("m", "-", x)
         , y = gsub("m", "-", y)
  ) %>%
  st_as_sf(coords = c("x", "y"), crs = 4326
           , remove = F) # keep co-ords
head(results)
plot(results %>% st_geometry())

xGroup <- results %>%
  group_by(res, continent, pairs, es)
xKG <- group_keys(xGroup)
head(xKG)

conts <- unique(xKG$continent)
allES <- unique(xKG$es)

#### loop for all continents ####
ppGrpah <- 0; ppList <- list(); ppMaxList <- list()
patch.list <- list(); patch.save <- 0
# pblapply(conts, function(iCont) {
for(iCont in conts){
  
  cat("\nCurrent cont:", iCont, "| ")
  #### get template rasters ####
  patt <- paste0(".*", iCont, "100.*")
  
  pattTemp <- list.files(file.path("data", "intermediate_outputs",
                                   "continents", "rasters")
                         , pattern = patt)
  
  template100 <- rast(file.path("data", "intermediate_outputs",
                                "continents", "rasters", 
                                pattTemp))
  plot(template100, main = pattTemp)
  
  #### loop for all ES ####
  for(iES in allES){
    cat("ES:", iES, "| ")
    
    # filter original results
    results.filterES <- results %>%
      filter(es == iES) %>%
      # filter by cont
      filter(continent == iCont)
    
    ## round x and y, then create them as character
    sr2 <- results.filterES %>%
      mutate(rcFid = paste(x, y 
                           , sep = "_")) %>%
      st_drop_geometry()
    head(sr2)
    ### make the df long
    ppp <- sr2 %>%
      # remove duplicate rows
      distinct()  %>%
      group_by(rcFid, pairs
               , perc_70, perc_80, perc_90) %>%
      summarise(finals = sum(count, na.rm = T))
    head(ppp)
    
    BDmax <- 11
    if(iES == "Carb"){
      ESmax <- 11
      esName <- "Carbon"
    } else if(iES == "Wate"){
      ESmax <- 8
      esName <- "Water"
    } else if(iES == "Recr"){
      ESmax <- 4
      esName <- "Recreation"
    } else {
      stop("ES not identified")
    }
    
    # start pdf for each combo of ES and continent - if it does not already exist
    xFile <- file.path("results", "figures", "maps"
                       , paste0("Maps_", iES, "_SR_", iCont, ".pdf")
                       
    )
    for(i in 1:ESmax){
      x1 <- factorial(ESmax) / (factorial(i) * factorial(ESmax - i))
      x2 <- factorial(BDmax) / (factorial(i) * factorial(BDmax - i))
      # calculate possible max amount
      fxx <- x1 * x2
      cat(paste0("p:", i, " [mxC: ", fxx , "]"), "|") # max combos
      
      # make longer, to sum the pixels that were hotspots for different percentiles
      pppFil <- ppp %>%
        # filter by pairs
        filter(pairs == i) %>%
        pivot_longer(!c(rcFid, pairs, finals)) %>%
        group_by(rcFid, pairs, value, name) %>%
        summarise(finalSum = sum(finals)) %>%
        filter(!is.na(value)) %>%
        rename(perc = value) %>%
        # convert x and y back
        separate_wider_delim(rcFid, "_", names = c('x', 'y')) %>%
        # create as point data
        st_as_sf(coords = c("x", "y")
                 , crs = 4326)
      # head(pppFil)
      # plot(pppFil[1])
      
      # increase, to indicate position
      ppGrpah <- ppGrpah + 1
      
      # for saving as a graph later, get all the results that equal max possible
      pppMax <- pppFil %>%
        filter(finalSum == fxx)
      head(pppMax)
      ## save, for graphing, using count to get total pixel number
      ppMaxList[[ppGrpah]] <- pppMax %>% st_drop_geometry() %>%
        group_by(pairs, perc) %>%
        summarise(total_100pc_pix = n()) %>%
        mutate(eco_sys = iES
               , continent = iCont)
      
      # convert back to the original raster (convert sf object to terra raster)
      ## for each percentile
      stop()
      
      normGraphsPC <- lapply(c("70", "80", "90"), function(percx){
        # percx <- "70"
        cat(" Current perc:", percx, "\n")
        
        # just current perc
        pppFilPerc <- pppFil %>%
          filter(perc == percx)
        # head(pppFilPerc)
        
        # rasterise
        pppRaster <- terra::rasterize(vect(pppFilPerc)
                                      , template100
                                      , field = "finalSum")
        pppRaster[is.na(pppRaster)] <- 0
        
        # plot(pppRaster, main = perc)
        # mask 
        pppRaster <- mask(pppRaster, template100)
        # plot(pppRaster, main = percx)
        # table(as.vector(pppRaster))
        
        # normalise
        pppRasterNorm <- normalise(pppRaster, mx = fxx)
        # plot(pppRasterNorm, main = percx)
        
        # create graph
        ppGrp <- ggplot() +
          geom_spatraster(data = pppRasterNorm, aes(fill = last)) +
          ggtitle(paste0(esName, " & Species richness in ", iCont, "| P", i, " | max: ", fxx, " | ", percx)) +
          scale_fill_gradientn(
            colours = c("grey", "darkblue", "lightblue", "peru"), # Add grey to the beginning
            values = scales::rescale(c(0, 0.01, 0.99, 1)), # Ensure 'grey' is mapped to 0, 
            # others are scaled accordingly
            limits = c(0, 1),  # Set the min and max limits of the scale
            na.value = "white"
          ) +
          geom_text() +
          annotate("text", label = percx, x = Inf, y = Inf, hjust = 1
                   , vjust = 1.2, size = 6, colour = "red") +
          theme_bw() +
          theme(
            axis.title = element_blank(), # Remove axis titles
            axis.text = element_blank(),  # Remove axis text
            axis.ticks = element_blank()  # Remove axis ticks
          )
        # plot(ppGrp)
        
        return(ppGrp)
      }) # percentile
      
      # combine the output
      patchOut <- (normGraphsPC[[1]] + guides(fill = "none") + theme(title = element_blank()) | 
                     normGraphsPC[[2]] + guides(fill = "none") + 
                     ggtitle("80") + theme(title = element_blank())   | 
                     normGraphsPC[[3]] + guides(fill = "none") + 
                     ggtitle("90") + theme(title = element_blank())   |
                     # get legend
                     get_legend(normGraphsPC[[1]]  + labs(fill="Normalised\nagreement"))
      ) +
        plot_layout(widths = c(5, 5, 5), heights = c(5, 5, 5)) +
        plot_layout(ncol = 2, nrow = 2) +
        plot_annotation(title = paste0(esName, " & Species richness in ", iCont, " (ensemble size: ", i, " | max: ", fxx, ")"),
                        theme = theme(plot.title = element_text(size = 14, vjust = 0, hjust = 0.5))
                        , caption = Sys.Date())
      plot(patchOut)
      
      # increase the number and then save patchout 
      patch.save <- patch.save + 1
      cat("Current patch.save number is", patch.save, "\n")
      patch.list[[patch.save]] <- patchOut
      
    } # number in ensemble
  } # "ES"
} # "cont

#### 2 - convert spatial info to a graph ####
# use 'ppMaxList' that was saved from part 1
allPerks <- bind_rows(ppMaxList) %>%
  rename(es = eco_sys)
head(allPerks)

dfx <- allPerks %>%
  mutate(perc = paste0(perc, "th"))
head(dfx)

coly <- "total_100pc_pix"
# create the dataset, selecting the columns of interest
dfnew <- dfx %>%
  dplyr::select(pairs, all_of(coly), perc, es, continent) %>%
  group_by(perc, pairs, es, continent) %>%
  # summarise the data
  summarise(median_result = median(get(coly), na.rm = TRUE)) %>%
  ungroup()
head(dfnew)

# create broken-stick regressions
## Create a data frame - one each for the different percentiles
dfnewSplit <- dfnew %>% ungroup() %>%
  group_split(es, continent)

# make grouped data frame
groupKeys <- dfnew %>%
  group_by(es, continent)
gk <- group_keys(groupKeys)

# get all percentiles used
udfStart <- unique(dfnew$perc)

dashed.list <- list(); dashed.save <- 0
for(xall in seq_along(dfnewSplit)){
  
  if(xall == 1){
    ## open a pdf to save all graphs
    pdf(file.path("results", "figures", "maps"
                  , paste0("spatial_agreement", Sys.Date(), ".pdf")))
    
    ### create a save list, to be matched with the graph outputs from 'uncert_graph.R'
    matchGraphList2 <- list()
  }
  
  cat("\nxall =", xall, "\n\n")
  # get df
  inDf <- dfnewSplit[[xall]]
  head(inDf)
  print(inDf[1, ])
  Sys.sleep(0.5)
  
  # get names
  nms <- gk[xall, ] %>% unlist()
  
  x <- "70th"
  udfx <- lapply(udfStart, function(x){
    cat("udf =", x, "\n")
    # loop one percentile at a time
    df2 <- inDf %>%
      filter(perc == x)
    return(df2)
  })
  
  # Create a colour mapping for lines
  colorMapping <- c("70th" = "#999999"
                    , "80th" = "#E69F00"
                    , "90th" = "#56B4E9"
                    , "95th" = "#04AF70")
  
  dfxToGraph <- dfx %>%
    filter(es == nms[1]
           , continent == nms[2])  %>%
    filter(!(es == "Carb" & pairs == 11))
  
  # get names
  if(nms[1] == "Carb"){
    cct <- "Aboveground carbon storage"
  } else if(nms[1] == "Wate"){
    cct <- "Water supply"
  } else if(nms[1] == "Recr"){
    cct <- "Recreation"
  } else if(nms[1] == "sppR"){
    cct <- "Species richness"
  } else {
    stop("ES not identified")
  }
  
  # create boxplots based on the non-summarised data
  ## create title for it
  gtit <- paste0("Spatial agreement for hotspots "
                 , nms[2], "\n", nms[1], " and species richness at 100 km2")
  gtit2 <- paste0(nms[2], " | ", paste0(nms[1], " & SppRich"))
  gtit2 <- paste0(nms[2], " | ", cct, " & vascular plant richness")
  
  ## and y lab
  yllab <- paste0("Spatial agreement of ES / BD hotspots (pixels)")
  dfbox <- ggplot() + 
    geom_point(data = dfxToGraph
               , aes(x=pairs, y=get(coly), col = perc)) +
    geom_line(data = dfxToGraph
              , aes(x = pairs, y = get(coly), color = perc)) +
    theme_minimal() +
    labs(x = "Models in ensemble"
         , y = yllab
         , colour = "Percentile") +
    ggtitle(gtit2) + 
    scale_color_manual(values = colorMapping) +
    scale_x_continuous(limits = c(1, max(dfxToGraph$pairs))
                       , breaks = seq(1, max(dfxToGraph$pairs), 1))
  
  dashed.save <- dashed.save + 1
  cat("Current dashed.save number is", dashed.save, "\n")
  dashed.list[[dashed.save]] <- dfbox
  
  # barchart
  dfbar <- ggplot(data = dfxToGraph 
                  , aes(x = pairs, y = get(coly), fill = perc)) +
    geom_bar(stat = "identity"
             , position = position_dodge(preserve = "single")
    ) +
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#04AF70")) +
    theme_minimal() +
    labs(x = "Models in ensemble"
         , y = yllab
         , fill = "Percentile") +
    ggtitle(gtit2) + 
    scale_color_manual(values = colorMapping) +
    scale_x_continuous(limits = c(0.5, max(dfxToGraph$pairs)+0.5)
                       , breaks = seq(1, max(dfxToGraph$pairs), 1)
    )
  
  # pdf(file.path("results", "figures", "maps"
  #               , paste0("spatial_agreement", nms[1], nms[2], ".pdf")))
  # # save graph
  print(dfbox)
  print(dfbar)
  # dev.off()
  
  # save the results, for easy navigation later
  matchGraphList2[[xall]] <- list(nms[1]
                              , nms[2]
                              , dfbar)
  
  if(xall == length(dfnewSplit)){
    dev.off()
    
    #### save matchGraphList
    save(matchGraphList2, file = "matchGraphList2.RData")
    cat("---------- matchGraphList saved ----------\n")
    
  }
}

#### 3 - save pdf that contains maps and graphs ####
for(xall in seq_along(dashed.list)){
  
  if(xall == 1){
    ## open a pdf to save all graphs
    pdf(file.path("results", "figures", "maps"
                  , paste0("spatAgree_maps_graphs", Sys.Date(), ".pdf")))
    cat("---------- PDF opened ----------\n")
    
    ### create a page of text explaining the results
    txt <- strwrap(
      paste(
        "Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13"
        , "\nProduced by 'map_creation_final.R'"
        , paste("\nDate: ", Sys.Date(), sep = "")  # Adding the date
        , "\nThe maps and graphs in this PDF were generated by overlaying ensembles of ecosystem services (ES) and biodiversity (BD) and observing where 'hotspots' of both spatially coincided . Model outputs related to ES and BD across different continents were obtained, with ES outputs representing carbon storage (Carb), water supply (Wate), and recreation (Recr), while BD outputs always represented species richness (SppRich).

\nFor each ES or BD, ensembles were created by overlaying different combinations of model outputs and calculating the resulting median value for each pixel. For example, three different model outputs for the same ES (e.g., ES1_mod1, ES1_mod2, ES1_mod3) were overlaid to create an ensemble, and the ensemble (i.e., median) value was calculated. All possible combinations of model outputs were then used to generate the full set of ensembles. The number of model outputs used in an ensemble is referred to as 'ensemble size,' which ranged from 1 to the maximum number of model outputs available for that ES or BD measure. Each ensemble was normalised, resulting in values between 0 and 1. Then, ES and BD ensembles that used the same number of models (i.e., ensemble size) were overlaid to see where ES and BD hotspots overlapped. 

\n'Hotspots' were identified as areas (pixels) where an ensemble exceeded specific thresholds, with the 70th, 80th, and 90th percentiles used as cut-offs.

\nThe results of all the same ensemble size category were summed and then normalised, creating the resulting maps below. '1' was equivalent to 100% of ensemble size ensembles ES and BD comparisons completely overlapping, this is known as 'spatial agreement'. '1' is indicated as brown on the maps, while '0' (i.e., never overlapped) is grey. The red number in the topright corner of each map indicates the percentile that was used as the cut-off. 

\nThe graphs that follow each set of maps shows the spatial agreement (i.e., only 100% agreement between hotspots) of each combination of ES and BD as the ensemble size increases. The dashed vertical lines show a statistically significant breakpoint in the data, where the gradients of the line before and after differ. In this case, it tends to indicate the point at which certainty in results starts increasing. " 
        , sep = "\n")
      , width = 105)
    
    plot.new()
    # readme txt
    text(x = .1, y = 1
         , "README"
         , font=1, cex=2, col="#F48024")
    # info txt
    text(x = -0.05, y = .3 # first 2 numbers are xy-coordinates within [0, 1]
         , paste(txt, collapse = "\n")
         , adj = c(0.5, 1) # 'centred'
         , cex = .65
         , pos = 4)  
  }
  
  # get the ES identifier from 'gk'
  esNow <- gk[xall, "es"] %>% as.character()
  cat("The current ES is", esNow, "\n")
  ## set the number of map figures to include per section (rhis changes per ES)
  useFigs <- ifelse(esNow == "Carb", 11
                    , ifelse(esNow == "Recr", 4
                             , ifelse(esNow == "Wate", 8, NA)))
  ## select the next figures from the 'patch.list' 
  if(xall == 1){
    useRefs <- 1:useFigs
    useRefStart <- max(useRefs) + 1
  } else {
    useRefs <- useRefStart: (useRefStart + useFigs - 1)
    useRefStart <- max(useRefs) + 1
    stopifnot(length(useRefs) == useFigs)
  }
  Sys.sleep(1)
  cat("The current useFigs is", useFigs, "\n")
  cat("The current useRefStart is", useRefStart, "\n")
  cat("The current useRefs is", useRefs, "\n")
  cat("\n")
  
  # print all of the maps for this combination
  for(xPrint in useRefs){
    print(patch.list[[xPrint]])
  }
  # ## then print the graph
  # print(dashed.list[[xall]])
  
  if(xall == length(dashed.list)){
    dev.off()
    cat("---------- PDF closed ----------\n")
  }
}

#### x - write readme ####
readmePath <- file.path("results", "figures", "maps", "readme.md")
# Create and open the file connection
fileConn <- file(readmePath)
writeLines(c(
  "Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13"
  , "Date created: 2024-08-20"
  , paste0("Last update:",  format(Sys.Date()))
  , "Produced by 'map_creation_final.R'

Description of files in directory:
  There are two types of pdf saved in this directory: 'Maps_[ES]_SR_[continent].pdf' and 'spatial_agreement[ES][continent].pdf' where:
    - [ES] = ecosystem service. It relates to the service being analysed, with 'Carb' = carbon storage, 'Wate' = water supply, and 'Recr' = recreation.
    - [continent] = continent of analysis. 
    - The 'SR' indicates that species richness was also included.
  
  The files with the 'Maps' prefix show where there was spatial agreement between different numbers of ensembled models. Spatial agreement is where the hotspots of different ensembled models overlap.
    For these maps, results have been normalised with 1 being the max (as shown in the map header)
  The files with the 'spatial_agreement' prefix represent the '1's from the maps (i.e., maximum agreement)")
  , fileConn)
close(fileConn)
