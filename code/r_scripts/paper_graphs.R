## ---------------------------
##
## Script name: paper_graphs.R
##
## Purpose of script: To create the final graphs that will feature in the
##                    paper for the project.
##                    
## ------------ Notes --------------  ##
## This script also contains code for adding some certain continents and 
## ES combinations to include in presentations.
## ------------ ----- --------------  ##
##
## Run after: map_creation_final.R
##
## Specific numbered tasks:
## 1 - Create figures for the main manuscript
## 2 - create graphs for presentations - 'singlePercGraphs'
## 3 - create final tables and PDFs for the SI
##
## list of final outputs:
##    results/tables/final_results.csv
##    results/tables/summExitResults.csv
##    all files in results/figures
##
## Author: Dr. Paul M. Evans
##
## Date Created: 2024-09-25
##
## Copyright (c) Paul M. Evans, 2024
## Email: paueva@ceh.ac.uk
## Git: https://github.com/pevans13

## ---------------------------
options(scipen = 6, digits = 4) # for non-scientific notation
## ---------------------------
# Create a colour mapping for lines
colorMapping <- c("70th" = "#999999"
                  , "80th" = "#E69F00"
                  , "90th" = "#56B4E9"
                  , "95th" = "#04AF70")

#### 0 - load libraries ####
## automatic install of packages if they are not installed already
list.of.packages <- c(
  "terra", "raster", "dplyr", "tidyr"
  , "ggplot2", "tidyterra", "gridExtra"
  , "patchwork"
  , "sf", "stars", "fasterize", "exactextractr"
  , "tictoc", "data.table", "parallel", "pbapply"
  , "matrixStats" # for 'weightedMedian' in dplyr
  # for use of the 'Winsorize' function [https://search.r-project.org/CRAN/refmans/DescTools/html/Winsorize.html]
  , "DescTools" 
  , "beepr" # for sounds when things are finished
  , "rnaturalearth" # for shapefiles
  , "foreach", "doParallel" # for parallel procesing
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

# Multiple plot function
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

### create a function to plot boxplots for two percentiles ###
# dfx = howManyResults[[1]]
# coly = "Overlap"
# amt = "total"

brokenBoxPlots <- function(dfx, coly, amt
                           , bdx = "species richness"
                           , yLab = paste0("TAO")
                           # , yLab = paste0("Overlapping ", amt, " hotspots (pixels)")
                           , xxlab = "Ensemble size"){
  # dfx <- results.combined[[1]]; xall = 9; x = "70th"; coly = "Overlap"; amt = "total"
  # dfx <- results.single2[[1]]; xall = 6; x = "90th";  coly = "median_HSs"; amt = "total"; bdx = "NA"
  dfx2 <- dfx %>%
    mutate(perc = paste0(perc, "th"))
  head(dfx2)
  
  # create the dataset, selecting the columns of interest
  dfnew <- dfx2 %>%
    dplyr::select(pairs, all_of(coly), perc, es, res, continent) %>%
    group_by(perc, pairs, es, res, continent) %>%
    # summarise the data
    summarise(median_result = median(get(coly), na.rm = TRUE)
              , .groups="keep") %>%
    ungroup()
  head(dfnew)
  
  # create broken-stick regressions
  ## Create a data frame - one each for the different percentiles
  dfnewSplit <- dfnew %>% ungroup() %>%
    group_split(es, res, continent)
  
  # make grouped data frame
  groupKeys <- dfnew %>%
    group_by(es, res, continent)
  gk <- group_keys(groupKeys)
  
  udf <- unique(dfnew$perc)
  
  pblapply(seq_along(dfnewSplit), function(xall) {
    
    cat("\nxall =", xall, "\n\n")
    # get df
    inDf <- dfnewSplit[[xall]]
    head(inDf)
    print(inDf[1, ])
    Sys.sleep(1)
    
    # get names
    nms <- gk[xall, ] %>% unlist()
    
    udfx <- lapply(udf, function(x){
      cat("udf =", x, "\n")
      # loop one percentile at a time
      df2 <- inDf %>%
        filter(perc == x)
      
      dataf <- data.frame(y = df2$median_result
                          , x = df2$pairs %>% as.numeric()
                          , perc = df2$perc)
      
      return(list(dataf))
    })
    
    dfxToGraph <- dfx %>%
      filter(es == nms[1]
             , res == nms[2]
             , continent == nms[3])
    
    # create boxplots based on the non-summarised data
    ## create title for it
    if(bdx == "species richness"){
      gtit <- paste0("Overlapping ", amt, " hotspots for "
                     , nms[3], "\n", nms[1], " and species richness at ", nms[2] , " km2")
      gtit2 <- paste0(nms[3], " | ", nms[1], " & SppRich")
    } else {
      gtit <- paste0("Overlapping ", amt, " hotspots for "
                     , nms[3], "\n", nms[1], " at ", nms[2] , " km2")
      gtit2 <- paste0(nms[3], " | ", nms[1])
    }
    gtit2 <- nms[3]
    
    # remove 11 esNo if carbon, due to only having one result
    if(es == "Carb"){
      dfxToGraph <- dfxToGraph %>%
        filter(pairs != 11)
    }
    
    ## and y lab
    yllab <- yLab
    dfbox <- ggplot() + 
      geom_boxplot(data = dfxToGraph
                   , aes(x=as.factor(pairs), y=get(coly), fill = perc)
                   , position=position_dodge(1)
                   , notch=TRUE
                   , outliers = FALSE
                   , outlier.shape = NA) +
      scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#04AF70")) +
      theme_minimal() +
      labs(x = xxlab
           , y = yllab
           , fill = "Percentile") +
      ggtitle(gtit2)
    
    dfbox <- dfbox + 
      scale_color_manual(values = colorMapping) +
      expand_limits(y = 0) +
      scale_y_continuous(expand = c(0, 0)) +
      guides(colour="none")
    return(dfbox)
  })
}

### create a function to plot boxplots for two percentiles ###
# and add the broken-stick analysis in addition
brokenBoxPlotsCountry <- function(dfx, coly, amt
                                  , bdx = "species richness"
                                  , yLab = paste0("TAO")
                                  , xxlab = "Ensemble size"){
  dfx2 <- dfx %>%
    mutate(perc = paste0(perc, "th"))
  head(dfx2)
  
  # create the dataset, selecting the columns of interest
  dfnew <- dfx2 %>%
    dplyr::select(pairs, all_of(coly), perc, es, res, continent, country) %>%
    group_by(perc, pairs, es, res, continent, country) %>%
    # summarise the data
    summarise(median_result = median(get(coly), na.rm = TRUE)
              , .groups="keep") %>%
    ungroup()
  head(dfnew)
  
  # create broken-stick regressions
  ## Create a data frame - one each for the different percentiles
  dfnewSplit <- dfnew %>% ungroup() %>%
    group_split(es, res, continent, country)
  
  # make grouped data frame
  groupKeys <- dfnew %>%
    group_by(es, res, continent, country)
  gk <- group_keys(groupKeys)
  
  udf <- unique(dfnew$perc)
  
  pblapply(seq_along(dfnewSplit), function(xall) {
    
    cat("\nxall =", xall, "\n\n")
    # get df
    inDf <- dfnewSplit[[xall]]
    head(inDf)
    print(inDf[1, ])
    Sys.sleep(1)
    
    # get names
    nms <- gk[xall, ] %>% unlist()
    
    udfx <- lapply(udf, function(x){
      cat("udf =", x, "\n")
      # loop one percentile at a time
      df2 <- inDf %>%
        filter(perc == x)
      
      dataf <- data.frame(y = df2$median_result
                          , x = df2$pairs %>% as.numeric()
                          , perc = df2$perc)
      
      return(list(dataf))
    })
    
    dfxToGraph <- dfx %>%
      filter(es == nms[1]
             , res == nms[2]
             , continent == nms[3]
             , country == nms[4])
    
    # create boxplots based on the non-summarised data
    ## title
    gtit1 <- paste0(nms[4], " [", nms[3], "] | ", cct, " & ", "species richness")
    gtit2 <- paste0(nms[4], " [", nms[3], "]")
    ## and y lab
    yllab <- yLab
    dfbox <- ggplot() + 
      geom_boxplot(data = dfxToGraph
                   , aes(x=as.factor(pairs), y=get(coly), fill = perc)
                   , position=position_dodge(1)
                   , notch=TRUE
                   , outliers = FALSE
                   , outlier.shape = NA) +
      scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#04AF70")) +
      theme_minimal() +
      labs(x = xxlab
           , y = yllab
           , fill = "Percentile") +
      ggtitle(gtit1)
    
    dfbox <- dfbox + 
      scale_color_manual(values = colorMapping) +
      expand_limits(y = 0) +
      scale_y_continuous(expand = c(0, 0)) +
      guides(colour="none")
    return(dfbox)
  })
}

### create a function to plot important information from the boxplots ###
# and add the broken-stick analysis in addition
brokenBoxDetails <- function(dfx, coly){
  # dfx <- howManyResults[[1]]; xall = 1; x = "80th"; coly = "Overlap"; amt = "total"
  dfx2 <- dfx %>%
    mutate(perc = paste0(perc, "th"))
  head(dfx2)
  
  # create the dataset, selecting the columns of interest
  dfnew <- dfx2 %>%
    dplyr::select(pairs, all_of(coly), perc, es, res, continent) %>%
    group_by(perc, pairs, es, res, continent) %>%
    # summarise the data
    summarise(median_result = median(get(coly), na.rm = TRUE)
              , .groups="keep") %>%
    ungroup()
  head(dfnew)
  
  # create broken-stick regressions
  ## Create a data frame - one each for the different percentiles
  dfnewSplit <- dfnew %>% ungroup() %>%
    group_split(es, res, continent)
  
  # make grouped data frame
  groupKeys <- dfnew %>%
    group_by(es, res, continent)
  gk <- group_keys(groupKeys)
  
  udf <- unique(dfnew$perc)
  
  pblapply(seq_along(dfnewSplit), function(xall) {
    
    cat("\nxall =", xall, "\n\n")
    # get df
    inDf <- dfnewSplit[[xall]]
    head(inDf)
    print(inDf[1, ])
    Sys.sleep(1)
    
    # get names
    nms <- gk[xall, ] %>% unlist()
    
    dfxToGraph <- dfx %>%
      filter(es == nms[[1]]
             , res == nms[[2]]
             , continent == nms[[3]])
    head(dfxToGraph)
    
    # remove 11 esNo if carbon, due to only having one result
    if(es == "Carb"){
      dfxToGraph <- dfxToGraph %>%
        filter(pairs != 11)
    }
    
    # summarise the data that would be dsplayed on a boxplot
    dfxToGraphSumm <- dfxToGraph %>%
      group_by(es, res, continent, perc, pairs) %>%
      summarise(
        median_value = median(Overlap, na.rm = TRUE),
        lower_quartile = quantile(Overlap, 0.25, na.rm = TRUE),
        upper_quartile = quantile(Overlap, 0.75, na.rm = TRUE),
        mean_value = mean(Overlap, na.rm = TRUE)
      )
    head(dfxToGraphSumm)
    return(dfxToGraphSumm)
  })
}

brokenBoxDetailsCountry <- function(dfx, coly){
  dfx2 <- dfx %>%
    mutate(perc = paste0(perc, "th"))
  head(dfx2)
  
  # create the dataset, selecting the columns of interest
  dfnew <- dfx2 %>%
    dplyr::select(pairs, all_of(coly), perc, es, res, continent, country) %>%
    group_by(perc, pairs, es, res, continent, country) %>%
    # summarise the data
    summarise(median_result = median(get(coly), na.rm = TRUE)
              , .groups="keep") %>%
    ungroup()
  head(dfnew)
  
  # create broken-stick regressions
  ## Create a data frame - one each for the different percentiles
  dfnewSplit <- dfnew %>% ungroup() %>%
    group_split(es, res, continent, country)
  
  # make grouped data frame
  groupKeys <- dfnew %>%
    group_by(es, res, continent, country)
  gk <- group_keys(groupKeys)
  
  udf <- unique(dfnew$perc)
  
  pblapply(seq_along(dfnewSplit), function(xall) {
    
    cat("\nxall =", xall, "\n\n")
    # get df
    inDf <- dfnewSplit[[xall]]
    head(inDf)
    print(inDf[1, ])
    Sys.sleep(1)
    
    # get names
    nms <- gk[xall, ] %>% unlist()
    
    dfxToGraph <- dfx %>%
      filter(es == nms[[1]]
             , res == nms[[2]]
             , continent == nms[[3]]
             , country == nms[[4]]
             , perc == "80")
    head(dfxToGraph)
    
    # summarise the data that would be dsplayed on a boxplot
    dfxToGraphSumm <- dfxToGraph %>%
      group_by(es, res, continent, country, perc, pairs) %>%
      summarise(
        median_value = median(Overlap, na.rm = TRUE),
        lower_quartile = quantile(Overlap, 0.25, na.rm = TRUE),
        upper_quartile = quantile(Overlap, 0.75, na.rm = TRUE),
        mean_value = mean(Overlap, na.rm = TRUE)
      )
    head(dfxToGraphSumm)
    return(dfxToGraphSumm)
  })
}

##### summarise the df #####
# create a function to load in all tables, and summarise that data, getting 
# mean and standard deviations for each group combinations.
summariseDfFunc <- function(inList # list of files to load in and combine
                            , coly = "Overlap" # column to analyse
                            , colSelt = c("es", "res", "continent", "perc", "pairs"
                                          , "Overlap", "one_hotpot", "no_hotspots")){
  
  # read in all results, extracting first four characters
  results <- lapply(inList, function(x) {
    cat(x, "\n")
    fread(x) %>%
      mutate(es = substring(basename(x), 1, 4))
  }) %>%
    bind_rows() %>%
    # get cols of interest, and remove 95th percentile
    dplyr::select(all_of(colSelt)) %>% filter(perc != "95")
  results$perc <- as.factor(results$perc)
  print(head(results))
  
  # Convert coly to a symbol for tidy evaluation
  coly_sym <- sym(coly)
  
  # Get averages - mean and std
  results.summary <- results %>%
    group_by(es, res, continent, perc, pairs) %>%
    summarise(
      mean_overlap = mean(!!coly_sym, na.rm = TRUE),
      sum_overlap = sum(!!coly_sym, na.rm = TRUE),
      sd_overlap = sd(!!coly_sym, na.rm = TRUE),
      count = n()
    )
  
  # make grouped data frame - to be able to keep track of the different combinations
  groupKeys <- results %>%
    group_by(es, res, continent)
  gk <- group_keys(groupKeys)
  
  # output both dfs
  return(list(results, results.summary, gk))
}

#### 0 - choose which bits to run ####
# mainFigures <- F
# SI <- T
# singlePercGraphs <- T

#### 0 - create empty lists to save ####
whereMapsComboList <- list(); whereMapsCombo <- 0 

#### 1 - Create figures for the main manuscript ####
# get dfs for each ES
if(mainFigures){
  # load in results that have 'END' in 
  resultsList <- list.files(file.path("results", "tables")
                            , pattern = "END.csv"
                            , full.names = T
                            , recursive = T)
  ## remove single-item
  resultsList <- resultsList[!grepl("/ES/|/BD/", resultsList)]
  ## keep only those with '100' in
  resultsList <- resultsList[grepl("100", resultsList)]
  resultsList
  
  ## get all the unique results in terms of ES
  esx <- gsub(".*/(.+)100.*", "\\1", resultsList) %>% unique()
  ## remove 'sppRich'
  esx <- esx[esx != "sppRich"]
  ### separate files into mapped ('where') and not ('howMany')
  howManyList <- resultsList[!grepl("map", resultsList, ignore.case = T)]
  whereList <- resultsList[grepl("map", resultsList, ignore.case = T)]
  
  # deal with each ES in turn
  for(es in esx){
    
    cat("-----------------")
    cat("starting", es)
    cat("-----------------\n")
    
    # get max combinations
    BDmax <- 11
    if(es == "Carb"){
      ESmax <- 11
      cct <- "Aboveground carbon storage"
      esColour.map <- "darkgreen"
    } else if(es == "Wate"){
      ESmax <- 8
      cct <- "Water supply"
      esColour.map <- "darkblue"
    } else if(es == "Recr"){
      ESmax <- 4
      cct <- "Recreation"
      esColour.map <- "#E69F00"
    } else {
      stop("ES not identified")
    }
    
    # reduce lists to current ES
    howManyES <- howManyList[grepl(es, howManyList)]
    whereES <- whereList[grepl(es, whereList)]
    
    ##### create summarised dfs for 'how many' graphs #####
    # load in data, and summarise
    howManyResults <- summariseDfFunc(inList = howManyES
                                      , colSelt = c("es", "res", "continent", "perc", "pairs"
                                                    , "Overlap"))
    head(howManyResults)
    
    ##### create boxplots for 'how many' graphs #####
    howManyBoxplot <- brokenBoxPlots(howManyResults[[1]]
                                     , coly = "Overlap"
                                     , amt = "total")
    ## get summary information
    howManyDetails <- brokenBoxDetails(howManyResults[[1]], coly = "Overlap") %>%
      bind_rows()
    head(howManyDetails)
    
    ##### create summarised dfs for 'where' graphs #####
    whereResults <- lapply(whereES, function(x) { fread(x) %>%
        mutate(es = substring(basename(x), 1, 4))}) %>%
      bind_rows() %>%
      ### convert back to spatial
      separate(xy, c("x", "y"), "_") %>%
      mutate(x = gsub("m", "-", x)
             , y = gsub("m", "-", y)
      ) %>%
      st_as_sf(coords = c("x", "y"), crs = 4326
               , remove = F) # keep co-ords
    head(whereResults)
    head(howManyResults)
    
    # use howManyResults[[3]] to ensure the continents are done in the same order
    conts <- howManyResults[[3]]$continent %>% unique()
    
    cnt <- 0
    togetherGraphs <- list(); togetherGraphsCount <- 0
    for(contx in conts){
      cnt <- cnt + 1
      # stop("stop, cnt")
      cat("\nCurrent cont:", contx, "| ")
      # get template rasters 
      patt <- paste0(".*", contx, "100.*")
      
      pattTemp <- list.files(file.path("data", "intermediate_outputs",
                                       "continents", "rasters")
                             , pattern = patt)
      
      template100 <- rast(file.path("data", "intermediate_outputs",
                                    "continents", "rasters", 
                                    pattTemp))
      plot(template100)
      cat("Pixels:"
          , table(template100 %>% as.vector()) %>%
            as.data.frame() %>% dplyr::select(Freq) %>% as.numeric()
          , "\n")
      # Sys.sleep(10)
      
      # filter original results
      whereResultsCont <- whereResults %>%
        # filter by cont
        filter(continent == contx)
      head(whereResultsCont)
      
      ## round x and y, then create them as character
      sr2 <- whereResultsCont %>%
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
      
      # start pdf for each combo of ES and continent - if it does not already exist
      xFile <- file.path("results", "figures", "maps"
                         , paste0("Maps_", es, "_SR_", contx, ".pdf"))
      cat("save file:", xFile
          , paste0("[", es, " | ", contx, "]")
          , "\n")
      
      # get number of pixels of max coverage
      ppMaxList <- lapply(1:ESmax, function(i){
        x1 <- factorial(ESmax) / (factorial(i) * factorial(ESmax - i))
        x2 <- factorial(BDmax) / (factorial(i) * factorial(BDmax - i))
        # calculate possible max amount
        fxx <- x1 * x2
        cat(paste0("p:", i, " [mxC: ", fxx , "]"), "\n") # max combos
        
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
          
          # graph label
          gLabel <- paste0("Ensemble size: ", i
                           # , " [Max = ", fxx, "]"
          )
          
          # set text position
          if(contx == "Africa"){
            hj <- 1.8
            vj <- 5
          } else if(contx == "South America"){
            hj <- 1.3
            vj <- 5
          } else {
            hj <- 0
            vj <- 5
          }
          
          # stop("before maps with blue-brown")
          
          # create graph
          ppGrp <- ggplot() +
            geom_spatraster(data = pppRasterNorm, aes(fill = last)) +
            # ggtitle(paste0(cct, " & vascular plant richness in ", contx)) +
            ggtitle(gLabel) +
            scale_fill_gradientn(
              colours = c("grey", esColour.map), # Add grey to the beginning
              values = scales::rescale(c(0, 1)), # Ensure 'grey' is mapped to 0, 
              # others are scaled accordingly
              limits = c(0, 1),  # Set the min and max limits of the scale
              na.value = "white"
            ) +
            geom_text() +
            theme_bw() +
            theme(
              axis.title = element_blank(), # Remove axis titles
              axis.text = element_blank(),  # Remove axis text
              axis.ticks = element_blank()  # Remove axis ticks
              , plot.title = element_text(size=14) # change title size
            )
          # plot(ppGrp)
          
          # create graph with only the ones
          ## limit to 1
          pppRay <- ifel(pppRasterNorm == 1, 1, 0)
          plot(pppRay)
          ppGrpRay <- ggplot() +
            geom_spatraster(data = pppRay, aes(fill = last)) +
            ggtitle(gLabel) +
            scale_fill_gradientn(
              colours = c("grey", esColour.map), # Add grey to the beginning
              values = scales::rescale(c(0, 1)), # Ensure 'grey' is mapped to 0,
              breaks = c(1)
              , labels = c("Spatial Agreement")
              # others are scaled accordingly
              , limits = c(0, 1),  # Set the min and max limits of the scale
              na.value = "white"
            ) +
            guides(
              fill = guide_legend(
                title = "Spatial Agreement", 
                override.aes = list(fill = esColour.map),   # Only show the legend for category '3'
                nrow = 1, ncol = 1,                        # Keep legend in a single row
                label.position = "bottom",                # Position labels if needed
                keyheight = unit(1, "cm"),                # Adjust the size of the legend square
                keywidth = unit(1, "cm")                  # Adjust the size of the legend square
              )) +
            geom_text() +
            theme_bw() +
            theme(
              axis.title = element_blank(), # Remove axis titles
              axis.text = element_blank(),  # Remove axis text
              axis.ticks = element_blank()  # Remove axis ticks
              , plot.title = element_text(size=14) # change title size
              , legend.title = element_blank()
            )
          return(ppGrpRay)
        }) # percentile
        
        # for saving as a graph later, get all the results that equal max possible
        pppMax <- pppFil %>%
          filter(finalSum == fxx)
        head(pppMax)
        ## save, for graphing, using count to get total pixel number
        ppM <- pppMax %>% st_drop_geometry() %>%
          group_by(pairs, perc) %>%
          summarise(total_100pc_pix = n()) %>%
          mutate(eco_sys = es, continent = contx)      
        return(list(normGraphsPC, ppM))
      })
      
      # str(ppMaxList)
      if(es == "Wate" & contx == "South America"){
        mapsaveWater <- lapply(c(1, 5, 8), function(pp){
          ppMaxList[[pp]][[1]][[2]]
          # stop("maps and water")
        })
      }
      
      if(es == "Carb" & contx == "Africa"){
        mapsaveCarb <- lapply(c(1, 6, 10), function(pp){
          ppMaxList[[pp]][[1]][[2]]
        })
      }
      
      # get max info
      whereDetails <- lapply(c(1:ESmax), function(pp){
        ppMaxList[[pp]][[2]]
      }) %>% bind_rows() %>%
        filter(perc == 80)
      # get max info
      whereDetailsAll <- lapply(c(1:ESmax), function(pp){
        ppMaxList[[pp]][[2]]
      }) %>% bind_rows() 
      stopifnot(length(unique(whereDetailsAll$total_100pc_pix)) > 2)
      if(contx == conts[[1]]){
        whereDetailsAllComb <- whereDetailsAll
      } else {
        whereDetailsAllComb <- bind_rows(whereDetailsAllComb
                                         , whereDetailsAll)
      }
      
      ##### create graphs for 'where' graphs #####
      # get all the maxes
      ppMaxDf <- lapply(1:ESmax, function(pp){
        ppMaxList[[pp]][[2]]
      }) %>% bind_rows()  %>%
        rename(es = eco_sys)
      head(ppMaxDf)
      # convert perc to factor
      whereGraphdf <- ppMaxDf %>%
        mutate(perc = paste0(perc, "th")) %>%
        mutate(esNo = as.factor(pairs))
      head(whereGraphdf)
      
      # remove 11 esNo if carbon, due to only having one result
      if(es == "Carb"){
        whereGraphdf <- whereGraphdf %>%
          filter(esNo != 11)
      }
      
      # create barcharts based on the non-summarised data
      ## create title for it
      gtit <- paste0("Spatial agreement for hotspots "
                     , contx, "\n", es, " and species richness at 100 km2")
      gtit2 <- paste0(contx, " | ", paste0(es, " & SppRich"))
      gtit2 <- paste0(contx)
      ## and y lab
      # yllab <- paste0("Spatial agreement of ES / BD hotspots (pixels)")
      yllab <- paste0("SpAg")
      # barchart
      dfbar <- ggplot(data = whereGraphdf
                      , aes(x = esNo, y = total_100pc_pix, fill = perc)) +
        geom_bar(stat = "identity"
                 , position = position_dodge(preserve = "single")) +
        scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#04AF70")) +
        theme_minimal() +
        labs(x = "Ensemble size"
             , y = yllab
             , fill = "Percentile") +
        ggtitle(gtit2)
      
      if(contx == conts[[1]]){
        togetherGraphsCount <- 1
        togetherGraphs[[togetherGraphsCount]] <- howManyBoxplot[[cnt]]
        togetherGraphsCount <- togetherGraphsCount + 1
        togetherGraphs[[togetherGraphsCount]] <- dfbar
      } else {
        togetherGraphsCount <- togetherGraphsCount + 1
        togetherGraphs[[togetherGraphsCount]] <- howManyBoxplot[[cnt]]
        togetherGraphsCount <- togetherGraphsCount + 1
        togetherGraphs[[togetherGraphsCount]] <- dfbar
      }
      
      print(length(togetherGraphs))
      # stop("end cont")
    } # end of 'contx'
    
    togetherGraphs[[7]]
    togetherGraphs[[8]]
    
    patchOut <- (
      togetherGraphs[[1]] + guides(fill = "none") + togetherGraphs[[2]] + guides(fill = "none") +
        togetherGraphs[[3]] + guides(fill = "none") + togetherGraphs[[4]] + guides(fill = "none") +
        togetherGraphs[[5]] + guides(fill = "none") + togetherGraphs[[6]] + guides(fill = "none") + 
        togetherGraphs[[7]] + guides(fill = "none") + togetherGraphs[[8]] +
        togetherGraphs[[9]] + guides(fill = "none") + togetherGraphs[[10]] + guides(fill = "none") +
        togetherGraphs[[11]] + guides(fill = "none") + togetherGraphs[[12]] + guides(fill = "none")
    ) + plot_layout(ncol = 4, nrow = 3, widths = rep(1, 2), heights = rep(1, 3), axis_titles = "collect") +
      plot_annotation(title = paste0(cct, " & vascular plant richness")
                      , theme = theme(plot.title = element_text(size = 16, vjust = 0, hjust = 0.5))
                      , caption = Sys.Date())
    patchOut
    
    # save the image
    cat("---------- patchOut saved ----------\n")
    thatsnumberwang <- ifelse(es == "Recr", 3
                              , ifelse(es == "Wate", 2
                                       , ifelse(es == "Carb", 1)
                              ))
    png(file.path("results", "figures",
                  paste0("Fig", thatsnumberwang, "_ES_", es, ".png"))
        , res = 200
        , height = 2000, width = 2500)
    print(patchOut)
    dev.off()
    
    # get the summary information for all ES
    if(es == esx[[1]]){
      summInfo <- howManyDetails
      si2 <- whereDetailsAllComb
    } else {
      summInfo <- bind_rows(summInfo
                            , howManyDetails)
      si2 <- bind_rows(si2, whereDetailsAllComb)
    }
    stopifnot("70" %in% unique(si2$perc))
    stopifnot("70" %in% unique(summInfo$perc))
    stopifnot(length(unique(si2$total_100pc_pix)) > 2)
    
  } # end of 'es'
  
  # stop("911 - es end")
  
  # Create the first patchwork of carbon maps
  patchOutMaps1 <- (
    mapsaveCarb[[1]] + theme(plot.title = element_text(size = 16)) +
      mapsaveCarb[[2]] + theme(plot.title = element_text(size = 16)) +
      mapsaveCarb[[3]]  + theme(plot.title = element_text(size = 16))
  ) + 
    plot_annotation(
      title = paste0("Aboveground carbon storage & vascular plant richness in Africa"),
      theme = theme(plot.title = element_text(size = 18, vjust = 0, hjust = 0.5))
    ) +
    plot_layout(guides = "collect") &
    theme(legend.position="bottom",
          legend.text = element_text(size = 16),
          legend.spacing.x = unit(0, 'cm')) &
    guides(fill = guide_legend(label.position = "bottom"))
  
  patchOutMaps1
  # save as a R data in case changes are required
  save(patchOutMaps1, file = file.path("results", "figures",
                                       paste0("Fig", 4, "_carbon_africa", ".RData")))
  
  png(file.path("results", "figures",
                paste0("Fig", 4, "_carbon_africa", ".png"))
      , res = 200
      , height = 1000, width = 2000)
  print(patchOutMaps1)
  dev.off()
  
  # Create the second patchwork of water maps
  patchOutMaps2 <- (
    mapsaveWater[[1]] +
      mapsaveWater[[2]] +
      mapsaveWater[[3]]
  ) + 
    plot_annotation(
      title = paste0("Water supply & vascular plant richness in South America"),
      theme = theme(plot.title = element_text(size = 16, vjust = 0, hjust = 0.5))
    ) +
    plot_layout(guides = "collect") &
    theme(legend.position="bottom",
          legend.text = element_text(size = 16),
          legend.spacing.x = unit(0, 'cm')) &
    guides(fill = guide_legend(label.position = "bottom"))
  
  patchOutMaps2
  # save as a R data in case changes are required
  save(patchOutMaps2, file = file.path("results", "figures",
                                       paste0("Fig", 5, "_water_Sout", ".RData")))
  
  png(file.path("results", "figures",
                paste0("Fig", 5, "_water_Sout", ".png"))
      , res = 200
      , height = 1000, width = 1500)
  print(patchOutMaps2)
  dev.off()
  patchOutMaps2
  
  multiplot(mapsaveWater[[1]] + theme(title = element_blank()) + guides(fill = "none")
            , mapsaveCarb[[1]] + theme(title = element_blank()) + guides(fill = "none")
            , mapsaveWater[[2]] + theme(title = element_blank()) + guides(fill = "none")
            , mapsaveCarb[[2]] + theme(title = element_blank()) + guides(fill = "none")
            , mapsaveWater[[3]] + theme(title = element_blank()) + guides(fill = "none")
            , mapsaveCarb[[3]] + theme(title = element_blank()) + guides(fill = "none")
            , cols = 3)
  
  # merge the table information together
  summExit <- merge(summInfo, si2 %>% rename(es = eco_sys)
                    , all = T) %>%
    rename(Percentile = perc
           , "Ecosystem service" = es
           , Contient = continent
           , SpAg = total_100pc_pix) %>%
    mutate(SpAg = ifelse(is.na(SpAg), 0, SpAg))
  head(summExit)
  # save the table
  fwrite(summExit, file.path("results", "tables", "final_results.csv"), row.names = F)
  
  # stop("before main table")
  
  ## summarise the table, determining per cent change with the addition of each next stage
  summExitGroups <- summExit %>%
    rename(es = "Ecosystem service") %>%
    group_split(es, Contient, Percentile)
  summExitResults <- pblapply(summExitGroups, function(x){
    x2 <- x %>% dplyr::select(-res) %>%
      mutate(across(everything(), unname)) %>%
      # get absolute differences between each esNo for medians
      mutate(abs_dif_med = median_value - lag(median_value)
             # get differences between each esNo for medians (in %)
             , pc_dif_med = (median_value - lag(median_value)) / lag(median_value) * 100
             # get cumulative differences between each esNo for medians (in %)
             , pc_cum_med = (median_value - first(median_value)) / first(median_value) * 100
      ) %>%
      ## repeat above but for spag
      mutate(abs_dif_SpAg = SpAg - lag(SpAg)
             # get differences between each esNo for medians (in %)
             , pc_dif_SpAg = (SpAg - lag(SpAg)) / lag(SpAg) * 100
             # get cumulative differences between each esNo for medians (in %)
             , pc_cum_SpAg = (SpAg - first(SpAg)) / first(SpAg) * 100
      ) %>%
      # get interquartile range
      mutate(IQR = upper_quartile - lower_quartile
             , pc_dif_IQR = (IQR - lag(IQR)) / lag(IQR) * 100
             , pc_cum_IQR = (IQR - first(IQR)) / first(IQR) * 100
      )
  }) %>% bind_rows() %>%
    rowwise() %>%
    # calculate SpAg as a percentage of country pixels
    mutate(
      pc_SpAg_ctry = case_when(
        Contient == "South America" ~ (SpAg / 2379) * 100,
        Contient == "North America" ~ (SpAg / 6145) * 100,
        Contient == "Oceania" ~ (SpAg / 1282) * 100,
        Contient == "Europe" ~ (SpAg / 2448) * 100,
        Contient == "Asia" ~ (SpAg / 8526) * 100,
        Contient == "Africa" ~ (SpAg / 3901) * 100,
        TRUE ~ NA_real_
      )
    ) %>%
    # organise by type
    relocate(es, Contient, Percentile, pairs
             , median_value, abs_dif_med, pc_dif_med, pc_cum_med
             , lower_quartile, upper_quartile, IQR, pc_dif_IQR, pc_cum_IQR) %>%
    # only to two decimal places
    mutate(across(5:last_col(), \(x) round(x, 2))) %>%
    filter(pairs != 11)
  head(summExitResults)
  # save the table
  fwrite(summExitResults, file.path("results", "tables", "summExitResults.csv"), row.names = F)
  
  # get results per continent
  ExitResultsShort <- summExitResults %>%
    group_by(es, Percentile, pairs) %>%
    mutate(across(everything(), unname)) %>%
    summarise(across(starts_with(c("median", "abs", "IQR", "pc")), list(mean = ~ mean(.x, na.rm = TRUE)))) %>%
    # and get per cent change
    mutate(pc_med = (median_value_mean - lag(median_value_mean)) / lag(median_value_mean) * 100
           , pc_iqr = (IQR_mean - lag(IQR_mean)) / lag(IQR_mean) * 100
    ) %>%
    filter(pairs != 11) %>%
    ungroup() %>%
    # only to two decimal places
    mutate(across(where(is.numeric), \(x) round(x, 2)))
  head(ExitResultsShort)
  # save the table
  fwrite(ExitResultsShort, file.path("results", "tables", "ExitResultsShort.csv"), row.names = F)
} # end 'mainFigures'

#### 2 - create graphs for presentations - 'singlePercGraphs' ####
if(singlePercGraphs){
  
  # stop("singlePercGraphs")
  
  # create graphs for only the specified percentiles
  pcRuns <- which(spg) # identifies which ones are true
  for(percentx in pcRuns){
    
    currentPerc <- ifelse(percentx == 1, 70
                          , ifelse(percentx == 2, 80
                                   , ifelse(percentx == 3, 90, NA)))
    
    # load in results that have 'END' in 
    resultsList <- list.files(file.path("results", "tables")
                              , pattern = "END.csv"
                              , full.names = T
                              , recursive = T)
    ## remove single-item
    resultsList <- resultsList[!grepl("/ES/|/BD/", resultsList)]
    ## keep only those with '100' in
    resultsList <- resultsList[grepl("100", resultsList)]
    resultsList
    
    ## get all the unique results in terms of ES
    esx <- gsub(".*/(.+)100.*", "\\1", resultsList) %>% unique()
    ## remove 'sppRich'
    esx <- esx[esx != "sppRich"]
    ### separate files into mapped ('where') and not ('howMany')
    howManyList <- resultsList[!grepl("map", resultsList, ignore.case = T)]
    whereList <- resultsList[grepl("map", resultsList, ignore.case = T)]
    
    cnt <- 0
    togetherGraphs <- list()
    areaGraphs <- list()
    
    # deal with each ES in turn
    for(es in esx){
      
      cat("-----------------")
      cat("starting", es)
      cat("-----------------\n")
      
      # get max combinations
      # get title (cct)
      # get colour for separate ES
      BDmax <- 11
      if(es == "Carb"){
        ESmax <- 11
        cct <- "Aboveground carbon storage"
        esColour.box <- "lightgreen"
        esColour.bar <- "darkgreen"
      } else if(es == "Wate"){
        ESmax <- 8
        cct <- "Water supply"
        esColour.box <- "lightblue"
        esColour.bar <- "darkblue"
      } else if(es == "Recr"){
        ESmax <- 4
        cct <- "Recreation"
        esColour.box <- "sandybrown"
        esColour.bar <- "#E69F00"
      } else {
        stop("ES not identified")
      }
      
      # reduce lists to current ES
      howManyES <- howManyList[grepl(es, howManyList)]
      whereES <- whereList[grepl(es, whereList)]
      
      ##### create summarised dfs for 'how many' graphs #####
      # load in data, and summarise
      howManyResults <- summariseDfFunc(inList = howManyES
                                        , colSelt = c("es", "res", "continent", "perc", "pairs"
                                                      , "Overlap"))
      head(howManyResults)
      
      # use howManyResults[[3]] to ensure the continents are done in the same order
      conts <- howManyResults[[3]]$continent %>% unique()
      conts
      
      for(contx in conts){
        cnt <- cnt + 1
        # stop("stop, cnt")
        cat("\nCurrent cont:", contx, "| ")
        # get template rasters 
        patt <- paste0(".*", contx, "100.*")
        
        pattTemp <- list.files(file.path("data", "intermediate_outputs",
                                         "continents", "rasters")
                               , pattern = patt)
        
        template100 <- rast(file.path("data", "intermediate_outputs",
                                      "continents", "rasters", 
                                      pattTemp))
        plot(template100, col = esColour.bar)
        cat("Pixels:"
            , table(template100 %>% as.vector()) %>%
              as.data.frame() %>% dplyr::select(Freq) %>% as.numeric()
            , "\n")
        # Sys.sleep(10)
        
        ##### create boxplots for 'how many' graphs #####
        dfx <- howManyResults[[1]] %>%
          filter(perc == currentPerc) %>%
          filter(continent == contx) %>%
          rename(esnw = es) %>%
          filter(esnw == es) %>%
          rename(es = esnw)
        
        dfx2 <- dfx %>%
          mutate(perc = paste0(perc, "th"))
        head(dfx2)
        
        # create the dataset, selecting the columns of interest
        dfnew <- dfx2 %>%
          dplyr::select(pairs, Overlap, perc, es, res, continent) %>%
          group_by(perc, pairs, es, res, continent) %>%
          # summarise the data
          summarise(median_result = median(Overlap, na.rm = TRUE)
                    , .groups="keep") %>%
          ungroup()
        head(dfnew)
        
        cat("\nxall =", 1, "\n\n")
        # get df
        inDf <- dfnew
        head(inDf)
        print(inDf[1, ])
        Sys.sleep(1)
        
        # get names
        dataf <- data.frame(y = inDf$median_result
                            , x = inDf$pairs %>% as.numeric()
                            , perc = inDf$perc)
        dfxToGraph <- dfx 
        
        # create boxplots based on the non-summarised data
        ## create title for it
        gtit2 <- paste0(contx, " | ", cct, " & vascular plant richness")
        
        # remove 11 esNo if carbon, due to only having one result
        if(es == "Carb"){
          dfxToGraph <- dfxToGraph %>%
            filter(pairs != 11)
        }
        
        ##### boxplot #####
        ## and y lab
        yllab <- expression(paste("Overlap Area (pixels)"))
        dfbox <- ggplot() + 
          geom_boxplot(data = dfxToGraph
                       , aes(x=as.factor(pairs), y=Overlap, fill = perc)
                       , position=position_dodge(1)
                       , notch=TRUE
                       , outliers = FALSE
                       , outlier.shape = NA) +
          scale_fill_manual(values = c("70" = "#999999", "80" = esColour.box
                                       , "90" = "#56B4E9", "95" = "#04AF70")) +
          theme_minimal() +
          labs(x = "Ensemble size"
               , y = yllab
               , fill = "Percentile") +
          ggtitle(gtit2)
        
        dfbox <- dfbox + 
          # expand_limits(y = 0) +
          scale_y_continuous(expand = c(0, 0)) +
          guides(colour="none")
        dfbox
        
        ## get summary information
        howManyDetails <- brokenBoxDetails(howManyResults[[1]], coly = "Overlap") %>%
          bind_rows()
        head(howManyDetails)
        
        ##### create summarised dfs for 'where' graphs #####
        whereResults <- lapply(whereES, function(x) { fread(x) %>%
            mutate(es = substring(basename(x), 1, 4))}) %>%
          bind_rows() %>%
          ### convert back to spatial
          separate(xy, c("x", "y"), "_") %>%
          mutate(x = gsub("m", "-", x)
                 , y = gsub("m", "-", y)
          ) %>%
          st_as_sf(coords = c("x", "y"), crs = 4326
                   , remove = F) # keep co-ords
        head(whereResults)
        head(howManyResults)
        # stop("howManyResults")
        
        # filter original results
        ## determine which two to not include
        notInc <- c("perc_70", "perc_80", "perc_90")
        Inc <- paste0("perc_", currentPerc)
        notInc <- notInc[-which(notInc %in% Inc)] # removes the one to keep
        whereResultsCont <- whereResults %>%
          # filter by cont
          filter(continent == contx) %>%
          dplyr::select(-c(all_of(notInc), perc_95)) %>%
          # filter out NAs in perc of interest
          filter(!is.na(.data[[Inc]]))
        head(whereResultsCont)
        
        ## round x and y, then create them as character
        sr2 <- whereResultsCont %>%
          mutate(rcFid = paste(x, y 
                               , sep = "_")) %>%
          st_drop_geometry()
        head(sr2)
        ### make the df long
        ppp <- sr2 %>%
          # remove duplicate rows
          distinct()  %>%
          group_by(rcFid, pairs, all_of(Inc)) %>%
          summarise(finals = sum(count, na.rm = T))
        head(ppp)
        
        # start pdf for each combo of ES and continent - if it does not already exist
        xFile <- file.path("results", "figures", "maps"
                           , paste0(es, "_", currentPerc, "_", contx
                                    , Sys.Date(), ".pdf"))
        cat("save file:", xFile
            , paste0("[", es, " | ", contx, "]")
            , "\n")
        
        # get number of pixels of max coverage
        ppMaxList <- lapply(1:ESmax, function(i){
          x1 <- factorial(ESmax) / (factorial(i) * factorial(ESmax - i))
          x2 <- factorial(BDmax) / (factorial(i) * factorial(BDmax - i))
          # calculate possible max amount
          fxx <- x1 * x2
          cat(paste0("p:", i, " [mxC: ", fxx , "]"), "\n") # max combos
          
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
          
          normGraphsPC <- lapply(currentPerc, function(percx){
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
            
            # graph label
            gLabel <- paste0("esNo: ", i, " [", fxx, "]")
            
            # set text position
            if(contx == "Africa"){
              hj <- 1.8
              vj <- 5
            } else if(contx == "South America"){
              hj <- 1.3
              vj <- 5
            } else {
              hj <- 0
              vj <- 5
            }
            
            ##### maps #####
            # create graph
            ppGrp <- ggplot() +
              geom_spatraster(data = pppRasterNorm, aes(fill = last)) +
              # ggtitle(paste0(cct, " & vascular plant richness in ", contx)) +
              ggtitle(gLabel) +
              scale_fill_gradientn(
                colours = c("grey", "darkblue", "lightblue", "peru"), # Add grey to the beginning
                values = scales::rescale(c(0, 0.01, 0.99, 1)), # Ensure 'grey' is mapped to 0, 
                # others are scaled accordingly
                limits = c(0, 1),  # Set the min and max limits of the scale
                na.value = "white"
              ) +
              geom_text() +
              # annotate("text", label = gLabel, x = Inf, y = Inf
              #          , hjust = hj
              #          # , vjust = "bottom"
              #          , vjust = vj
              #          , size = 4, colour = "red") +
              theme_bw() +
              theme(
                axis.title = element_blank(), # Remove axis titles
                axis.text = element_blank(),  # Remove axis text
                axis.ticks = element_blank()  # Remove axis ticks
              )
            # plot(ppGrp)
            
            return(ppGrp)
          }) # percentile
          
          # for saving as a graph later, get all the results that equal max possible
          pppMax <- pppFil %>%
            filter(finalSum == fxx)
          head(pppMax)
          ## save, for graphing, using count to get total pixel number
          ppM <- pppMax %>% st_drop_geometry() %>%
            group_by(pairs, perc) %>%
            summarise(total_100pc_pix = n()) %>%
            mutate(eco_sys = es, continent = contx)      
          return(list(normGraphsPC, ppM))
        })
        
        # get max info
        whereDetails <- lapply(c(1:ESmax), function(pp){
          ppMaxList[[pp]][[2]]
        }) %>% bind_rows() %>%
          filter(perc == currentPerc)
        # get max info
        whereDetailsAll <- lapply(c(1:ESmax), function(pp){
          ppMaxList[[pp]][[2]]
        }) %>% bind_rows() 
        stopifnot(length(unique(whereDetailsAll$total_100pc_pix)) > 2)
        if(contx == conts[[1]]){
          whereDetailsAllComb <- whereDetailsAll
        } else {
          whereDetailsAllComb <- bind_rows(whereDetailsAllComb
                                           , whereDetailsAll)
        }
        
        ##### create graphs for 'where' graphs #####
        # get all the maxes
        ppMaxDf <- lapply(1:ESmax, function(pp){
          ppMaxList[[pp]][[2]]
        }) %>% bind_rows()  %>%
          rename(es = eco_sys)
        head(ppMaxDf)
        # convert perc to factor
        whereGraphdf <- ppMaxDf %>%
          mutate(perc = paste0(perc, "th")) %>%
          mutate(esNo = as.factor(pairs))
        head(whereGraphdf)
        
        # remove 11 esNo if carbon, due to only having one result
        if(es == "Carb"){
          whereGraphdf <- whereGraphdf %>%
            filter(esNo != 11)
        }
        
        # Add the new levels with value 0 (if not all are present)
        ## see which are current
        cnn <- as.numeric(as.character(whereGraphdf$esNo))
        cFac <- which(!c(1:max(cnn)) %in% cnn)
        whereGraphdf <- bind_rows(
          data.frame(esNo = factor(c(cFac), levels = c(1:max(cnn))), value = c(rep(0, length(cFac)))),
          whereGraphdf
        ) %>%
          mutate(perc = paste0("perc_", currentPerc, "th"))
        
        # create barcharts based on the non-summarised data
        ## create title for it
        gtit2 <- paste0(contx, " | ", cct, " & vascular plant richness")
        ## and y lab
        # yllab <- paste0("Spatial agreement of ES / BD hotspots (pixels)")
        yllab <- expression(paste("Spatial Agreement (pixels)"))
        
        ##### barchart #####
        # barchart
        dfbar <- ggplot(data = whereGraphdf
                        , aes(x = esNo, y = total_100pc_pix, fill = perc)) +
          geom_bar(stat = "identity"
                   , position = position_dodge(preserve = "single")) +
          scale_fill_manual(values = c("perc_70th" = "#999999", "perc_80th" = esColour.bar
                                       , "perc_90th" = "#56B4E9", "perc_95th" = "#04AF70")
                            , labels = c(currentPerc)) +
          # scale_x_discrete(expand = c(0, 0), limits = c(0, max(cnn))) +
          theme_minimal() +
          labs(x = "Ensemble size"
               , y = yllab
               , fill = "Percentile") +
          ggtitle(gtit2)
        dfbar
        
        togetherGraphs[[cnt]] <- dfbar
        areaGraphs[[cnt]] <- dfbox
        
        print(length(togetherGraphs))
        # stop("end cont")
      } # end of 'contx'
      
      # save the image
      thatsnumberwang <- ifelse(es == "Recr", 3
                                , ifelse(es == "Wate", 2
                                         , ifelse(es == "Carb", 1)
                                ))
    }
    
    # start pdf for each combo of ES and continent - if it does not already exist
    xFile <- file.path("results", "figures"
                       , paste0("singlePercGraphs_", currentPerc, "_", Sys.Date(), ".pdf"))
    cat("save file:", xFile
        , "\n")
    pdf(xFile)
    for(i in seq_along(togetherGraphs)){
      
      plot(togetherGraphs[[i]])
      plot(areaGraphs[[i]])
      
    }
    dev.off()
    cat("---------- patchOut saved ----------\n")
    
    # save all graphs together
    ## first, Extract the title for all entries in the list togetherGraphs
    titles <- strsplit(sapply(togetherGraphs, function(graph) graph$labels$title), " | ",  fixed=TRUE)
    ### ensure it matches total area
    titlesAO <- strsplit(sapply(areaGraphs, function(graph) graph$labels$title), " | ",  fixed=TRUE)
    stopifnot(identical(titles, titlesAO))
    titlesp1 <- sapply(titles, function(x) x[[1]]); titlesp1
    titlesp2 <- sapply(titles, function(x) x[[2]]); titlesp2
    utitle <- unique(gsub("^(.+) &.*", "\\1", titlesp2))
    
    xxList <- list() # list to save the 6 continent graphs for each ES
    for(x in utitle){
      cat(x, "| ")
      xChoice <- which(grepl(x, titlesp2)) # extract appropriate graphs
      
      for(xx in xChoice){
        cat(xx, "\n")
        
        if(!xx %in% c(2, 5, 8, 11, 14, 17)){
          
          xxList[[xx]] <- (areaGraphs[[xx]] + 
                             theme(axis.title.y = element_blank()) + 
                             theme(plot.title = element_blank()
                                   , axis.title.x = element_text(size = 14)) + guides(fill = "none") 
                           | 
                             togetherGraphs[[xx]] +
                             theme(axis.title.y = element_blank()) + 
                             theme(plot.title = element_blank()
                                   , axis.title.x = element_text(size = 14)) + guides(fill = "none")) +
            plot_annotation(title = titlesp1[[xx]],
                            tag_levels = 'A'
                            , theme = theme(plot.title = element_text(hjust = 0.5)  # Centre the title horizontally
                            ))
        } else {
          
          xxList[[xx]] <- (areaGraphs[[xx]] + 
                             theme(plot.title = element_blank()
                                   , axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14)
                             ) + guides(fill = "none")
                           | 
                             togetherGraphs[[xx]] + 
                             theme(plot.title = element_blank()
                                   , axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14)) + 
                             guides(fill = "none")) +
            plot_annotation(title = titlesp1[[xx]],
                            tag_levels = 'A'
                            , theme = theme(plot.title = element_text(hjust = 0.5)  # Centre the title horizontally
                            ))
        }
      }
      
      patchOut <- (wrap_elements(xxList[[xChoice[1]]]) / wrap_elements(xxList[[xChoice[2]]]) / wrap_elements(xxList[[xChoice[3]]])
                   |
                     wrap_elements(xxList[[xChoice[4]]]) / wrap_elements(xxList[[xChoice[5]]]) / wrap_elements(xxList[[xChoice[6]]])
      ) +
        plot_annotation(title = titlesp2[[xChoice[1]]]
                        , theme = theme(plot.title = element_text(size = 16, vjust = 0, hjust = 0.5))
                        , caption = Sys.Date())
      patchOut
      # stop("after patch in single-perc")
      
      # save the image
      cat("---------- patchOut saving", x, "----------\n")
      es <- ifelse(grepl("carbon", x, ignore.case = T), "Carb"
                   , ifelse(grepl("recreation", x, ignore.case = T), "Recr"
                            , ifelse(grepl("wate", x, ignore.case = T), "Wate")))
      
      thatsnumberwang <- ifelse(es == "Recr", 3
                                , ifelse(es == "Wate", 2
                                         , ifelse(es == "Carb", 1)))
      
      png(file.path("results", "figures",
                    paste0("Fig", thatsnumberwang, "_", currentPerc, "_", es, Sys.Date(), ".png"))
          , res = 200
          , height = 2000, width = 2500)
      print(patchOut)
      dev.off()
      
      # save as a R data in case changes are required
      save(patchOut, file = file.path("results", "figures",
                                      paste0("Fig", thatsnumberwang, "_", currentPerc, "_", es, ".RData")))
      
      # stop("after saving bar and box plots")
      
    } # different title 'utitle'
  } # percentx in pcRuns
} # end 'singlePercGraphs'

#### 3 - create final tables and PDFs for the SI ####
if(SI){
  
  # start pdf for each combo of ES and continent - if it does not already exist
  xFile <- file.path("results", "figures"
                     , paste0("SI_maps_SpatAgree", ".pdf")
                     
  ); xFile
  
  ##### 3.1 - continents and Spatial Agreement #####
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
  # plot(results %>% st_geometry())
  
  xGroup <- results %>%
    group_by(res, continent, pairs, es)
  xKG <- group_keys(xGroup)
  head(xKG)
  
  # get unique features
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
        cct <- "Aboveground carbon storage"
        esColour.box <- "lightgreen"
        esColour.bar <- "darkgreen"
      } else if(iES == "Wate"){
        ESmax <- 8
        esName <- "Water"
        cct <- "Water supply"
        esColour.box <- "lightblue"
        esColour.bar <- "darkblue"
      } else if(iES == "Recr"){
        ESmax <- 4
        esName <- "Recreation"
        cct <- "Recreation"
        esColour.box <- "sandybrown"
        esColour.bar <- "#E69F00"
      } else {
        stop("ES not identified")
      }
      esColour.map <- esColour.bar
      
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
          pppRasterNorm[pppRasterNorm != 1] <- 0  # make 1/0
          
          gLabel <- paste0(esName, " & Species richness in ", iCont, " | Ens", i, " | max: ", fxx, " | ", percx)
          # create graph
          ppGrp <- ggplot() +
            geom_spatraster(data = pppRasterNorm, aes(fill = last)) +
            ggtitle(gLabel) +
            scale_fill_gradientn(
              colours = c("grey", esColour.map), # Add grey to the beginning
              values = scales::rescale(c(0, 1)), # Ensure 'grey' is mapped to 0,
              breaks = c(1)
              , labels = c("Spatial Agreement")
              # others are scaled accordingly
              , limits = c(0, 1),  # Set the min and max limits of the scale
              na.value = "white"
            ) +
            guides(
              fill = guide_legend(
                title = "Spatial Agreement", 
                override.aes = list(fill = esColour.map),   # Only show the legend for category '3'
                nrow = 1, ncol = 1,                        # Keep legend in a single row
                label.position = "bottom",                # Position labels if needed
                keyheight = unit(1, "cm"),                # Adjust the size of the legend square
                keywidth = unit(1, "cm")                  # Adjust the size of the legend square
              )) +
            geom_text() +
            annotate("text", label = percx, x = Inf, y = Inf, hjust = 1
                     , vjust = 1.2, size = 6, colour = "red") +
            theme_bw()
          theme_bw() +
            theme(
              axis.title = element_blank(), # Remove axis titles
              axis.text = element_blank(),  # Remove axis text
              axis.ticks = element_blank()  # Remove axis ticks
              , plot.title = element_text(size=14) # change title size
              , legend.title = element_blank()
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
          plot_annotation(title = paste0(esName, " & Species richness in ", iCont, " (esNo: ", i, " | max: ", fxx, ")"),
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
  
  ### save pdf that contains maps and graphs 
  for(xall in seq_along(patch.list)){
    
    if(xall == 1){
      ## open a pdf to save all graphs
      pdf(file.path("results", "figures"
                    , paste0("Spatial_Agreement_figures_", Sys.Date(), ".pdf")))
      cat("---------- PDF opened ----------\n")
      
      ### create a page of text explaining the results
      txt <- strwrap(
        paste(
          "Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13"
          , "\nProduced by 'paper_graphs.R'"
          , paste("\nDate: ", Sys.Date(), sep = "")  # Adding the date
          , "\nThe maps and graphs in this PDF were generated by overlaying ensembles of ecosystem services (ES) and biodiversity (BD) and observing where 'hotspots' of both spatially coincided. Model outputs related to ES and BD across different continents were obtained, with ES outputs representing carbon storage (Carb), water supply (Wate), and recreation (Recr), while BD outputs always represented species richness (SppRich).

\nFor each ES or BD, ensembles were created by overlaying different combinations of model outputs and calculating the resulting median value for each pixel. For example, three different model outputs for the same ES (e.g., ES1_mod1, ES1_mod2, ES1_mod3) were overlaid to create an ensemble, and the ensemble (i.e., median) value was calculated. All possible combinations of model outputs were then used to generate the full set of ensembles. The number of model outputs used in an ensemble is referred to as 'esNo,' which ranged from 1 to the maximum number of model outputs available for that ES or BD measure. Each ensemble was normalised, resulting in values between 0 and 1. Then, ES and BD ensembles that used the same number of models (i.e., esNo) were overlaid to see where ES and BD hotspots overlapped. 

\n'Hotspots' were identified as areas (pixels) where an ensemble exceeded specific thresholds, with the 70th, 80th, and 90th percentiles used as cut-offs.

\nThe results of all the same esNo category were summed and then normalised, creating the resulting maps below. '1' was equivalent to 100% of esNO ensembles ES and BD comparisons completely overlapping, this is known as 'spatial agreement'. '1' is indicated as brown on the maps, while '0' (i.e., never overlapped) is grey. The red number in the topright corner of each map indicates the percentile that was used as the cut-off. 

\nThe graphs that follow each set of maps shows the spatial agreement (i.e., only 100% agreement between hotspots) of each combination of ES and BD as the esNo increases. The dashed vertical lines show a statistically significant breakpoint in the data, where the gradients of the line before and after differ. In this case, it tends to indicate the point at which certainty in results starts increasing. " 
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
    
    cat("xall =", xall, "\n")
    patch.list[[xall]]
    # print all of the maps for this combination
    print(patch.list[[xall]])
    
    if(xall == length(patch.list)){
      dev.off()
      cat("---------- PDF closed ----------\n")
    }
  }
  
  
  
  ##### 3.2 - countries and Spatial Agreement / Overlap #####    
  ppMaxList <- list()
  
  ## ------------ Notes --------------  ##
  ## These tables consist of continental- and country-level data
  ## ------------ ----- --------------  ##
  # load in results that have 'END' in 
  countriesList <- list.files(file.path("results", "tables")
                              , pattern = "END.csv"
                              , full.names = T
                              , recursive = T)
  ## remove single-item
  countriesList <- countriesList[!grepl("/ES/|/BD/", countriesList)]
  ## make 2 lists - continent and country
  ### keep those with '100' in (i.e., continent data)
  continentList <- countriesList[grepl("100", countriesList)]
  ### remove those with '100' in (i.e., continent data)
  countriesList <- countriesList[!grepl("100", countriesList)]
  countriesList
  continentList
  
  ## get all the unique results in terms of ES
  esx <- gsub(".*/(.+?)_.*", "\\1", countriesList) %>% unique()
  ### separate files into mapped ('where') and not ('howMany')
  howManyList <- countriesList[!grepl("map", countriesList, ignore.case = T)]
  whereList <- countriesList[grepl("map", countriesList, ignore.case = T)]
  
  ## ------------ Notes --------------  ##
  ## Circle through each final result, saving a row of summary results,
  ## and saving the PDF of the final graphs
  
  ## Do this per ES, to save computer memory
  ## ------------ ----- --------------  ##
  
  for(esi in c("Wate", "Recr", "Carb")[2]){
    
    # get colour for separate ES
    BDmax <- 11
    if(esi == "Carb"){
      ESmax <- 11
      cct <- "Aboveground carbon storage"
      esColour.box <- "lightgreen"
      esColour.bar <- "darkgreen"
    } else if(esi == "Wate"){
      ESmax <- 8
      cct <- "Water supply"
      esColour.box <- "lightblue"
      esColour.bar <- "darkblue"
    } else if(esi == "Recr"){
      ESmax <- 4
      cct <- "Recreation"
      esColour.box <- "sandybrown"
      esColour.bar <- "#E69F00"
    } else {
      stop("ES not identified")
    }
    
    # get hw many for countries
    howManyListES <- howManyList[grepl(esi, howManyList)]
    
    # start pdf
    pdf(paste0("results/figures/SI_graphs", esi, Sys.Date(), ".pdf")
        , paper="USr"
        , height = 8, width = 10)
    fcCount <- 0
    for(fr in howManyListES){
      fcCount <- fcCount + 1
      
      # load in
      cat("Reading in", fr, "\n")
      frIn <- fread(fr) %>%
        mutate(perc = paste0(perc, "th"))
      head(frIn)
      
      # get names
      nameCountry <- frIn$country[1]
      nameContinent <- frIn$continent[1]
      nameContinent <- ifelse(nameCountry == "Russia", "Asia", nameContinent)
      nameEsx <- gsub(".*/(.+?)_.*", "\\1", fr)
      nameRun <- gsub(".*/(.+)_END.*", "\\1", fr)
      
      cat("\n")
      cat(nameEsx, "for", nameCountry)
      cat(" started at", format(Sys.time(), "%Y-%m-%d %H:%M"), "| number:", fcCount)
      cat("\n")
      
      # create the dataset, selecting the columns of interest
      frIn2 <- frIn %>%
        mutate(es = nameEsx) %>%
        dplyr::select(pairs, Overlap, perc, es, res, country, continent) %>%
        group_by(perc, pairs, es, res, country, continent) %>%
        # summarise the data
        summarise(median_result = median(Overlap, na.rm = TRUE)
                  , percentile_25 = quantile(Overlap, probs = 0.25, na.rm = TRUE)
                  , percentile_75 = quantile(Overlap, probs = 0.75, na.rm = TRUE)
                  # see how many combinations were used
                  , combinations = n()
                  , .groups="keep") %>%
        ungroup()
      head(frIn2)
      # get unique combinations - this is to check that some combinations were not possible
      uniCom <- frIn2 %>% dplyr::select(pairs, combinations) %>% distinct()
      
      # get max combinations
      BDmax <- 11
      if(nameEsx == "Carb"){
        ESmax <- 11
        cct <- "Aboveground carbon storage"
      } else if(nameEsx == "Wate"){
        ESmax <- 8
        cct <- "Water supply"
      } else if(nameEsx == "Recr"){
        ESmax <- 4
        cct <- "Recreation"
      } else {
        stop("ES not identified")
      }
      
      # load in equivalent df produced from maps
      currentMap <- whereList[grepl(nameRun, whereList)]
      stopifnot(length(currentMap) == 1)
      ## load in 
      currentMapIn <- fread(currentMap)
      head(currentMapIn)
      ## get top pairs
      topEsNo <- max(currentMapIn$pairs, na.rm = T)
      ### limit to ten if carbon
      if(nameEsx == "Carb"){
        topEsNo <- ifelse(topEsNo > 10, 10, topEsNo)
      }
      ## round x and y, then create them as character
      sr2 <- currentMapIn %>%
        separate(xy, c("x", "y"), "_") %>%
        mutate(x = gsub("m", "-", x)
               , y = gsub("m", "-", y)) %>%
        mutate(rcFid = paste(x, y 
                             , sep = "_")) 
      head(sr2)
      ### make the df long
      ppp <- sr2 %>%
        pivot_longer(c(perc_70, perc_80, perc_90)) %>%
        rename(finals = value
               , perc = name)
      head(ppp, 10)
      
      for(i in 1:ESmax){
        # make longer, to sum the pixels that were hotspots for different percentiles
        pppFil <- ppp %>%
          # filter by pairs
          filter(pairs == i) %>%
          group_by(rcFid, pairs, perc) %>%
          summarise(finalSum = sum(finals)) %>%
          filter(!is.na(finalSum)) 
        head(pppFil)
        
        fxx <- uniCom[uniCom$pairs == i, 2] %>% as.numeric()
        # for saving as a graph later, get all the results that equal max possible
        pppMax <- pppFil %>%
          filter(finalSum == uniCom[uniCom$pairs == i, 2] %>% as.numeric())
        head(pppMax)
        ## save, for graphing, using count to get total pixel number
        pp2 <- pppMax %>% st_drop_geometry() %>%
          group_by(pairs, perc) %>%
          summarise(total_100pc_pix = n()) 
        if(i == 1){
          ppOut <- pp2
        } else {
          ppOut <- bind_rows(ppOut, pp2)
        }
        
        # if 'i' is the highest, create maps
        if(i == topEsNo){
          
          # get reference raster based on continent
          patt <- paste0(".*", nameContinent, "100.*")
          pattTemp <- list.files(file.path("data", "intermediate_outputs",
                                           "continents", "rasters")
                                 , pattern = patt)
          template100 <- rast(file.path("data", "intermediate_outputs",
                                        "continents", "rasters", 
                                        pattTemp))
          ## get country, to be able to mask 
          pattTemp <- list.files(file.path("data", "intermediate_outputs",
                                           "countries")
                                 , full.names = T
                                 , pattern = paste0("^", nameCountry, "_final.gpkg$"))
          pattTemp
          ### get rasterised gpkg
          cGpkg <- st_read(pattTemp)
          
          # just current perc
          pppFilPerc <- pppFil %>%
            # convert x and y back
            separate_wider_delim(rcFid, "_", names = c('x', 'y')) %>%
            # create as point data
            st_as_sf(coords = c("x", "y")
                     , crs = 4326) %>%
            filter(perc == "perc_80")
          ## rasterise df
          pt2 <- crop(template100, cGpkg) %>%
            mask(., cGpkg)
          # plot(pt2)
          
          pppRaster <- terra::rasterize(vect(pppFilPerc)
                                        , pt2
                                        , field = "finalSum")
          pppRaster[is.na(pppRaster)] <- 0
          ## mask
          pppRaster <- mask(pppRaster, cGpkg)
          # plot(pppRaster)
          ## normalise
          pppRasterNorm <- normalise(pppRaster, mx = fxx)
          # plot(pppRasterNorm)
          
          # set title
          gLabel <- paste0(cct, " & species richness in ", nameCountry, "\nesNo: ", i, " | max: ", fxx, " | 80th")
          
          # create graph
          ppGrp <- ggplot() +
            geom_spatraster(data = pppRasterNorm, aes(fill = last)) +
            ggtitle(gLabel) +
            scale_fill_gradientn(
              colours = c("grey", esColour.map), # Add grey to the beginning
              values = scales::rescale(c(0, 1)), # Ensure 'grey' is mapped to 0,
              breaks = c(1)
              , labels = c("Spatial Agreement")
              # others are scaled accordingly
              , limits = c(0, 1),  # Set the min and max limits of the scale
              na.value = "white"
            ) +
            guides(
              fill = guide_legend(
                title = "Spatial Agreement", 
                override.aes = list(fill = esColour.map),   # Only show the legend for category '3'
                nrow = 1, ncol = 1,                        # Keep legend in a single row
                label.position = "bottom",                # Position labels if needed
                keyheight = unit(1, "cm"),                # Adjust the size of the legend square
                keywidth = unit(1, "cm")                  # Adjust the size of the legend square
              )) +
            geom_text() +
            theme_bw() +
            theme(
              axis.title = element_blank(), # Remove axis titles
              axis.text = element_blank(),  # Remove axis text
              axis.ticks = element_blank()  # Remove axis ticks
              , plot.title = element_text(size=14) # change title size
              , legend.title = element_blank()
            )
        }
      }
      if(nrow(ppOut)>0){
        head(ppOut)
        ppOut <- ppOut %>%
          mutate(perc = gsub("perc_", "", perc)
                 , perc = paste0(perc, "th"))
        head(ppOut)
        # Create a complete grid of all levels
        complete_grid <- expand.grid(pairs = 1:max(ppOut$pairs, na.rm = T)
                                     , perc = c("70th", "80th", "90th"))
        ppOut <- ppOut %>% merge(., complete_grid, all = T)  %>%
          mutate(total_100pc_pix = replace_na(total_100pc_pix, 0))
        head(ppOut)
        
        ## ------------ Notes --------------  ##
        ## Print (and combine) the boxplot of 'how many' and barchart
        ## of 'where' so that they are both on the same page
        ## ------------ ----- --------------  ##
        comboSummary <- frIn2 %>%
          merge(ppOut, all = T) %>% dplyr::select(-combinations) %>%
          rename(spag = total_100pc_pix)
        head(comboSummary)
        
        # boxplot
        dfbox <- ggplot() + 
          geom_boxplot(data = frIn
                       , aes(x=as.factor(pairs), y=Overlap, fill = perc)
                       , position=position_dodge(1)
                       , notch=TRUE
                       , outlier.shape = NA) +
          scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#04AF70")) +
          theme_minimal() +
          labs(x = "Ensemble size"
               , y = "TAO"
               , fill = "Percentile") +
          ggtitle(paste0(nameCountry, " [in ", nameContinent, "] | ", cct)) +
          scale_color_manual(values = colorMapping) +
          expand_limits(y = 0) +
          scale_y_continuous(expand = c(0, 0)) +
          guides(colour="none")
        
        # bar chart
        dfbar <- ggplot(data = ppOut
                        , aes(x = as.factor(pairs), y = total_100pc_pix, fill = perc)) +
          geom_bar(stat = "identity"
                   , position = position_dodge(preserve = "single")) +
          scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#04AF70")) +
          theme_minimal() +
          labs(x = "Ensemble size"
               , y = "SpAg"
               , fill = "Percentile")
        
        ## create a combined plot to print
        combined_plot <- 
          (dfbox + 
             guides(fill = "none") +
             theme(plot.title = element_blank(),
                   axis.title.x = element_blank())) /  # Remove x-axis title from the top plot
          (dfbar +
             theme(plot.title = element_blank(), legend.position = "bottom") +
             guides(fill = guide_legend(nrow = 1))) +
          plot_annotation(title = paste0(nameCountry, " [in ", nameContinent, "] | ", cct),
                          theme = theme(plot.title = element_text(size = 16, vjust = 0, hjust = 0.5)),
                          caption = Sys.Date())
        
        print(combined_plot)
        print(ppGrp)
        
        if(fcCount == 1){
          finTable <- comboSummary
        } else {
          finTable <- finTable %>% bind_rows(., comboSummary)
        }
        
      }
      cat("Finished at", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
      # if(fcCount == 83){
      #   stop("stop 83")
      # }
    }
    dev.off()
    
    finTable <- finTable %>%
      relocate(es, continent, country) %>%
      filter(perc == "80th") %>%
      arrange(es, continent, country, pairs)
    
    # save the final table
    fwrite(finTable
           , file.path("results/tables", paste0("SI_table_", esi, Sys.Date(), ".csv"))
           , row.names = F)
    
  } # end of ES
} # end 'SI'
