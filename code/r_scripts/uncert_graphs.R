## ---------------------------
##
## Script name: uncert_graphs.R
##
## Purpose of script: create the final graphs for the uncertainty project. The final graphs
##                    will consist of boxplots for the 'how much area' graphs, and
##                    barcharts for the 'where' graphs
##                    
## ------------ Notes --------------  ##
## As of 24/09/24, the graphs should be split into three panels, based on the ES.
## Each panel will have 4 x 3 graphs, with the graphs of each continent being adjacent,
## box on the left, bar on the right.
## ------------ ----- --------------  ##
##
## Run after: hotspot_analysis.R
##
## Run before: map_creation_final.R
##
## list of final outputs:
##    'singleES_[es][cont].pdf' where:
##        'es' = the ES or BD (which is SppRich)
##        'cont' = the continent assessed"
##
## Author: Dr. Paul M. Evans
##
## Date Created: 2024-07-26
##
## Copyright (c) Paul M. Evans, 2024
## Email: paueva@ceh.ac.uk
## Git: https://github.com/pevans13

cat('\n\n####################################'
    , '\n\n uncert_graphs.R'
    , '\n\n####################################\n\n')

#### 0 - load functions ####
##### summarise the df #####
# create a function to load in all tables, and summarise that data, getting 
# mean and standard deviations for each group combinations.
summariseDfFunc <- function(inList # list of files to load in and combine
                            , coly = "Overlap" # column to analyse
                            , colSelt = c("es", "res", "continent", "perc", "pairs"
                                          , "Overlap", "one_hotpot", "no_hotspots")
){
  
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

##### create a function to plot boxplots for two percentiles #####
# and add the broken-stick analysis in addition
brokenBoxPlots <- function(dfx, coly, amt
                           , bdx = "species richness"
                           , yLab = paste0("Overlapping ", amt, " hotspots (pixels)")
                           , xxlab = "Models in ensemble"){
  # dfx <- results.combined[[1]]; xall = 9; x = "70th"; coly = "Overlap"; amt = "total"
  # dfx <- results.single2[[1]]; xall = 6; x = "90th";  coly = "median_HSs"; amt = "total"; bdx = "NA"
  dfx2 <- dfx %>%
    filter(!(es == "Carb" & pairs == 11)) %>%
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
    
    # Create a colour mapping for lines
    colorMapping <- c("70th" = "#999999"
                      , "80th" = "#E69F00"
                      , "90th" = "#56B4E9"
                      , "95th" = "#04AF70")
    
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
      gtit2 <- paste0(nms[3], " | ", cct, " & vascular plant richness")
    } else {
      gtit <- paste0("Overlapping ", amt, " hotspots for "
                     , nms[3], "\n", nms[1], " at ", nms[2] , " km2")
      gtit2 <- paste0(nms[3], " | ", nms[1])
      gtit2 <- paste0(nms[3], " | ", cct, " & vascular plant richness")
    }
    ## and y lab
    yllab <- yLab
    dfbox <- ggplot() + 
      geom_boxplot(data = dfxToGraph %>%
                     filter(!(es == "Carb" & pairs == 11))
                   , aes(x=as.factor(pairs), y=get(coly), fill = perc)
                   , position=position_dodge(1)
                   , notch=TRUE
                   , outliers = FALSE) +
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

#### 1 - load ####
##### combined dfs #####
# load in results that have 'END' in 
resultsList <- list.files(file.path("results", "tables")
                          , pattern = "END.csv"
                          , full.names = T)
## remove maps
resultsListCSV <- resultsList[!grepl("map", tolower(resultsList))]
### keep with '100'
resultsListCSV <- resultsListCSV[grepl("100", resultsListCSV)]
resultsListCSV

### get final dfs
## ------------ Notes --------------  ##
## summariseDfFunc produces a list 3 in length:
## 1: original data loaded in
## 2: mean and standard deviation of each combination
## 3: the groups that were analysed (consisting of ES, resolution, continent)
## ------------ ----- --------------  ##
results.combined <- summariseDfFunc(inList = resultsListCSV)
head(results.combined)

##### single-item dfs #####
# load in results that have 'END' in - for only ES or BD
resultsList2 <- list.files(file.path("results", "tables")
                           , pattern = "END.csv"
                           , recursive = T
                           , full.names = T)
## remove maps
resultsListCSV2 <- resultsList2[!grepl("map", tolower(resultsList2))]
## keep with '100'
resultsListCSV2 <- resultsListCSV2[grepl("100", resultsListCSV2)]
## keep only the ones
resultsListCSV2 <- resultsListCSV2[grepl("/es/|/bd/", resultsListCSV2, ignore.case = TRUE)]
## and original (i.e., non-'ens' END)
resultsListCSV_end <- resultsListCSV2[!grepl("_ens", resultsListCSV2, ignore.case = TRUE)]
## and 'ens' (i.e., 'ens' END)
resultsListCSV_ens <- resultsListCSV2[grepl("_ens", resultsListCSV2, ignore.case = TRUE)]
resultsListCSV_end; resultsListCSV_ens
### get final dfs
## ------------ Notes --------------  ##
## summariseDfFunc produces a list 3 in length:
## 1: original data loaded in
## 2: mean and standard deviation of each combination
## 3: the groups that were analysed (consisting of ES, resolution, continent)
## ------------ ----- --------------  ##
results.single1 <- summariseDfFunc(inList = resultsListCSV_end
                                   , colSelt = c("es", "res", "continent", "perc", "pairs"
                                                 , "Overlap"))
results.single2 <- summariseDfFunc(inList = resultsListCSV_ens
                                   , coly = "median_HSs" # column to analyse
                                   , colSelt = c("es", "res", "continent", "perc", "pairs"
                                                 , "median_HSs"))
### merge the results of 1 and 2
names(results.single1[[1]]); dim(results.single1[[1]])
names(results.single2[[1]]); dim(results.single2[[1]])
results.single <- bind_rows(results.single1[[1]]
                            , results.single2[[1]])

#### 2 - box plots - BD/ES ####
# for both ES and BD
## total overlap
completebox.combined <- brokenBoxPlots(dfx = results.combined[[1]]
                                       , coly = "Overlap"
                                       , amt = "total"
                                       )

### loop through each combination
for(i in seq_len(nrow(results.combined[[3]]))){
  
  esx <- results.combined[[3]][i, "es"] %>% as.character()
  contx <- results.combined[[3]][i, "continent"] %>% as.character()
  
  if(i == 1){
    ## open a pdf to save all graphs
    pdf(file.path("results/figures", "combined"
                  , paste0("combinedES", Sys.Date(), ".pdf")))
    cat("---------- PDF opened ----------\n")
    
    ### create a save list, to be matched with the graph outputs from 'map_creation_final.R'
    matchGraphList <- list()
    
    ### create a page of text explaining the results
    txt <- strwrap(
      paste(
        "Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13"
        , "\nProduced by 'uncert_graphs.R'"
        , paste("\nDate: ", Sys.Date(), sep = "")  # Adding the date
        , "\nThe diagrams in this PDF were generated by overlaying ensembles of ecosystem services (ES) and biodiversity (BD). Model outputs related to ES and BD across different continents were obtained, with ES outputs representing carbon storage (Carb), water supply (Wate), and recreation (Recr), while BD outputs always represented species richness (SppRich).

\nFor each ES or BD, ensembles were created by overlaying different combinations of model outputs and calculating the resulting median value for each pixel. For example, three different model outputs for the same ES (e.g., ES1_mod1, ES1_mod2, ES1_mod3) were overlaid to create an ensemble, and the ensemble (i.e., median) value was calculated. All possible combinations of model outputs were then used to generate the full set of ensembles. The number of model outputs used in an ensemble is referred to as 'esNo,' which ranged from 1 to the maximum number of model outputs available for that ES or BD measure. Each ensemble was normalized, resulting in values between 0 and 1. Then, ES and BD ensembles that used the same number of models (i.e., esNo) were overlaid to see where ES and BD hotspots overlapped. 

\n'Hotspots' were identified as areas (pixels) where an ensemble exceeded specific thresholds, with the 70th, 80th, and 90th percentiles used as cut-offs.

\nThe resulting graphs show the number of overlapping hotspots as esNo increases, highlighting how closely important ecosystem services and biodiversity align." 
        , sep = "\n")
      , width = 90)
    
    plot.new()
    # readme txt
    text(x = .1, y = .9
         , "README"
         , font=1, cex=2, col="#F48024")
    # info txt
    text(x = -0.05, y = .3 # first 2 numbers are xy-coordinates within [0, 1]
         , paste(txt, collapse = "\n")
         , adj = c(0.5, 1) # 'centred'
         , cex = .65
         , pos = 4)  
  }
  cat(i, esx, contx, "\n")
  
  # print to total pdf
  xga <- grid.arrange(completebox.combined[[i]])
  print(xga)
  
  # # print single pdf
  # png(file.path("results/figures", "combined"
  #               , paste0("combinedES", esx, "_", contx, ".png"))
  #     , width = 1500, height = 1000, res = 200)
  # grid.arrange(completebox.combined[[i]] 
  #              # + guides(fill = "none")
  #              # , completePerN.combined[[i]]
  # )
  # dev.off()
  # 
  
  ##### uncertainty quantification #####
  # Below, we test the variance of information within each category (i.e., per 'esNo')
  ## set df up
  analyseNow <- results.combined[[1]] %>%
    filter(es == esx
           , continent == contx
           , perc == 80) %>%
    mutate(esNo = as.factor(pairs))
  esNo <- unique(analyseNow$pairs) %>% as.character() # difference esNos
  
  # save the results, for easy navigation later
  matchGraphList[[i]] <- list(esx, contx
                              , xga)

  if(i == nrow(results.combined[[3]])){
  # if(i == 1){
    # dev.off()
    dev.off()
    # stop()
    cat("---------- PDF closed ----------\n")
    
    #### save matchGraphList
    save(matchGraphList, file = "matchGraphList.RData")
    cat("---------- matchGraphList saved ----------\n")
  }
}

#### 3 - box plots - one ####
## ------------ Steps for this part --------------  ##
## For each combination of ES and continent:
##  - combine all data output for single-item focus (i.e., either BD or an ES)
##  - reduce down to dataframe with just one of each element
##  - get medians, in order to be able to perform broken-stick analysis on them
##  - do broken-stick analysis on each percentile within the median df
##  - create a boxplot graph including brokwn-stick analysis
## ------------ ------------------- --------------  ##
# for either ES and BD
## total overlap - for ensemble output
completebox.singEns <- brokenBoxPlots(results.single2[[1]]
                                      , bdx = "NA"
                                      , coly = "median_HSs"
                                      , amt = "total"
                                      , yLab = "Hotspots amount (pixels)")

## total overlap - for comparison between models output
completebox.singComp <- brokenBoxPlots(results.single1[[1]]
                                       , bdx = "NA"
                                       , coly = "Overlap"
                                       , amt = "total")

for(i in seq_len(nrow(results.single2[[3]]))){
  
  esx <- results.single2[[3]][i, "es"] %>% as.character()
  contx <- results.single2[[3]][i, "continent"] %>% as.character()
  
  if(i == 1){
    ## open a pdf to save all graphs
    pdf(file.path("results/figures", "single_item"
                  , paste0("singleES", Sys.Date(), ".pdf")))
    cat("---------- PDF opened ----------\n")
    
    ### create a page of text explaining the results
    txt <- strwrap(
      paste(
        "Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13"
        , "\nProduced by 'uncert_graphs.R'"
        , paste("\nDate: ", Sys.Date(), sep = "")  # Adding the date
        , "\nThe diagrams in this PDF were generated by assessing ensembles of ecosystem services (ES) or biodiversity (BD). Model outputs related to ES and BD across different continents were obtained, with ES outputs representing carbon storage (Carb), water supply (Wate), and recreation (Recr), while BD outputs always represented species richness (SppRich).

\nFor each a single metric of ES or BD, ensembles were created by overlaying different combinations of model outputs and calculating the resulting median value for each pixel. For example, three different model outputs for the same ES (e.g., ES1_mod1, ES1_mod2, ES1_mod3) were overlaid to create an ensemble, and the ensemble (i.e., median) value was calculated. All possible combinations of model outputs were then used to generate the full set of ensembles. The number of model outputs used in an ensemble is referred to as 'esNo,' which ranged from 1 to the maximum number of model outputs available for that ES or BD measure. Each ensemble was normalized, resulting in values between 0 and 1.

\n'Hotspots' were identified as areas (pixels) where an ensemble equalled or exceeded specific thresholds, with the 70th, 80th, and 90th percentiles used as cut-offs.

\nThe resulting graphs illustrate the number of hotspots as the esNo increases, highlighting the stability of hotspot numbers when using different model combinations within the ensembles." 
        , sep = "\n")
      , width = 90)
    
    plot.new()
    # readme txt
    text(x = .1, y = .9
         , "README"
         , font=1, cex=2, col="#F48024")
    # info txt
    text(x = -0.05, y = .3 # first 2 numbers are xy-coordinates within [0, 1]
         , paste(txt, collapse = "\n")
         , adj = c(0.5, 1) # 'centred'
         , cex = .65
         , pos = 4)  
  }
  cat(i, esx, contx, "\n")
  
  # print to total pdf
  grid.arrange(completebox.singEns[[i]]
               # + guides(fill = "none")
               # , completePerN.single[[i]]
  )
  
  # print single pdf
  png(file.path("results/figures", "single_item"
                , paste0("singleES", esx, "_", contx, ".png"))
      , width = 1500, height = 1000, res = 200)
  grid.arrange(completebox.singEns[[i]]
               # + guides(fill = "none")
               # , completePerN.single[[i]]
  )
  dev.off()
  
  if(i == nrow(results.single2[[3]])){
    dev.off()
    cat("---------- PDF closed ----------\n")
  }
}

#### x - write readme ####
readmePath <- file.path("results/figures", "combined", "readme_combinedES.md")
# Create and open the file connection
fileConn <- file(readmePath)
writeLines(c(
  "Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13"
  , "Date created: 2024-08-20"
  , paste0("Last update :",  format(Sys.Date()))
  , "Produced by 'uncert_graphs.R'

Description of files in directory:
    The files in this directory contain graphs of hotspot comparisons between different ecosystem services (ES) and biodversity, and continents that were used as part of the uncertainty project.
    The graphs were created by assessing the 'hotspots' of different model outputs. Specifically, for the top graphs:
      - All ES or BD outputs were normalised, creating 0 - 1 values for each model output
        - The ES assessed were carbon (which had 11 model outputs), water supply (8 outputs), and recreation (4)
        - The BD, which was always species richness, had 11 model outputs
      - Ensembles were formed by combining all possible outputs for each ES and all possible combinations for BD. For instance, an ES with 8 model outputs had a total of 255 possible ensembles, created by combining 1, 2, 3, and so on, of the model outputs.
        - the number of model inputs that go into each ensemble were known as 'esNo'
      - The ensembles were created by using the median values for all input ensemble models for each pixel.
      - Ensembles for BD and a single ES were overlaid when the number of model inputs that went in to each ensemble was the same (e.g., all carbon ensembles produced from 3 carbon models outputs were compared with all 3 BD model output combinations)
        - The calculation for the final ES amount was: factorial(max model outputs for ES) / (factorial(i) * factorial(max model outputs for ES - i)) where current 'i' is the number of model in the ensembles (i.e., 'esNo')
        - The calculation for the final BD amount was: factorial(11 [max BD outputs]) / (factorial(i) * factorial(11 - i)) where current 'i' is the number of model in the ensembles (i.e., 'esNo')
        - for each esNo, these need to be mulitplied together to get all possible combinations comparisons
          - e.g., at 6 esNo for water and BD comparisons: factorial(8) / (factorial(6) * factorial(8 - 6)) * factorial(11) / (factorial(6) * factorial(11 - 6)) = 12,936 combinations. 
      - Different percentiles were chosen to indicate above which value would be considered a 'hotspot'. We chose 70th, 80th, and 90th percentiles.
      - The individual model output pixels were classified as being a hotspot (≥ percentile) or not 
      - All combinations of models were assessed for overlapping hotspot pixels, whereby the total number of hotspot pixels that always agreed were calculated
      - each set of comparisons are on the x-axis of each graph, titled '# esNo'
      - y-axis is 'Overlapping total hotspots (pixels)'
      
      - For the lower graphs on each page, a per-n was extracted for the y-axis. 
      - The 'per-n' value was calculated by first summing all the overlapping hotspot pixels within each unique combination of ES, continent, pairs, and percentage group. Then, this sum was divided by the number of combinations in that group.
      
      - For both graph, the bold lines indicated the point trajectory as determined by breakpoint analysis, while the vertical dashed lines indicate the location identified as most likely a breakpoint by the analysis.  
      
File names:
    'combined_[es][cont].pdf' where:
      'es' = the ES or BD (which is SppRich)
      'cont' = the continent assessed")
  , fileConn)
close(fileConn)

mx <- 11
for(i in 1:11){
  x1 <- factorial(mx) / (factorial(i) * factorial(mx - i))
  # calculate possible max amount
  fxx <- x1 
  cat(paste0("p:", i, " [mxC: ", fxx , "]"), "\n") # max combos
}

#### x - write readme ####
readmePath <- file.path("results/figures", "single_item", "readme_singleES.md")
# Create and open the file connection
fileConn <- file(readmePath)
writeLines(c(
  "Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13"
  , "Date created: 2024-08-02"
  , paste0("Last update: ",  format(Sys.Date()))
  , "Produced by 'uncert_graphs.R'

Description of files in directory:
    The files in this directory contain graphs of different ecosystem services (ES), biodversity (BD), and continents that were used as part of the uncertainty project.
    The graphs were created by assessing the 'hotspots' of different model outputs. Specifically, 
      - All ES or BD outputs were normalised, creating 0 - 1 values
      - Different percentiles were chosen to indicate above which value would be considered a 'hotspot'. We chose 70th, 80th, and 90th percentiles.
      - The individual model output pixels were classified as being a hotspot (≥ percentile) or not 
      - All combinations of models were assessed for overlapping hotspot pixels, whereby the total number of hotspot pixels that always agreed were calculated
      - The number of model comparisons depended on the number of model outputs, with all comparisons being used. 
        - For example, if an ES had 11 model outputs, 2,047 total comparisons were used: 11 1 x 1 comparisons; 55 2 x 2 comparisons; 165 3 x 3 comparisons; 330 4 x 4 comparisons ... 1 11 x 11 comparisons. 
      - each set of comparisons are on the x-axis of each graph, titled '# in summed output'
      - y-axis is 'Overlapping total hotspots (pixels)'
      
File names:
    'singleES_[es][cont].pdf' where:
      'es' = the ES or BD (which is SppRich)
      'cont' = the continent assessed")
  , fileConn)
close(fileConn)
