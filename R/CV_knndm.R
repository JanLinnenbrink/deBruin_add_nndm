# *****************************************************************************
# R Script implementing conventional random f-fold cross-validation.  
# Related to the paper "Dealing with clustered samples for assessing map 
# accuracy by cross-validation".
# Contact: Sytze de Bruin, Wageningen University, Laboratory of Geo-information
# Science and Remote Sensing, email: sytze.debruin@wur.nl
# May 3, 2022
# *****************************************************************************

# sleep for random time so multiple jobs dont take the same job out of the csv
# Sys.sleep(round(runif(1, min = 1, max = 240)))

# ****** load required library *******
.libPaths("/home/j/jlinnenb/r_packages/")
library(ranger)
library(CAST)
library(sf)
library(raster)
library(caret)
library(parallel)

# ************ GLOBALS ***************
infolder <- "/scratch/tmp/jlinnenb/deBruin_add_nndm/samples"
outfolder <- "/scratch/tmp/jlinnenb/deBruin_add_nndm/CVresults/knndm"
datafolder <- "/scratch/tmp/jlinnenb/deBruin_add_nndm/data"


n <- 700 # usually 5000
n_CV <- 3

# create outfolders if they don't exist
if(!dir.exists(outfolder))
  dir.create(outfolder, recursive = TRUE)

# ************ FUNCTIONS ***************

sumSquares <- function(ref, pred){
  muref <- mean(ref, na.rm=T)
  SSR <- sum((ref - pred)^2)
  SST <- sum((ref - muref)^2)
  return(c(SSR, SST))
}

err_fu <- function(obs, pred){
  rmse <- sqrt(mean((obs-pred)^2))
  muref <- mean(obs)
  SSR <- sum((obs - pred)^2)
  SST <- sum((obs - muref)^2)
  mec <- 1 - SSR/SST
  me <- mean(obs - pred)
  list(me = me, rmse = rmse, mec = mec)
}

nndmCV <- function(smpl, number, variate) {
  
  fname <- paste0(variate, "data", sprintf("%03d", number), ".Rdata")
  f_in <- file.path(infolder,smpl,fname)
  load(f_in)
  
  RMSE <- numeric(n_CV)
  
  # load sample file containing coordinates
  load(file.path(infolder,smpl,paste0(sprintf("%03d", number), "_coords", ".Rdata")))
  if(class(pts)[1] != "sf") {
    sample_df <- as.data.frame(pts)
    sample_sf <- st_as_sf(sample_df, coords = c("x", "y"))
  } else {
    sample_sf <- pts
  }
  
  agb_raster <- raster::raster(file.path(datafolder, "agb_resampled.tif")) # load agb raster
  
  folds <- list()
  
  for(i_CV in 1:n_CV) {
    
    sample_subset <- st_cast(st_sample(sample_sf, n), to = "POINT") # subset sample
    st_crs(sample_subset) <- st_crs(agb_raster) # set crs
    
    #st_as_sf(raster::rasterToPoints(agb_raster[[1]], spatial = TRUE))
    knndm <- knndm(sample_subset, agb_raster)
    
    # Evaluate RF model using NDM CV
    trainControl_knndm <- trainControl(method = "cv",
                                       index=knndm$indx_train,
                                       indexOut=knndm$indx_test,
                                       savePredictions = "final") # save predictions final to avoid writing CV myself
    
    paramGrid <-  data.frame(mtry = 2, min.node.size = 5, splitrule = "variance")
    
    if (variate == "AGB") {
      training_data <- AGBdata; rf_form <- agb~.
    } else {
      training_data <- OCSdata; rf_form <- ocs~.
    }
    
    # next: execute this and see where it goes from there
    mod_knndm <- train(rf_form,
                       method = "ranger",
                       trControl = trainControl_knndm,
                       tuneGrid = paramGrid, 
                       data = training_data)
    
    RMSE[i_CV] <- global_validation(mod_knndm)["RMSE"][[1]]
    folds <- append(folds, knndm)
    
  }
  
  fname  <-  paste0(variate, "_", smpl, sprintf("%03d", number), ".Rdata")
  f_out  <- file.path(outfolder,fname)
  save(RMSE, folds, file=f_out)
}

# ************ CALL THE FUNCTIONS ************ 

samples <- c("clusterMedium", "clusterStrong", "clusterGapped", "regular", 
             "simpleRandom")
n_samp <- 100
cores <- 20

mclapply(seq(n_samp), function(i) {
  for(smpl in samples) {
    nndmCV(smpl = smpl, number = i, variate = "AGB")
    # nndmCV(smpl, i, "OCS")
  }
}, mc.cores = cores)
