# *****************************************************************************
# R Script implementing conventional random f-fold cross-validation.  
# Related to the paper "Dealing with clustered samples for assessing map 
# accuracy by cross-validation".
# Contact: Sytze de Bruin, Wageningen University, Laboratory of Geo-information
# Science and Remote Sensing, email: sytze.debruin@wur.nl
# May 3, 2022
# *****************************************************************************

# ****** load required library *******
#.libPaths("/home/j/jlinnenb/r_packages/")
library(ranger)
library(CAST)
library(sf)
library(raster)
library(caret)
library(parallel)
library(mlr3)
library(mlr3learners)
source("./R/global_validation.R")

# ************ GLOBALS ***************
#setwd("/scratch/tmp/jlinnenb/deBruin_add_nndm/")

samples   <- c("simpleRandom")

infolder <- "./samples"
outfolder <- "./CVresults"
startseed <- 1234567
n_CV   <- 5  # number of cross validation replications
n_samp <- 3  # number of sample replicates (for each design)
cores <- 20

# create outfolders if they don't exist
if(!dir.exists(outfolder))
  dir.create(outfolder)

if(!dir.exists(paste0(outfolder, "/random_mlr")))
  dir.create(paste0(outfolder, "/random_mlr"))

# ************ FUNCTIONS ***************


knndmCV <- function(smpl, number, variate, seed){
  
  fname  <-  paste0(variate, "_", smpl, sprintf("%03d", number), ".Rdata")
  f_out  <- file.path(outfolder,"random_mlr", fname)
  
  if(!file.exists(f_out)) {
    fname <- paste0(variate, "data", sprintf("%03d", number), ".Rdata")
    f_in <- file.path(infolder,smpl,fname)
    load(f_in)
    
    # load ppoints
    load(file.path(infolder, "ppoints.Rdata"))
    
    MEC=RMSE=time=WS=time_mod <- numeric(n_CV)
    
    if(variate == "AGB"){
      pts_df <- data.frame(x=AGBdata$xcoord * 1000, y=AGBdata$ycoord * 1000)
      pts_sf <- st_as_sf(pts_df, coords = c("x","y"), crs = st_crs(ppoints))
    } else {
      pts_df <- data.frame(x=OCSdata$xcoord * 1000, y=OCSdata$ycoord * 1000)
      pts_sf <- st_as_sf(pts_df, coords = c("x","y"), crs = st_crs(ppoints))
    }
    
    n <- length(pts_df$x)
    
    for(i_CV in 1:n_CV){
      
      set.seed(seed)
      knndm <- knndm(pts_sf, ppoints = ppoints, k = 10, maxp = 0.8)
      WS[i_CV] <- knndm$W
      fold <- knndm$clusters
      
      set.seed(seed)
      if(variate == "AGB"){
        
        task_regr = as_task_regr(AGBdata, target = "agb")
        learner = lrn("regr.ranger",
                      mtry=floor(sqrt(ncol(AGBdata[,-1]))),
                      num.trees=500)
        
      } else{
        
        task_regr = as_task_regr(OCSdata, target = "ocs")
        learner = lrn("regr.ranger",
                      mtry=floor(sqrt(ncol(OCSdata[,-1]))))
      }
      
      custom_cv = rsmp("cv", folds=10)
      rr = resample(task_regr, learner, custom_cv,store_backends=TRUE)
      
      RMSE[i_CV] <- global_validation(rr)[[1]]
      
      seed <- seed + 1
    } # loop over i_CV
    
    save(RMSE, WS, file=f_out)
  } 
}


# ************ CALL THE FUNCTIONS ************ 
library(doParallel)
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)
pbapply::pblapply(seq(n_samp), function(i) {
  for(smpl in samples) {
    knndmCV(smpl, i, "AGB", startseed)
    knndmCV(smpl, i, "OCS", startseed)
  }
})
stopCluster(cl)
rm("cl")

mclapply(seq(n_samp), function(i) {
  for(smpl in samples) {
    knndmCV(smpl, i, "AGB", startseed)
    knndmCV(smpl, i, "OCS", startseed)
  }
}, mc.cores = cores)

