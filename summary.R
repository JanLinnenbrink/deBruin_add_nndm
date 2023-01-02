
## ---- echo=F, warning=F, message=F-----------------------------------------------------------------------------------------------------------------------------------------
library(sf)
library(raster)
library(ggplot2)
library(ggpubr)
library(openxlsx)



## ---- echo=F---------------------------------------------------------------------------------------------------------------------------------------------------------------

agb <- raster::raster("./data/agb_resampled.tif")
sample_plots <- lapply(list.files("./samples"), function(dir) {
  r <- floor(runif(1, 1, 100))
  # print(paste0(dir, ": plotting the ", r, "th sampling design."))
  samp_name <- paste0("./samples/", dir, "/", sprintf("%03d", r), "_coords", ".Rdata")
  load(samp_name)
  if(class(pts)[1] != "sf") {
    sample_df <- as.data.frame(pts)
    sample_sf <- st_as_sf(sample_df, coords = c("x", "y"))
    st_crs(sample_sf) <- st_crs(agb) # set crs EPSG:3035
  } else {
    sample_sf <- pts
  }
  
  agb_df <- as.data.frame(agb, xy=TRUE)
  
  m1 <- ggplot() + 
    geom_raster(data=agb_df, aes(x=x, y=y, fill=agb_resampled)) + 
    geom_sf(data=sample_sf, shape = 2, size = 0.1) + 
    theme_light() +
    scale_fill_gradientn(colours=terrain.colors(10), na.value="transparent") + #, na.value = "transparent") +
    ggtitle(dir)
  
  return(m1)
}) 


## ---- eval=T---------------------------------------------------------------------------------------------------------------------------------------------------------------


infolder  <- "./CVresults"
outfolder <- "./material"

mets <- c("exhaustive", "random", "spatial", "intensity",
          "modelbased", "heteroscedastic", "nndm", "knndm", "bnndm", "bnndm_10", "bnndm_10_wo")

colnms <- c("method", "variate", "design", "number", "RMSE", "MEC")
outtab <- data.frame(matrix(NA, 0, 10))
names(outtab) <- colnms

for(m in mets){
  p <- file.path(infolder, m)
  f_ins <- list.files(p, glob2rx("???_*.Rdata"))
  for(f_in in f_ins){
    lchar <- nchar(f_in)
    variate <- substr(f_in, 1, 3)
    design <- substr(f_in, 5, lchar-9)
    number <- as.numeric(substr(f_in, lchar-8, lchar-6))
    load(file.path(p, f_in))
    if(m == "modelbased" | m == "heteroscedastic"){
      MEC    <- mean(MECs)
      RMSE   <- mean(RMSEs)
    } else{
      if(length(MEC) > 1){
        MEC  <- mean(MEC)
        RMSE <- mean(RMSE)
      }
    }
    
    newrow <- data.frame(method = m, variate = variate, design = design,
                         number = number, RMSE = RMSE, MEC = MEC)
    outtab <- rbind(outtab, newrow)
  }
}

write.xlsx(outtab, file.path(outfolder, "outtab100.xlsx"), overwrite = T)


outtab <- read.xlsx(file.path(outfolder, "outtab100.xlsx"))
outtab$methodID <- with(outtab, ifelse(method == "random", 3, 
                                       ifelse(method == "spatial", 8,
                                              ifelse(method == "intensity", 13,
                                                     ifelse(method == "modelbased", 18,
                                                            ifelse(method == "heteroscedastic", 23,
                                                                   ifelse(method == "nndm", 28,
                                                                          ifelse(method == "knndm", 33,
                                                                                 ifelse(method == "bnndm", 38,
                                                                                        ifelse(method=="bnndm_10",43,
                                                                                               ifelse(method=="bnndm_10_wo",48,
                                                                                        0)))))))))))

# relative RMSE & MEC
outtab$rRMSE <- NA
outtab$rMEC  <- NA

# some rows missing in exhaustive method for OCS data
numbers <- outtab[outtab$variate == "OCS" & outtab$method=="exhaustive" & outtab$design=="clusterGapped",]$number
outtab <- outtab[outtab$number %in% numbers, ]


for(variate in c("AGB", "OCS")){
  for(design in unique(outtab$design)){
    for(number in numbers){
      idx1 <- which(outtab$design == design & outtab$variate == variate & 
                      outtab$number == number)
      idx2 <- which(outtab$design == design & outtab$variate == variate & 
                      outtab$methodID == 0 & outtab$number == number)
      
      outtab$rRMSE[idx1] <- 100 * (outtab$RMSE[idx1] - outtab$RMSE[idx2])/
        outtab$RMSE[idx2]
      outtab$rMEC[idx1] <- 100 * (outtab$MEC[idx1] - outtab$MEC[idx2])/
        outtab$MEC[idx2]
    }
  }
}

outtab$method <- factor(outtab$method, levels=c("exhaustive","random", "spatial", "intensity", "modelbased", "heteroscedastic", "nndm", "knndm", "bnndm","bnndm_10", "bnndm_10_wo"))
outtab$design <- factor(outtab$design, levels=c("simpleRandom", "regular", "clusterMedium", "clusterStrong", "clusterGapped"))
xlabs <- c("conventional", "spatial", "weighted", "hom.sced.","het.sced.","NNDM","kNNDM","bNNDM")
collabs <- c("SRS","syst","clustMed","clustStr","clustGap")
cols <- c("white", "brown", "pink", "orange", "lightgreen")


(rmse <- ggplot(outtab[outtab$method!="exhaustive"&outtab$method!="bnndm"&outtab$method!="bnndm_10",], 
               aes(x=method, y=rRMSE, fill=design)) +
  geom_boxplot(linetype = "dashed", outlier.shape = NA) +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..)) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..)) +
    stat_summary(fun.y=mean, geom="point", shape=20, size=3, colour="orange",
                 position = position_dodge2(width = 0.75,   
                                            preserve = "single")) +
  geom_hline(yintercept = 0, linetype="solid", alpha=0.6) +
  scale_fill_manual(values=cols, labels=collabs) +
  scale_y_continuous(limits=c(-60,60)) +
  scale_x_discrete(labels=xlabs) +
  ylab("relative RMSE [%]") + 
  xlab("") +
  theme_classic() +
  theme(strip.text=element_text(hjust = 0.55,  # hjust = 0.5 centers the title
                                size = 14,
                                face = "bold"),
        panel.border = element_rect(linetype = "solid",
                                    colour = "black", fill = "NA", size = 0.5),
        legend.position = NaN) +
  facet_wrap(~variate)  )

for(groups in unique(outtab$method)) {
  print(median(outtab[outtab$method==groups,]$rMEC))
}

mec <- ggplot(outtab[outtab$method!="exhaustive"&outtab$method!="bnndm"&outtab$method!="bnndm_10",], aes(x=method, y=rMEC, fill=design)) +
  geom_boxplot(linetype = "dashed", outlier.shape = NA) +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..)) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..)) +
  stat_summary(fun.y=mean, geom="point", shape=20, size=3, colour="orange",
               position = position_dodge2(width = 0.75,   
                                          preserve = "single")) +
  geom_vline(xintercept=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5),  color="grey60", linetype="dashed", alpha=0.6) +
  geom_hline(yintercept = 0, linetype="solid", alpha=0.6) +
  scale_fill_manual(values=cols, labels=collabs) +
  scale_y_continuous(limits=c(-40,40)) +
  scale_x_discrete(labels=xlabs) + 
  xlab("") +
  ylab("relative MEC [%]") + 
  theme_classic() +
  theme(strip.text=element_text(hjust = 0.5,  # hjust = 0.5 centers the title
                                size = 14,
                                face = "bold"),
        panel.border = element_rect(linetype = "solid",
                                    colour = "black", fill = "NA", size = 0.5),legend.box.just = "left",
        legend.justification = c(1,0),
        legend.position = c(1,0),
        legend.text.align = 0,
        legend.title = element_blank(),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
  facet_wrap(~variate)

bpl <- gridExtra::grid.arrange(rmse, mec)

ggsave(file.path(outfolder, "plot_both6.pdf"), bpl, height = unit(8, "cm"), width=unit(12, "cm"))

## ---- eval=T---------------------------------------------------------------------------------------------------------------------------------------------------------------
mytab <- outtab[outtab$variate == "OCS",]
# mytab <- mytab[mytab$number < 6,]
mytab$rRMSE <- NA
mytab$rMEC  <- NA

for (design in unique(mytab$design)) {
  for (number in c(1:85, 94:100)) {
    measurement <- which(mytab$design == design & mytab$number == number)
    validation <- which(mytab$design == design & mytab$number == number & mytab$method == "exhaustive")
    mytab$rRMSE[measurement] <- 100 * (mytab$RMSE[measurement] - mytab$RMSE[validation])/
        mytab$RMSE[validation]
    mytab$rMEC[measurement] <- 100 * (mytab$MEC[measurement] - mytab$MEC[validation])/
        mytab$MEC[validation]
  }
}

# boxplot(rRMSE~method, data=mytab)

mytab$design[mytab$design == "clusterGapped"] <- "e_clusterGapped"
mytab$design[mytab$design == "simpleRandom"] <- "a_simpleRandom"
mytab$design[mytab$design == "regular"] <- "b_regular"
mytab$design[mytab$design == "clusterMedium"] <- "c_clusterMedium"
mytab$design[mytab$design == "clusterStrong"] <- "d_clusterStrong"

mytab$method[mytab$method == "exhaustive"] <- "a_exhaustive"
mytab$method[mytab$method == "random"] <- "b_random"
mytab$method[mytab$method == "spatial"] <- "c_spatial"
mytab$method[mytab$method == "intensity"] <- "d_intensity"
mytab$method[mytab$method == "modelbased"] <- "e_modelbased"
mytab$method[mytab$method == "heteroscedastic"] <- "f_heteroscedastic"
mytab$method[mytab$method == "nndm"] <- "g_nndm"
mytab$method[mytab$method == "knndm"] <- "h_knndm"



levels(as.factor(mytab$design))

ggplot(data=mytab[mytab$method != "a_exhaustive",]) +
  geom_boxplot(aes(x=method, y=rRMSE, color=design)) +
  ggtitle("OCS: relative RMSE (%) by CV Method")

