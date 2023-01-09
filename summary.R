
## ---- echo=F, warning=F, message=F-----------------------------------------------------------------------------------------------------------------------------------------
library(sf)
library(raster)
library(ggplot2)
library(ggpubr)
library(openxlsx)

setwd("C:/0_Msc_Loek/Z_Palma/deBruin_add_nndm/")

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
          "modelbased", "heteroscedastic", "nndm_test", "knndm","bnndm_5", "bnndm_def")

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
    
    if(m %in% c("nndm_test", "knndm", "bnndm_5")) {
      time = time
      time_mod = time_mod
      WS = WS
    } else WS = 0
    
    newrow <- data.frame(method = m, variate = variate, design = design,
                         number = number, RMSE = RMSE, MEC = MEC, WS=WS)
    outtab <- rbind(outtab, newrow)
  }
}

write.xlsx(outtab, file.path(outfolder, "outtab100.xlsx"), overwrite = T)


outtab <- read.xlsx(file.path(outfolder, "outtab100.xlsx"))
outtab$methodID <- with(outtab, ifelse(method=="exhaustive",0,1))
                                       
                                       
# relative RMSE & MEC
outtab$rRMSE <- NA
outtab$rMEC  <- NA

# some rows missing in exhaustive method for OCS data
numbers <- outtab[outtab$variate == "OCS" & outtab$method=="exhaustive" & outtab$design=="clusterGapped",]$number
outtab <- outtab[outtab$number %in% numbers, ]

n <- 1:100

n[!n %in% numbers]
86
88
93

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

outtab$method <- factor(outtab$method, levels=c("exhaustive","random", "spatial", "intensity", "modelbased", "heteroscedastic", 
                                                "nndm_test", "knndm", "bnndm_5", "bnndm_def"))
outtab$design <- factor(outtab$design, levels=c("simpleRandom", "regular", "clusterMedium", "clusterStrong", "clusterGapped"))
outtab$variate <- as.factor(outtab$variate)
xlabs <- c("conventional", "spatial", "weighted", "hom.sced.","het.sced.","NNDM","kNNDM","bNNDM","bNNDM_def")
collabs <- c("SRS","syst","clustMed","clustStr","clustGap")
cols <- c("white", "brown", "pink", "orange", "lightgreen")


(rmse <- ggplot(outtab[outtab$methodID!=0 & outtab$method != "bnndm_def",], aes(x=method, y=rRMSE, fill=design)) +
  geom_boxplot(linetype = "dashed", outlier.shape = NA) +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..)) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..)) +
  stat_summary(fun.y=mean, geom="point", shape=20, size=2, colour="black",
               position = position_dodge2(width = 0.75,   
                                          preserve = "single")) +
    geom_vline(xintercept=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5),  color="grey60", linetype="dashed", alpha=0.6) +
  geom_hline(yintercept = 0, linetype="solid", alpha=0.6) +
  scale_fill_manual(values=cols, labels=collabs) +
  scale_x_discrete(labels=xlabs) + 
  scale_y_continuous(limits=c(-60,60), breaks=seq(-60,60,20)) +
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

(mec <- ggplot(outtab[outtab$methodID!=0 & outtab$method != "bnndm_def",], aes(x=method, y=rMEC, fill=design)) +
  geom_boxplot(linetype = "dashed", outlier.shape = NA) +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..)) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..)) +
  stat_summary(fun.y=mean, geom="point", shape=20, size=2, colour="black",
               position = position_dodge2(width = 0.75,   
                                          preserve = "single")) +
  geom_vline(xintercept=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5),  color="grey60", linetype="dashed", alpha=0.6) +
  geom_hline(yintercept = 0, linetype="solid", alpha=0.6) +
  scale_fill_manual(values=cols, labels=collabs) +
  scale_y_continuous(limits=c(-40,40), breaks=seq(-40,40,20)) +
  scale_x_discrete(labels=xlabs) + 
  xlab("") +
  ylab("relative MEC [%]") + 
  theme_classic() +
  theme(strip.text=element_text(hjust = 0.5,  # hjust = 0.5 centers the title
                                size = 14,
                                face = "bold"),
        panel.border = element_rect(linetype = "solid",
                                    colour = "black", fill = "NA", size = 0.5),legend.box.just = "left",
        legend.justification = c(0,0),
        legend.position = c(0,0),
        legend.text.align = 0,
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
  facet_wrap(~variate))

bpl <- gridExtra::grid.arrange(rmse, mec)

ggsave(file.path(outfolder, "plot_try2.pdf"), bpl, height = unit(8, "cm"), width=unit(12, "cm"))

(pl <- ggplot(data=outtab[outtab$method %in% c("random", "spatial", "knndm", "bnndm_5"),], aes(fill=method)) +
    geom_density(aes(y=rRMSE), alpha=0.5) +
    geom_hline(yintercept = 0) +
    scale_y_continuous(limits=c(-60,60)) +
    xlab("") +
    ylab("relative RMSE [%]") + 
    theme_bw() +
    facet_wrap(~design))
ggsave(paste0(outfolder, "rrmse_density.pdf"),pl)

(pl_WS <- ggplot(outtab[outtab$method%in%c("knndm","bnndm_5"),], aes(x=abs(rRMSE),y=WS,color=design)) +
    geom_point() +
    facet_wrap(~method))
(pl_WS <- ggplot(outtab, aes(x=method,y=WS,fill=design)) +
    geom_boxplot(linetype = "dashed", outlier.shape = NA) +
    stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = NA) +
    stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..)) +
    stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..)) +
    stat_summary(fun.y=mean, geom="point", shape=20, size=2, colour="black",
                 position = position_dodge2(width = 0.75,   
                                            preserve = "single")) +
    geom_vline(xintercept=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5),  color="grey60", linetype="dashed", alpha=0.6) +
    geom_hline(yintercept = 0, linetype="solid", alpha=0.6) +
    scale_fill_manual(values=cols, labels=collabs) +
    facet_wrap(~variate))
ggsave(paste0(outfolder, "WS.pdf"),pl_WS)


# time
infolder  <- "./CVresults/"
outfolder <- "./material/"
ktab <- data.frame(matrix(NA, 0, 10))

mets <- c("nndm_test", "knndm","bnndm_5")

colnms <- c("method", "variate", "design", "WS",  "time")
ktab <- data.frame(matrix(NA, 0, 11))
names(ktab) <- colnms
for(m in mets){
  p <- file.path(infolder, m)
  f_ins <- list.files(p, glob2rx("???_*.Rdata"))
  for(f_in in f_ins){
    lchar <- nchar(f_in)
    variate <- substr(f_in, 1, 3)
    design <- substr(f_in, 5, lchar-9)
    number <- as.numeric(substr(f_in, lchar-8, lchar-6))
    load(file.path(p, f_in))
    time <- time
    time_mod <- time_mod
    
    newrow <- data.frame(method = m, variate = variate, design = design, time=time, time_mod=time_mod)
    ktab <- rbind(ktab, newrow)
  }
}


ktab$design <- factor(ktab$design, levels=c("simpleRandom", "regular", "clusterMedium", "clusterStrong", "clusterGapped"))
ktab$variate <- as.factor(ktab$variate)
ktab$time_f <- ktab$time
ktab$time <- ktab$time_f + ktab$time_mod
collabs <- c("SRS","syst","clustMed","clustStr","clustGap")
cols <- c("white", "brown", "pink", "orange", "lightgreen")


(pl <- ggplot(data=ktab, aes(x=method, y=time, fill=design)) +
    geom_violin() +
    scale_y_log10() +
    ylab("time [s]"))


(pl_ws <- ggplot(data=ktab, aes(fill=variate, alpha=0.5)) +
    geom_density(aes(x=WS)) +
    geom_jitter(data = ktab, aes(x=WS,y=0.000001),size=0.4, height=0.0000001, alpha=0.5) +
    theme_bw() +
    facet_wrap(~design, scales="free"))
ggsave(paste0(outfolder, "WS_design.pdf"),pl_ws)

(pl_rmse <- ggplot(data=outtab[outtab$method=="knndm",], aes(fill=variate, alpha=0.5)) +
    geom_density(aes(x=rRMSE)) +
    geom_jitter(data = outtab, aes(x=rRMSE,y=0.000),size=0.4, height=0.0001, alpha=0.5) +
    geom_vline(xintercept = 0)+
    theme_bw() +
    facet_wrap(~design, scales="free_y"))
ggsave(paste0(outfolder, "rmse_design.pdf"),pl_rmse)


(pl <- ggplot(data=ktab, aes(fill=method)) +
    geom_density(aes(x=time)) +
    geom_jitter(aes(x=time,y=0.0001),size=0.4, height=0.0001, alpha=0.5) +
    theme_bw()+
    facet_wrap(~design, scales="free_y"))
ggsave(paste0(outfolder, "knndm_time.pdf"),pl)

(pl <- ggplot(data=outtab[outtab$method %in% c("random","knndm", "bnndm_5"),], aes(x=method, y=rRMSE, fill=design)) +
    geom_violin(linetype = "dashed", outlier.shape = NA) +
    geom_point(position=position_jitterdodge(dodge.width = 0.9), shape=1) +
    geom_vline(xintercept=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5),  color="grey60", linetype="dashed", alpha=0.6) +
    geom_hline(yintercept = 0, linetype="solid", alpha=0.6) +
    scale_fill_manual(values=cols, labels=collabs) +
    scale_y_continuous(limits=c(-60,60)) +
    ylab("relative RMSE [%]") + 
    xlab("") +
    theme_classic() +
    theme(strip.text=element_text(hjust = 0.55,  # hjust = 0.5 centers the title
                                  size = 14,
                                  face = "bold"),
          panel.border = element_rect(linetype = "solid",
                                      colour = "black", fill = "NA", size = 0.5),
          legend.position = c(0.05,0.15)) +
    facet_wrap(~variate))

ggsave(paste0(outfolder, "/rmse_jitter.pdf"),pl)


