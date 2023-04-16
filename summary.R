library(ggplot2)
infolder  <- "./CVresults/"
outfolder <- "./material"
mets <- c("exhaustive_test", "random_caret_new", "spatial", "intensity","knndm_caret_new")
colnms <- c("method", "variate", "design", "number", "RMSE")
outtab <- data.frame(matrix(NA, 0, 5))
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
      RMSE   <- mean(RMSE)
    
    newrow <- data.frame(method = m, variate = variate, design = design,
                         number = number, RMSE = RMSE)
    outtab <- rbind(outtab, newrow)
  }
}

write.csv(outtab, file.path(outfolder, "outtab100.csv"))


outtab <- read.csv(file.path(outfolder, "outtab100.csv"))
outtab$methodID <- with(outtab, ifelse(method=="exhaustive_test",0,1))
                                       
                                       
# relative RMSE & MEC
outtab$rRMSE <- NA
#outtab$rMEC  <- NA
library(dplyr)
outtab <- outtab |> 
  group_by(method, design, variate) |> 
  slice(1:50)


numbers=19:50
#18
 

for(variate in c("AGB", "OCS")){
  for(design in unique(outtab$design)){
    for(number in numbers){
      idx1 <- which(outtab$design == design & outtab$variate == variate & 
                      outtab$number == number)
      idx2 <- which(outtab$design == design & outtab$variate == variate & 
                      outtab$methodID == 0 & outtab$number == number)
      
      outtab$rRMSE[idx1] <- 100 * (outtab$RMSE[idx1] - outtab$RMSE[idx2])/
        outtab$RMSE[idx2]
      #outtab$rMEC[idx1] <- 100 * (outtab$MEC[idx1] - outtab$MEC[idx2])/
      #  outtab$MEC[idx2]
    }
  }
}

outtab$method <- factor(outtab$method, levels=c("exhaustive_test", "random_caret_new", "spatial", "intensity","knndm_caret_new"))
outtab$design <- factor(outtab$design, levels=c("simpleRandom", "regular", "clusterMedium", 
                                                "clusterStrong", "clusterGapped"))
outtab$variate <- as.factor(outtab$variate)
xlabs <- c("conventional", "spatial", "weighted", "hom.sced.","het.sced.","NNDM","knndm_test")
collabs <- c("SRS","syst","clustMed","clustStr","clustGap")
cols <- c("white", "brown", "pink", "orange", "lightgreen")
outtab[!outtab$method%in%c("exhaustive_test")&outtab$variate=="OCS",]

(rmse <- ggplot(outtab[!outtab$method%in%c("exhaustive_test"),], 
                aes(x=method, y=rRMSE, fill=design)) +
  geom_boxplot(linetype = "dashed", outlier.shape = NA) +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..)) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..)) +
  stat_summary(fun.y=mean, geom="point", shape=1, size=1.5, stroke = 1, colour="black",
               position = position_dodge2(width = 0.75,   
                                          preserve = "single")) +
    geom_vline(xintercept=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5),  color="grey60", linetype="dashed", alpha=0.6) +
  geom_hline(yintercept = 0, linetype="solid", alpha=0.6) +
  scale_fill_manual(values=cols, labels=collabs) +
 #   scale_x_discrete(labels=xlabs) +
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


(mec <- ggplot(outtab[outtab$methodID!=0,], aes(x=method, y=rMEC, fill=design)) +
  geom_boxplot(linetype = "dashed", outlier.shape = NA) +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..)) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..)) +
  stat_summary(fun.y=mean, geom="point", shape=1, size=1.5, colour="black",
               position = position_dodge2(width = 0.75,   
                                          preserve = "single")) +
  geom_vline(xintercept=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5),  color="grey60", linetype="dashed", alpha=0.6) +
  geom_hline(yintercept = 0, linetype="solid", alpha=0.6) +
  scale_fill_manual(values=cols, labels=collabs) +
    scale_x_discrete(labels=xlabs) +
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

ggsave(file.path(outfolder, "comp_global.pdf"), bpl, height = unit(8, "cm"), width=unit(12, "cm"))


## density
# global validation results in under/overestimation of map error in random - regualar cases
(pl <- ggplot(data=outtab[outtab$method %in% c("random", "spatial", "knndm_sample"),], aes(fill=method)) +
    geom_density(aes(y=rRMSE), alpha=0.5) +
    geom_hline(yintercept = 0) +
    scale_y_continuous(limits=c(-60,60)) +
    xlab("") +
    ylab("relative RMSE [%]") + 
    theme_bw() +
    facet_wrap(~design))
ggsave(paste0(outfolder, "rrmse_density10.pdf"),pl)


# WS
(pl_WS <- ggplot(outtab[outtab$method=="knndm",], aes(x=abs(rRMSE),y=WS,color=design)) +
    geom_point() +
    facet_wrap(~method))
(pl_WS <- ggplot(outtab[outtab$method%in%c("knndm","pnndm_10", "nndm"),], aes(x=method,y=WS,fill=design)) +
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


#### test #####
library(ggplot2)
infolder  <- "./CVresults"
outfolder <- "./material"
mets <- c("exhaustive_test","knndm_sample_global_fold2indx_test")
colnms <- c("method", "variate", "design", "number", "RMSE_rand_caret","RMSE_rand_manual",
            "RMSE_indx","RMSE_kndm","RMSE_manual", "RMSE_mlr","RMSE_mlr_global")
outtab <- data.frame(matrix(NA, 0, length(colnms)))
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
    
    if (m == "knndm_sample_global_fold2indx_test") {
      newrow <- data.frame(method = m, variate = variate, design = design,
                           number = number, RMSE_rand_caret = RMSE_rand_caret, RMSE_rand_manual=RMSE_rand_manual,
                           RMSE_mlr=RMSE_mlr, RMSE_mlr_global = RMSE_mlr_global,
                           RMSE_kndm_f2id=RMSE_indx,RMSE_kndm_caret=RMSE_kndm, RMSE_kndm_manual=RMSE_manual)
    } else {
      newrow <- data.frame(method = m, variate = variate, design = design,
                           number = number, RMSE_rand_caret = RMSE, RMSE_rand_manual=RMSE,
                           RMSE_mlr=RMSE, RMSE_mlr_global = RMSE,
                           RMSE_kndm_f2id=RMSE,RMSE_kndm_caret=RMSE, RMSE_kndm_manual=RMSE)
    }
    
    outtab <- rbind(outtab, newrow)
  }
}


outtab$methodID <- with(outtab, ifelse(method=="exhaustive_test",0,1))
outtab$rRMSE_rand_caret=outtab$rRMSE_rand_manual=outtab$rRMSE_kndm_f2id=outtab$rRMSE_kndm_caret=outtab$rRMSE_kndm_manual  <- NA
outtab$rRMSE_kndm_mlr=outtab$rRMSE_kndm_mlr_global <- NA

# some rows missing in exhaustive method for OCS data
numbers <- outtab[outtab$variate == "OCS" & outtab$method=="exhaustive_test" & outtab$design=="clusterGapped",]$number
outtab <- outtab[outtab$number %in% numbers, ]


for(variate in c("AGB", "OCS")){
  for(design in unique(outtab$design)){
    for(number in numbers){
      idx1 <- which(outtab$design == design & outtab$variate == variate & 
                      outtab$number == number)
      idx2 <- which(outtab$design == design & outtab$variate == variate & 
                      outtab$methodID == 0 & outtab$number == number)
      
      outtab$rRMSE_rand_caret[idx1] <- 100 * (outtab$RMSE_rand_caret[idx1] - outtab$RMSE_rand_caret[idx2])/
        outtab$RMSE_rand_caret[idx2]
      outtab$rRMSE_rand_manual[idx1] <- 100 * (outtab$RMSE_rand_manual[idx1] - outtab$RMSE_rand_manual[idx2])/
        outtab$RMSE_rand_manual[idx2]
      outtab$rRMSE_kndm_f2id[idx1] <- 100 * (outtab$RMSE_kndm_f2id[idx1] - outtab$RMSE_kndm_f2id[idx2])/
        outtab$RMSE_kndm_f2id[idx2]
      outtab$rRMSE_kndm_caret[idx1] <- 100 * (outtab$RMSE_kndm_caret[idx1] - outtab$RMSE_kndm_caret[idx2])/
        outtab$RMSE_kndm_caret[idx2]
      outtab$rRMSE_kndm_manual[idx1] <- 100 * (outtab$RMSE_kndm_manual[idx1] - outtab$RMSE_kndm_manual[idx2])/
        outtab$RMSE_kndm_manual[idx2]
      outtab$rRMSE_kndm_mlr[idx1] <- 100 * (outtab$RMSE_mlr[idx1] - outtab$RMSE_mlr[idx2])/
        outtab$RMSE_mlr[idx2]
      outtab$rRMSE_kndm_mlr_global[idx1] <- 100 * (outtab$RMSE_mlr_global[idx1] - outtab$RMSE_mlr_global[idx2])/
        outtab$RMSE_mlr_global[idx2]
    }
  }
}

#outtab$method <- factor(outtab$method, levels=c("exhaustive","kndm"))
#outtab$design <- factor(outtab$design, levels=c("exhaustive","knndm"))
outtab$variate <- as.factor(outtab$variate)
#xlabs <- c("RMSE_kndm", "RMSE_indx")
collabs <- c("SRS","syst","clustMed","clustStr","clustGap")
cols <- c("white", "brown", "pink", "orange", "lightgreen")



library(dplyr)
library(tidyverse)
out_long <- outtab |>
  dplyr::filter(method!="exhaustive") |> 
  dplyr::select(!c(method, number, methodID)) |>
  dplyr::select(!starts_with("RMSE")) |> 
  pivot_longer(!c(variate,design))
out_long$design <- fct_relevel(out_long$design, c("simpleRandom", "regular", "clusterMedium", "clusterStrong", "clusterGapped"))

write.csv(out_long, file.path(outfolder, "outtab10_test.csv"))


(gp <- ggplot(out_long[out_long$name %in% c("rRMSE_rand_caret", "rRMSE_rand_manual"),],
              aes(x=design, fill=name)) +
  geom_boxplot(aes(y=value, alpha=0.9)) +
  facet_wrap(~variate) +
  geom_hline(yintercept = 0))

ggsave("material/test_plot10.pdf", height = unit(4, "cm"), width=unit(10, "cm"))


#### test #####
library(ggplot2)
infolder  <- "./CVresults"
outfolder <- "./material"
mets <- c("exhaustive","random_caret10")
colnms <- c("method", "variate", "design", "number", "RMSE")
outtab <- data.frame(matrix(NA, 0, length(colnms)))
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
      newrow <- data.frame(method = m, variate = variate, design = design,
                           number = number, RMSE=RMSE)
    
    outtab <- rbind(outtab, newrow)
  }
}

unique(outtab$design)

outtab$methodID <- with(outtab, ifelse(method=="exhaustive",0,1))
outtab$rRMSE <- NA

# some rows missing in exhaustive method for OCS data
outtab <- outtab[outtab$number %in% numbers, ]
#outtab <- outtab[outtab$design=="simpleRandom",]

for(variate in c("AGB", "OCS")){
  for(design in unique(outtab$design)){
    for(number in 1:10){
      idx1 <- which(outtab$design == design & outtab$variate == variate & 
                      outtab$number == number)
      idx2 <- which(outtab$design == design & outtab$variate == variate & 
                      outtab$methodID == 0 & outtab$number == number)
      
      outtab$rRMSE[idx1] <- 100 * (outtab$RMSE[idx1] - outtab$RMSE[idx2])/
        outtab$RMSE[idx2]
    }
  }
}

#outtab$method <- factor(outtab$method, levels=c("exhaustive","kndm"))
#outtab$design <- factor(outtab$design, levels=c("exhaustive","knndm"))
outtab$variate <- as.factor(outtab$variate)
#xlabs <- c("RMSE_kndm", "RMSE_indx")
collabs <- c("SRS","syst","clustMed","clustStr","clustGap")
cols <- c("white", "brown", "pink", "orange", "lightgreen")


library(tidyverse)
out_long <- outtab |>
  dplyr::filter(method!="exhaustive") |> 
  dplyr::select(!c(method, number, methodID)) |>
  dplyr::select(!starts_with("RMSE")) |> 
  pivot_longer(!c(variate,design))
#out_long$design <- fct_relevel(out_long$design, c("simpleRandom", "regular", "clusterMedium", "clusterStrong", "clusterGapped"))

#write.csv(out_long, file.path(outfolder, "outtab10_test.csv"))


(gp <- ggplot(out_long,
              aes(x=design, fill=name)) +
    geom_boxplot(aes(y=value, alpha=0.9)) +
    scale_y_continuous(limits=c(-50,30)) +
    facet_wrap(~variate) +
    geom_hline(yintercept = 0))

#ggsave("material/test_plot10.pdf", height = unit(4, "cm"), width=unit(10, "cm"))

