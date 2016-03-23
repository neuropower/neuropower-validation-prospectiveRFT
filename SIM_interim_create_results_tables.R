library(reshape)

subs <- 46
sims <- 600

effs <- rep(c(0.5,1,1.5,2),each=4)
acts <- rep(c(2,4,6,8),4)
effn <- rep(c("half","one","onehalf","two"),each=4)

RESDIR <- "/Users/Joke/Documents/Onderzoek/Studie_4_propow/InterimPower_Results/interim/"

##################################
## READ IN NON-ADAPTIVE RESULTS ##
##################################

# results model table
res.nonad <- data.frame()
for(p in 1:sims){
  file <- paste(RESDIR,"nonadaptive/estimation_sim_",p,".csv",sep="")
  if(!file.exists(file)){next}
  res <- read.csv(file,header=FALSE,na.strings=c("nan"))
  res$simulation <- p
  res$condition <- 1:16
  res <- res[,c(11,12,1:10)]
  names(res) <- c("simulation","condition","eff","p","pi1e","pi1t","effe","efft","effex","sde","sdt","a")
  res.nonad <- rbind(res.nonad,res)
}
write.table(res.nonad,paste(RESDIR,"tables/results_nonadaptive.txt",sep=""))

# predicted power
power.pre.nonad <- data.frame()
for(p in 1:sims){
  print(p)
  for (c in 1:16){
    file <- paste(RESDIR,"nonadaptive/powpre_sim_",p,"_w_",acts[c],"_e_",effn[c],".csv",sep="")
    if(!file.exists(file)){next}
    res <- read.table(file,sep=",",dec=".",header=TRUE)
    if(is.null(res$BH)){res$BH <- NA}
    res <- data.frame(p,c,15:60,res$UN,res$BH,res$BF,res$RFT)
    power.pre.nonad <- rbind(power.pre.nonad,res)
  }
}
power.pre.nonad <- reshape(power.pre.nonad,
  varying=c("res.UN","res.BH","res.BF","res.RFT"),
  direction="long",
  sep=".",
  timevar="method"
)
names(power.pre.nonad) <- c("simulation","condition","subjects","mcp","power","id")
write.table(power.pre.nonad,paste(RESDIR,"tables/power_predicted_nonadaptive.txt",sep=""))


# observed power
power.obs.nonad <- data.frame()
for(p in 1:sims){
  print(p)
  for (c in 1:16){
    file <- paste(RESDIR,"nonadaptive/powtru_sim_",p,"_w_",acts[c],"_e_",effn[c],".csv",sep="")
    if(!file.exists(file)){next}
    res <- read.table(file,sep=",",dec=".",header=TRUE)
    resBF <- data.frame(p,1,c,res$BF_TP,res$BF_FP,res$BF_TN,res$BF_FN)
    resBH <- data.frame(p,2,c,res$BH_TP,res$BH_FP,res$BH_TN,res$BH_FN)
    resRFT <- data.frame(p,3,c,res$RFT_TP,res$RFT_FP,res$RFT_TN,res$RFT_FN)
    resUN <- data.frame(p,4,c,res$UN_TP,res$UN_FP,res$UN_TN,res$UN_FN)
    resTOT <- data.frame(mapply(c,resBF,resBH,resRFT,resUN))
    names(resTOT) <- c("simulation","mcp","condition","TP","FP","TN","FN")
    resTOT$mcp <- ifelse(resTOT$mcp==1,"BF",ifelse(resTOT$mcp==2,"BH",ifelse(resTOT$mcp==3,"RFT","UN")))
    resTOT$subs <- 15:60
    power.obs.nonad <- rbind(power.obs.nonad,resTOT)
  }
}
power.obs.nonad$TPR = power.obs.nonad$TP/(power.obs.nonad$TP+power.obs.nonad$FN)
power.obs.nonad$FPR = power.obs.nonad$FP/(power.obs.nonad$FP+power.obs.nonad$TN)
power.obs.nonad$FDR = power.obs.nonad$FP/(power.obs.nonad$FP+power.obs.nonad$TP)
power.obs.nonad$FWER = ifelse(power.obs.nonad$FP>0,1,0)
write.table(power.obs.nonad,paste(RESDIR,"tables/power_observed_nonadaptive.txt",sep=""))

##############################
## READ IN ADAPTIVE RESULTS ##
##############################

# results model table
res.ad <- data.frame()
for(p in 1:sims){
  file <- paste(RESDIR,"adaptive/estimation_sim_",p,".csv",sep="")
  if(!file.exists(file)){next}
  res <- read.csv(file,header=FALSE,na.strings=c("nan"))
  res$simulation <- p
  res$condition <- 1:(dim(res)[1])
  res <- res[,c(11,12,1:10)]
  names(res) <- c("simulation","condition","eff","p","pi1e","pi1t","effe","efft","effex","sde","sdt","a")
  res.ad <- rbind(res.ad,res)
}
write.table(res.ad,paste(RESDIR,"tables/results_adaptive.txt",sep=""))

# predicted power
power.pre.ad <- data.frame()
for(p in 1:sims){
  print(p)
  for (c in 1:16){
    file <- paste(RESDIR,"adaptive/powpre_sim_",p,"_w_",acts[c],"_e_",effn[c],".csv",sep="")
    if(!file.exists(file)){next}
    res <- read.table(file,sep=",",dec=".",header=TRUE)
    if(is.null(res$BH)){res$BH <- NA}
    res <- data.frame(p,c,15:60,res$UN,res$BH,res$BF,res$RFT)
    power.pre.ad <- rbind(power.pre.ad,res)
  }
}
power.pre.ad <- reshape(power.pre.ad,
  varying=c("res.UN","res.BH","res.BF","res.RFT"),
  direction="long",
  sep=".",
  timevar="method"
)
names(power.pre.ad) <- c("simulation","condition","subjects","mcp","power","id")
write.table(power.pre.ad,paste(RESDIR,"tables/power_predicted_adaptive.txt",sep=""))


# observed power
power.obs.ad <- data.frame()
for(p in 1:sims){
  print(p)
  for (c in 1:16){
    file <- paste(RESDIR,"adaptive/powtru_sim_",p,"_w_",acts[c],"_e_",effn[c],".csv",sep="")
    if(!file.exists(file)){next}
    res <- read.table(file,sep=",",dec=".",header=TRUE)
    resBF <- data.frame(p,1,c,res$BF_TP,res$BF_FP,res$BF_TN,res$BF_FN)
    resBH <- data.frame(p,2,c,res$BH_TP,res$BH_FP,res$BH_TN,res$BH_FN)
    resRFT <- data.frame(p,3,c,res$RFT_TP,res$RFT_FP,res$RFT_TN,res$RFT_FN)
    resUN <- data.frame(p,4,c,res$UN_TP,res$UN_FP,res$UN_TN,res$UN_FN)
    resTOT <- data.frame(mapply(c,resBF,resBH,resRFT,resUN))
    names(resTOT) <- c("simulation","mcp","condition","TP","FP","TN","FN")
    resTOT$mcp <- ifelse(resTOT$mcp==1,"BF",ifelse(resTOT$mcp==2,"BH",ifelse(resTOT$mcp==3,"RFT","UN")))
    resTOT$subs <- 15:60
    power.obs.ad <- rbind(power.obs.ad,resTOT)
  }
}
power.obs.ad$TPR = power.obs.ad$TP/(power.obs.ad$TP+power.obs.ad$FN)
power.obs.ad$FPR = power.obs.ad$FP/(power.obs.ad$FP+power.obs.ad$TN)
power.obs.ad$FDR = power.obs.ad$FP/(power.obs.ad$FP+power.obs.ad$TP)
power.obs.ad$FWER = ifelse(power.obs.ad$FP>0,1,0)
write.table(power.obs.ad,paste(RESDIR,"tables/power_observed_adaptive.txt",sep=""))

##########################
## READ IN BIAS RESULTS ##
##########################

# results model table
res.bias <- data.frame()
for(p in 1:sims){
  for(c in 1:2){
    str <- ifelse(c == 1,"_cor","_uncor")
    file <- paste(RESDIR,"bias/estimation_sim_",p,str,".csv",sep="")
    if(!file.exists(file)){next}
    res <- read.csv(file,header=FALSE,na.strings=c("nan"))
    names(res) <- c("effect","width","pi1e","pi1t","effe","efft","expeff","sde","sdt","buma",
    "SS_UN","SS_BH","SS_RFT","SS_BF","discrepancy",
    "UN_TP","UN_FP","UN_FN","UN_TN",
    "BH_TP","BH_FP","BH_FN","BH_TN",
    "RFT_TP","RFT_FP","RFT_FN","RFT_TN",
    "BF_TP","BF_FP","BF_FN","BF_TN")
    ressh <- reshape(res,
      varying=list(c("SS_UN","SS_BH","SS_RFT","SS_BF"),
      c("UN_TP","BH_TP","BF_TP","RFT_TP"),
      c("UN_FP","BH_FP","BF_FP","RFT_FP"),
      c("UN_FN","BH_FN","BF_FN","RFT_FN"),
      c("UN_TN","BH_TN","BF_TN","RFT_TN")),
      v.names=c("SS","TP","FP","FN","TN"),
      direction="long",
      sep="_",
      timevar="method"
      )
    ressh$method <- ifelse(ressh$method==1,"UN",ifelse(ressh$method==2,"BH",ifelse(ressh$method==3,"RFT","BF")))
    ressh$simulation <- p
    ressh$cor <- ifelse(c==1,0,1)
    ressh <- ressh[,c(19,18,1,2,12,20,3:11,13:17)]
    names(ressh) <- c("simulation","condition","effect","width","method","correction","pi1e","pi1t","effe","efft","expeff","sde","sdt","buma","discrepancy","SS","TP","FP","FN","TN")
    ressh$TPR = ressh$TP/(ressh$TP+ressh$FN)
    ressh$FPR = ressh$FP/(ressh$FP+ressh$TN)
    ressh$FDR = ressh$FP/(ressh$FP+ressh$TP)
    ressh$FWER = ifelse(ressh$FP>0,1,0)
    res.bias <- rbind(res.bias,ressh)
  }
}
write.table(res.bias,paste(RESDIR,"tables/results_bias.txt",sep=""))
