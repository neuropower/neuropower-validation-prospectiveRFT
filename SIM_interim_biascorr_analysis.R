setwd("~/Documents/Onderzoek/Studie_4_propow/InterimPower_Results/interim_biascorr/")

results <- array(NA,dim=c(1040,16,13))
for(i in 1:1040){
  simres <- array(NA,dim=c(16,13))
  file <- paste("interim_SIM_biascorr_15/estimation_sim_",i,".csv",sep="")
   if(file.exists(file)){
     helparr <- as.matrix(read.table(file,sep=",",dec="."))
     simres[1:(dim(helparr)[1]),] <- helparr
     results[i,,] <- simres        
  }
}

result_names = (
  "es",
  "width",
  "pi1e",
  "pi1t",
  "effe",
  "efft",
  "sde",
  "sdt",
  "buma",
  "estss",
  "TPR",
  "FP"
  )

results[,,13] <- ifelse(results[,,13]>0,1,0)
avres <- apply(results,c(2,3),mean,na.rm=TRUE)
varrest <- apply(results,c(2,3),var,na.rm=TRUE)
