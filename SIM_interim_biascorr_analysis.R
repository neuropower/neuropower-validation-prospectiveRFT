setwd("~/Documents/Onderzoek/Studie_4_propow/InterimPower_Results/interim_SIM_biascorr_15/")
sims <- 
results <- array(NA,dim=c(2000,16,16))
for(i in 1:1040){
  simres <- array(NA,dim=c(16,16))
  file <- paste("estimation_sim_",i,".csv",sep="")
   if(file.exists(file)){
     helparr <- as.matrix(read.table(file,sep=",",dec="."))
     simres[1:(dim(helparr)[1]),] <- helparr
     results[i,,] <- simres        
  }
}

result_names = c(
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
  "TPR_u",
  "FPR_u",
  "TPR_c",
  "FPR_c",
  "discr"
  )

avres <- apply(results,c(2,3),mean,na.rm=TRUE)
varrest <- apply(results,c(2,3),var,na.rm=TRUE)



