library(RColorBrewer)
library(reshape)
library(plotrix)
library(scales)
library(lattice)
library(latticeExtra)
library(colorRamps)
library(gplots)
library(gridExtra)
library(vioplot)
library(data.table)
library(ggplot2)

HOMEDIR <- "~/Documents/Onderzoek/Studie_4_propow/ProspectivePowerValidation/"
FIGDIR <- "~/Documents/Onderzoek/Studie_4_propow/ProspectivePower-Paper/Studie_4_v1.4/Figures/"

subs <- 45
sims <- 100

effs <- rep(c(0.5,1,1.5,2),each=4)
acts <- rep(c(2,4,6,8),4)
effn <- rep(c("half","one","onehalf","two"),each=4)

##################################
## READ IN NON-ADAPTIVE RESULTS ##
##################################

RESDIR <- "/Users/Joke/Documents/Onderzoek/Studie_4_propow/InterimPower_Results/power_SIM/"

# results model table
res.nonad <- data.frame()
for(p in 1:sims){
  file <- paste(RESDIR,"estimation_sim_",p,".csv",sep="")
  if(!file.exists(file)){next}
  res <- read.csv(file,header=FALSE,na.strings=c("nan"))
  res$simulation <- p
  res <- res[,c(11,1:10)]
  res.nonad <- rbind(res.nonad,res)
  names(res.nonad) <- c("simulation","eff","p","pi1e","pi1t","effe","efft","effex","sde","sdt","a")
}

# predicted power
power.pre.nonad <- data.frame()
for(p in ((1:sims))){
  for (c in 1:16){
    file <- paste(RESDIR,"powpre_sim_",p,"_w_",acts[c],"_e_",effn[c],".csv",sep="")
    if(!file.exists(file)){next}
    res <- read.table(file,sep=",",dec=".",header=TRUE)
    if(is.null(res$BH)){res$BH <- NA}
    res <- data.frame(p,c,1:59,res$UN,res$BH,res$BF,res$RFT)
    power.pre.nonad <- rbind(power.pre.nonad,res)
  }
}
    
# observed power
power.obs.nonad <- data.frame()
for(p in 1:sims){
  for (c in 1:16){
    file <- paste(RESDIR,"powtru_sim_",p,"_w_",acts[c],"_e_",effn[c],".csv",sep="")
    if(!file.exists(file)){next}
    res <- read.table(file,sep=",",dec=".",header=TRUE)
    resBF <- data.frame(p,1,c,res$BF_TP,res$BF_FP,res$BF_TN,res$BF_FN)
    resBH <- data.frame(p,2,c,res$BH_TP,res$BH_FP,res$BH_TN,res$BH_FN)
    resRFT <- data.frame(p,3,c,res$RFT_TP,res$RFT_FP,res$RFT_TN,res$RFT_FN)    
    resUN <- data.frame(p,4,c,res$UN_TP,res$UN_FP,res$UN_TN,res$UN_FN)
    resTOT <- data.frame(mapply(c,resBF,resBH,resRFT,resUN))
    names(resTOT) <- c("simulation","mcp","condition","TP","FP","TN","FN")
    resTOT$mcp <- ifelse(resTOT$mcp==1,"BF",ifelse(resTOT$mcp==2,"BH",ifelse(resTOT$mcp==3,"RFT","UN")))
    resTOT$subs <- 15:59
    power.obs.nonad <- rbind(power.obs.nonad,resTOT)
  }
}
power.obs.nonad$TPR = power.obs.nonad$TP/(power.obs.nonad$TP+power.obs.nonad$FN)
power.obs.nonad$FPR = power.obs.nonad$FP/(power.obs.nonad$FP+power.obs.nonad$TN)
power.obs.nonad$FDR = power.obs.nonad$FP/(power.obs.nonad$FP+power.obs.nonad$TP)
power.obs.nonad$FWER = ifelse(power.obs.nonad$FP>0,1,0)






# read in results: model estimation
estimation <- array(NA,dim=c(sims,16,10))
range <- (1:sims)[-c(91)]
for (p in range){
    
    
  }
  pred <- read.table,sep=",",dec=".",header=TRUE)
  
  file <- paste(RESDIR,"estimation_sim_",p,".csv",sep="")
  res <- read.csv(file,header=FALSE,na.strings=c("nan"))
  names(res) <- c("eff","p","pi1e","pi1t","effe","efft","effex","sde","sdt","a")
  estimation[p,,] <- data.matrix(res)
}


# read in results: power estimation
powpred <- array(NA,dim=c(sims,4,4,subs,16)) 
powtrue <- array(NA,dim=c(sims,4,4,subs,16)) 

for(s in range){
  print(s)
  for(p in 1:4){
    for(e in 1:4){
      per <- c(2,4,6,8)[p]
      eff <- c("half","one","onehalf","two")[e]
      
      pred <- read.table(paste(RESDIR,"powpre_sim_",s,"_w_",per,"_e_",eff,".csv",sep=""),sep=",",dec=".",header=TRUE)
      if(is.null(pred$BH)){pred$BH <- 0}
      pred <- data.frame(pred$BF,pred$UN,pred$RFT,pred$BH)   
      powpred[s,p,e,,] <-  data.matrix(pred)
      
      true <- read.table(paste(RESDIR,"powtru_sim_",s,"_w_",per,"_e_",eff,".csv",sep=""),sep=",",dec=".",header=TRUE)
      true <- data.frame(true$BF,true$UN,true$RFT,true$BH)
      true[true=="NaN"] <- 0
      powtrue[s,p,e,,] <-  data.matrix(true)
    }
  }
}


powpred.av <- apply(powpred,c(2,3,4,5),mean,na.rm=TRUE)
powtrue.av <- apply(powtrue,c(2,3,4,5),mean,na.rm=TRUE)

powpred3D <- array(NA,dim=c(16,subs,4))
powtrue3D <- array(NA,dim=c(16,subs,4))
k <- 0
for(p in 1:4){
  for(e in 1:4){
    k <- k+1
    powpred3D[k,,] <- powpred.av[p,e,,]
    powtrue3D[k,,] <- powtrue.av[p,e,,]
  }
}


###############################
## EVALUATE ESTIMATION MODEL ##
###############################

Greens <- brewer.pal(9,"Greens")[c(3,5,7,9)]
Blues <- brewer.pal(9,"Blues")[c(3,5,7,9)]
Greys <- brewer.pal(9,"Greys")[c(3,5,7,9)]
Reds <- brewer.pal(9,"Reds")[c(3,5,7,9)]

cols.d <- c()
for(i in 1:4){
  cols.d <- c(cols.d,c(Greens[i],Blues[i],Greys[i],Reds[i]))
}

transp <- 0.2
cxp <- 0.5


pdf(paste(FIGDIR,"FIG_SIM_modelestimation_ss10.pdf",sep=""),width=8,height=5)

par(mar=c(4,4,1.5,1),oma=c(0,0,0,0))
layout(matrix(c(1,2,3,3),2,2,byrow=TRUE),widths=c(0.5,0.5),heights=c(0.7,0.3))

plot(seq(0,1,length=10),seq(0,1,length=10),col=NA,xlab="True prevalence of activation",ylab="Estimated prevalence of activation",axes=FALSE,main="Prevalence of activation")
abline(0,1,lwd=1,col="grey50")
box()
axis(1)
axis(2)
for(i in 16:1){
a <- data.frame(estimation[,i,3],estimation[,i,4])
points(a[,2],a[,1],col=alpha(cols.d[i],transp),pch=16,cex=cxp)
}

plot(seq(2.3,6,length=10),seq(2.3,6,length=10),col=NA,xlab="True expected effect size",ylab="Estimated expected effect size",axes=FALSE,main="Scatter plot of effect size")
abline(0,1,lwd=1,col="grey50")
box()
axis(1)
axis(2)
for(i in 16:1){
  a <- data.frame(estimation[,i,5],estimation[,i,6])
  points(a[,2],a[,1],col=alpha(cols.d[i],transp),pch=16,cex=cxp)
}


par(mar=c(0,0,0,0))
plot(1:3, 1:3, col=NA,axes=FALSE,xlab="",ylab="")
y <- rep(seq(1.2,2.2,length=4),each=4)
x <- rep(seq(2,2.2,length=4),times=4)
points(x,y,pch=15,col=cols.d,cex=2.8)
text(mean(x),2.8,"% of brain active",font=2)
for(i in 1:4){text(x[i],2.5,c(2,4,6,8)[i])}
text(1.8,mean(y),"Effect size",font=2,srt=90)
for(i in 1:4){text(1.9,y[i*4],c(0.5,1,1.5,2)[i])}
dev.off()

#################################
## EVALUATION POWER ESTIMATION ##
#################################

colsmatrix <- matrix(c(Greens,Blues,Greys,Reds),nrow=4,byrow=TRUE)

cons_list <- paste("EFFECT:",effs," - ","ACTIVATION SIZE:",acts,sep="")
methname <- c("Bonferroni","Uncorrected","Random Field Theory","False Discovery Rate")

bias <- powpred3D - powtrue3D
method = 1

col.pow <- colorRampPalette(brewer.pal(9,"YlOrRd"))(100)
at.pow <- seq(0, 1, length.out=length(col.pow)+1)
ckey.pow <- list(at=at.pow, col=col.pow)

col.bias <- colorRampPalette(brewer.pal(11,"RdBu")[c(1:5,6,6,7:11)])(50)
at.bias <- seq(-1, 1, length.out=length(col.bias)-1)
ckey.bias <- list(at=at.bias, col=col.bias)

level.scales1 <- list(y=list(labels=cons_list,at=1:16,cex=0.7),
					 x=list(labels=c(10,15,20,25,30),at=c(1,6,11,16,21,26,30)))
level.scales2 <- list(y=list(labels=rep("",47),at=1:16,cex=0.7),
					 x=list(labels=c(10,15,20,25,30),at=c(1,6,11,16,21,26,30)))

plots <- list(0)
for (method in 1:4){
t1 <- levelplot(t(powpred3D[,,method]),
                scales=level.scales1,
                colorkey=TRUE,
                at=at.pow,
                col.regions=col.pow,
                xlab="",
                par.settings=list(layout.heights=list(top.padding=-3)),                                
                ylab=methname[method],
                aspect="fill")
t2 <- levelplot(t(powtrue3D[,,method]),
                scales=level.scales2,
                colorkey=FALSE,
                at=at.pow,
                col.regions=col.pow,
                xlab="",
                par.settings=list(layout.heights=list(top.padding=-3)),
                ylab=methname[method],
                aspect="fill")
t3 <- levelplot(t(bias[,,method]),
                scales=level.scales2,
                colorkey=TRUE,
                at=at.bias,
                col.regions=col.bias,
                xlab="",
                par.settings=list(layout.heights=list(top.padding=-3)),
                ylab="",
                aspect="fill")
t1t2 <- c(t1,t2,layout=c(2,1),merge.legends=TRUE,x.same=TRUE,y.same=TRUE)
plots[[method]] <- list(t1t2,t3)
}


pdf(paste(FIGDIR,"FIG_SIM_power_ss10.pdf",sep=""),width=15,height=17)
grid.arrange(plots[[1]][[1]],plots[[1]][[2]],
             plots[[2]][[1]],plots[[2]][[2]],
             plots[[3]][[1]],plots[[3]][[2]],             
             plots[[4]][[1]],plots[[4]][[2]],
             nrow=4,widths=c(7.7,4))
dev.off()

############################
## MISMATCH DISTRIBUTIONS ##
############################

# what are the average estimated beta-distributions
ava <- apply(estimation,c(2,3),mean,na.rm=TRUE)[,10]
xr <- seq(from=0,to=1,length=1000)
y_bet <- array(NA,dim=c(16,1000))
k <- 0
for(k in 1:16){y_bet[k,] <- ava[k]*xr^(ava[k]-1)
}


# what are the average effect size distributions
y_exp <- array(NA,dim=c(16,299))
mns <- apply(estimation,c(2,3),mean,na.rm=TRUE)[,6]
sds <- apply(estimation,c(2,3),sd,na.rm=TRUE)[,9]
exc <- 2.3

for(c in 1:16){
  print(c)
  tr <- rnorm(100000000,mns[c],sds[c])
  tr <- tr[tr>exc]
  pr <- e^(-exc*(tr-exc))
  a <- hist(pr,breaks=seq(0,1,length=300),plot=FALSE)
  y_exp[c,] <- a$density
}



pdf(paste(FIGDIR,"FIG_SIM_distribution_mismatch.pdf",sep=""),width=8,height=5)

par(mar=c(4,4,4.5,1),oma=c(0,0,0,0))

layout(matrix(c(1,2,3,3),2,2,byrow=TRUE),widths=c(0.5,0.5),heights=c(0.7,0.3))

plot(xr,y_bet[1,],type="l",col=cols.d[1],ylim=c(0,6),lwd=2,xlab="p-values",ylab="density",main="Beta-distribution")
for(k in 2:16){lines(xr,y_bet[k,],col=cols.d[k],lwd=2)}

plot(seq(from=0,to=1,length=10),seq(from=0,to=6,length=10),type="l",pch="n",col="NA",xlab="p-values",ylab="density",main="Random Field Theory empirical distribution \n of peak p-values under activation")
for(c in 1:16){lines(a$mids,y_exp[c,],col=cols.d[c],lwd=2)}

par(mar=c(0,0,0,0))
plot(1:3, 1:3, col=NA,axes=FALSE,xlab="",ylab="")

y <- rep(seq(1.2,2.2,length=4),each=4)
x <- rep(seq(2,2.2,length=4),times=4)
points(x,y,pch=15,col=cols.d,cex=2.8)

text(mean(x),2.8,"% of brain active",font=2)
for(i in 1:4){text(x[i],2.5,c(2,4,6,8)[i])}
text(1.8,mean(y),"Effect size",font=2,srt=90)
for(i in 1:4){text(1.9,y[i*4],c(0.5,1,1.5,2)[i])}
dev.off()


######################################
## REQUIRED VS OBTAINED SAMPLE SIZE ##
######################################

pilot_sub <- 15
reqsub <- array(NA,dim=c(sims,4,4,4))
obsub <- array(NA,dim=c(sims,4,4,4))

for(i in 1:sims){
  print(i)
  for(j in 1:4){
    for(k in 1:4){
      for(l in 1:4){
        ind <- which(powpred[i,j,k,,l]>0.6)
        minind <- ifelse(sum(ind)>0,min(ind),NA)
        reqsub[i,j,k,l] <- ifelse(sum(ind)>0,(pilot_sub:(pilot_sub+45))[minind],NA)
        ind <- which(powtrue[i,j,k,,l]>0.6)
        minind <- ifelse(sum(ind)>0,min(ind),NA)
        obsub[i,j,k,l] <- ifelse(sum(ind)>0,(pilot_sub:(pilot_sub+45))[minind],NA)
      }
    }
  }
}


vals <- c()
truest <- c()
cond <- c()
MCP <- c()

k <- 0
for(i in 1:4){
  for(j in 1:4){
    k <- k+1
    for(l in 1:4){
      vals <- c(vals,reqsub[,i,j,l],obsub[,i,j,l])
      truest <- c(truest,rep("estimated",sims),rep("observed",sims))
      MCP <- c(MCP,rep(l,sims*2))
      cond <- c(cond,rep(k,sims*2))
    }  
  }
}

res <- data.frame(vals,truest,MCP,cond)
names(res) <- c("samplesize","truest","MCP","condition")
res$newcon = ifelse(res$truest=="estimated",res$condition*2,res$condition*2-1)

bias <- res$samplesize[res$truest=="estimated"]-res$samplesize[res$truest=="observed"]
bias <- data.frame(bias,res$MCP[res$truest=="estimated"],c(res$newcon[res$truest=="estimated"]))
names(bias) <- c("samplesize","MCP","newcon")

resn <- res[$MCP==3,]
p <- ggplot(bias[bias$MCP,],aes(factor(bias$newcon),bias$samplesize))
p+geom_violin(adjust=3,width=2,trim=FALSE)
p+geom_boxplot() + geom_jitte




FDR <- power.obs.nonad[power.obs.nonad$mcp=="BH",]
p <- ggplot(FDR,aes(factor(condition),FDR))
p + geom_violin(trim=FALSE,aes(fill=factor(subs)))

nana <- ddply(FDR,~condition+subs,summarise,mean=mean(FDR))
plot()

p <- ggplot(data=nana, aes(x=subs, y=mean, group = condition, colour = condition))
p+geom_line()

enter image description here

FWE <- power.obs.nonad[power.obs.nonad$mcp=="RFT",]
nana <- ddply(FWE,~condition+subs,summarise,mean=mean(FWER),na.rm=TRUE)
UN <- power.obs.nonad[power.obs.nonad$mcp=="UN",]
nana <- ddply(UN,~condition+subs,summarise,mean=mean(FPR),na.rm=TRUE)
linmod <- lm(FPR~condition+subs,data=power.obs.nonad[power.obs.nonad$mcp=="UN",])


