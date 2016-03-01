library(RColorBrewer)
library(reshape)
library(plotrix)
library(scales)
library(lattice)
library(latticeExtra)
library(colorRamps)
library(gplots)
library(gridExtra)


RESDIR <- "/Users/Joke/Documents/Onderzoek/Studie_4_propow/ProspectivePowerValidation_Results/power_peak_SIM/"
HOMEDIR <- "~/Documents/Onderzoek/Studie_4_propow/ProspectivePowerValidation/"
FIGDIR <- "~/Documents/Onderzoek/Studie_4_propow/ProspectivePowerValidation_Figures/SIM/"

setwd(FIGDIR)

effs <- rep(c(0.5,1,1.5,2),each=4)
acts <- rep(c(2,4,6,8),4)
effn <- rep(c("half","one","onehalf","two"),each=4)
cons <- data.frame(effs,acts)

# read in results: model estimation
estimation <- array(NA,dim=c(500,16,7))
range <- (1:500)[-c(95,252,351,371,377,403,472,475)]
for (p in range){
  file <- paste(RESDIR,"estimation_sim_",p,".csv",sep="")
  res <- read.csv(file,header=FALSE,na.strings=c("nan"))
  names(res) <- c("eff","p","pi1e","pi1t","effe","efft","effec")
  estimation[p,,] <- as.numeric(as.matrix(res))
}

# read in results: power estimation
powtru <- array(NA,dim=c(500,16,20,4))
powpre <- array(NA,dim=c(500,16,20,4))
for (p in range){
  for(c in 1:16){
    prefile <- paste(RESDIR,"powpre_sim_",p-1,"_w_",acts[c],"_e_",effn[c],".csv",sep="")
    if(file.exists(prefile)){
    pre <- read.csv(prefile,na.strings=c("nan"),header=TRUE)
    if(dim(pre)[2]==3){pre$BH=NA}
    pre <- data.frame(pre$BF,pre$UN,pre$RFT,pre$BH)
    powpre[p,c,,] <- as.numeric(as.matrix(pre))
    }
    
    trufile <- paste(RESDIR,"powtru_sim_",p-1,"_w_",acts[c],"_e_",effn[c],".csv",sep="")
    if(file.exists(trufile)){
      tru <- read.csv(trufile,na.strings=c("nan"),header=TRUE)
      if(dim(tru)[2]==3){tru$BH=NA}
      pre <- data.frame(tru$BF,tru$UN,tru$RFT,tru$BH)
      powtru[p,c,,] <- as.numeric(as.matrix(tru))   
    }
  }
}


## EVALUATE ESTIMATION PI0

Greens <- brewer.pal(9,"Greens")[c(3,5,7,9)]
Blues <- brewer.pal(9,"Blues")[c(3,5,7,9)]
Greys <- brewer.pal(9,"Greys")[c(3,5,7,9)]
Reds <- brewer.pal(9,"Reds")[c(3,5,7,9)]

cols.d <- c()
for(i in 1:4){
  cols.d <- c(cols.d,c(Greens[i],Blues[i],Greys[i],Reds[i]))
}

transp <- 0.4


pdf(paste(FIGDIR,"FIG_SIM_modelestimation.pdf",sep=""),width=8,height=5)

par(mar=c(4,4,1.5,1),oma=c(0,0,0,0))

layout(matrix(c(1,2,3,3),2,2,byrow=TRUE),widths=c(0.5,0.5),heights=c(0.75,0.25))
plot(seq(0,1,length=10),seq(0,1,length=10),col=NA,xlab="True prevalence of activation",ylab="Estimated prevalence of activation",axes=FALSE,main="Prevalence of activation")
abline(0,1,lwd=1,col="grey50")
box()
axis(1)
axis(2)

for(i in 16:1){
a <- data.frame(estimation[,i,3],estimation[,i,4])
points(a[,2],a[,1],col=alpha(cols.d[i],transp),pch=16,cex=0.75)
}


plot(seq(2.3,6,length=10),seq(2.3,6,length=10),col=NA,xlab="True expected effect size",ylab="Estimated expected effect size",axes=FALSE,main="Scatter plot of effect size")
abline(0,1,lwd=1,col="grey50")
box()
axis(1)
axis(2)

for(i in 16:1){
  a <- data.frame(estimation[,i,7],estimation[,i,6])
  points(a[,2],a[,1],col=alpha(cols.d[i],transp),pch=16,cex=0.75)
}


par(mar=c(0,0,0,0))
plot(1:3, 1:3, col=NA,axes=FALSE,xlab="",ylab="")

y <- rep(seq(1.8,2.4,length=4),each=4)
x <- rep(seq(2,2.12,length=4),times=4)
points(x,y,pch=15,col=cols.d,cex=2.8)

text(2.05,2.8,"% of brain active",font=2)
for(i in 1:4){text(x[i],2.6,(1:4)[i])}
text(1.85,2.1,"Effect size",font=2,srt=90)
for(i in 1:4){text(1.9,y[i*4],c(0.5,1,1.5,2)[i])}


dev.off()

################################################

powpred <- array(NA,dim=c(500,4,4,20,4)) 
powtrue <- array(NA,dim=c(500,4,4,20,4)) 

colsmatrix <- matrix(c(Greens,Blues,Greys,Reds),nrow=4,byrow=TRUE)

for(s in (1:500)[-c(34)]){
  print(s)
  for(p in 1:4){
    for(e in 1:4){
      per <- c(2,4,6,8)[p]
      eff <- c("half","one","onehalf","two")[e]
      
      pred <- read.table(paste(RESDIR,"powpre_sim_",s,"_w_",per,"_e_",eff,".csv",sep=""),sep=",",dec=".",header=TRUE)
      if(is.null(pred$BH)){pred$BH <- 0}
      pred <- data.frame(pred$BF,pred$UN,pred$RFT,pred$BH)   
      powpred[s,p,e,,] <-  as.numeric(as.matrix(pred))
      
      true <- read.table(paste(RESDIR,"powtru_sim_",s,"_w_",per,"_e_",eff,".csv",sep=""),sep=",",dec=".",header=TRUE)
      true <- data.frame(true$BF,true$UN,true$RFT,true$BH)
      true[true=="NaN"] <- 0
      powtrue[s,p,e,,] <-  as.numeric(as.matrix(true))
    }
  }
}

powpred.av <- apply(powpred,c(2,3,4,5),mean,na.rm=TRUE)
powtrue.av <- apply(powtrue,c(2,3,4,5),mean,na.rm=TRUE)


# FIGURE AS TOM WANTS
par(mar=c(4,4,1.5,1),oma=c(0,0,0,0),mfrow=c(1,4))

for(m in 1:4){
plot(seq(10,29,length=10),seq(0,1,length=10),col=NA,xlab="Subjects",ylab="Power",axes=FALSE,main="Power estimation")
abline(0,1,lwd=1,col="grey50")
box()
axis(1)
axis(2)

for(p in 1:4){
  for(e in 1:4){
  lines(10:29,powpred.av[p,e,,m],col=colsmatrix[e,p],lwd=3)  
  lines(10:29,powtrue.av[p,e,,m],col=colsmatrix[e,p],lty=2,lwd=3)
  }
}
}

# FIGURE AS IT WAS

powpred3D <- array(NA,dim=c(16,20,4))
powtrue3D <- array(NA,dim=c(16,20,4))

k <- 0
for(p in 1:4){
  for(e in 1:4){
    k <- k+1
    powpred3D[k,,] <- powpred.av[p,e,,]
    powtrue3D[k,,] <- powtrue.av[p,e,,]
  }
}


cons_list <- paste("EFFECT:",effs," - ","ACTIVATION SIZE:",acts,sep="")
methname <- c("Uncorrected","Bonferroni","Random Field Theory","False Discovery Rate")

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


pdf(paste(FIGDIR,"FIG_SIM_power.pdf",sep=""),width=15,height=17)
grid.arrange(plots[[1]][[1]],plots[[1]][[2]],
             plots[[2]][[1]],plots[[2]][[2]],
             plots[[3]][[1]],plots[[3]][[2]],             
             plots[[4]][[1]],plots[[4]][[2]],
             nrow=4,widths=c(7.7,4))
dev.off()