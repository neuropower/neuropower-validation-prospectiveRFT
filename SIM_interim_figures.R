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
RESDIR <- "/Users/Joke/Documents/Onderzoek/Studie_4_propow/InterimPower_Results/interim/"

effs <- rep(c(0.5,1,1.5,2),each=4)
acts <- rep(c(2,4,6,8),4)

# read in results

res.nonad <- read.table(paste(RESDIR,"tables/results_nonadaptive.txt",sep=""))
pre.nonad <- read.table(paste(RESDIR,"tables/power_predicted_nonadaptive.txt",sep=""))
obs.nonad <- read.table(paste(RESDIR,"tables/power_observed_nonadaptive.txt",sep=""))
res.ad <- read.table(paste(RESDIR,"tables/results_adaptive.txt",sep=""))
pre.ad <- read.table(paste(RESDIR,"tables/power_predicted_adaptive.txt",sep=""))
obs.ad <- read.table(paste(RESDIR,"tables/power_observed_adaptive.txt",sep=""))

res.bias <- read.table(paste(RESDIR,"tables/results_bias.txt",sep=""))





# powpred.av <- apply(powpred,c(2,3,4,5),mean,na.rm=TRUE)
# powtrue.av <- apply(powtrue,c(2,3,4,5),mean,na.rm=TRUE)
# 
# powpred3D <- array(NA,dim=c(16,subs,4))
# powtrue3D <- array(NA,dim=c(16,subs,4))
# k <- 0
# for(p in 1:4){
#   for(e in 1:4){
#     k <- k+1
#     powpred3D[k,,] <- powpred.av[p,e,,]
#     powtrue3D[k,,] <- powtrue.av[p,e,,]
#   }
# }


########################################
## EVALUATE ESTIMATION MODEL ADAPTIVE ##
########################################

Greens <- brewer.pal(9,"Greens")[c(3,5,7,9)]
Blues <- brewer.pal(9,"Blues")[c(3,5,7,9)]
Greys <- brewer.pal(9,"Greys")[c(3,5,7,9)]
Reds <- brewer.pal(9,"Reds")[c(3,5,7,9)]
cols <- c()
for(i in 1:4){cols <- c(cols,c(Greens[i],Blues[i],Greys[i],Reds[i]))}
transp <- 0.2
cxp <- 0.3

pdf(paste(FIGDIR,"FIG_SIM_modelestimation_ss10.pdf",sep=""),width=8,height=5)

# plot layout

par(mar=c(4,4,1.5,1),oma=c(0,0,0,0))
layout(
  matrix(
    c(
      1,2,
      3,3
      ),
    2,2,byrow=TRUE),
  widths=c(0.5,0.5),
  heights=c(0.7,0.3)
  )

# plot 1: pi1 estimation

plot(seq(0,1,length=10),
     seq(0,1,length=10),
     col=NA,
     xlab="True prevalence of activation",
     ylab="Estimated prevalence of activation",
     axes=FALSE,main="Prevalence of activation"
     )

abline(0,1,lwd=1,col="grey50")
box();axis(1);axis(2)
points(res.ad$pi1t,
       res.ad$pi1e,
       col=alpha(cols[res$condition]),
       pch=16,cex=cxp
       )

# plot 2: model estimation

plot(seq(2.3,6,length=10),
     seq(2.3,6,length=10),
     col=NA,
     xlab="True expected effect size",
     ylab="Estimated expected effect size",
     axes=FALSE,
     main="Scatter plot of effect size")

abline(0,1,lwd=1,col="grey50")
box();axis(1);axis(2)

points(res.nonad$efft,
       res.nonad$effex,
       col=alpha(cols[res$condition]),
       pch=16,cex=cxp
)

# legend

par(mar=c(0,0,0,0))
plot(1:3, 1:3, col=NA,axes=FALSE,xlab="",ylab="")
x <- rep(seq(2,2.2,length=4),times=4)

y <- rep(seq(1.2,2.2,length=4),each=4)
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

########################
## BIAS WHEN ADAPTIVE ##
########################

head(res.bias)



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


