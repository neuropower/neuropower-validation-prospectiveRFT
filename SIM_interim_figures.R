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
library(plyr)
library(gridExtra)


HOMEDIR <- "~/Documents/Onderzoek/Studie_4_propow/ProspectivePowerValidation/"
RESDIR <- "/Users/Joke/Documents/Onderzoek/Studie_4_propow/InterimPower_Results/interim/"

effs <- rep(c(0.5,1,1.5,2),each=4)
acts <- rep(c(2,4,6,8),4)
cons_list <- paste("EFFECT:",effs," - ","ACTIVATION SIZE:",acts,sep="")
methods_n <- c("Uncorrected","Benjamini-Hochberg","Bonferroni","Random Field Theory")

## FUNCTIONS

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

# read in results

res.nonad <- read.table(paste(RESDIR,"tables/results_nonadaptive.txt",sep=""))
pre.nonad <- read.table(paste(RESDIR,"tables/power_predicted_nonadaptive.txt",sep=""))
obs.nonad <- read.table(paste(RESDIR,"tables/power_observed_nonadaptive.txt",sep=""))
res.ad <- read.table(paste(RESDIR,"tables/results_adaptive.txt",sep=""))
pre.ad <- read.table(paste(RESDIR,"tables/power_predicted_adaptive.txt",sep=""))
obs.ad <- read.table(paste(RESDIR,"tables/power_observed_adaptive.txt",sep=""))
res.bias <- read.table(paste(RESDIR,"tables/results_bias.txt",sep=""))

# take mean of power predictions over simulations

pre.ad.mn <- ddply(pre.ad,
                   ~subjects+condition+mcp,
                   summarise,
                   TPR=mean(power)
                   )
obs.ad.mn <- ddply(obs.ad,
                   ~subs+condition+mcp,
                   summarise,
                   TPR=mean(TPR),
                   FPR=mean(FPR),
                   FDR=mean(FDR),
                   FWER=mean(FWER)
                   )
obs.ad.mn$condition <- factor(obs.ad.mn$condition)

pre.nonad.mn <- ddply(pre.nonad,
                   ~subjects+condition+mcp,
                   summarise,
                   TPR=mean(power)
)
obs.nonad.mn <- ddply(obs.nonad,
                   ~subs+condition+mcp,
                   summarise,
                   TPR=mean(TPR),
                   FPR=mean(FPR),
                   FDR=mean(FDR),
                   FWER=mean(FWER)
)
obs.nonad.mn$condition <- factor(obs.nonad.mn$condition)

res.bias.mn <- ddply(res.bias,
                     ~condition+effect+width+method+correction,
                     summarise,
                     SS = mean(SS),
                     TPR=mean(TPR),
                     FPR=mean(FPR),
                     FDR=mean(FDR),
                     FWER=mean(FWER))
bias.ad.mn <- pre.ad.mn[,1:3]
bias.ad.mn$bias <- pre.ad.mn$TPR - obs.ad.mn$TPR
res.bias.mn$correction <- factor(res.bias.mn$correction)
res.bias.mn$effect <- factor(res.bias.mn$effect)
res.bias.mn$width <- factor(res.bias.mn$width)




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

pdf(paste(RESDIR,"figures/modelestimation.pdf",sep=""),width=8,height=5)

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
       col=alpha(cols[res.ad$condition]),
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

points(res.ad$efft,
       res.ad$effex,
       col=alpha(cols[res.ad$condition]),
       pch=16,cex=cxp
)

# legend

par(mar=c(0,0,0,0))
plot(1:3, 1:3, col=NA,axes=FALSE,xlab="",ylab="")
x <- rep(seq(2,2.2,length=4),times=4)

y <- rep(seq(1.2,2.2,length=4),each=4)
points(x,y,pch=15,col=cols,cex=2.8)
text(mean(x),2.8,"% of brain active",font=2)
for(i in 1:4){text(x[i],2.5,c(2,4,6,8)[i])}
text(1.8,mean(y),"Effect size",font=2,srt=90)
for(i in 1:4){text(1.9,y[i*4],c(0.5,1,1.5,2)[i])}
dev.off()

#################################
## EVALUATION POWER ESTIMATION ##
#################################

# parameters for text
methname <- c("Uncorrected","False Discovery Rate","Bonferroni","Random Field Theory")
methods <- c("UN","BH","BF","RFT")

# color ramps for power and bias
col.pow <- colorRampPalette(brewer.pal(9,"YlOrRd"))(100)
at.pow <- seq(0, 1, length.out=length(col.pow)+1)
ckey.pow <- list(at=at.pow, col=col.pow)

col.bias <- colorRampPalette(brewer.pal(11,"RdBu")[c(1:5,6,6,7:11)])(50)
at.bias <- seq(-1, 1, length.out=length(col.bias)-1)
ckey.bias <- list(at=at.bias, col=col.bias)

# define scales for x and y axis
level.scales1 <- list(y=list(labels=cons_list,at=1:16,cex=0.7),
					 x=list(labels=seq(15,60,by=5),at=seq(1,46,by=5)))
level.scales2 <- list(y=list(labels=rep("",16),at=1:16,cex=0.7),
					 x=list(labels=seq(15,60,by=5),at=seq(1,46,by=5)))

plots <- list(0)
for (method in methods){
t1 <- levelplot(t(array(pre.ad.mn[pre.ad.mn$mcp==method,]$TPR,dim=c(16,46))),
                scales=level.scales1,
                colorkey=TRUE,
                at=at.pow,
                col.regions=col.pow,
                xlab="",
                par.settings=list(layout.heights=list(top.padding=-3)),                                
                ylab=method,
                aspect="fill")
t2 <- levelplot(t(array(obs.ad.mn[obs.ad.mn$mcp==method,]$TPR,dim=c(16,46))),
                scales=level.scales2,
                colorkey=FALSE,
                at=at.pow,
                col.regions=col.pow,
                xlab="",
                par.settings=list(layout.heights=list(top.padding=-3)),
                ylab=methname[method],
                aspect="fill")
t3 <- levelplot(t(array(bias.ad.mn[bias.ad.mn$mcp==method,]$bias,dim=c(16,46))),
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


pdf(paste(RESDIR,"figures/powerpredictions.pdf",sep=""),width=15,height=17)
grid.arrange(plots[["UN"]][[1]],plots[["UN"]][[2]],
             plots[["BH"]][[1]],plots[["BH"]][[2]],
             plots[["BF"]][[1]],plots[["BF"]][[2]],             
             plots[["RFT"]][[1]],plots[["RFT"]][[2]],
             nrow=4,widths=c(7.7,4))
dev.off()

###########################
## BIAS WHEN NONADAPTIVE ##
###########################

# make graphs

p <- list()
p[[1]] <- p[[2]] <- p[[3]] <- p[[4]] <- list()
p[[1]][[1]] <- ggplot(obs.ad.mn[obs.ad.mn$mcp=="UN",],aes(x=subs,y=FPR,group=condition,colour=condition))
p[[1]][[2]] <- ggplot(obs.nonad.mn[obs.nonad.mn$mcp=="UN",],aes(x=subs,y=FPR,group=condition,colour=condition))
p[[2]][[1]] <- ggplot(obs.ad.mn[obs.ad.mn$mcp=="BH",],aes(x=subs,y=FDR,group=condition,colour=condition))
p[[2]][[2]] <- ggplot(obs.nonad.mn[obs.nonad.mn$mcp=="BH",],aes(x=subs,y=FDR,group=condition,colour=condition))
p[[3]][[1]] <- ggplot(obs.ad.mn[obs.ad.mn$mcp=="BF",],aes(x=subs,y=FWER,group=condition,colour=condition))
p[[3]][[2]] <- ggplot(obs.nonad.mn[obs.nonad.mn$mcp=="BF",],aes(x=subs,y=FWER,group=condition,colour=condition))
p[[4]][[1]] <- ggplot(obs.ad.mn[obs.ad.mn$mcp=="RFT",],aes(x=subs,y=FWER,group=condition,colour=condition))
p[[4]][[2]] <- ggplot(obs.nonad.mn[obs.nonad.mn$mcp=="RFT",],aes(x=subs,y=FWER,group=condition,colour=condition))

# create figures

fig <- list()
fig[[1]] <- fig[[2]] <- fig[[3]] <- fig[[4]] <- list()
titles <- c("FPR control (uncorrected) \n adaptive design","FPR control (uncorrected) \n fixed design", 
            "FDR control (Benjamini-Hochberg) \n adaptive design","FDR control (Benjamini-Hocbherg) \n fixed design",
            "FWER control (Bonferroni) \n adaptive design","FWER control (Bonferroni)  \n fixed design",
            "FWER control (RFT) adaptive design","FWER control (RFT) fixed design")
control <- rep(c("FPR","FDR","FWER","FWER") ,each=2)

k <- 0
for(m in 1:4){
  for(n in 1:2){
    k <- k+1
    fig[[m]][[n]] <- p[[m]][[n]] +
      geom_point(size=3) + 
      geom_line() +
      geom_hline(aes(yintercept=0.05)) +
      theme(panel.background = element_rect(fill = NA, colour = "white"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(fill=NA,colour="black"),
            legend.position="none"
      ) +
      scale_colour_manual(labels=rep("",16), 
                          values = cols) +
      ylim(c(0,0.052)) +
      labs(x="Sample size",
           y = control[k],
           title=titles[k])    
  }
}

# create legend

empty <- ggplot(data.frame()) + geom_point() + xlim(0, 10) + ylim(0, 100) + theme_minimal() + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank())

dat <- data.frame(a = factor(acts), b = factor(effs),c=cons_list)
leg <- ggplotGrob(
  ggplot(unique(subset(dat, select = a:b)), 
         aes(a, b, colour=cons_list)) + 
    geom_point(size=10) +
    theme(panel.background = element_rect(fill = NA, colour = "white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.ticks = element_blank()) +
    scale_colour_manual(labels=rep("",16), 
                        values = cols) +
    labs(x = "effect size",
         y = "percentage activation")
  )


pdf(paste(RESDIR,"figures/FPR_adaptive_vs_nonadaptive.pdf",sep=""),width=15,height=17)
grid.arrange(
  arrangeGrob(fig[[1]][[1]],
              fig[[2]][[1]],
              fig[[3]][[1]],
              fig[[4]][[1]],
              nrow=4),
  arrangeGrob(fig[[1]][[2]],
              fig[[2]][[2]],
              fig[[3]][[2]],
              fig[[4]][[2]],
              nrow=4),
  arrangeGrob(empty,leg,empty,nrow=3),
  ncol=3,
  widths = c(2/5,2/5,1/5),
  main="Comparison of adaptive and non-adaptive designs \nin the coltrol of false positives"
)
dev.off()

########################
## BIAS WHEN ADAPTIVE ##
########################

#False positive rate

p <- list()
p[[1]] <- ggplot(res.bias.mn[res.bias.mn$method=="UN",],aes(x=width,y=FPR,group=interaction(correction,effect),linetype=correction,colour=effect))
p[[2]] <- ggplot(res.bias.mn[res.bias.mn$method=="BH",],aes(x=width,y=FPR,group=interaction(correction,effect),linetype=correction,colour=effect))
p[[3]] <- ggplot(res.bias.mn[res.bias.mn$method=="BF",],aes(x=width,y=FPR,group=interaction(correction,effect),linetype=correction,colour=effect))
p[[4]] <- ggplot(res.bias.mn[res.bias.mn$method=="RFT",],aes(x=width,y=FPR,group=interaction(correction,effect),linetype=correction,colour=effect))

control <- c("FPR","FDR","FWER","FWER")
method <- c("Uncorrected","Benjamini-Hochberg","Bonferroni","Random Field Theory")

fig <- list()
for(m in 1:4){
fig[[m]] <- p[[m]] +
  geom_point(size=3) + 
  geom_line() +
  geom_hline(aes(yintercept=0.05)) +
  theme(panel.background = element_rect(fill = NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
        ) +
  scale_colour_manual(labels=rep(c(0.5,1,1.5,2),4), 
                      values = cols[9:12]) +
    ylim(c(0,0.06)) +
  labs(x="Percentage of brain activated",
       y = control[m],
       title=method[m])
}

legend <- g_legend(p[[1]] + geom_point(size=3) + geom_line() + theme(legend.key = element_rect(fill = NA, colour = NA)) + scale_colour_manual(labels=rep(c(0.5,1,1.5,2),4), 
                                         values = cols[9:12]) )


pdf(paste(RESDIR,"figures/FPR_adaptive_with_and_without_correction.pdf",sep=""),width=15,height=17)
grid.arrange(
  arrangeGrob(fig[[1]],fig[[2]],fig[[3]],fig[[4]],ncol=2),
  legend,
  ncol=2,
  widths=c(3/4,1/4))
dev.off()

######################################
## REQUIRED VS OBTAINED SAMPLE SIZE ##
######################################

mcps <- c("UN","BH","BF","RFT")

# compute minimal sample size

SS.pre <- SS.obs <- eff <- act <- method <- condition <- c()
k <- 0
for(c in 1:16){
  for(m in 1:4){
    k <- k+1
    PRE <- pre.ad.mn[pre.ad.mn$condition==c & pre.ad.mn$mcp==mcps[m],]
    SS.pre[k] <- (15:60)[min(which(PRE$TPR>0.8))]
    OBS <- obs.ad.mn[obs.ad.mn$condition==c & pre.ad.mn$mcp==mcps[m],]
    SS.obs[k] <- (15:60)[min(which(OBS$TPR>0.8))]
    eff[k]<- effs[c]
    act[k] <- acts[c]
    condition[k] <- c
    method[k] <- m
  }
}  

results <- data.frame(SS.pre,SS.obs,eff,act,condition,method)

# reshape data to long format

ressh <- reshape(results,
                 varying=c("SS.pre","SS.obs"),
                 direction="long",
                 v.names="samplesize",
                 timevar="predicted",
                 times=c("predicted","observed"),
                 sep="."
                 )
names(ressh) <- c("eff","act","condition","method","predicted","samplesize","id")
    
ressh$eff <- factor(ressh$eff)
ressh$act <- factor(ressh$act)

# make plots

p <- list()
for(m in 1:4){
im <- ggplot(ressh[ressh$method==m,],
            aes(x=act,
                y=samplesize,
                group=interaction(predicted,eff),
                linetype=predicted,
                colour=eff,
                )
            )
p[[m]] <- im + 
  geom_line() + 
  geom_point() + 
  theme(panel.background = element_rect(fill = NA, colour = 'grey'),
        legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
        ) +
  scale_colour_manual(labels=rep(c(0.5,1,1.5,2),4), 
                      values = cols[9:12]) +
  ylim(c(15,60))+
  labs(x="Percentage of brain active",
       y = "Sample size",
       title=method[m])
}

legend <- g_legend(im+geom_line()+geom_point() + theme(legend.key = element_rect(fill = NA, colour = NA)) +scale_colour_manual(labels=rep(c(0.5,1,1.5,2),4), values = cols[9:12]))

pdf(paste(RESDIR,"figures/samplesize.pdf",sep=""),width=15,height=17)
grid.arrange(
  arrangeGrob(p[[1]],p[[3]],nrow=2),
  arrangeGrob(p[[2]],p[[4]],nrow=2),
  legend,
  ncol=3,
  widths = c(2/5,2/5,1/5)
    )
dev.off()
  

