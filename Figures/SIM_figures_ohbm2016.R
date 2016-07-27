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


RESDIR <- "/Users/Joke/Documents/Onderzoek/Presentations/2016_HBM_Geneve/InterimPoster/"


effs <- rep(c(0.5,1,1.5,2),each=4)
acts <- rep(c(2,4,6,8),4)
cons_list <- paste("EFFECT:",effs," - ","ACTIVATION SIZE:",acts,sep="")
methods_n <- c("Uncorrected","Benjamini-Hochberg","Bonferroni","Random Field Theory")
u <- 3.2

## FUNCTIONS

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

# result files name

res.nonad.file = paste(RESDIR,"estimation_sim_predictive_",u,"_CS.csv",sep="")
pre.nonad.file = paste(RESDIR,"prediction_sim_predictive_",u,"_CS.csv",sep="")
true.nonad.file = paste(RESDIR,"true_sim_predictive_",u,"_CS.csv",sep="")
cond.nonad.file = paste(RESDIR,"conditional_sim_predictive_",u,"_CS.csv",sep="")
res.ad.file = paste(RESDIR,"estimation_sim_adaptive_",u,"_CS.csv",sep="")
pre.ad.file = paste(RESDIR,"prediction_sim_adaptive_",u,"_CS.csv",sep="")
true.ad.file = paste(RESDIR,"true_sim_adaptive_",u,"_CS.csv",sep="")
cond.ad.file = paste(RESDIR,"conditional_sim_adaptive_",u,"_CS.csv",sep="")

# read results

res.nonad <- read.table(res.nonad.file,header=TRUE,sep=",")
pre.nonad <- read.table(pre.nonad.file,header=TRUE,sep=",")
obs.nonad <- read.table(true.nonad.file,header=TRUE,sep=",")
cond.nonad <- read.table(cond.nonad.file,header=TRUE,sep=",")
cond.nonad$adaptive <- "predictive"
res.ad <- read.table(res.ad.file,header=TRUE,sep=",")
pre.ad <- read.table(pre.ad.file,header=TRUE,sep=",")
obs.ad <- read.table(true.ad.file,header=TRUE,sep=",")
cond.ad <- read.table(cond.ad.file,header=TRUE,sep=",")
cond.ad$adaptive <- "adaptive"

cntr <- 0
for (eff in c(0.5,1,1.5,2)){
  for (act in c(2,4,6,8)){
    cntr <- cntr+1
    res.nonad$condition[res.nonad$es == eff & res.nonad$activation == act] = cntr
    res.ad$condition[res.ad$es == eff & res.ad$activation == act] = cntr
    
  }
}


cond.tot <- rbind(cond.nonad,cond.ad)
cond.tot$condition = factor(cond.tot$condition)

# average over simulations
# nonadaptive

pre.nonad.mn <- ddply(pre.nonad,
                   ~subjects+condition+mcp,
                   summarise,
                   TPR=mean(power)
)

obs.nonad.mn <- ddply(obs.nonad,
                   ~subjects+condition+mcp,
                   summarise,
                   TPR=mean(power))
obs.nonad.mn$condition <- factor(obs.nonad.mn$condition)

bias.nonad.mn = pre.nonad.mn[,1:3]
bias.nonad.mn$bias = pre.nonad.mn[,4]-obs.nonad.mn[,4]

cond.nonad.mn <- ddply(cond.nonad,
                    ~condition + mcp,
                    summarise,
                    power = mean(power),
                    predicted = mean(predicted),
                    true = mean(true),
                    FWE = mean(FWE,na.rm=TRUE),
                    FDR = mean(FDR,na.rm=TRUE),
                    FPR = mean(FPR,na.rm=TRUE)
                    )

# average over simulations
# adaptive

pre.ad.mn <- ddply(pre.ad,
                      ~subjects+condition+mcp,
                      summarise,
                      TPR=mean(power)
)
pre.ad.mn$true <- "predicted"
pre.ad.mn$condition <- factor(pre.ad.mn$condition)


obs.ad.mn <- ddply(obs.ad,
                      ~subjects+condition+mcp,
                      summarise,
                      TPR=mean(power))
obs.ad.mn$true <- "true"
obs.ad.mn$condition <- factor(pre.ad.mn$condition)

ad.mn <- rbind(obs.ad.mn,pre.ad.mn)


bias.ad.mn = pre.ad.mn[,1:3]
bias.ad.mn$bias = pre.ad.mn[,4]-obs.ad.mn[,4]

cond.ad.mn <- ddply(cond.ad,
                       ~condition + mcp,
                       summarise,
                       power = mean(power),
                       predicted = mean(predicted),
                       true = mean(true),
                      FWE = mean(FWE,na.rm=TRUE),
                      FDR = mean(FDR,na.rm=TRUE),
                      FPR = mean(FPR,na.rm=TRUE)                    
)

# average over simulations
# conditional total

cond.mn <- ddply(cond.tot,
                       ~condition + mcp + adaptive,
                       summarise,
                       power = mean(power),
                       predicted = mean(predicted),
                       FPR = mean(FPR,na.rm=TRUE),
                       FDR = mean(FDR,na.rm=TRUE),
                       FWE = mean(FWE,na.rm=TRUE),
                       true = mean(true,na.rm=TRUE)
)


########################################
## EVALUATE ESTIMATION MODEL ADAPTIVE ##
########################################

Greens <- c('#F5BCBC','#E86464','#C81E1E','#6F1111')
Blues <- c('#D7D8DA','#A2A5A9','#6E7277','#3D3F42')
Greys <- c('#E7DECB','#C7B285','#9D8248','#574828')
Reds <- c('#C1F0EA','#6FDCCE','#2DB9A6','#19675C')
cols <- c()
for(i in 1:4){cols <- c(cols,c(Greens[i],Blues[i],Greys[i],Reds[i]))}
cols <- c(Greens,Blues,Greys,Reds)
transp <- 0.2
cxp <- 0.6

pdf(paste(RESDIR,"figures/modelestimation.pdf",sep=""),width=3,height=3)

ggplot(res.ad,aes(x=est,y=ese,colour=factor(condition))) +
  geom_abline(intercept = 0, slope = 1) + 
  geom_point() +
  theme(panel.background = element_rect(fill = NA, colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.ticks = element_blank(), 
        legend.title=element_blank(),
        legend.position="none",     
        panel.border = element_rect(fill=NA,colour="black")
  ) +
  scale_colour_manual(values=cols) + 
  ylim(c(2,6)) +
  xlim(c(2,6)) +
  labs(x = "True effect size",
       y = "Estimated effect size",
       title="") 

dev.off()

############
## LEGEND ##
############


pdf(paste(RESDIR,"figures/legend.pdf",sep=""),width=2,height=2)

par(mar=c(0,0,0,0))
plot(1:3, 1:3, col=NA,axes=FALSE,xlab="",ylab="")
x <- rep(seq(1.8,2.8,length=4),times=4)

y <- rep(seq(1.2,2.2,length=4),each=4)
points(x,y,pch=15,col=cols,cex=2.8)
text(mean(x),2.8,"% of brain active",font=2)
for(i in 1:4){text(x[i],2.5,c(2,4,6,8)[i])}
text(1,mean(y),"Effect size",font=2,srt=90)
for(i in 1:4){text(1.3,y[i*4],c(0.5,1,1.5,2)[i])}

dev.off()



##################
## POWER CURVES ##
##################


pdf(paste(RESDIR,"figures/powercurves.pdf",sep=""),width=4,height=4)

ggplot(ad.mn[ad.mn$mcp=="RF",],aes(x=subjects,y=TPR,colour=factor(condition),linetype=factor(true))) +
  geom_line(size=0.5) +
  theme(panel.background = element_rect(fill = NA, colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.background = element_rect(fill = "transparent",colour = NA),
        #plot.background = element_rect(fill = "transparent",colour = NA),
        axis.ticks = element_blank(), 
        legend.title=element_blank(),
        legend.position="bottom",    
        legend.key = element_rect(fill="transparent",colour="transparent"),
        panel.border = element_rect(fill=NA,colour="black")
  ) +
  scale_colour_manual(values=cols) + 
  ylim(c(0,1)) +
  labs(x = "Sample size",
       y = "Power",
       title="")  +
  guides(colour=FALSE)

dev.off()

###########################
## BIAS WHEN NONADAPTIVE ##
###########################



pUN <- ggplot(data=cond.tot[cond.tot$mcp=="UN",],aes(interaction(factor(adaptive),factor(condition)),FPR,fill=condition,alpha=adaptive)) +
  geom_violin(trim=TRUE) +
  geom_hline(aes(yintercept=0.05)) +
  theme(panel.background = element_rect(fill = NA, colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        legend.title=element_blank(),
        legend.position="bottom",     
        panel.border = element_rect(fill=NA,colour="black")
  ) +
  scale_fill_manual(values = cols) +
  scale_colour_manual(values=c("black",NA)) + 
  ylim(c(0,0.08)) +
  labs(x = "Condition",
       y = "FPR",
       title="False Positive Rate Control") +
  guides(fill=FALSE)


pRFT <- ggplot(cond.mn[cond.mn$mcp=="RF",],aes(x=factor(condition),y=FWE,colour=condition,shape=adaptive)) +
  geom_point(size=3) +
  geom_hline(aes(yintercept=0.05)) +
    theme(panel.background = element_rect(fill = NA, colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #legend.title=element_blank(),
        legend.position="none",
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        legend.position="bottom",     
        panel.border = element_rect(fill=NA,colour="black")
  ) +
  scale_colour_manual(values=cols) +
  ylim(c(0.0,0.05)) +
  labs(x = "Condition",
       y = "FWE",
       title="Familywise Error Rate Control")    



pdf(paste(RESDIR,"figures/FPR_adaptive_vs_nonadaptive.pdf",sep=""),width=4,height=5)
grid.arrange(pUN,pRFT,
  nrow=2,
  main="",
  heights=c(3/5,2/5)
)
dev.off()


###########################
## FINAL POWER ##
###########################

pUN <- ggplot(cond.ad[cond.ad$mcp=="UN",],aes(x=factor(condition),y=power,fill=factor(condition))) +
  geom_hline(aes(yintercept=0.8)) +
  geom_boxplot() + 
  theme(panel.background = element_rect(fill = NA, colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        legend.position="none",     
        panel.border = element_rect(fill=NA,colour="black")
  ) +
  scale_fill_manual(values = cols[-1]) +
  scale_colour_manual(values=c("black",NA)) + 
  ylim(c(0.6,1)) +
  labs(x = "",
       y = "Power",
       title="False Positive Rate Control")    


pRFT <- ggplot(cond.ad[cond.ad$mcp=="RF",],aes(factor(condition),power,fill=factor(condition))) +
  geom_boxplot() + 
  #geom_violin(trim=TRUE) +
  geom_hline(aes(yintercept=0.8)) +
   theme(panel.background = element_rect(fill = NA, colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        legend.position="none",     
        panel.border = element_rect(fill=NA,colour="black")
  ) +
  scale_fill_manual(values=cols[-(1:4)]) +
  ylim(c(0.6,1)) +
  labs(x = "",
       y = "Power",
       title="Familywise Error Rate Control")    



pdf(paste(RESDIR,"figures/ResultingPower.pdf",sep=""),width=4,height=5)
grid.arrange(pUN,pRFT,
             nrow=2,
             main=""
)

dev.off()


