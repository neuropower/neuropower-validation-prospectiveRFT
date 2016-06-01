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


RESDIR <- "/Users/Joke/Documents/Onderzoek/ProjectsOngoing/Power/ValidationResults/"


effs <- rep(c(0.5,1,1.5,2),each=4)
acts <- rep(c(2,4,6,8),4)
cons_list <- paste("EFFECT:",effs," - ","ACTIVATION SIZE:",acts,sep="")
methods_n <- c("Uncorrected","Benjamini-Hochberg","Bonferroni","Random Field Theory")
u <- 2

## FUNCTIONS

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

# read in results

res.nonad.file = paste(RESDIR,"estimation_sim_predictive_2.3_CS.csv",sep="")
pre.nonad.file = paste(RESDIR,"prediction_sim_predictive_2.3_CS.csv",sep="")
true.nonad.file = paste(RESDIR,"true_sim_predictive_2.3_CS.csv",sep="")
cond.nonad.file = paste(RESDIR,"conditional_sim_predictive_2.3_CS.csv",sep="")

res.nonad <- read.table(res.nonad.file,header=TRUE,sep=",")
pre.nonad <- read.table(pre.nonad.file,header=TRUE,sep=",")
obs.nonad <- read.table(true.nonad.file,header=TRUE,sep=",")
cond.nonad <- read.table(cond.nonad.file,header=TRUE,sep=",")

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

cond.mn <- ddply(cond.nonad,
                    ~condition + mcp,
                    summarise,
                    power = mean(power),
                    predicted = mean(predicted),
                    true = mean(true)
                    )

########################################
## EVALUATE ESTIMATION MODEL ADAPTIVE ##
########################################

Greens <- brewer.pal(9,"Greens")[c(3,5,7,9)]
Blues <- brewer.pal(9,"Blues")[c(3,5,7,9)]
Greys <- brewer.pal(9,"Greys")[c(3,5,7,9)]
Reds <- brewer.pal(9,"Reds")[c(3,5,7,9)]
cols <- c()
for(i in 1:4){cols <- c(cols,c(Greens[i],Blues[i],Greys[i],Reds[i]))}
cols <- c(Greens,Blues,Greys,Reds)
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
points(res.nonad$pi1t,
       res.nonad$pi1e,
       col=alpha(cols[res.nonad$condition]),
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

points(res.nonad$est,
       res.nonad$ese,
       col=alpha(cols[res.nonad$condition]),
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
methods <- c("UN","BH","BF","RF")

# color ramps for power and bias
col.pow <- colorRampPalette(brewer.pal(9,"YlOrRd"))(100)
at.pow <- seq(0, 1, length.out=length(col.pow)+1)
ckey.pow <- list(at=at.pow, col=col.pow)

col.bias <- colorRampPalette(brewer.pal(11,"RdBu")[c(1:5,6,6,7:11)])(50)
at.bias <- seq(-1, 1, length.out=length(col.bias)-1)
ckey.bias <- list(at=at.bias, col=col.bias)

# define scales for x and y axis
level.scales1 <- list(y=list(labels=cons_list,at=1:16,cex=0.7),
					 x=list(labels=seq(15,40,by=5),at=seq(1,26,by=5)))
level.scales2 <- list(y=list(labels=rep("",16),at=1:16,cex=0.7),
					 x=list(labels=seq(15,40,by=5),at=seq(1,26,by=5)))

plots <- list(0)
for (method in methods){
t1 <- levelplot(t(array(pre.nonad.mn[pre.nonad.mn$mcp==method,]$TPR,dim=c(16,26))),
                scales=level.scales1,
                colorkey=TRUE,
                at=at.pow,
                col.regions=col.pow,
                xlab="",
                par.settings=list(layout.heights=list(top.padding=-3)),                                
                ylab=method,
                aspect="fill")
t2 <- levelplot(t(array(obs.nonad.mn[obs.nonad.mn$mcp==method,]$TPR,dim=c(16,26))),
                scales=level.scales2,
                colorkey=FALSE,
                at=at.pow,
                col.regions=col.pow,
                xlab="",
                par.settings=list(layout.heights=list(top.padding=-3)),
                ylab=method,
                aspect="fill")
t3 <- levelplot(t(array(bias.nonad.mn[bias.nonad.mn$mcp==method,]$bias,dim=c(16,26))),
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
             plots[["RF"]][[1]],plots[["RF"]][[2]],
             nrow=4,widths=c(7.7,4))
dev.off()

###########################
## BIAS WHEN NONADAPTIVE ##
###########################

pUN <- ggplot(cond.tot[cond.tot$mcp=="UN",],aes(interaction(factor(adaptive),factor(condition)),FPR,fill=condition,alpha=adaptive)) +
  geom_violin(trim=TRUE) +
  geom_hline(aes(yintercept=0.05)) +
  theme(panel.background = element_rect(fill = NA, colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        #legend.position="none",     
        panel.border = element_rect(fill=NA,colour="black")
  ) +
  scale_fill_manual(values = cols) +
  scale_colour_manual(values=c("black",NA)) + 
  ylim(c(0.015,0.08)) +
  labs(x = "",
       y = "False Positive Rate",
       title="Uncorrected")    

pBH <- ggplot(cond.tot[cond.tot$mcp=="BH",],aes(interaction(factor(adaptive),factor(condition)),FDR,fill=condition,alpha=adaptive)) +
  geom_violin(trim=TRUE) +
  geom_hline(aes(yintercept=0.05)) +
    theme(panel.background = element_rect(fill = NA, colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        legend.position="none",     
        panel.border = element_rect(fill=NA,colour="black")
  ) +
  scale_fill_manual(values = cols) +
  scale_colour_manual(values=c("black",NA)) + 
  ylim(c(0.0,0.05)) +
  labs(x = "",
       y = "False Discovery Rate",
       title="FDR control")    

pRFT <- ggplot(cond.mn[cond.mn$mcp=="RF",],aes(x=factor(condition),y=FWER,colour=adaptive,group=adaptive)) +
  geom_line() +
  geom_hline(aes(yintercept=0.05)) +
    theme(panel.background = element_rect(fill = NA, colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none",     
        panel.border = element_rect(fill=NA,colour="black")
  ) +
  scale_colour_manual(values=brewer.pal(2,"Set2")) +
  ylim(c(0.0,0.05)) +
  labs(x = "Condition",
       y = "Familywise error rate",
       title="RFT control")    

pBF <- ggplot(cond.mn[cond.mn$mcp=="BF",],aes(x=factor(condition),y=FWER,colour=adaptive,group=adaptive)) +
  geom_line() +
  geom_hline(aes(yintercept=0.05)) +
    theme(panel.background = element_rect(fill = NA, colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none",     
        panel.border = element_rect(fill=NA,colour="black")
  ) +
  scale_colour_manual(values=brewer.pal(3,"Set2")[1:2]) +
  ylim(c(0.0,0.05)) +
  labs(x = "Condition",
       y = "Familywise error rate",
       title="Bonferroni control")    

empty <- ggplot(data.frame()) + geom_point() + xlim(0, 10) + ylim(0, 100) + theme_minimal() + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank())

dat <- data.frame(a = factor(acts), b = factor(effs),c=cons_list)
leg <- ggplotGrob(
  ggplot(unique(subset(dat, select = a:b)), 
         aes(a, b, colour=cons_list)) + 
    geom_point(size=10, shape=15) +
    theme(panel.background = element_rect(fill = NA, colour = "white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.ticks = element_blank()) +
    scale_colour_manual(labels=rep("",16), 
                        values = cols) +
    labs(x = "percentage active",
         y = "effect size")
)

legend <- g_legend(ggplot(cond.mn[cond.mn$mcp=="BF",],aes(x=factor(condition),y=FWER,colour=adaptive,group=adaptive)) + geom_line() + theme_bw())



pdf(paste(RESDIR,"figures/FPR_adaptive_vs_nonadaptive.pdf",sep=""),width=15,height=17)
grid.arrange(
  arrangeGrob(pUN,pBH,
              arrangeGrob(empty,leg,empty,nrow=3,heights=c(1/5,3/5,1/5)),
              ncol=3,
              widths=c(3/5,3/5,1/5)),
  arrangeGrob(pBF,pRFT,
              arrangeGrob(empty,legend,empty,nrow=3,heights=c(3/7,1/7,3/7)),
              ncol=3,
              widths=c(3/5,3/5,1/5)),
  nrow=2,
  main="Comparison of adaptive and non-adaptive designs \nin the coltrol of false positives"
)
dev.off()

######################################
## REQUIRED VS OBTAINED SAMPLE SIZE ##
######################################

transp = 0.2
koeleurtjes <- alpha(cols,transp)

minmap <- 15
maxmap <- 60
subbias <- 5
x <- c(minmap+subbias,minmap,minmap,   maxmap-subbias, maxmap,maxmap )
y <- c(minmap   ,minmap,minmap+subbias,maxmap,    maxmap,maxmap-subbias )
polygon <- data.frame(x,y)


p <- list()

for(m in 1:4){
  method <- c("BF","UN","RF","BH")[m] 
  p[[m]] <- ggplot() + 
    geom_polygon(data=polygon, mapping=aes(x=x, y=y),fill="gray90") +
    geom_abline(colour="grey50") +
    geom_jitter(data=cond.ad[cond.ad$mcp==method,],
                aes(x=true,
                    y=predicted,
                    group=factor(condition),
                    colour=factor(condition)
                ),
                size=3,
                alpha=1/50
    ) +
    geom_jitter(data=cond.ad.mn[cond.ad.mn$mcp==method,],
                aes(x=true,
                    y=predicted,
                    group=factor(condition),
                    colour=factor(condition)
                ),
                size=5,
                alpha=1,
                shape=18
    ) +
    theme(panel.background = element_rect(fill = NA, colour = 'grey'),
          legend.position="none",        
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    ) +
    scale_colour_manual(values=cols.d) +
    xlim(c(15,60))+  
    ylim(c(15,60))+
    #coord_map(xlim = c(15, 50),ylim = c(15, 50)) + 
    labs(x="True required sample size",
         y = "Predicted required sample size",
         title=method
    )
}


cons_list <- paste("EFFECT:",effs," - ","ACTIVATION SIZE:",acts,sep="")
dat <- data.frame(a = factor(acts), b = factor(effs),c=cons_list)

empty <- ggplot(data.frame()) + geom_point() + xlim(0, 10) + ylim(0, 100) + theme_minimal() + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank())


leg <- ggplotGrob(
  ggplot(unique(subset(dat, select = a:b)), 
         aes(a, b, colour=cons_list)) + 
    geom_point(size=10, shape=15) +
    theme(panel.background = element_rect(fill = NA, colour = "white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.ticks = element_blank()) +
    scale_colour_manual(labels=rep("",16), 
                        values = cols) +
    labs(x = "percentage active",
         y = "effect size")
)

pdf(paste(RESDIR,"figures/sstruereq.pdf",sep=""),width=10,height=9)

grid.arrange(
  arrangeGrob(p[[1]],p[[3]],nrow=2),
  arrangeGrob(p[[2]],p[[4]],nrow=2),
  arrangeGrob(empty,leg,empty,nrow=3,heights=c(2/5,1/5,2/5)),
  ncol=3,
  widths = c(2/5,2/5,1/5))
dev.off()

