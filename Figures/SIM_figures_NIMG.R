library(RColorBrewer)
library(reshape)
library(plotrix)
library(scales)
library(lattice)
library(latticeExtra)
library(colorRamps)
library(gridExtra)
library(data.table)
library(ggplot2)
library(plyr)
library(gridExtra)

args <- commandArgs(TRUE)
RESDIR <- args[1]
HOMEDIR <- args[2]
FIGDIR <- args[3]
u <- args[4]
print(args)

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

# result files name
res.nonad.file = paste(RESDIR,"/estimation_SIM_predictive_",u,"_RFT.csv",sep="")
pre.nonad.file = paste(RESDIR,"/prediction_SIM_predictive_",u,"_RFT.csv",sep="")
true.nonad.file = paste(RESDIR,"/true_SIM_predictive_",u,"_RFT.csv",sep="")
cond.nonad.file = paste(RESDIR,"/conditional_SIM_predictive_",u,"_RFT.csv",sep="")

# read results

res.nonad <- read.table(res.nonad.file,header=TRUE,sep=",",dec=".")
names(res.nonad)[1] <- "condition"
res.nonad$condition <- res.nonad$condition+1
pre.nonad <- read.table(pre.nonad.file,header=TRUE,sep=",")
obs.nonad <- read.table(true.nonad.file,header=TRUE,sep=",")
cond.nonad <- read.table(cond.nonad.file,header=TRUE,sep=",")
cond.nonad$adaptive <- "predictive"
names(cond.nonad)[1] <- "ind"


# average over simulations

res.nonad.mn <- ddply(res.nonad,
                      ~condition,
                      summarise,
                      pi1e = mean(pi1e,na.rm=TRUE),
                      pi1t = mean(pi1t,na.rm=TRUE),
                      ese = mean(ese,na.rm=TRUE),
                      est = mean(est,na.rm=TRUE),
                       esexp = mean(esexp,na.rm=TRUE),
                      sde = mean(sde,na.rm=TRUE),
                      sdt = mean(sdt,na.rm=TRUE),
                      bumpar = mean(bumpar,na.rm=TRUE)
)


pre.nonad.mn <- ddply(pre.nonad,
                   ~subjects+condition+mcp,
                   summarise,
                   TPR=mean(power,na.rm=TRUE)
)

obs.nonad.mn <- ddply(obs.nonad,
                   ~subjects+condition+mcp,
                   summarise,
                   TPR=mean(power,na.rm=TRUE))
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

pdf(paste(FIGDIR,"FIG_SIM_modelestimation.pdf",sep=""),width=8,height=5)

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

plot(seq(2,6,length=10),
     seq(2,6,length=10),
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

min <- 15
max <- 60
size <- max-min+1
# parameters for text
methname <- c("UN","FDR","FWER-BF","FWER-RFT")
methods <- c("UN","BH","BF","RF")

# color ramps for power and bias
col.pow <- colorRampPalette(brewer.pal(9,"YlOrRd"))(100)
at.pow <- seq(0, 1, length.out=length(col.pow)+1)
ckey.pow <- list(at=at.pow, col=col.pow)

col.bias <- colorRampPalette(brewer.pal(11,"RdBu")[c(1:5,6,6,7:11)])(50)
at.bias <- seq(-1,1, length.out=length(col.bias)-1)
ckey.bias <- list(at=at.bias, col=col.bias)

# define scales for x and y axis
level.scales1 <- list(y=list(labels=cons_list,at=1:16,cex=0.7),
                      x=list(labels=seq(min,max,by=5),at=seq(1,size,by=5)))
level.scales2 <- list(y=list(labels=rep("",16),at=1:16,cex=0.7),
                      x=list(labels=seq(min,max,by=5),at=seq(1,size,by=5)))

plots <- list(0)
x <- 0
for (method in methods){
  x <- x+1
  methn <- methname[x]

  t1 <- levelplot(t(array(pre.nonad.mn[pre.nonad.mn$mcp==method,]$TPR,dim=c(16,size))),
                  scales=level.scales1,
                  colorkey=TRUE,
                  at=at.pow,
                  col.regions=col.pow,
                  xlab="",
                  par.settings=list(layout.heights=list(top.padding=-3)),
                  ylab=methn,
                  aspect="fill")
  t2 <- levelplot(t(array(obs.nonad.mn[obs.nonad.mn$mcp==method,]$TPR,dim=c(16,size))),
                  scales=level.scales2,
                  colorkey=FALSE,
                  at=at.pow,
                  col.regions=col.pow,
                  xlab="",
                  par.settings=list(layout.heights=list(top.padding=-3)),
                  ylab=methn,
                  aspect="fill")
  t3 <- levelplot(t(array(bias.nonad.mn[bias.nonad.mn$mcp==method,]$bias,dim=c(16,size))),
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


pdf(paste(FIGDIR,"FIG_SIM_power_ss15.pdf",sep=""),width=15,height=17)
grid.arrange(plots[["UN"]][[1]],plots[["UN"]][[2]],
             plots[["BH"]][[1]],plots[["BH"]][[2]],
             plots[["BF"]][[1]],plots[["BF"]][[2]],
             plots[["RF"]][[1]],plots[["RF"]][[2]],
             nrow=4,widths=c(7.7,4))
dev.off()


######################################
## REQUIRED VS OBTAINED SAMPLE SIZE ##
######################################

transp = 0.5
koeleurtjes <- alpha(cols,transp)

minmap <- 15
maxmap <- 60
subbias <- 5
x <- c(minmap+subbias,minmap,minmap,   maxmap-subbias, maxmap,maxmap )
y <- c(minmap   ,minmap,minmap+subbias,maxmap,    maxmap,maxmap-subbias )
polygon <- data.frame(x,y)


p <- list()

for(m in 1:4){
  methname <- c("UN","FDR","FWER-BF","FWER-RFT")
  method <- c("UN","BH","BF","RF")[m]
  p[[m]] <- ggplot() +
    geom_polygon(data=polygon, mapping=aes(x=x, y=y),fill="gray90") +
    geom_abline(colour="grey50") +
    geom_jitter(data=cond.nonad[cond.nonad$mcp==method,],
                aes(x=true,
                    y=predicted,
                    group=factor(condition),
                    colour=factor(condition)
                ),
                size=3,
                alpha=1/50
    ) +
    geom_jitter(data=cond.nonad.mn[cond.nonad.mn$mcp==method,],
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
    scale_colour_manual(values=cols) +
    xlim(c(minmap,maxmap))+
    ylim(c(minmap,maxmap))+
    #coord_map(xlim = c(15, 50),ylim = c(15, 50)) +
    labs(x="True required sample size",
         y = "Predicted required sample size",
         title=methname[m]
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
    labs(x = "signal extent",
         y = "signal height")
)

pdf(paste(FIGDIR,"FIG_SIM_sscalc.pdf",sep=""),width=10,height=9)

grid.arrange(
  arrangeGrob(p[[1]],p[[3]],nrow=2),
  arrangeGrob(p[[2]],p[[4]],nrow=2),
  arrangeGrob(empty,leg,empty,nrow=3,heights=c(2/5,1/5,2/5)),
  ncol=3,
  widths = c(2/5,2/5,1/5))
dev.off()
