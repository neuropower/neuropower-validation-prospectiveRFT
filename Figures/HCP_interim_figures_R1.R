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
HOMEDIR <- "/Users/Joke/Documents/Onderzoek/ProjectsOngoing/Power/Validation/"
FIGDIR <- "/Users/Joke/Documents/Onderzoek/ProjectsOngoing/Power/ValidationFigures/"

pars <- read.table(paste(HOMEDIR,"HcpFiles/HCP_paradigms.txt",sep=""))$V1
cons <- read.table(paste(HOMEDIR,"HcpFiles/HCP_contrasts.txt",sep=""))$V1
contrast <- c(
  "Random",
  "Theory of Mind",
  "ToM - Random",
  "Faces",
  "Shapes",
  "Faces - Shapes",
  "Match",
  "Relational",
  "Match-Relational",
  "Math",
  "Story",
  "Story-Math",
  "Punish",
  "Reward",
  "Punish-Reward",
  "2back-body",
  "2back-face",
  "2back-place",
  "2back-tool",
  "0back-body",
  "0back-face",
  "0back-place",
  "0back-tool",
  "2back",
  "0back",
  "2back-0back",
  "body",
  "face",
  "place",
  "tool",
  "body-average",
  "face-average",
  "place-average",
  "tool-average",
  "Cue",
  "Left foot",
  "Left hand",
  "Right foot",
  "Right hand",
  "Toes",
  "Average",
  "Cue - average",
  "Left foot - average",
  "Left hand - average",
  "Right foot - average",
  "Right hand - average",
  "Toes - average"
)
confull <- paste(pars,":",contrast)

u <- 2.3

## FUNCTIONS

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

# read in results

res.nonad <- read.table(paste(RESDIR,"estimation_HCP_predictive_",u,"_RFT.csv",sep=""),header=TRUE,sep=",")
pre.nonad <- read.table(paste(RESDIR,"prediction_HCP_predictive_",u,"_RFT.csv",sep=""),header=TRUE,sep=",")
obs.nonad <- read.table(paste(RESDIR,"true_HCP_predictive_",u,"_RFT.csv",sep=""),header=TRUE,sep=",")
cond.nonad <- read.table(paste(RESDIR,"conditional_HCP_predictive_",u,"_RFT.csv",sep=""),header=TRUE,sep=",")

names(cond.nonad) <- names(cond.ad) <- c("ind","mcp","predicted","power","true","simulation","condition")
cond.nonad$condition <- factor(cond.nonad$condition)

# take mean of power predictions over simulations

res.nonad.mn <- ddply(res.nonad,
                      ~contrast,
                      summarise,
                      pi1e = mean(pi1e),
                      pi1t = mean(pi1t),
                      pi1c = mean(pi1c),
                      ese = mean(ese),
                      est = mean(est),
                      esc = mean(esc),
                      esexp = mean(esexp),
                      sde = mean(sde),
                      sdt = mean(sdt),
                      bumpar = mean(bumpar)
)


pre.nonad.mn <- ddply(pre.nonad,
                   ~subjects+condition+mcp,
                   summarise,
                   TPR=mean(power)
)
obs.nonad.mn <- ddply(obs.nonad,
                   ~subjects+condition+mcp,
                   summarise,
                   TPR=mean(power)
)
obs.nonad.mn$condition <- factor(obs.nonad.mn$condition)

bias.nonad.mn <- pre.nonad.mn[,1:3]
bias.nonad.mn$bias <- pre.nonad.mn$TPR - obs.nonad.mn$TPR


cond.mn <- ddply(cond.nonad,
                    ~condition + mcp,
                    summarise,
                    TPR = mean(power,na.rm=TRUE),
                    predicted = mean(predicted,na.rm=TRUE),
                    true = mean(true,na.rm=TRUE)
                    )


cond.nonad$par <- cond.nonad$con <- NA
cond.mn$par <- cond.mn$con <- NA
for (i in 1:47){
  cond.nonad$par[cond.nonad$condition==i] <- pars[i]
  cond.nonad$con[cond.nonad$condition==i] <- cons[i]
  cond.mn$par[cond.mn$condition==i] <- pars[i]
  cond.mn$con[cond.mn$condition==i] <- cons[i]
}

cond.nonad$con <- factor(cond.nonad$con)
cond.mn$con <- factor(cond.mn$con)

########################################
## EVALUATE ESTIMATION MODEL ADAPTIVE ##
########################################

cols.b <- c(brewer.pal(10,"Paired")[seq(1,10,2)],brewer.pal(9,"Greys")[4],brewer.pal(11,"PiYG")[4])
cols.f <- c(brewer.pal(10,"Paired")[seq(2,10,2)],brewer.pal(9,"Greys")[7],brewer.pal(11,"PiYG")[2])


pchs <- c(16,17,18,15,1,2,3,4,5,6,21,22,23,24,25,7,8,9,10,11,12,13,14)
transp <- 0.1
cx <- 0.6
cxav <- 1
cxp <- 0.5
cxtx <- 0.5


pdf(paste(FIGDIR,"hcp_modelestimation.pdf",sep=""),width=9,height=5)

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
points(res.nonad$pi1c,
       res.nonad$pi1e,
       col=alpha(cols.b[pars],transp),
       pch=pchs[cons],
       cex=cxp
       )


points(res.nonad.mn$pi1c,
       res.nonad.mn$pi1e,
       col=cols.f[pars],
       pch=pchs[cons],
       cex=cxav)


# plot 2: model estimation

plot(seq(2.3,5,length=10),
     seq(2.3,5,length=10),
     col=NA,
     xlab="True expected effect size",
     ylab="Estimated expected effect size",
     axes=FALSE,
     main="Scatter plot of effect size")

abline(0,1,lwd=1,col="grey50")
box();axis(1);axis(2)

points(res.nonad$esc,
       res.nonad$ese,
       col=alpha(cols.b[pars],transp),
       pch=pchs[cons],
       cex=cxp
)


points(res.nonad.mn$esc,
       res.nonad.mn$ese,
       col=cols.f[pars],
       pch=pchs[cons],
       cex=cxav)



# legend

par(mar=c(0,0,0,0))
plot(1:3, 1:3, col=NA,axes=FALSE,xlab="",ylab="")
legend(1,3,paste(pars[1:10],":",contrast[1:10]),col=cols.f[pars[1:10]],pch=pchs[cons[1:10]],bty="n",cex=cxtx)
legend(1.4,3,paste(pars[11:20],":",contrast[11:20]),col=cols.f[pars[11:20]],pch=pchs[cons[11:20]],bty="n",cex=cxtx)
legend(1.8,3,paste(pars[21:30],":",contrast[21:30]),col=cols.f[pars[21:30]],pch=pchs[cons[21:30]],bty="n",cex=cxtx)
legend(2.2,3,paste(pars[31:40],":",contrast[31:40]),col=cols.f[pars[31:40]],pch=pchs[cons[31:40]],bty="n",cex=cxtx)
legend(2.6,3,paste(pars[41:47],":",contrast[41:47]),col=cols.f[pars[41:47]],pch=pchs[cons[41:47]],bty="n",cex=cxtx)

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
level.scales1 <- list(y=list(labels=confull,at=1:47,cex=0.5),
                      x=list(labels=seq(15,64,10),at=seq(1,49,10)))
level.scales2 <- list(y=list(labels=rep("",47),at=1:47,cex=0.5),
                      x=list(labels=seq(15,64,10),at=seq(1,49,10)))


plots <- list(0)
for (method in methods){
  
  effectorder <- order(res.nonad.mn$est)
  
  powpred3D <- t(array(pre.nonad.mn[pre.nonad.mn$mcp==method,]$TPR,dim=c(46,50)))
  mm <- ifelse(method=="RF","RFT",method)
  powtrue3D <- t(array(obs.nonad.mn[obs.nonad.mn$mcp==method,]$TPR,dim=c(46,50)))
  
  powpred3D <- powpred3D[,]
  powtrue3D <- powtrue3D[,]
  confull <- confull[]
  bias3D <- powpred3D-powtrue3D
  
  
t1 <- levelplot(powpred3D,
                scales=level.scales1,
                colorkey=TRUE,
                at=at.pow,
                col.regions=col.pow,
                xlab="",
                par.settings=list(layout.heights=list(top.padding=-3)),                                
                ylab=method,
                aspect="fill")
t2 <- levelplot(powtrue3D,
                scales=level.scales2,
                colorkey=FALSE,
                at=at.pow,
                col.regions=col.pow,
                xlab="",
                par.settings=list(layout.heights=list(top.padding=-3)),
                ylab=method,
                aspect="fill")
t3 <- levelplot(bias3D,
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


pdf(paste(FIGDIR,"powerpredictions_HCP.pdf",sep=""),width=15,height=17)
grid.arrange(plots[["UN"]][[1]],plots[["UN"]][[2]],
             plots[["BH"]][[1]],plots[["BH"]][[2]],
             plots[["BF"]][[1]],plots[["BF"]][[2]],             
             plots[["RF"]][[1]],plots[["RF"]][[2]],
             nrow=4,widths=c(7.7,4))
dev.off()


######################################
## REQUIRED VS OBTAINED SAMPLE SIZE ##
######################################

transp = 0.2
koeleurtjes <- alpha(cols.f,transp)

minmap <- 15
maxmap <- 80
subbias <- 5
x <- c(minmap+subbias,minmap,minmap,   maxmap-subbias, maxmap,maxmap )
y <- c(minmap   ,minmap,minmap+subbias,maxmap,    maxmap,maxmap-subbias )
polygon <- data.frame(x,y)

subbias <- 10
x <- c(minmap+subbias,minmap,minmap,   maxmap-subbias, maxmap,maxmap )
y <- c(minmap   ,minmap,minmap+subbias,maxmap,    maxmap,maxmap-subbias )
polygon10 <- data.frame(x,y)




cond.nonad$par <- factor(cond.nonad$par)
cond.mn$par <- factor(cond.mn$par)

p <- list()

for(m in 1:4){
  method <- c("BF","UN","RF","BH")[m] 
  p[[m]] <- ggplot() + 
    geom_polygon(data=polygon10, mapping=aes(x=x, y=y),fill="gray95") +
    geom_polygon(data=polygon, mapping=aes(x=x, y=y),fill="gray90") +
    geom_abline(colour="grey50") +
    geom_jitter(data=cond.nonad[cond.nonad$mcp==method,],
                aes(x=true,
                    y=predicted,
                    group=interaction(par,con),
                    colour=par,
                    shape=con
                ),
                 size=1,
                alpha=1/10
    ) +
    geom_jitter(data=cond.mn[cond.mn$mcp==method,],
                aes(x=true,
                    y=predicted,
                    group=interaction(par,con),
                    colour=par,
                    shape=con
                ),
                size=2,
                alpha=1,
    ) +
    
    theme(panel.background = element_rect(fill = NA, colour = 'grey'),
          legend.position="none",        
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    ) +
    scale_colour_manual(values=cols.f) +
    scale_shape_manual(values=1:100) +
    xlim(c(minmap,maxmap))+  
    ylim(c(minmap,maxmap))+
    #coord_map(xlim = c(15, 50),ylim = c(15, 50)) + 
    labs(x="True required sample size",
         y = "Predicted required sample size",
         title=method
    )
}


legcol <- c()
legsh <- c()
for(i in 1:47){
  legcol[i] <- cols.f[pars[i]]
  legsh[i] <- (1:100)[cons[i]]
}

dat <- data.frame(a=1,b=confull,c=confull,legcol=factor(legcol),legsh=factor(legsh))
dat$b <- dat$c <- factor(dat$b,levels=dat$b[order(1:47)])

leg <- ggplotGrob(  
  ggplot(dat,aes(a,b,group=c,colour=c,shape=c)) + 
    geom_point() +
    scale_colour_manual(values=legcol,labels=confull) +
    scale_shape_manual(values=legsh) +
    theme(panel.background = element_rect(fill = NA, colour = "white"),
          legend.position="none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank()
    ) +
    labs(x="",y="")
)



pdf(paste(FIGDIR,"sstruereq_hcp.pdf",sep=""),width=10,height=9)
grid.arrange(
  arrangeGrob(p[[1]],p[[3]],nrow=2),
  arrangeGrob(p[[2]],p[[4]],nrow=2),
  arrangeGrob(leg),
  ncol=3,
  widths = c(1/3,1/3,1/3))
dev.off()

