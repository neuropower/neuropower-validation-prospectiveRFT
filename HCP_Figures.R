library(RColorBrewer)
library(reshape)
library(plotrix)
library(scales)
library(lattice)
library(latticeExtra)
library(colorRamps)
library(gplots)
library(gridExtra)


RESDIR <- "/Users/Joke/Documents/Onderzoek/Studie_4_propow/ProspectivePowerValidation_Results/power_peak_HCP/power_peak_HCP_15/"
HOMEDIR <- "~/Documents/Onderzoek/Studie_4_propow/ProspectivePower-Validation/"
FIGDIR <- "~/Documents/Onderzoek/Studie_4_propow/ProspectivePower-Paper/Studie_4_v1.4/Figures/"


sims <- 500

setwd(FIGDIR)

pars <- read.table(paste(HOMEDIR,"HCP_paradigms.txt",sep=""))$V1
cons <- read.table(paste(HOMEDIR,"HCP_contrasts.txt",sep=""))$V1
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

# read in results: model estimation
estimation <- array(NA,dim=c(sims,47,6))
range <- (1:sims)
for (p in range){
  file <- paste(RESDIR,"estimation_hcp_",p,".csv",sep="")
  res <- read.csv(file,header=FALSE,na.strings=c("nan"))
  names(res) <- c("c","pi1e","pi1t","pi1c","effe","efft","effc","effec","sde","sdt")
  estimation[p,,] <- data.matrix(res[,c(2,4,7,8,9,10)])
}

# read in results: power estimation
powtru <- array(NA,dim=c(sims,47,50,4))
powpre <- array(NA,dim=c(sims,47,50,4))
for (p in range){
  print(p)
  for(c in 1:47){
    prefile <- paste(RESDIR,"powpre_hcp_",p,"_contrast_",c-1,".csv",sep="")
    if(file.exists(prefile)){
    pre <- read.csv(prefile,na.strings=c("nan"),header=TRUE)
    if(dim(pre)[2]==3){pre$BH=NA}
    pre <- data.frame(pre$UN,pre$BF,pre$RFT,pre$BH)
    powpre[p,c,,] <- data.matrix(pre)
    }
    
    trufile <- paste(RESDIR,"powtru_hcp_",p,"_contrast_",c-1,".csv",sep="")
    if(file.exists(trufile)){
      tru <- read.csv(trufile,na.strings=c("nan"),header=TRUE)
      if(dim(tru)[2]==3){tru$BH=NA}
      tru <- data.frame(tru$UN,tru$BF,tru$RFT,tru$BH)
      powtru[p,c,,] <- as.numeric(as.matrix(tru))   
    }
  }
}


## EVALUATE ESTIMATION PI0

cols.b <- c(brewer.pal(10,"Paired")[seq(1,10,2)],brewer.pal(9,"Greys")[4],brewer.pal(11,"PiYG")[4])
cols.f <- c(brewer.pal(10,"Paired")[seq(2,10,2)],brewer.pal(9,"Greys")[7],brewer.pal(11,"PiYG")[2])


pchs <- c(16,17,18,15,1,2,3,4,5,6,21,22,23,24,25,7,8,9,10,11,12,13,14)
transp <- 0.1
cx <- 0.6
cxav <- 1
cxtx <- 0.5
estimation.av <- apply(estimation,c(2,3),mean,na.rm=TRUE)

pdf(paste(FIGDIR,"FIG_HCP_modelestimation_ss10.pdf",sep=""),width=8,height=5)

par(mar=c(4,4,3,1),oma=c(0,0,0,0))
layout(matrix(c(1,2,3,3),2,2,byrow=TRUE),widths=c(0.5,0.5),heights=c(0.8,0.2))


plot(seq(0,0.8,length=10),seq(0,0.8,length=10),col=NA,xlab="True prevalence of activation",ylab="Estimated prevalence of activation",axes=FALSE,main="Prevalence of activation")
abline(0,1,lwd=1,col="grey50")
box()
axis(1)
axis(2)

for(i in 47:1){
a <- data.frame(estimation[,i,2],estimation[,i,1])
points(a[,1],a[,2],col=alpha(cols.b[pars[i]],transp),pch=pchs[cons[i]],cex=cx)
}

points(estimation.av[,2],estimation.av[,1],col=cols.f[pars],pch=pchs[cons],cex=cxav)




plot(seq(2.7,5,length=10),seq(2.7,5,length=10),col=NA,xlab="True expected effect size",ylab="Estimated expected effect size",axes=FALSE,main="Expected effect size: \n all samples")
abline(0,1,lwd=1,col="grey50")
box()
axis(1)
axis(2)

for(i in 47:1){
  a <- data.frame(estimation[,i,4],estimation[,i,3])
  points(a[,1],a[,2],col=alpha(cols.b[pars[i]],transp),pch=pchs[cons[i]],cex=cx)
}

points(estimation.av[,4],estimation.av[,3],col=cols.f[pars],pch=pchs[cons],cex=cxav)

par(mar=c(0,0,0,0))
plot(1:3, 1:3, col=NA,axes=FALSE,xlab="",ylab="")
legend(1,3,paste(pars[1:10],":",contrast[1:10]),col=cols.f[pars[1:10]],pch=pchs[cons[1:10]],bty="n",cex=cxtx)
legend(1.4,3,paste(pars[11:20],":",contrast[11:20]),col=cols.f[pars[11:20]],pch=pchs[cons[11:20]],bty="n",cex=cxtx)
legend(1.8,3,paste(pars[21:30],":",contrast[21:30]),col=cols.f[pars[21:30]],pch=pchs[cons[21:30]],bty="n",cex=cxtx)
legend(2.2,3,paste(pars[31:40],":",contrast[31:40]),col=cols.f[pars[31:40]],pch=pchs[cons[31:40]],bty="n",cex=cxtx)
legend(2.6,3,paste(pars[41:47],":",contrast[41:47]),col=cols.f[pars[41:47]],pch=pchs[cons[41:47]],bty="n",cex=cxtx)




dev.off()


################################################

effectorder <- order(estimation.av[,4])
  
powpred3D <- apply(powpre,c(2,3,4),mean,na.rm=TRUE)
powtrue3D <- apply(powtru,c(2,3,4),mean,na.rm=TRUE)

powpred3D <- powpred3D[effectorder,,]
powtrue3D <- powtrue3D[effectorder,,]
confull <- confull[effectorder]

methname <- c("Uncorrected","Bonferroni","Random Field Theory","False Discovery Rate")

bias <- powpred3D - powtrue3D
method = 1

col.pow <- colorRampPalette(brewer.pal(9,"YlOrRd"))(100)
at.pow <- seq(0, 1, length.out=length(col.pow)+1)
ckey.pow <- list(at=at.pow, col=col.pow)

col.bias <- colorRampPalette(brewer.pal(11,"RdBu")[c(1:5,6,6,7:11)])(50)
at.bias <- seq(-1, 1, length.out=length(col.bias)-1)
ckey.bias <- list(at=at.bias, col=col.bias)

level.scales1 <- list(y=list(labels=confull,at=1:47,cex=0.5),
					 x=list(labels=seq(15,65,10),at=seq(1,50,10)))
level.scales2 <- list(y=list(labels=rep("",47),at=1:47,cex=0.5),
					 x=list(labels=seq(15,65,10),at=seq(1,50,10)))

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


pdf(paste(FIGDIR,"FIG_HCP_power.pdf",sep=""),width=15,height=17)
grid.arrange(plots[[1]][[1]],plots[[1]][[2]],
             plots[[2]][[1]],plots[[2]][[2]],
             plots[[3]][[1]],plots[[3]][[2]],             
             plots[[4]][[1]],plots[[4]][[2]],
             nrow=4,widths=c(7.7,4))
dev.off()


##########################################
## HOW LARGE SHOULD THE SAMPLE SIZE BE? ##
##########################################

nec <- data.frame(array(NA,dim=c((500*47*4),6)))
names(nec) <- c("pars","con","sim","mcp","pred","true")
k <- 0

range <- 1:sims
for (p in range){
  print(p)
  for(c in 1:47){
    k <- k+1
    rg <- ((k-1)*4+1):((k-1)*4+4)
    nec$pars[rg] <- pars[c]
    nec$con[rg] <- cons[c]
    nec$sim[rg] <- s
    nec$mcp[rg] <- c("UN","BF","RFT","BH")
    
    
    prefile <- paste(RESDIR,"powpre_hcp_",p,"_contrast_",c-1,".csv",sep="")
    if(file.exists(prefile)){
      pre <- read.csv(prefile,na.strings=c("nan"),header=TRUE)
      if(dim(pre)[2]==3){pre$BH=NA}
      pre <- data.frame(pre$UN,pre$BF,pre$RFT,pre$BH)
      ind <- data.frame(which(pre>0.7,arr.ind=TRUE))
      if(nrow(ind) !=0){
        nectab <- ddply(ind,~col,summarise,min=min(row))
        necl <- c(ifelse(length(nectab$min[nectab$col==1])!=0,nectab$min[nectab$col==1],NA),
                  ifelse(length(nectab$min[nectab$col==2])!=0,nectab$min[nectab$col==2],NA),
                  ifelse(length(nectab$min[nectab$col==3])!=0,nectab$min[nectab$col==3],NA),                  
                  ifelse(length(nectab$min[nectab$col==4])!=0,nectab$min[nectab$col==4],NA)
        )
        nec$pred[rg] <- (PILOT:61)[necl]       
      }
    }
    
    trufile <- paste(RESDIR,"powtru_hcp_",p,"_contrast_",c-1,".csv",sep="")
    if(file.exists(trufile)){
      tru <- read.csv(trufile,na.strings=c("nan"),header=TRUE)
      if(dim(tru)[2]==3){tru$BH=NA}
      tru <- data.frame(tru$UN,tru$BF,tru$RFT,tru$BH)
      ind <- data.frame(which(tru>0.7,arr.ind=TRUE))
      if(nrow(ind) !=0){
        nectab <- ddply(ind,~col,summarise,min=min(row))
        necl <- c(ifelse(length(nectab$min[nectab$col==1])!=0,nectab$min[nectab$col==1],NA),
                  ifelse(length(nectab$min[nectab$col==2])!=0,nectab$min[nectab$col==2],NA),
                  ifelse(length(nectab$min[nectab$col==3])!=0,nectab$min[nectab$col==3],NA),                  
                  ifelse(length(nectab$min[nectab$col==4])!=0,nectab$min[nectab$col==4],NA)
        )
        nec$true[rg] <- (PILOT:61)[necl]       
      }
      
    }
  }
}
nec$con <- factor(nec$con)
nec$pars <- factor(nec$pars)

meanss <- ddply(nec,~pars+con+mcp,summarise,pred=mean(pred,na.rm=TRUE),true=mean(true,na.rm=TRUE))
meanss <- meanss[1:188,]

cols.f <- c(brewer.pal(10,"Paired")[seq(2,10,2)],brewer.pal(9,"Greys")[7],brewer.pal(11,"PiYG")[2])


minmap <- 15
maxmap <- 60
subbias <- 5
x <- c(minmap+subbias,minmap,minmap,   maxmap-subbias, maxmap,maxmap )
y <- c(minmap   ,minmap,minmap+subbias,maxmap,    maxmap,maxmap-subbias )
polygon5 <- data.frame(x,y)

subbias <- 10
x <- c(minmap+subbias,minmap,minmap,   maxmap-subbias, maxmap,maxmap )
y <- c(minmap   ,minmap,minmap+subbias,maxmap,    maxmap,maxmap-subbias )
polygon10 <- data.frame(x,y)



p <- list()

for(m in 1:4){
  method <- c("BF","UN","RFT","BH")[m] 
  p[[m]] <- ggplot() + 
    geom_polygon(data=polygon10, mapping=aes(x=x, y=y),fill="gray95") +
    geom_polygon(data=polygon5, mapping=aes(x=x, y=y),fill="gray90") +
    geom_abline(colour="grey50") +
    geom_jitter(data=nec[nec$mcp==method,],
                aes(x=true,
                    y=pred,
                    group=interaction(pars,con),
                    colour=pars,
                    shape=con
                ),
                size=1,
                alpha=1/50
    ) +
    geom_jitter(data=meanss[meanss$mcp==method,],
                aes(x=true,
                    y=pred,
                    group=interaction(pars,con),
                    colour=pars,
                    shape=con
                ),
                size=3,
                alpha=1,
    ) +
    theme(panel.background = element_rect(fill = NA, colour = 'grey'),
          legend.position="none",        
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    ) +
    scale_colour_manual(values=cols.f) +
    scale_shape_manual(values=1:100) +
    xlim(c(15,maxmap))+  
    ylim(c(15,maxmap))+
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


pdf(paste(FIGDIR,"FIG_HCP_sscalc.pdf",sep=""),width=10,height=9)

grid.arrange(
  arrangeGrob(p[[1]],p[[3]],nrow=2),
  arrangeGrob(p[[2]],p[[4]],nrow=2),
  arrangeGrob(leg),
  ncol=3,
  widths = c(2/6,2/6,2/6)
)
dev.off()

