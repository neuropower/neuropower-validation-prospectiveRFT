library(RColorBrewer)
library(reshape)
library(plotrix)
library(scales)
library(lattice)
library(latticeExtra)
library(colorRamps)
library(gplots)
library(gridExtra)

PILOT <- 15
subrange <- 60-PILOT

RESDIR <- paste("/Users/Joke/Documents/Onderzoek/Studie_4_propow/ProspectivePowerValidation_Results/power_peak_SIM/SIM_",PILOT,"/",sep="")
HOMEDIR <- "~/Documents/Onderzoek/Studie_4_propow/ProspectivePowerValidation/"
FIGDIR <- "~/Documents/Onderzoek/Studie_4_propow/ProspectivePower-Paper/Studie_4_v1.4/Figures/"

setwd(FIGDIR)

effs <- rep(c(0.5,1,1.5,2),each=4)
acts <- rep(c(2,4,6,8),4)
effn <- rep(c("half","one","onehalf","two"),each=4)
cons <- data.frame(effs,acts)

# read in results: model estimation
estimation <- array(NA,dim=c(500,16,10))
range <- (1:500)[-c(91)]
for (p in range){
  file <- paste(RESDIR,"estimation_sim_",p,".csv",sep="")
  res <- read.csv(file,header=FALSE,na.strings=c("nan"))
  names(res) <- c("eff","p","pi1e","pi1t","effe","efft","effex","sde","sdt","a")
  estimation[p,,] <- data.matrix(res)
}


# read in results: power estimation
powpred <- array(NA,dim=c(500,4,4,subrange,4)) 
powtrue <- array(NA,dim=c(500,4,4,subrange,4)) 

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
      true$RFT <- true$RFT_TP/(true$RFT_TP+true$RFT_FN)
      true$BH <- true$BH_TP/(true$BH_TP+true$BH_FN)
      true$BF <- true$BF_TP/(true$BF_TP+true$BF_FN)      
      true$UN <- true$UN_TP/(true$UN_TP+true$UN_FN)
      true <- data.frame(true$BF,true$UN,true$RFT,true$BH)
      true[true=="NaN"] <- 0
      powtrue[s,p,e,,] <-  data.matrix(true)
    }
  }
}


powpred.av <- apply(powpred,c(2,3,4,5),mean,na.rm=TRUE)
powtrue.av <- apply(powtrue,c(2,3,4,5),mean,na.rm=TRUE)

powpred3D <- array(NA,dim=c(16,subrange,4))
powtrue3D <- array(NA,dim=c(16,subrange,4))
k <- 0
for(e in 1:4){
  for(p in 1:4){
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


##########################################
## HOW LARGE SHOULD THE SAMPLE SIZE BE? ##
##########################################

PILOT <- 15

  nec <- data.frame(array(NA,dim=c(32000,6)))
  names(nec) <- c("eff","ac","sim","mcp","pred","true")
  k <- 0
  for(s in range){
    print(s)
    for(p in 1:4){
      for(e in 1:4){
        per <- c(2,4,6,8)[p]
        eff <- c("half","one","onehalf","two")[e]
        
        k <- k+1
        rg <- ((k-1)*4+1):((k-1)*4+4)
        nec$eff[rg] <- e
        nec$ac[rg] <- p
        nec$sim[rg] <- s
        nec$mcp[rg] <- c("BF","UN","RFT","BH")
        
        pred <- read.table(paste(RESDIR,"powpre_sim_",s,"_w_",per,"_e_",eff,".csv",sep=""),sep=",",dec=".",header=TRUE)
        if(is.null(pred$BH)){pred$BH <- 0}
        pred <- data.frame(pred$BF,pred$UN,pred$RFT,pred$BH)           
        ind <- data.frame(which(pred>0.7,arr.ind=TRUE))
        if(nrow(ind) !=0){
        nectab <- ddply(ind,~col,summarise,min=min(row))
        necl <- c(ifelse(length(nectab$min[nectab$col==1])!=0,nectab$min[nectab$col==1],NA),
                  ifelse(length(nectab$min[nectab$col==2])!=0,nectab$min[nectab$col==2],NA),
                  ifelse(length(nectab$min[nectab$col==3])!=0,nectab$min[nectab$col==3],NA),                  
                  ifelse(length(nectab$min[nectab$col==4])!=0,nectab$min[nectab$col==4],NA)
        )
        nec$pred[rg] <- (PILOT:60)[necl]
        
        }
        
        true <- read.table(paste(RESDIR,"powtru_sim_",s,"_w_",per,"_e_",eff,".csv",sep=""),sep=",",dec=".",header=TRUE)
        RFT <- true$RFT_TP/(true$RFT_TP+true$RFT_FN)
        BH <- true$BH_TP/(true$BH_TP+true$BH_FN)
        BF <- true$BF_TP/(true$BF_TP+true$BF_FN)      
        UN <- true$UN_TP/(true$UN_TP+true$UN_FN)
        true <- data.frame(BF,UN,RFT,BH)
        ind <- data.frame(which(true>0.7,arr.ind=TRUE))
        if(nrow(ind) !=0){
          nectab <- ddply(ind,~col,summarise,min=min(row))
          necl <- c(ifelse(length(nectab$min[nectab$col==1])!=0,nectab$min[nectab$col==1],NA),
                    ifelse(length(nectab$min[nectab$col==2])!=0,nectab$min[nectab$col==2],NA),
                    ifelse(length(nectab$min[nectab$col==3])!=0,nectab$min[nectab$col==3],NA),                  
                    ifelse(length(nectab$min[nectab$col==4])!=0,nectab$min[nectab$col==4],NA)
          )
          
          nec$true[rg] <- (PILOT:60)[necl]
        }
        
      }
    }
  }
  
}


meanss <- ddply(nec,~ac+eff+mcp,summarise,pred=mean(pred),true=mean(true))
meanss <- meanss[1:64,]

koeleurtjes <- alpha(cols.d,transp)

nec$eff <- factor(nec$eff)
nec$ac <- factor(nec$ac)

minmap <- 15
maxmap <- 60
subbias <- 5
x <- c(minmap+subbias,minmap,minmap,   maxmap-subbias, maxmap,maxmap )
y <- c(minmap   ,minmap,minmap+subbias,maxmap,    maxmap,maxmap-subbias )
polygon <- data.frame(x,y)


p <- list()

for(m in 1:4){
  method <- c("BF","UN","RFT","BH")[m] 
 p[[m]] <- ggplot() + 
  geom_polygon(data=polygon, mapping=aes(x=x, y=y),fill="gray90") +
  geom_abline(colour="grey50") +
  geom_jitter(data=nec[nec$mcp==method,],
             aes(x=true,
                 y=pred,
                 group=interaction(ac,eff),
                 colour=interaction(ac,eff)
             ),
             size=3,
             alpha=1/50
             ) +
  geom_jitter(data=meanss[meanss$mcp==method,],
             aes(x=true,
                 y=pred,
                 group=interaction(ac,eff),
                 colour=interaction(ac,eff)
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

pdf(paste(FIGDIR,"FIG_SIM_sscalc.pdf",sep=""),width=10,height=9)

grid.arrange(
    arrangeGrob(p[[1]],p[[3]],nrow=2),
    arrangeGrob(p[[2]],p[[4]],nrow=2),
  arrangeGrob(empty,leg,empty,nrow=3,heights=c(2/5,1/5,2/5)),
  ncol=3,
  widths = c(2/5,2/5,1/5))
dev.off()
