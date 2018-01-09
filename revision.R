require(plyr)
require(ggplot2)
require(reshape2)
require(dplyr)
require(tidyr)
require(gridExtra)
require(mvtnorm)
disease.color <- "#CC6666"
nondisease.color <- "#9999CC"

## a. visit schedule
dat <- read.csv('final.csv')
dat <- subset(dat,select=c('ID','DOV','DOB'))
dat$DOV <- strptime(dat$DOV,format='%m/%d/%Y')
dat$DOV <- as.POSIXct(dat$DOV)
dat$DOB <- strptime(dat$DOB,format='%m/%d/%Y')
dat$DOB <- as.POSIXct(dat$DOB)
## plot(y=rep(1,50),x=dat$DOV[dat$ID==1])
## ggplot(dat,aes(y=ID,x=DOV))+geom_point()
dat2 <- dat %>% group_by(ID) %>% summarize(n.visits=n(),mean.bw=round(mean(diff(DOV))),sd.bw=round(sd(diff(DOV))),first.visit=min(DOV)) %>% as.data.frame
## sink('1026.txt')
## head(subset(dat2,dat2$n.visits>20),10)
## sink()
dat <- left_join(x=dat,by='ID',y=dat2)
plt <- ggplot(subset(dat,n.visits>20),aes(y=ID,x=DOV))+geom_point(size=.5)+geom_point(aes(x=DOB,y=ID),shape=3)+theme_classic()+scale_shape_identity()+labs(x='date of visit',y='')
plt
## ggsave('081817a.png',plt)
## ggsave('revision_a.png',plt)

dat <- read.csv('final.csv')
dat <- subset(dat,select=c('ID','DOV','DOB'))
dat$DOV <- strptime(dat$DOV,format='%m/%d/%Y')
dat$DOV <- as.POSIXct(dat$DOV)
dat$DOB <- strptime(dat$DOB,format='%m/%d/%Y')
dat$DOB <- as.POSIXct(dat$DOB)
## plot(y=rep(1,50),x=dat$DOV[dat$ID==1])
## ggplot(dat,aes(y=ID,x=DOV))+geom_point()
dat2 <- dat %>% group_by(ID) %>% summarize(n.visits=n(),mean.bw=round(mean(diff(DOV))),sd.bw=round(sd(diff(DOV))),first.visit=min(DOV)) %>% as.data.frame
split.data <- split(dat,dat$ID)
diffs <- sapply(split.data,function(df)as.numeric(diff(df$DOV),units='days'))
diffs <- unlist(diffs)
summary(as.numeric(diffs))
visits <- sapply(split.data,function(df)df$n.visits[1])
summary(visits)

## b. correlation heatmap
dat <- read.csv('final.csv')
dat$blip <- as.numeric(dat$VL>1e3)
dat <- subset(dat,select=c('ID','Blip_YN50','CD4','CD4P','Visits','blip','Sex','Transmission','trtm','AgeDiag','TimeDiag','AGE','NumberVisits'))
dat <- subset(dat,select=c('ID','CD4','CD4P','Visits','blip','Blip_YN50','NumberVisits'))
dat <- rename(dat,Visit=Visits)
dat$CD4[dat$CD4>6000] <- NA
dat$CD4P[dat$CD4P>=99] <- NA
dat$Visit <- unlist(sapply(rle(dat$ID)$lengths,function(x)seq_len(x)))
dat$blip <- factor(dat$blip)
dat$ID <- factor(dat$ID)
dat$Blip_YN50 <- factor(as.numeric(dat$Blip_YN50)-1)
dat <- droplevels(dat)

group_by(dat,ID) %>% summarize(blip=(function(x)x[1])(Blip_YN50)) %>% as.data.frame %>% select(blip) %>% table()

## ggplot(dat,aes(x=CD4P,group=Blip_YN50))+geom_histogram(aes(fill=Blip_YN50),alpha=.5,color='black',position='identity')
## ggsave('011517a.pdf')


## plt <- ggplot(subset(dat,NumberVisits>0),aes(x=Visit,y=CD4,group=factor(ID),color=Blip_YN50))+geom_line()+stat_smooth(method='lm',se=F)+labs(title='CD4 over time, including regression lines')
## plt

## dat$Blip_YN50 <- as.numeric(dat$Blip_YN50)-1


## sample correlations
## by.visit <- subset(dat,NumberVisits>=10 & Visit<=10)
## by.visit$Visit <- paste0('vis',by.visit$Visit)
## by.visit <- dcast(by.visit,ID+Blip_YN50 ~ Visit,value.var='CD4P')
## by.visit <- na.omit(by.visit)

by.visit <- subset(dat,NumberVisits>=40 & Visit<=40)
#by.visit$Visit <- paste0('vis',by.visit$Visit)
by.visit.CD4 <- dcast(by.visit, ID ~ Visit,value.var='CD4')
by.visit.CD4P <- dcast(by.visit, ID ~ Visit,value.var='CD4P')
by.visit.CD4 <- na.omit(by.visit.CD4)[,-1]
by.visit.CD4P <- na.omit(by.visit.CD4P)[,-1]
ggplot(data.frame(test='CD4',melt(abs(cor(by.visit.CD4)))), aes(x=Var1,y=Var2,fill=value))+geom_tile()+theme_classic()+scale_fill_gradient2(low='white',high='black')+labs(x='',y='')
## plt
ggsave('081817b.png')



## c. faceted plots of individual cd4s
set.seed(2)
min.visits <- 25
dat <- read.csv('final.csv')
dat <- subset(dat,select=c('ID','blip','CD4','CD4P','Visits'))
dat <- dplyr::rename(dat,Visit=Visits)
dat$CD4[dat$CD4>6000] <- NA
dat$CD4P[dat$CD4P>=99] <- NA
dat <- droplevels(dat)
dat <- na.omit(dat)
##dat$ID <- factor(dat$ID)
n.visits <- group_by(dat,ID) %>% summarize(n.visits=n()) %>% filter(n.visits>=min.visits)
dat <- right_join(dat,n.visits,'ID')

dat$Visit <- unlist(sapply(rle(dat$ID)$lengths,function(x)seq_len(x)))

dat.CD4 <- subset(dat,Visit>=n.visits-min.visits+1,c('ID','Visit','CD4'))
dat.CD4$Visit <- unlist(lapply(rle(dat.CD4$ID)$lengths,function(x)seq_len(x)))
dat.blip <- subset(dat,Visit>=n.visits-min.visits+1,c('ID','Visit','blip'))
dat.blip$Visit <- unlist(lapply(rle(dat.blip$ID)$lengths,function(x)seq_len(x)))
dat.CD4 <- spread(dat.CD4,Visit,CD4)
dat.blip <- spread(dat.blip,Visit,blip)



gg.df <- gather(dat.CD4,visit,CD4,-ID) %>% arrange(ID)
gg.df <- left_join(gg.df,gather(dat.blip,visit,blip,-ID),by=c('ID','visit'))
gg.df$visit <- as.numeric(gg.df$visit)
gg.df$blip <- factor(gg.df$blip,labels=c('negative','positive'))
ggplot(filter(gg.df,ID%in%sample(unique(gg.df$ID),9)), aes(x=visit,y=CD4))+geom_point(aes(shape=blip),size=3)+geom_line()+facet_wrap(~ID)+theme_classic()+scale_shape_manual(values=c(1,3),name='blip status')



## ggsave('revision_c.png')


## d. misspecified model--plot of roc curves -- see 081817d.R

## e. misspecified model--table of coverages/biaes -- see 081817e.R

## f. synthetic data -- combined plots

source('ragon.R')
load('simulations/f.RData')
gg.df <- rbind(attr(plt.roc1,'gg.df'),attr(plt.roc3,'gg.df'))
gg.df$group <- rep(c('roc1','roc3'),c(nrow(attr(plt.roc1,'gg.df')),nrow(attr(plt.roc3,'gg.df'))))
plts.roc1.roc3 <- ggplot(gg.df,aes(x=fpr,y=tpr,group=group))+geom_path(aes(linetype=group))+geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.5,fill='lightgray') + geom_abline(slope=1,col='lightgray',linetype=2) +
    theme_classic() + labs(x='FPR',y='TPR') + scale_linetype_discrete(name = 'ROC type',labels=c(expression(ROC[1]),expression(ROC[3])))
plts.roc2 <- do.call(grid.arrange,list(grobs=plts.roc2))
## png('revision_f.png',width=1920,height=1200)
grid.arrange(plts.roc1.roc3,plts.roc2,nrow=1)
dev.off()
## ggsave('revision.k.png',arrangeGrob(plts.roc1.roc3,plts.roc2,nrow=1))



## g. synthetic data -- ROC 3

## h. synthetic data -- ROC 1

load('simulations/g/g.RData')
table3 <- table
load('simulations/h/h.RData')
table1 <- table

table <- abind(table1,table3,along=3)
dimnames(table) <- list(fpr=sprintf('%.*f',2,fprs),visit=as.character(ts),c('$ROC_1$','$ROC_3$'))

cover.table <- ftable(table,row.vars=c(2,3))
names(attr(cover.table,'row.vars')) <- c('\\textbf{visit}','\\textbf{ROC}')
names(attr(cover.table,'col.vars')) <- '\\textbf{FPR}'
require(xtable)
caption0 <- 'Nominal 95\\% CI coverage (CVL), bias (BS), mean standard error (MSE), and empirical standard error (ESE) of $ROC_1$ and $ROC_3$ for FPRs 10\\%, 25\\%, 50\\%, and 75\\% at visits 5, 15, 25, and 35 (synthetic data using hyperparameters estimated from the pediatric HIV data, $N=80$ patients).'

out <- print(xtableFtable(cover.table,label='tab:1',caption=caption0),include.rownames=F,include.colnames=F,quote=F,sanitize.rownames.function=function(x)x,sanitize.colnames.function=function(x)x)
out <- strsplit(out,split='\\\\\\\\')[[1]]
out[1] <- gsub('multicolumn\\{1\\}\\{l\\}','multicolumn\\{1\\}\\{c\\}',out[1])
out[1] <- sub('rrrr','llll',out[1])
out <- append(out,paste(c('& & ',rep('& CVL (BS,MSE,ESE)',4)),collapse=''),2)
out <- paste(out,collapse='\\\\')
## sink('figures/revision_g2.tex')
cat(out)
sink()




## i. HIV data -- roc 1 and roc 2 -- i.R

## j. HIV data -- roc 3 -- j.R

## k. HIV data -- roc 4 -- k.R
source('ragon.R')
load('simulations/i/i.RData')
load('simulations/k/k.RData')
gg.df <- rbind(attr(plt.roc1,'gg.df'),attr(plt.roc4,'gg.df'))
gg.df$group <- rep(c('roc1','roc4'),c(nrow(attr(plt.roc1,'gg.df')),nrow(attr(plt.roc4,'gg.df'))))
plts.roc1.roc3 <- ggplot(gg.df,aes(x=fpr,y=tpr,group=group))+geom_path(aes(linetype=group))+geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.5,fill='lightgray') + geom_abline(slope=1,col='lightgray',linetype=2) +
    theme_classic() + labs(x='FPR',y='TPR') + scale_linetype_discrete(name = 'ROC type',labels=c(expression(ROC[1]),expression(ROC[4])))
plts.roc2 <- do.call(grid.arrange,list(grobs=plts))
## png('revision_k.png',width=1920,height=1200)
grid.arrange(plts.roc1.roc3,plts.roc2,nrow=1)
dev.off()
## ggsave('revision.k.png',arrangeGrob(plts.roc1.roc3,plts.roc2,nrow=1))
## attr(plt.roc4,'gg.df')


## hiv.roc.1.plt <- plot(roc=plts.roc1[[1]]$observed, bootstrap.rocs=plts.roc1[[1]]$bootstrap, alpha=.05)
## gg.df.1 <- attr(hiv.roc.1.plt,'gg.df')
## load('simulations/k/k.RData')
## source('ragon.R')
## hiv.roc.4.plt <- plot(hiv.roc.4,bootstrap.rocs=hiv.roc.4.bs,alpha=.05)
## gg.df.4 <- attr(hiv.roc.4.plt,'gg.df')
## gg.df <- rbind(gg.df.1,gg.df.4)
## gg.df$group <- gl(2,nrow(gg.df.1),labels=c('roc.1','roc.4'))
## ggplot(gg.df, aes(x=fpr,y=tpr,group=group))+geom_line()+geom_ribbon(aes(ymax=upper,ymin=lower,alpha=.3))


## m. comparison of proposed method with glm/fitted hats method
## load('simulations/k/k.RData')
## load('simulations/m/m.RData')
## source('ragon.R')
## gg.df <- rbind(attr(plt.roc4,'gg.df'),attr(plt.glmer,'gg.df'))
## gg.df$group <- rep(c('roc4','glmer'),c(nrow(attr(plt.roc4,'gg.df')),nrow(attr(plt.glmer,'gg.df'))))
## plt <- ggplot(gg.df,aes(x=fpr,y=tpr,group=group))+geom_path(aes(linetype=group))+geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.5,fill='lightgray') + geom_abline(slope=1,col='lightgray',linetype=2) +
##     theme_classic() + labs(x='FPR',y='TPR') + scale_linetype_discrete(name = 'ROC type',labels=c(expression(ROC[4]),expression(GLM)))
## plt
## ## png('revision_m.png')
## plt
## dev.off()

load('simulations/i/i.RData')
load('simulations/k/k.RData')
load('simulations/m/m.RData')
source('ragon.R')
gg.df <- rbind(attr(plt.roc1,'gg.df'),attr(plt.roc4,'gg.df'),attr(plt.glmer,'gg.df'))
gg.df$group <- rep(c('roc1','roc4','glmer'),c(nrow(attr(plt.roc1,'gg.df')),nrow(attr(plt.roc4,'gg.df')),nrow(attr(plt.glmer,'gg.df'))))
plt <- ggplot(gg.df,aes(x=fpr,y=tpr,group=group))+geom_path(aes(linetype=group))+geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.5,fill='lightgray') + geom_abline(slope=1,col='lightgray',linetype=2) +
    theme_classic() + labs(x='FPR',y='TPR') + scale_linetype_discrete(name = 'ROC type',labels=c(expression(GLMM),expression(ROC[1]),expression(ROC[4])))
plt
png('figures/revision_m2.png')
plt
dev.off()

## n. new time series plots and histograms without using 50% dividing line


dat <- read.csv('final.csv')
dat$blip <- as.numeric(dat$VL>1e3)
dat <- subset(dat,select=c('ID','Blip_YN50','CD4','CD4P','Visits','blip','Sex','Transmission','trtm','AgeDiag','TimeDiag','AGE','DOV','TIME'))
dat <- subset(dat,select=c('ID','CD4','CD4P','Visits','blip','Blip_YN50','TIME'))
dat <- dplyr::rename(dat,Visit=Visits)
dat$CD4[dat$CD4>6000] <- NA
dat$CD4P[dat$CD4P>=99] <- NA
dat$Visit <- unlist(sapply(rle(dat$ID)$lengths,function(x)seq_len(x)))
dat$blip <- factor(dat$blip)
dat$ID <- factor(dat$ID)
dat <- droplevels(dat)

disease.color <- "#CC6666"
nondisease.color <- "#9999CC"
ggplot(na.omit(dat),aes(x=Visit,y=CD4,group=blip,shape=blip))+geom_point()+geom_smooth(color='gray')+theme_classic()+scale_shape_manual(values=c(1,3))
ggsave('revision_n1.png')


ggplot(na.omit(dat),aes(x=TIME,y=CD4,group=blip,shape=blip))+geom_point()+geom_smooth(color='gray')+theme_classic()+scale_shape_manual(values=c(1,3))+labs(x='years since diagnosis')
ggsave('revision_n2.png')

require(reshape2)
gg.df <- melt(dat,id.vars=c('ID','blip'),measure.vars=c('CD4','CD4P'),variable.name='test',value.name='marker')
gg.df <- gg.df[gg.df$test!='CD4P',]
gg.df$blip <- factor(gg.df$blip,labels=c('negative','positive'))
ggplot(na.omit(gg.df),aes(x=marker,fill=blip))+geom_histogram(alpha=.4,position='identity',aes(y=..density..))+theme_classic()+scale_fill_manual(values=c('grey20','grey50'),name='blip status')
ggsave('revision_n3.png')


## o. CIs for effect estimates
source('ragon.R')
data <- get.hiv.data()
N <- length(data)
hiv.params <- get.rho.sigma(data)
eta <- pop.fit(data,rho=hiv.params$rho)
est <- c(unlist(hiv.params),ar=eta[1,],mc=eta[2,])
est <- c(est,ar.dff=(est['ar.theta1'] - est['ar.theta0'])/est['sigma'])
est <- c(est,mc.dff=(est['mc.theta1'] - est['mc.theta0']))
est.bs <- replicate(3e3, {
    data.star <- data[sample(1:N,N,replace=TRUE)]
    rho.sigma.star <- unlist(get.rho.sigma(data.star))
    eta.star <- pop.fit(data.star,rho=hiv.params$rho)
    c(rho.sigma.star,ar=eta.star[1,],mc=eta.star[2,])
},simplify=FALSE)
## save.image('o.RData')
est.bs <- lapply(est.bs,function(l)unlist(l))
est.bs <- do.call(rbind,est.bs)
est.bs <- as.data.frame(est.bs)
est.bs$ar.diff <- with(est.bs, (ar.theta1 - ar.theta0) / sigma)
est.bs$mc.diff <- with(est.bs, (mc.theta1 - mc.theta0))
alpha <- .05
lower <- 2*est - apply(est.bs,2,quantile,1-alpha/2)
upper <- 2*est - apply(est.bs,2,quantile,alpha/2)


## p. overfitting in GLMM ROC -- see p/p.R
