## setwd('/dropbox/ragon')
.libPaths(c('/afs/ir/users/h/a/habnice/Desktop/R/lib/R/library',"/usr/lib/R/library"))
source('../g/ragon.old.R')

## ## simulation data based on HIV data hyperparameters--large number of patients

## dat <- read.csv('../../final.csv')
## dat$blip <- as.numeric(dat$VL>1e3)
## dat <- subset(dat,select=c('ID','Blip_YN50','CD4','CD4P','Visits','blip','Sex','Transmission','trtm','AgeDiag','TimeDiag','AGE'))
## dat <- subset(dat,select=c('ID','CD4','CD4P','Visits','blip','Blip_YN50'))
## dat <- rename(dat,Visit=Visits)
## dat$CD4[dat$CD4>6000] <- NA
## dat$CD4P[dat$CD4P>=99] <- NA
## dat$Visit <- unlist(sapply(rle(dat$ID)$lengths,function(x)seq_len(x)))
## dat$blip <- factor(dat$blip)
## dat$ID <- factor(dat$ID)
## dat <- droplevels(dat)

## dat <- split(dat,dat$ID)
## dat <- lapply(dat,function(dat.i)list(X=dat.i$CD4P,S=2-as.numeric(dat.i$blip)))


## ans <- sapply(dat,function(dat.i) {
##     runs <- rle(dat.i$S)
##     ars <- split(dat.i$X,rep(1:length(runs$lengths),runs$lengths))
##     ars.pos <- ars[runs$values==1]; ars.neg <- ars[runs$values==0]
##     ols.pos <- sapply(ars.pos,function(ar){
##         if(length(ar)>=20) {
##             ar0 <- arima(ar,order=c(1,0,0),method='ML')#,demean=F,intercept=T)
##             ## if(ar0$order==0) ar0$ar <- 0
##             ## c(rho=as.numeric(ar0$coef['ar1']),sigma2=ar0$sigma2,rho.se=as.numeric(ar0$var.coef['ar1','ar1'])
##             coefs <- c(ar0$coef,diag(ar0$var.coef),ar0$sigma2) %>% as.numeric()
##             names(coefs) <- c('rho','intercept','rho.se','intercept.se','sigma2')
##             coefs
##         } else {c('rho'=NA,'intercept'=NA,'rho.se'=NA,'intercept.se'=NA,'sigma2'=NA)}
##     })
##     ols.neg <- sapply(ars.neg,function(ar){
##         if(length(ar)>=20) {
##             ar0 <- arima(ar,order=c(1,0,0),method='ML')#,demean=F,intercept=T)
##             ## if(ar0$order==0) ar0$ar <- 0
##             ## c(rho=as.numeric(ar0$coef['ar1']),sigma2=ar0$sigma2,rho.se=as.numeric(ar0$var.coef['ar1','ar1'])
##             coefs <- c(ar0$coef,diag(ar0$var.coef),ar0$sigma2) %>% as.numeric()
##             names(coefs) <- c('rho','intercept','rho.se','intercept.se','sigma2')
##             coefs
##         } else {c('rho'=NA,'intercept'=NA,'rho.se'=NA,'intercept.se'=NA,'sigma2'=NA)}
##     })
##     ## ols.neg <- sapply(ars.neg,function(ar){
##     ##     if(length(ar)>=2) {
##     ##         ar0 <- ar.ols(ar,order=1,demean=F,intercept=T)
##     ##         if(ar0$order==0) ar0$ar <- 0
##     ##         c(rho=ar0$ar,sigma2=var(ar0$resid,na.rm=T))
##     ##     } else {c(rho=NA,sigma2=NA)}
##     ## })
##     if(length(ols.pos)==0) ols.pos <- as.matrix(c('rho'=NA,'intercept'=NA,'rho.se'=NA,'intercept.se'=NA,'sigma2'=NA),nrow=5)
##     if(length(ols.neg)==0) ols.neg <- as.matrix(c('rho'=NA,'intercept'=NA,'rho.se'=NA,'intercept.se'=NA,'sigma2'=NA),nrow=5)
##     ols.pos['rho.se',] <- sqrt(ols.pos['rho.se',])
##     ols.neg['rho.se',] <- sqrt(ols.neg['rho.se',])
##     pos.est <- c(rho=sum(ols.pos['rho',]/ols.pos['rho.se',],na.rm=T)/sum(1/ols.pos['rho.se',],na.rm=T),
##                  intercept=sum(ols.pos['intercept',]/ols.pos['intercept.se',],na.rm=T)/sum(1/ols.pos['intercept.se',],na.rm=T),
##                  sigma2=mean(ols.pos['sigma2',],na.rm=T))
##     neg.est <- c(rho=sum(ols.neg['rho',]/ols.neg['rho.se',],na.rm=T)/sum(1/ols.neg['rho.se',],na.rm=T),
##                  intercept=sum(ols.neg['intercept',]/ols.neg['intercept.se',],na.rm=T)/sum(1/ols.neg['intercept.se',],na.rm=T),
##                  sigma2=mean(ols.neg['sigma2',],na.rm=T))
##     apply(rbind(pos.est,neg.est),2,mean,na.rm=T)
## })
## hiv.params <- list(rho=mean(ans['rho',],na.rm=T),
##   intercept=mean(ans['intercept',],na.rm=T),
##   sigma=sqrt(mean(ans['sigma2',],na.rm=T)))


## eta <- pop.fit(dat,hiv.params$rho)
## eta['mc',1:2] <- plogis(eta['mc',1:2])

## ## set.seed(2)
## N <- 5e2
## T <- 1e1
## t <- 7
## bootstrap.reps <- 5
## B <- 5e1
## mu.theta.prior <- eta['ar',1:2]; mu.p.prior <- eta['mc',1:2]
## sigma.theta.prior <- matrix(eta['ar',c(3,4,4,5)],nrow=2);
## sigma.p.prior <- matrix(eta['mc',c(3,4,4,5)],nrow=2)
## theta <- rmvnorm(N,mean=mu.theta.prior,sigma=sigma.theta.prior)
## p <- rmvnorm(N,mean=mu.p.prior,sigma=sigma.p.prior) %>% plogis
## data <- rpatient(hiv.params$rho,hiv.params$sigma,theta,p,T)



## eta <- pop.fit(data,hiv.params$rho)
## xi <- subj.fit(data,hiv.params$rho,hiv.params$sigma,eta,t-1)


## roc3
## xi <- subj.fit(data,hiv.params$rho,hiv.params$sigma,eta,t-1)


## save.image('f.init.RDatatmp')
load('f.init.RDatatmp')
set.seed(NULL)
bootstrap.reps <- 1

counter <- 0
sim.roc.3.bs <- replicate(bootstrap.reps,{
    cat('.')
    if(counter%%1e2==0)print(counter)
    counter <<- counter+1
    tryCatch({
        ## weights <- rexp(N)        
        ## eta.star <- pop.fit(data,rho=hiv.params$rho,weights=weights)
        ## rocs <- list()
        ## xi.star <- subj.fit(data,hiv.params$rho,hiv.params$sigma,eta.star,t-1)
        ## roc(data,xi.star,rho=hiv.params$rho,sigma=hiv.params$sigma,t)
        roc.3(data[sample(1:N,replace=T)],rho=hiv.params$rho,sigma=hiv.params$sigma,t=t)
    }, error=function(e){cat('*');NA})
},simplify=F)

sim.roc.3.bs <- sim.roc.3.bs[!is.na(sim.roc.3.bs)]
## save.image('revision_f3.RData')

## plt1 <- plot(sim.roc.3,bootstrap.rocs=sim.roc.3.bs,alpha=.05)+ggtitle(substitute(paste(ROC[3],'(u|t) at time t=',t),list(t=t)))+theme_classic()+labs(x='FPR',y='TPR')



## 2. roc1
## plot(hiv.roc.1)

sim.roc.1.bs <- replicate(bootstrap.reps,{
    cat('.')
    ## weights <- rexp(N)
    ## eta.star <- pop.fit(data,hiv.params$rho,weights=weights)
    eta.star <- pop.fit(data[sample(1:N,replace=T)],hiv.params$rho)
    tryCatch(roc.1(eta.star,hiv.params$rho,hiv.params$sigma,t,B=B,roc.hat.reps=1e2),error=function(e){cat('*');NA})           
},simplify=F)

sim.roc.1.bs <- sim.roc.1.bs[!is.na(sim.roc.1.bs)]

save.image(paste0('f-',round(runif(1)*1e8),'.RData'))


## roc2
patient.idx <- 1:15
roc.2.observed <- lapply(patient.idx,function(i)roc.2(xi=xi[i,],sigma=hiv.params$sigma,rho=hiv.params$rho))
do.call(grid.arrange,lapply(roc.2.observed,plot))
bootstrap.reps <- 1.5e2
sim.roc.2.bs <-
    replicate(bootstrap.reps, {
        cat('.')
        tryCatch({
            ## weights <- rexp(N)
            eta.star <- pop.fit(data[sample(1:N,replace=T)],rho=hiv.params$rho)
            ## eta.star <- pop.fit(data,rho=hiv.params$rho,weights=weights)
            ## eta.star['mc',1:2] <- plogis(eta.star['mc',1:2])
            xi <- subj.fit(data,hiv.params$rho,hiv.params$sigma,eta.star,t-1)
            lapply(patient.idx,function(i)roc.2(xi=xi[i,],hiv.params$sigma,hiv.params$rho))
            },error=function(e)NA)
    },simplify=F)


sim.roc.2.bs <- sim.roc.2.bs[!is.na(sim.roc.2.bs)]

sim.roc.2 <- lapply(1:length(patient.idx),function(i)list(roc=roc.2.observed[[i]],bootstrap.rocs=lapply(sim.roc.2.bs,function(x)x[[i]])))
## save.image('revision_f2.RData')
plts.roc2 <- lapply(c(1,13,11,7),function(i) {
    cat('.')
    plot(sim.roc.2[[i]]$roc,bootstrap.rocs=sim.roc.2[[i]]$bootstrap.rocs,alpha=.05,gaussian.CI=TRUE,N=1)
    })


## plts <- lapply(1:length(patient.idx),function(i)plts[[i]]+ggtitle(paste0('patient #',patient.idx[i])))
## plt3 <- do.call(grid.arrange,list(grobs=plts,top=grid.text(expression(paste('Limiting ROC curve ',ROC[i*infinity],' for selected patients')))))



source('../../ragon.R')
filelist <- dir()
filelist <- filelist[grep('f-[0-9]{8}.RData$',filelist)]
k <- 1

sim.roc.1.bs.full <- list()
sim.roc.3.bs.full <- list()
for(file in filelist) {
    print(k)
    k <- k+1
    load(file)

    sim.roc.1.bs.full <- c(sim.roc.1.bs.full,list(sim.roc.1.bs))

    sim.roc.3.bs.full <- c(sim.roc.3.bs.full,list(sim.roc.3.bs))
}
bootstrap.rocs.1 <- unlist(sim.roc.1.bs.full)
bootstrap.rocs.3 <- unlist(sim.roc.3.bs.full)
sim.roc.3 <- roc(data,xi,rho=hiv.params$rho,sigma=hiv.params$sigma,t)
sim.roc.1 <- roc.1(eta,hiv.params$rho,hiv.params$sigma,t,B=1e2,roc.hat.reps=4e2)
plot(sim.roc.1)

plt.roc1 <- plot(sim.roc.1,bootstrap.rocs=bootstrap.rocs.1,alpha=.05)
plt.roc1
plt.roc3 <- plot(sim.roc.3,bootstrap.rocs=bootstrap.rocs.3,alpha=.05)
plt.roc3

save.image('../f.RData')
