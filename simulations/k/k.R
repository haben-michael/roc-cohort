.libPaths(c('/afs/ir/users/h/a/habnice/Desktop/R/lib/R/library',"/usr/lib/R/library"))
source('../../ragon.R')

## 2. HIV data -- ROC3

dat <- read.csv('../../final.csv')
dat <- subset(dat,select=c('ID','Blip_YN50','CD4','CD4P','Visits','blip','Sex','Transmission','trtm','AgeDiag','TimeDiag','AGE'))
dat <- subset(dat,select=c('ID','CD4','CD4P','Visits','blip','Blip_YN50'))
dat <- rename(dat,Visit=Visits)
dat$CD4[dat$CD4>6000] <- NA
dat$CD4P[dat$CD4P>=99] <- NA
dat$Visit <- unlist(sapply(rle(dat$ID)$lengths,function(x)seq_len(x)))
dat$blip <- factor(dat$blip)
dat$ID <- factor(dat$ID)
dat <- droplevels(dat)

## ggplot(na.omit(subset(dat, (ID %in% unique(dat$ID[which(dat$Visit>50)])))),aes(x=Visit,y=CD4P,group=ID))+geom_point(aes(color=Blip_YN50))+geom_path()+facet_wrap(~ID)

## ggplot(na.omit(dat),aes(x=CD4P,fill=blip))+geom_histogram(position='identity',alpha=.5)

## group_by(dat,ID) %>% summarize(mean(blip,na.rm=T)) %>% as.data.frame %>% plot

dat <- split(dat,dat$ID)
dat <- lapply(dat,function(dat.i)list(X=dat.i$CD4,S=2-as.numeric(dat.i$blip)))


## dat.i <- dat[[1]]
ans <- sapply(dat,function(dat.i) {
    runs <- rle(dat.i$S)
    ars <- split(dat.i$X,rep(1:length(runs$lengths),runs$lengths))
    ars.pos <- ars[runs$values==1]; ars.neg <- ars[runs$values==0]
    ols.pos <- sapply(ars.pos,function(ar){
        if(length(ar)>=20) {
            ar0 <- arima(ar,order=c(1,0,0),method='ML')#,demean=F,intercept=T)
            ## if(ar0$order==0) ar0$ar <- 0
            ## c(rho=as.numeric(ar0$coef['ar1']),sigma2=ar0$sigma2,rho.se=as.numeric(ar0$var.coef['ar1','ar1'])
            coefs <- c(ar0$coef,diag(ar0$var.coef),ar0$sigma2) %>% as.numeric()
            names(coefs) <- c('rho','intercept','rho.se','intercept.se','sigma2')
            coefs
        } else {c('rho'=NA,'intercept'=NA,'rho.se'=NA,'intercept.se'=NA,'sigma2'=NA)}
    })
    ols.neg <- sapply(ars.neg,function(ar){
        if(length(ar)>=20) {
            ar0 <- arima(ar,order=c(1,0,0),method='ML')#,demean=F,intercept=T)
            ## if(ar0$order==0) ar0$ar <- 0
            ## c(rho=as.numeric(ar0$coef['ar1']),sigma2=ar0$sigma2,rho.se=as.numeric(ar0$var.coef['ar1','ar1'])
            coefs <- c(ar0$coef,diag(ar0$var.coef),ar0$sigma2) %>% as.numeric()
            names(coefs) <- c('rho','intercept','rho.se','intercept.se','sigma2')
            coefs
        } else {c('rho'=NA,'intercept'=NA,'rho.se'=NA,'intercept.se'=NA,'sigma2'=NA)}
    })
    ## ols.neg <- sapply(ars.neg,function(ar){
    ##     if(length(ar)>=2) {
    ##         ar0 <- ar.ols(ar,order=1,demean=F,intercept=T)
    ##         if(ar0$order==0) ar0$ar <- 0
    ##         c(rho=ar0$ar,sigma2=var(ar0$resid,na.rm=T))
    ##     } else {c(rho=NA,sigma2=NA)}
    ## })
    if(length(ols.pos)==0) ols.pos <- as.matrix(c('rho'=NA,'intercept'=NA,'rho.se'=NA,'intercept.se'=NA,'sigma2'=NA),nrow=5)
    if(length(ols.neg)==0) ols.neg <- as.matrix(c('rho'=NA,'intercept'=NA,'rho.se'=NA,'intercept.se'=NA,'sigma2'=NA),nrow=5)
    ols.pos['rho.se',] <- sqrt(ols.pos['rho.se',])
    ols.neg['rho.se',] <- sqrt(ols.neg['rho.se',])
    pos.est <- c(rho=sum(ols.pos['rho',]/ols.pos['rho.se',],na.rm=T)/sum(1/ols.pos['rho.se',],na.rm=T),
                 intercept=sum(ols.pos['intercept',]/ols.pos['intercept.se',],na.rm=T)/sum(1/ols.pos['intercept.se',],na.rm=T),
                 sigma2=mean(ols.pos['sigma2',],na.rm=T))
    neg.est <- c(rho=sum(ols.neg['rho',]/ols.neg['rho.se',],na.rm=T)/sum(1/ols.neg['rho.se',],na.rm=T),
                 intercept=sum(ols.neg['intercept',]/ols.neg['intercept.se',],na.rm=T)/sum(1/ols.neg['intercept.se',],na.rm=T),
                 sigma2=mean(ols.neg['sigma2',],na.rm=T))
    apply(rbind(pos.est,neg.est),2,mean,na.rm=T)
})
hiv.params <- list(rho=mean(ans['rho',],na.rm=T),
  intercept=mean(ans['intercept',],na.rm=T),
  sigma=sqrt(mean(ans['sigma2',],na.rm=T)))


dat <- lapply(dat,function(dat.i) {
    bad.idx <- unique(c(which(is.na(dat.i$X)),which(is.na(dat.i$S))))
    list(X=dat.i$X[-bad.idx],S=dat.i$S[-bad.idx])
})
## plot(ecdf(sapply(dat,function(d)length(d$X))))
T <- 20
## sum(sapply(dat,function(d)length(d$X))>T)
## for (T in 5:50) {
## dat <- dat0
dat <- dat[sapply(dat,function(d)length(d$X))>T]
dat <- lapply(dat,function(dat.i)list(X=dat.i$X[1:T],S=dat.i$S[1:T]))
## dat <- dat[sapply(dat,function(d)(1 %in% d$S) & (0 %in% d$S))]
N <- length(dat)
## cat(T,' ',N,'\n')
## }
attr(dat,'N') <- N
attr(dat,'T') <- T
attr(dat,'class') <- 'PatientClass'
attr(dat,'names') <- NULL
hiv.data <- dat

## set.seed(1)
B <- 1e2
bootstrap.reps <- 3e2




## roc4
## source('../../ragon.R')
## hiv.data <- get.hiv.data('../../final.csv')
## hiv.params <- get.rho.sigma(hiv.data)

## T <- attr(hiv.data,'T')
## N <- attr(hiv.data,'N')
## B <- 1e2
## bootstrap.reps <- 3e0

ts <- 8:10
ll <- lapply(5:10,function(ts){
rocs <- lapply(ts,function(t) {
    eta <- pop.fit(hiv.data,rho=hiv.params$rho,t=t-1)
    xi <- subj.fit(hiv.data,hiv.params$rho,hiv.params$sigma,eta=eta,t=t-1)
    roc(hiv.data,xi,rho=hiv.params$rho,sigma=hiv.params$sigma,t)
})
})
ll <- unlist(ll)
class(ll) <- c(class(ll),'PatientROC')
plot(ll)

hiv.roc.4 <- average.roc(rocs)
plot(hiv.roc.4)

hiv.roc.4.bs <- replicate(bootstrap.reps, {
    cat('.')
    hiv.data.star <- hiv.data[sample(1:N,replace=TRUE)]

    tryCatch({
        rocs <- lapply(ts,function(t) {
            eta <- pop.fit(hiv.data.star,rho=hiv.params$rho,t=t-1)
            xi <- subj.fit(hiv.data.star,hiv.params$rho,hiv.params$sigma,eta,t=t-1)
            roc(hiv.data.star,xi,rho=hiv.params$rho,sigma=hiv.params$sigma,t)
        })
        average.roc(rocs)
    },error=function(e)NA)
    
}, simplify=FALSE)

hiv.roc.4.bs <- hiv.roc.4.bs[!is.na(hiv.roc.4.bs)]
plt.roc4 <- plot(hiv.roc.4,bootstrap.rocs=hiv.roc.4.bs,alpha=.05,N=N,gaussian.CI=TRUE)
plt.roc4
## save.image('k.RData')
load('k.RData')
## save.image(paste0('j-',round(runif(1)*1e8),'.RData'))

## ## from corn
## source('../../ragon.R')
## filelist <- dir()
## filelist <- filelist[grep('j-[0-9]{8}.RData$',filelist)]
## roc.3.bs.full <- list()
## ## k <- 1
## for(file in filelist) {
##     ## print(k)
##     ## k <- k+1
##     load(file)
##     roc.3.bs.full <- c(roc.3.bs.full,hiv.roc.3.bs)
## }

## hiv.roc.3.bs <- roc.3.bs.full
## hiv.roc.3.bs <- hiv.roc.3.bs[!is.na(hiv.roc.3.bs)]
## plot(hiv.roc.3,bootstrap.rocs=hiv.roc.3.bs,alpha=.05,resolution=1e2)+ylim(0,1.1)
