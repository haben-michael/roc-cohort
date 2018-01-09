.libPaths(c('/afs/ir/users/h/a/habnice/Desktop/R/lib/R/library',"/usr/lib/R/library"))
source('../../ragon.R')

## 2. HIV data

## set.seed(1)
hiv.data <- get.hiv.data('../../final.csv')
hiv.params <- get.rho.sigma(hiv.data)

T <- attr(hiv.data,'T')
N <- attr(hiv.data,'N')
B <- 1e2
bootstrap.reps <- 3e0
ts <- c(5,15,25,35)
ts <- t <- 5

## 1. roc1

## plts.roc1 <- lapply(ts, function(t) {
hiv.params <- get.rho.sigma(hiv.data)
eta <- pop.fit(hiv.data,hiv.params$rho,t=t-1)   
## hiv.roc.1 <- roc.1(eta,hiv.params$rho,hiv.params$sigma,t,B=B,roc.hat.reps=3e2)
## plot(hiv.roc.1)

hiv.roc.1.bs <- replicate(bootstrap.reps,{
    ## ## weights <- rexp(N)
    ## ## eta.star <- pop.fit(hiv.data,hiv.params$rho,weights=weights)
    ## ## eta.star['mc',1:2] <- plogis(eta.star['mc',1:2])
    ## hiv.data.star <- hiv.data[sample(1:N,replace=TRUE)]
    ## hiv.params.star <- get.rho.sigma(hiv.data.star)
    ## eta.star <- pop.fit(hiv.data.star,hiv.params.star$rho,t=t)
    ## tryCatch(roc.1(eta.star,hiv.params.star$rho,hiv.params.star$sigma,t,B=B,roc.hat.reps=1e2),error=function(e)NA)
    tryCatch(roc.1(eta,hiv.params$rho,hiv.params$sigma,t,B=N,roc.hat.reps=3e2),error=function(e)NA)
   
})

## hiv.roc.1.bs <- hiv.roc.1.bs[!is.na(hiv.roc.1.bs)]

## list(observed=hiv.roc.1,bootstrap=hiv.roc.1.bs)
## plt.roc1 <- plot(hiv.roc.1,bootstrap.rocs=hiv.roc.1.bs,alpha=.05)+ggtitle(substitute(paste(ROC[1],'(u|',t,') at time t=',t),list(t=t)))
## plt1.roc1
## })

save.image(paste0('i-',round(runif(1)*1e8),'.RData'))
## plot(plts.roc1[[1]]$observed,bootstrap.rocs=plts.roc1[[1]]$bootstrap,alpha=.05)
## save.image('i.RData')

## ## save(plt1,file='020917_plt1.RData') ## to use in 011417_corn.R
## ## ggsave('011117a.png',plt1)
## t <- 4
## r1 <- roc.1(pop.fit(hiv.data,rho=hiv.params$rho,t=t),rho=hiv.params$rho,sigma=hiv.params$sigma,t=t,B=B,roc.hat.reps=1e2)
## plot(r1)
## hiv.data.star <- hiv.data[sample(1:N,replace=TRUE)]
## r2 <- roc.1(pop.fit(hiv.data.star,rho=hiv.params$rho,t=5),rho=hiv.params$rho,sigma=hiv.params$sigma,t=5,B=B,roc.hat.reps=1e2)
## ll <- list(r1,r2)
## class(ll) <- c(class(ll),'PatientROC')
## plot(ll)
## bootstrap.rocs <- replicate(5,{
##     cat('.')
##     hiv.data.star <- hiv.data[sample(1:N,replace=TRUE)]
##     roc.1(pop.fit(hiv.data.star,rho=hiv.params$rho,t=5+1),rho=hiv.params$rho,sigma=hiv.params$sigma,t=5,B=B,roc.hat.reps=1e2)
## },simplify=FALSE)
## ll <- c(r1,bootstrap.rocs[20:30])
## class(ll) <- c(class(ll),'PatientROC')
## plot(ll)

## dd
## ## r2 <- roc.1(pop.fit(hiv.data,rho=hiv.params$rho),hiv.params$rho,hiv.params$sigma,5,B=B,roc.hat.reps=1e2)
## ## ll <- list(r1,r2)
## ## class(ll) <- c(class(ll),'PatientROC')
## ## plot(ll)


## bootstrap.rocs <- replicate(5,roc.1(pop.fit(hiv.data[sample(1:N,replace=FALSE)],rho=hiv.params$rho),hiv.params$rho,hiv.params$sigma,5,B=B,roc.hat.reps=1e2),simplify=FALSE)
## bootstrap.rocs <- replicate(5,roc.1(pop.fit(hiv.data[sample(1:N,replace=FALSE)],rho=hiv.params$rho),hiv.params$rho,hiv.params$sigma,5,B=B,roc.hat.reps=1e2),simplify=FALSE)
## plot(r1,bootstrap.rocs=bootstrap.reps,alpha=.05)
## ll <- c(r1,bootstrap.rocs)
## class(ll) <- c(class(ll),'PatientROC')
## plot(ll)


## ll <- plts.roc1[[1]]$bootstrap[1:5]
## ll <- c(plts.roc1[[1]]$observed,ll)
## class(ll) <- c(class(ll),'PatientROC')
## plot(ll)

## ## 2. roc2
## t <- 5
## xi <- subj.fit(hiv.data,hiv.params$rho,hiv.params$sigma,eta,t)
## patient.idx <- c(6,14,17,61)
## hiv.roc.2 <- lapply(patient.idx,function(i)roc.2(xi=xi[i,],sigma=hiv.params$sigma,rho=hiv.params$rho))
## class(hiv.roc.2) <- c(class(hiv.roc.2),'PatientROC')
## plot(hiv.roc.2)

## ## repeat {
## hiv.roc.2.bs <-
##     replicate(bootstrap.reps, {
##         cat('.')
##         tryCatch({
##             ## weights <- rexp(N)
##             ## eta.star <- pop.fit(hiv.data,rho=hiv.params$rho,weights=weights)
##             ## eta.star['mc',1:2] <- plogis(eta.star['mc',1:2])
##             eta.star <- pop.fit(hiv.data[sample(1:N,replace=TRUE)],rho=hiv.params$rho)
##             xi <- subj.fit(hiv.data,hiv.params$rho,hiv.params$sigma,eta.star,t)
##             lapply(patient.idx,function(i)roc.2(xi=xi[i,],hiv.params$sigma,hiv.params$rho))
##             },error=function(e)NA)
##     },simplify=F)
## ##     if(!('error' %in% class(bootstrap.rocs))) break
## ## }

## hiv.roc.2.bs <- hiv.roc.2.bs[!is.na(hiv.roc.2.bs)]

## hiv.roc.2 <- lapply(1:length(patient.idx),function(i)list(roc=hiv.roc.2[[i]],bootstrap.rocs=lapply(hiv.roc.2.bs,function(x)x[[i]])))

## plts <- lapply(1:length(patient.idx),function(i)
##     plot(hiv.roc.2[[i]]$roc,bootstrap.rocs=hiv.roc.2[[i]]$bootstrap.rocs,alpha=.05)
##     )
## plts.roc2 <- lapply(1:length(patient.idx),function(i)plts[[i]]+ggtitle(paste0('patient #',patient.idx[i])))
## ## plts.roc2 <- do.call(grid.arrange,list(grobs=plts))
## ## plts.roc2 <- do.call(grid.arrange,list(grobs=plts,top=grid.text(expression(paste('Limiting ROC curve ',ROC[i*infinity],' for selected patients')))))
## ## ggsave('011117c.png',plt3)
## ## save.image('020917_hiv.RData')








## ##combined plots--roc3 here is bad--use combined roc3. roc1 and 4 roc2s are good.
## ## plt3 <- do.call(grid.arrange,list(grobs=plts,ncol=1,top=grid.text(expression(paste(ROC[i*infinity],', selected patients')))))+ylim(0,1.1)
## ## plt1 <- plt1+ggtitle(substitute(paste(ROC[1],'(u|t) at time t=',t),list(t=t)))
## ## plt12 <- grid.arrange(plt1,plt2,ncol=1)
## ## ggsave('020917.png',grid.arrange(plt12,plt3,ncol=2,top='ROC estimates for pediatric HIV data'))



##############################################################################
## from corn
source('../../ragon.R')
filelist <- dir()
filelist <- filelist[grep('i-[0-9]{8}.RData$',filelist)]
bootstrap.rocs <- list()
## k <- 1
for(file in filelist) {
    ## print(k)
    ## k <- k+1
    load(file)
    bootstrap.rocs <- c(bootstrap.rocs,hiv.roc.1.bs)
}
bootstrap.rocs <- bootstrap.rocs[!is.na(bootstrap.rocs)]
eta <- pop.fit(hiv.data,rho=hiv.params$rho,t=t-1)
hiv.roc.1 <- roc.1(eta,rho=hiv.params$rho,sigma=hiv.params$sigma,t,B=N,roc.hat.reps=3e2)
plt.roc1 <- plot(hiv.roc.1,bootstrap.rocs=bootstrap.rocs,alpha=.05,resolution=2e2)
## save.image('i.RData')
load('i.RData')
plt.roc1


## ll <- c(hiv.roc.1,bootstrap.rocs[1:10])
## class(ll) <- c(class(ll),'PatientROC')
## plot(ll)
