.libPaths(c("/usr/lib/R/library",'/home/habnice/R/lib/R/library'))
source('../../ragon.old.R')

## roc1 -- synthetic data 
B <- 1e2
N <- 8e1
T <- 4e1

ts <- c(5,15,25,35)
fprs <- c(.1,.25,.5,.75)
alpha <- .05
q <- qnorm(1-alpha/2)

## hiv.data <- get.hiv.data('../../final.csv')
## hiv.params <- get.rho.sigma(hiv.data)
## eta <- pop.fit(hiv.data,rho=hiv.params$rho)
## mu.theta.prior <- eta['ar',1:2]; mu.p.prior <- eta['mc',1:2]
## sigma.theta.prior <- matrix(eta['ar',c(3,4,4,5)],nrow=2);
## sigma.p.prior <- matrix(eta['mc',c(3,4,4,5)],nrow=2)
## theta <- rmvnorm(B,mean=mu.theta.prior,sigma=sigma.theta.prior)
## p <- rmvnorm(B,mean=mu.p.prior,sigma=sigma.p.prior) %>% plogis

rho <- .5
sigma <- 1
mu.theta.prior <- c(3,4)
sigma.theta.prior <- matrix(c(1,.1,.1,2),nrow=2)
mu.p.prior <- c(3,0)
sigma.p.prior <- matrix(c(2.7,-2.3,-2.3,3),nrow=2)
eta <- rbind(ar=c(mu.theta.prior,sigma.theta.prior[c(1,2,4)]),mc=c(mu.p.prior,sigma.p.prior[c(1,2,4)]))


## roc.1.true <- roc.1(eta,rho,sigma,ts,B,roc.hat.reps=5e2)
## plot(roc.1.true[[2]])
## tpr.true <- sapply(roc.1.true,function(r)r(fprs))
## require(MASS)
## write.matrix(tpr.true,'true.txt')


## xi <- c(theta,p)


##     roc.1.star <- roc.1(eta.star,rho,sigma,ts=ts,B,roc.hat.reps=5e1)

## ll <- structure(c(roc.1.true,roc.1.star),class=c('list','PatientROC'))

bootstrap.reps <-2e1
## reps <- 1
## roc.1.sim <- list()
## roc1s <- replicate(reps, {
##     cat('-')

theta <- rmvnorm(N,mean=mu.theta.prior,sigma=sigma.theta.prior)
p <- rmvnorm(N,mean=mu.p.prior,sigma=sigma.p.prior) %>% plogis
data <- rpatient(rho=rho,sigma=sigma,theta=theta,p=p,T=T)
eta.star <- pop.fit(data,rho=rho)

roc.1.star <- roc.1(eta.star,rho,sigma,ts=ts,B,roc.hat.reps=5e1)
## ll <- structure(list(roc.1.true[[1]],roc.1.star[[1]]),class=c('list','PatientROC'))
## plot(ll)

bootstrap.rocs <- replicate(bootstrap.reps,{
    ## system.time({
    cat('.')      
    eta.star <- pop.fit(data[sample(1:N,replace=TRUE)],rho)

    tryCatch(
        roc.1(eta.star,rho,sigma,ts=ts,B,roc.hat.reps=5e1),
        error=function(e){cat('*');NA}
    )
    
    ## tryCatch(
    ##     roc.1(eta.star,rho,sigma,ts=ts,B,roc.hat.reps=1e2),
    ##     error=function(e){cat('*');NA}
    ##     )
    ## })
},simplify=F)

tpr.true <- as.matrix(read.table('true.txt',as.is=TRUE))

tpr.observed <- sapply(roc.1.star,function(r)r(fprs))

if(sum(is.na(bootstrap.rocs))>0)cat('!')
bootstrap.rocs <- bootstrap.rocs[!is.na(bootstrap.rocs)]
tpr.bs <- lapply(bootstrap.rocs, function(r)
    sapply(r,function(rr)rr(fprs)))
tpr.bs <- abind(tpr.bs,along=3)
boot.sds <- apply(tpr.bs,1:2,sd)
rownames(boot.sds) <- fprs
colnames(boot.sds) <- ts

upper <- tpr.observed + q*boot.sds
lower <- tpr.observed - q*boot.sds
coverage.1 <- tpr.true <= upper & tpr.true >= lower

lower <- apply(tpr.bs,1:2,quantile,probs=alpha/2)
upper <- apply(tpr.bs,1:2,quantile,probs=1-alpha/2)
coverage.3 <- tpr.true <= upper & tpr.true >= lower 

bias <- tpr.true - tpr.observed
roc.1.sim <- list(coverage.1=coverage.1,coverage.3=coverage.3,bias=bias,boot.sds=boot.sds,mu.star=tpr.observed)   

save(roc.1.sim,file=paste0('h-',round(runif(1)*1e8),'.RData'))

## ## boot.sds <- sapply(rocs,function(r)r(fprs)) %>% apply(1,sd)
## ## CI.star <- rbind(upper=roc.1.star(fprs)+q*boot.sds,
## ##                  lower=roc.1.star(fprs)-q*boot.sds)
## ## colnames(CI.star) <- fprs
## ## cover <- roc.1.true(fprs)<CI.star['upper',] & roc.1.true(fprs)>CI.star['lower',]
## ## print(CI.star)
## ## print(cover)

## ## ## roc.1.sim[[i]] <- list(roc=roc.1.star,bootstrap.rocs=rocs)
## ## roc.1.sim[[i]] <- list(CI=CI.star,coverage=cover)

## boot.sds <- abind(lapply(rocs,function(r)
##     sapply(r,function(rr)rr(fprs))), along=3) %>% apply(1:2,sd)
## rownames(boot.sds) <- fprs; colnames(boot.sds) <- ts

## mu.star <- sapply(roc.1.star,function(r)r(fprs))
## upper <- mu.star+q*boot.sds
## lower <- mu.star-q*boot.sds
## coverage <- true<=upper & true>=lower

## print(coverage)

## ## roc.1.sim[[i]] <- list(roc=roc.1.star,bootstrap.rocs=rocs)
## roc.1.sim[[i]] <- coverage


## if(i%%10==0)save.image(paste0('011717_roc1.RData'))
## roc.1.sim <- list(roc.star=roc.1.star,bootstrap.rocs=rocs)


## save.image(paste0('h-',round(runif(1)*1e8),'.RData'))







## ## from barley
## source('../../ragon.R')

## filelist <- dir()
## filelist <- filelist[grep('h-[0-9]{8}.RData$',filelist)]
## roc.1.sim.full <- lapply(filelist[300:366], function(file) {
##     load(file)
##     ## bootstrap.rocs <- c(bootstrap.rocs,hiv.roc.1.bs)

##     list(coverage.1=roc.1.sim$coverage.1,coverage.3=roc.1.sim$coverage.3,bias=roc.1.sim$bias,boot.sds=roc.1.sim$boot.sds,mu.star=roc.1.sim$mu.star)   
## })
## roc.1.sim <- roc.1.sim.full

## coverages <- lapply(roc.1.sim,function(l)l$coverage.1)
## cc <- abind(coverages,along=3)
## cc1 <- apply(cc,1:2,mean)


## biases <- lapply(roc.1.sim,function(l)l$bias)
## bb <- abind(biases,along=3)
## bb1 <- apply(bb,1:2,mean)

## ss <- lapply(roc.1.sim,function(r)r$boot.sds)
## ss <- abind(ss,along=3)
## ss1 <- apply(ss,1:2,mean)

## mm <- lapply(roc.1.sim,function(r)r$mu.star)
## mm <- abind(mm,along=3)
## mm1 <- apply(mm,1:2,sd)

## table <- abind(cc1,bb1,ss1,mm1,along=3) %>%
##     apply(1:3,function(x)sprintf('%.*f',2,x)) %>%
##     apply(1:2,function(x)paste0(x[1],' (',x[2],',',x[3],',',x[4],')'))
## rownames(table) <- fprs;colnames(table) <- ts
## dimnames(table) <- list(fpr=sprintf('%.*f',2,fprs),visit=as.character(ts))
## table

## ## table <- abind(cc1,bb1,along=3) %>%
## ##     apply(1:2,function(x)paste0(x[1],' (',x[2],')'))
## ## rownames(table) <- fprs;colnames(table) <- ts
## ## dimnames(table) <- list(fpr=sprintf('%.*f',2,fprs),visit=as.character(ts))

## cover.table <- ftable(table,row.vars=c(2))
## names(attr(cover.table,'row.vars')) <- '\\textbf{t}'
## names(attr(cover.table,'col.vars')) <- '\\textbf{FPR}'
## save.image('h.RData')
## ## bootstrap.rocs <- lapply(bootstrap.rocs,function(l)l[[1]])
## ## plot(roc.observed[[1]],bootstrap.rocs=bootstrap.rocs,alpha=.05)
## ## ll <- structure(c(roc.observed[[2]],roc.1.true[[2]]),class=c('list','PatientROC'))
## ## plot(ll)
