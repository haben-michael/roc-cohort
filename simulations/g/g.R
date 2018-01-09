

######################################################################
.libPaths(c("/usr/lib/R/library",'/home/habnice/R/lib/R/library'))
source('ragon.old.R')
rho <- .5
sigma <- 2
mu.theta.prior <- c(0,0)
sigma.theta.prior <- matrix(c(2,.4,.4,2),nrow=2)
mu.p.prior <- c(0,0)
sigma.p.prior <- matrix(c(1,.4,.4,1),nrow=2)
## eta <- rbind(ar=c(mu.theta.prior,sigma.theta.prior[c(1,2,4)]),mc=c(mu.p.prior,sigma.p.prior[c(1,2,4)]))
## xi <- c(theta,p)

N <-8e1#3e2
T <- 5e1

ts <- c(5,15,25,35)
fprs <- c(.1,.25,.5,.75)
alpha <- .05
q <- qnorm(1-alpha/2)

eta <- rbind(ar=c(mu.theta.prior,sigma.theta.prior[c(1,2,4)]),mc=c(mu.p.prior,sigma.p.prior[c(1,2,4)]))

## roc.3.true <- lapply(ts,function(t)roc.3(eta,rho,sigma,t,B=5e3))
## plot(roc.3.true)
## tpr.true <- sapply(roc.3.true,function(r)r(fprs))
## tpr.true
## require(MASS)
## write.matrix(tpr.true,'true.txt')


## testing...
## uu <- replicate(1e2,{
##     cat('.')
##     roc.3.true <- roc.3(eta,rho,sigma,5,B=5e1)
##     roc.3.true(fprs)
## })
## hist(uu[1,])
## abline(v=tpr.true[1,1],col='red')
    
## tprs <- replicate(1e2,
## {theta <- rmvnorm(N,mean=mu.theta.prior,sigma=sigma.theta.prior)
##     cat('.')
## p <- rmvnorm(N,mean=mu.p.prior,sigma=sigma.p.prior) %>% plogis
## data <- rpatient(rho=rho,sigma=sigma,theta=theta,p=p,T=T)
## roc.observed <- roc.3(data,rho,sigma,ts=35)
## roc.observed(fprs)
## })
## hist(tprs[4,])
## abline(v=tpr.true[4,4],col='red')
## abline(v=mean(tprs[4,]),col='blue')




bootstrap.reps <- 1e2
reps <- 1
roc.3.sim <- list()

theta <- rmvnorm(N,mean=mu.theta.prior,sigma=sigma.theta.prior)
p <- rmvnorm(N,mean=mu.p.prior,sigma=sigma.p.prior) %>% plogis
data <- rpatient(rho=rho,sigma=sigma,theta=theta,p=p,T=T)
## eta.star <- pop.fit(data,rho=rho)

roc.observed <- roc.3(data,rho,sigma,ts=ts)

bootstrap.rocs <- replicate(bootstrap.reps,{
    cat('.')
    tryCatch({
        ## weights <- rexp(N)
        ## roc.3(data,rho=rho,sigma=sigma,t=t,weights=weights)
        roc.3(data[sample(1:N,replace=T)],rho=rho,sigma=sigma,ts=ts)

    },error=function(e)NA)
},simplify=F)

tpr.true <- as.matrix(read.table('true.txt',as.is=TRUE))

tpr.observed <- sapply(roc.observed,function(r)r(fprs))
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
roc.3.sim <- list(coverage.1=coverage.1,coverage.3=coverage.3,bias=bias,boot.sds=boot.sds,mu.star=tpr.observed)   

save(roc.3.sim,file=paste0('g-',round(runif(1)*1e8),'.RData'))


## from barley
## source('../../ragon.R')
## filelist <- dir()
## filelist <- filelist[grep('g-[0-9]{8}.RData$',filelist)]
## k <- 1
## roc.3.sim.full <- lapply(filelist, function(file) {
##     load(file)
##     roc.3.sim
## })
## roc.3.sim <- roc.3.sim.full
## coverages <- lapply(roc.3.sim,function(l)l$coverage.3)
## cc <- abind(coverages,along=3)
## cc3 <- apply(cc,1:2,mean)


## biases <- lapply(roc.3.sim,function(l)l$bias)
## bb <- abind(biases,along=3)
## bb3 <- apply(bb,1:2,mean)


## ss <- lapply(roc.3.sim,function(r)r$boot.sds)
## ss <- abind(ss,along=3)
## ss3 <- apply(ss,1:2,mean)

## mm <- lapply(roc.3.sim,function(r)r$mu.star)
## mm <- abind(mm,along=3)
## mm3 <- apply(mm,1:2,sd)

## table <- abind(cc3,bb3,ss3,mm3,along=3) %>%
##     apply(1:3,function(x)sprintf('%.*f',2,x)) %>%
##     apply(1:2,function(x)paste0(x[1],' (',x[2],',',x[3],',',x[4],')'))
## rownames(table) <- fprs;colnames(table) <- ts
## dimnames(table) <- list(fpr=sprintf('%.*f',2,fprs),visit=as.character(ts))
## table

## cover.table <- ftable(table,row.vars=c(2))
## names(attr(cover.table,'row.vars')) <- '\\textbf{t}'
## names(attr(cover.table,'col.vars')) <- '\\textbf{FPR}'
## ## save.image('g.RData')
## load('g.RData')
## require(xtable)
## caption0 <- 'Misspecified model: Nominal 95\\% CI coverage, bias, mean standard error, and empirical standard deviation of $ROC_3(u|t)$ for FPRs 10\\%, 25\\%, 50\\%, and 75\\% at visits 5, 15, 25, and 35 (random slope/intercept logistic model, $N=80$ patients).'
## ## sink('revision_e.tex')
## print(xtableFtable(cover.table,label='tab:2',caption=caption0),include.rownames=F,include.colnames=F,quote=F,sanitize.rownames.function=function(x)x,sanitize.colnames.function=function(x)x)
## sink()
 
 
