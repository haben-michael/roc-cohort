.libPaths(c('/afs/ir/users/h/a/habnice/Desktop/R/lib/R/library',"/usr/lib/R/library"))
## coverage/bias under misspecification
source('ragon.R')

N <- 1e2
T <- 1e2
t <- 5
rho <- .5
ts <- c(5,15,25,30)
fprs <- c(.1,.25,.5,.75)
q <- qnorm(1-.05/2)
sigma <- 1

gen.data <- function(N) {
    alpha <- rnorm(N,sd=2)
    beta <- rnorm(N,2,sd=1)
    x <- matrix(rnorm(N*T,2,4),ncol=T)
    p <- plogis(alpha + diag(beta)%*%x)

    data <- data.frame(p=as.numeric(t(p)),x=as.numeric(t(x)),id=gl(N,T),t=rep(1:T,N))
    data$y <- rbinom(N*T,1,prob=data$p)
    ## boxplot(x~y,data)
    data.list <- split(data,data$id)
    data.list <- lapply(data.list,function(df)list(X=df$x,S=df$y))
    attributes(data.list) <- list(class='PatientClass',N=N,T=T)
    return(data.list)
}


## data <- gen.data(1e3)
## roc.3.true <- roc.3(data,rho,sigma,ts=ts)
## true <- sapply(roc.3.true,function(r)r(fprs))
## require(MASS)
## write.matrix(true,'true.txt')


N <- 1e2
bootstrap.reps <- 3e2
reps <- 2
roc.3.sim <- list()
for (i in 1:reps) {
    data <- gen.data(N)
    roc.3.star <- tryCatch(roc.3(data,rho,sigma,ts=ts),error=function(e)NA)
    if (sum(is.na(roc.3.star))>0) next

    rocs <- replicate(bootstrap.reps,{
        cat('.')
        tryCatch({
            ## weights <- rexp(N)
            ## roc.3(data,rho=rho,sigma=sigma,t=t,weights=weights)
            data.star <- data[sample(1:N,replace=T)]
            roc.3(data.star,rho=rho,sigma=sigma,ts=ts)

            },error=function(e){traceback();NA})
    },simplify=F)
    rocs <- rocs[sapply(rocs,function(r)sum(is.na(r)))==0]

    true <- as.matrix(read.table('true.txt',as.is=TRUE))
    boot.sds <- abind(lapply(rocs,function(r)
        sapply(r,function(rr)rr(fprs))), along=3) %>% apply(1:2,sd)
    rownames(boot.sds) <- fprs; colnames(boot.sds) <- ts

    mu.star <- sapply(roc.3.star,function(r)r(fprs))
    upper <- mu.star+q*boot.sds
    lower <- mu.star-q*boot.sds
    coverage <- true<=upper & true>=lower
    bias <- abs(mu.star-true)

    roc.3.sim[[i]] <- list(coverage=coverage,bias=bias)
    ## roc.3.sim[[i]] <- rocs

}

save.image(paste0('081817e-',round(runif(1)*1e8),'.RData'))


## ## from corn
## filelist <- dir()
## filelist <- filelist[grep('081817e-[0-9]{8}.RData$',filelist)]
## roc.3.sim.full <- list()
## k <- 1
## for(file in filelist) {
##     print(k)
##     k <- k+1
##     load(file)
##     roc.3.sim.full <- c(roc.3.sim.full,roc.3.sim)
## }

## roc.3.sim <- roc.3.sim.full
## ## roc.3.sim <- roc.3.sim[sapply(roc.3.sim,length)>5]

## ## roc.3.sim <- roc.3.sim[!sapply(roc.3.sim,function(x)is.null(x))]
## ## length(roc.3.sim)

## cc <- lapply(roc.3.sim,function(r)r$coverage)
## cc <- abind(cc,along=3)
## cc3 <- apply(cc,1:2,mean)
## cc3


## bb <- lapply(roc.3.sim,function(r)r$bias)
## bb <- abind(bb,along=3)
## bb3 <- apply(bb,1:2,mean)

## table <- abind(cc3,bb3,along=3) %>%
##     apply(1:2,function(x)paste0(sprintf('%.*f',2,x[1]),' (',round(x[2],3),')'))
## rownames(table) <- fprs;colnames(table) <- ts
## dimnames(table) <- list(fpr=sprintf('%.*f',2,fprs),visit=as.character(ts))

## cover.table <- ftable(table,row.vars=c(2))
## names(attr(cover.table,'row.vars')) <- '\\textbf{t}'
## names(attr(cover.table,'col.vars')) <- '\\textbf{FPR}'
## require(xtable)
## caption0 <- 'Misspecified model: Nominal 95\\% CI coverage and bias of $ROC_3(u|t)$ for FPRs 10\\%, 25\\%, 50\\%, and 75\\% at visits 3, 5, and 9 (random slope/intercept logistic model, $N=300$ patients).'
## sink('081817b.tex')
## print(xtableFtable(cover.table,label='tab:2',caption=caption0),include.rownames=F,include.colnames=F,quote=F,sanitize.rownames.function=function(x)x,sanitize.colnames.function=function(x)x)
## sink()
 
