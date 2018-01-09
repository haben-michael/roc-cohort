## .libPaths(c('/afs/ir/users/h/a/habnice/Desktop/R/lib/R/library',"/usr/lib/R/library"))
source('../../ragon.R')
setwd('/dropbox/ragon/simulations/n')
N <- 1e2
T <- 1e2
t <- 5
sigma <- 1
sd.alpha <- 1
sd.beta <- 1
mean.x <- .5
sd.x <- 1

T <- 2e2
x <- matrix(rnorm(N*T),ncol=T)
alpha <- rnorm(N,sd=sd.alpha)
beta <- rnorm(N,2,sd=sd.beta)
x <- matrix(rnorm(N*T,mean.x,sd.x),ncol=T)
p <- plogis(alpha + diag(beta)%*%x)

data <- data.frame(p=as.numeric(t(p)),x=as.numeric(t(x)),id=gl(N,T),t=rep(1:T,N))
data$y <- rbinom(N*T,1,prob=data$p)
## boxplot(x~y,data)
data.list <- split(data,data$id)
data.list <- lapply(data.list,function(df)list(X=df$x,S=df$y))
attributes(data.list) <- list(class='PatientClass',N=N,T=T)
true.plt <- plot(roc(data$x,data$y,interpolator=TRUE))


T <- 1e1
x <- matrix(rnorm(N*T),ncol=T)
alpha <- rnorm(N,sd=sd.alpha)
beta <- rnorm(N,2,sd=sd.beta)
x <- matrix(rnorm(N*T,mean.x,sd.x),ncol=T)
p <- plogis(alpha + diag(beta)%*%x)

data <- data.frame(p=as.numeric(t(p)),x=as.numeric(t(x)),id=gl(N,T),t=rep(1:T,N))
data$y <- rbinom(N*T,1,prob=data$p)
## boxplot(x~y,data)
data.list <- split(data,data$id)
data.list <- lapply(data.list,function(df)list(X=df$x,S=df$y))
attributes(data.list) <- list(class='PatientClass',N=N,T=T)


glmer0 <- glmer(y ~ x|id, data=data, family=binomial)
fits <- fitted(glmer0)
roc.fitted <- roc(fits,data$y,interpolator=T)
bootstrap.rocs <- replicate(2e2, {
    tryCatch({
        bootstrap.ids <- sample(unique(data$id),replace=TRUE)
        data.list.bs <- data.list[bootstrap.ids]
        data.list.bs <- lapply(data.list.bs,as.data.frame)
        data.bs <- do.call(rbind,data.list.bs)
        data.bs$id <- rep(bootstrap.ids,each=T)
        glmer.bs <- glmer(S ~ X|id, data=data.bs, family=binomial)
        roc(fitted(glmer.bs),data.bs$S,interpolator=TRUE)
    },error=function(e)NA)
},simplify=FALSE)
bootstrap.rocs <- bootstrap.rocs[!is.na(bootstrap.rocs)]
length(bootstrap.rocs)
roc.fitted.plt <- plot(roc.fitted,bootstrap.rocs=bootstrap.rocs,alpha=.05)
gg.df <- attr(true.plt,'gg.df')
gg.df$upper <- gg.df$lower <- gg.df$tpr
gg.df <- rbind(attr(roc.fitted.plt,'gg.df'),gg.df)
gg.df$ROC <- rep(c('GLMM','true'),c(nrow(attr(true.plt,'gg.df')),nrow(attr(roc.fitted.plt,'gg.df'))))
plt <- ggplot(gg.df,aes(x=fpr,y=tpr,group=ROC))+geom_line(aes(linetype=ROC))+geom_abline(slope=1,col='gray',linetype=2)+theme_classic()+labs(x='FPR',y='TPR')+geom_ribbon(aes(ymax=upper,ymin=lower),alpha=.5,fill='lightgray')

png('../../figures/revision_p.png')
plt
dev.off()



table <- abind(cc3,bb3,ss3,mm3,along=3) %>%
    apply(1:3,function(x)sprintf('%.*f',2,x)) %>%
    apply(1:2,function(x)paste0(x[1],' (',x[2],',',x[3],',',x[4],')'))
rownames(table) <- fprs;colnames(table) <- ts
dimnames(table) <- list(fpr=sprintf('%.*f',2,fprs),visit=as.character(ts))
table
table <- cbind(table,rep('',nrow(table)))
table <- table[,c(ncol(table),1:(ncol(table)-1))]
## colnames(table)[1] <- ''
addtorow <- list()
addtorow$pos <- list(0,0)
addtorow$command <- c('&\\textbf{FPR}& \\multicolumn{1}{c}{5}&  \\multicolumn{1}{c}{15} &  \\multicolumn{1}{c}{25} &  \\multicolumn{1}{c}{35}\\\\\n',' \\textbf{visit} & & CVL (BS,MSE,ESE) & CVL (BS,MSE,ESE) &CVL (BS,MSE,ESE) & CVL (BS,MSE,ESE)\\\\\n')
align <- rep('p{.25\\textwidth}',6)
align[1] <- align[2] <- 'p{.4cm}'
sink('test.tex')
print(xtable(table,align=align),add.to.row=addtorow,include.colnames=FALSE)
sink()


## cover.table <- ftable(table,row.vars=c(2))
## names(attr(cover.table,'row.vars')) <- '\\textbf{t}'
## names(attr(cover.table,'col.vars')) <- '\\textbf{FPR}'
## ## save.image('g.RData')
## load('g.RData')
## require(xtable)
## caption0 <- 'Misspecified model: Nominal 95\\% CI coverage, bias, mean standard error, and empirical standard deviation of $ROC_3(u|t)$ for FPRs 10\\%, 25\\%, 50\\%, and 75\\% at visits 5, 15, 25, and 35 (random slope/intercept logistic model, $N=80$ patients).'
## sink('test.tex')
## print(xtableFtable(cover.table,label='tab:2',caption=caption0),include.rownames=F,include.colnames=F,quote=F,sanitize.rownames.function=function(x)x,sanitize.colnames.function=function(x)x)
## sink()


## table <- abind(cc3,bb3,ss3,mm3,along=3) %>%
##     apply(1:3,function(x)sprintf('%.*f',2,x)) %>%
##     apply(1:2,function(x)paste0(x[1],' (',x[2],',',x[3],',',x[4],')'))
## rownames(table) <- fprs;colnames(table) <- ts
## dimnames(table) <- list(fpr=sprintf('%.*f',2,fprs),visit=as.character(ts))
## names(attr(table,'row.vars')) <- '\\textbf{visit}'
## names(attr(table,'col.vars')) <- '\\textbf{FPR}'
## dimnames(table) <- list(fpr=fprs,visit=ts)
## print(xtable(table))

## table <- as.data.frame(ftable(table,row.vars=c(2)))
## ftable(table,row.vars = 3)

## names(attr(cover.table,'row.vars')) <- '\\textbf{visit}'
## names(attr(cover.table,'col.vars')) <- '\\textbf{FPR}'

## ## save.image('g.RData')
## ## load('g.RData')
## require(xtable)
## caption0 <- 'Misspecified model: Nominal 95\\% CI coverage, bias, mean standard error, and empirical standard deviation of $ROC_3(u|t)$ for FPRs 10\\%, 25\\%, 50\\%, and 75\\% at visits 5, 15, 25, and 35 (random slope/intercept logistic model, $N=80$ patients).'
## sink('test.tex')
## print(xtableFtable(cover.table,label='tab:2',caption=caption0),add.to.row=addtorow,include.rownames=F,include.colnames=FALSE,quote=F,sanitize.rownames.function=function(x)x,sanitize.colnames.function=function(x)x)
## sink()


## res <- replicate(20,{
## n <- 20
## A <- rbinom(n,1,1/2)
## Y <- rnorm(n)
## u <- (mean(A*Y)/mean(A) - mean(Y))/(1-mean(A))
## v <- mean(Y)/2-mean(A*Y)
## c(u,v)
## }
## )
