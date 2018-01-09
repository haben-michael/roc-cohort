.libPaths(c('/afs/ir/users/h/a/habnice/Desktop/R/lib/R/library',"/usr/lib/R/library"))
source('../../ragon.R')

N <- 1e2
T <- 5e3
t <- 5
ts <- c(5,15,25,35)
rho <- .5
sigma <- 1

## i <- 1
##     set.seed(i)
##     cat('.')
##     i <- i+1
alpha <- rnorm(N,sd=2)
beta <- rnorm(N,.5,sd=1)
x <- matrix(rnorm(N*T),ncol=T)
alpha <- rnorm(N,sd=2)
beta <- rnorm(N,1,sd=1)
x <- matrix(rnorm(N*T,2,1),ncol=T)
p <- plogis(alpha + diag(beta)%*%x)

data <- data.frame(p=as.numeric(t(p)),x=as.numeric(t(x)),id=gl(N,T),t=rep(1:T,N))
data$y <- rbinom(N*T,1,prob=data$p)
## boxplot(x~y,data)
data.list <- split(data,data$id)
data.list <- lapply(data.list,function(df)list(X=df$x,S=df$y))
attributes(data.list) <- list(class='PatientClass',N=N,T=T)
roc.3.true <- roc.3(data.list,rho,sigma,ts=ts)




## eta <- pop.fit(data.list,rho)
## ## sigma.theta.prior <- matrix(eta['ar',c(3,4,4,5)],nrow=2)
## ## xi <- tryCatch(subj.fit(data.list,rho,sigma,eta,t=50),error=function(e)e)
## xi <- subj.fit(data.list,rho,sigma,eta,t=50)

glmer0 <- glmer(y ~ x|id, data=data, family=binomial)
plot(fitted(glmer0) ~ data$p)
fits <- fitted(glmer0)[data$t<=t]
plt1 <- plot(roc(fits,data$y[data$t==t],interpolator=T))
## plt1
plt2 <- plot(roc(fitted(glmer0),data$y,interpolator=T))
roc.fitted <- roc(fits,data$y[data$t<=t],interpolator=T)
plot(roc.fitted)


## save.image('081817d.RData')


fpr <- seq(0,1,length=1e2)
tprs <- lapply(roc.3.true,function(r)r(fpr))
tprs <- lapply(1:length(ts),function(i)data.frame(fpr=fpr,tpr=tprs[[i]],t=ts[i]))
tprs.fitted <- data.frame(fpr=fpr,tpr=roc.fitted(fpr),t='fitted')
tprs.fitted <- rbind(data.frame(fpr=0,tpr=0,t='fitted'),tprs.fitted)
gg.df <- do.call(rbind,tprs)
gg.df <- rbind(gg.df,tprs.fitted)

## plt <-
##     ggplot(gg.df,aes(x=fpr,y=tpr,color=factor(t)))+geom_line()+geom_abline(slope=1,col='lightblue',linetype=2)+scale_color_manual('',breaks=c(ts,'fitted'),values=cm.colors(5),labels=c(expression(paste(ROC[3],'(u|5)')),expression(paste(ROC[3],'(u|15)')),expression(paste(ROC[3],'(u|25)')),expression(paste(ROC[3],'(u|30)')),'benchmark'))+ theme_classic()
## ggsave('081817d.png',plt)

plt <- ggplot(gg.df,aes(x=fpr,y=tpr,linetype=factor(t)))+geom_line()+geom_abline(slope=1,col='gray',linetype=2)+scale_linetype_manual('',breaks=c(ts,'fitted'),values=factor(5:1),labels=c(expression(paste(ROC[3],'(u|5)')),expression(paste(ROC[3],'(u|15)')),expression(paste(ROC[3],'(u|25)')),expression(paste(ROC[3],'(u|35)')),'benchmark'))+  theme_classic()+labs(x='FPR',y='TPR')
## ggsave('081817d.png',plt)
## ggsave('figures/revision_d.png,plt)


## c(expression(paste(ROC[3],'(u|5)')),expression(paste(ROC[3],'(u|15)')),expression(paste(ROC[3],'(u|25)')),'benchmark')
## dd <- sapply(ts,function(t)substitute(paste(ROC[3],'(u|',t,')'),list(t=t)))
plt
