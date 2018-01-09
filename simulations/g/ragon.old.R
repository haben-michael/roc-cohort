require(abind)
require(ggplot2)
require(plyr)
require(dplyr)
require(tidyr)
require(mvtnorm)
require(lme4)
require(stringr)
require(grid)
require(gridExtra)
require(RColorBrewer)
require(nnet)


rpatient <- function(rho,sigma,theta,p,T) {
    stopifnot(ncol(theta)==2 & ncol(p)==2)
    params <- data.frame(rho=rho,sigma=sigma,theta1=theta[,1],theta2=theta[,2],p.11=p[,1],p.22=p[,2],T=T)
    ans <- apply(params,1,function(param) {
        trans.mtx <- matrix(c(param['p.11'],1-param['p.22'],1-param['p.11'],param['p.22']),nrow=2)
        S <- sample(0:1,1)
        for(t in 2:as.numeric(param['T'])) S[t] <- sample(0:1,1,prob=trans.mtx[S[t-1]+1,])

        P <- diag(param['T'])
        P <- param['rho']^(row(P)-col(P))
        P[col(P)>row(P)] <- 0

        theta <- param[c('theta1','theta2')]
        Theta <- theta[S+1]
        attr(Theta,'dim' ) <- c(length(Theta),1)

        X <- rmvnorm(1,mean=P%*%Theta,sigma=P%*%t(P)*param['sigma']^2)%>%as.numeric

        return(list(X=X,S=S))
    })
    return(structure(ans,class='PatientClass',T=params$T,N=length(ans)))
}

## rename PatientClass -> Patient throughout
plot.PatientClass <- function(dat){
    gg.df <- lapply(1:length(dat),
                function(idx)data.frame(X=dat[[idx]]$X,S=as.factor(dat[[idx]]$S),patient=idx,time=1:length(dat[[idx]]$X)))
    gg.df <- do.call(rbind,gg.df)
    ## gg.df <- data.frame(X=dat$X,S=as.factor(dat$S),time=1:length(dat$X))
    plt <- ggplot(gg.df,aes(x=time,y=X,group=patient))+geom_point(aes(color=S))+geom_line()
    return(plt)
}

`[.PatientClass` <- function(x, n, ...) {
    attrs <- attributes(x)
    out <- unclass(x)
    out <- out[n]
    if(length(list(...))>0) {
        t <- as.numeric(...)
        out <- lapply(out,function(dat.i)list(X=dat.i$X[t],S=dat.i$S[t]))
    }
    attrs$names <- attrs$names[n]
    attributes(out) <- attrs
    attr(out,'N') <- length(out)
    attr(out,'T') <- sapply(out,function(df)length(df$X))
    out
}

surv.positive <- function(quantile,theta,sigma,rho,p,T,pi0=rep(1/2,2)) { ## make asymptotic
    alpha0 <- -diff(theta^2)/(2*sigma^2)+diff(log(p[1,]))
    alpha1 <- log(prod(diag(p))/(p[1,2]*p[2,1]))
    alpha2 <- -rho*diff(theta)/sigma^2
    alpha3 <- diff(theta)/sigma^2

    eigen0 <- eigen(p)
    prior <- eigen0$vectors %*% diag(eigen0$values^(T-1)) %*% solve(eigen0$vectors)
    prior <- pi0%*%prior

    tpr <- pnorm(((quantile-alpha0-alpha1*c(0,1))/alpha3 - theta[2])/sigma,lower.tail=F) %*% t(p[,2]*prior/sum(p[,2]*prior))
    fpr <- pnorm(((quantile-alpha0-alpha1*c(0,1))/alpha3 - theta[1])/sigma,lower.tail=F) %*% t(p[,1]*prior/sum(p[,1]*prior))

    return(c(tpr=tpr,fpr=fpr))
}

plot.PatientROC <- function(roc,resolution=1e2,alpha=0,bootstrap.rocs=list(),diag=TRUE) {
    grid <- seq(1/resolution,1-1/resolution,length.out=resolution)
    gg.df <- data.frame(fpr=grid) %>% mutate(tpr=roc(fpr))
    gg.df <- rbind(data.frame(fpr=0,tpr=0),gg.df)
    if(alpha==0) {
        ggplot(gg.df,aes(x=fpr,y=tpr))+geom_line()+geom_abline(slope=1,col='lightblue',linetype=2)+ylim(0,1)
    } else {
        tpr.bs <- sapply(bootstrap.rocs,function(bootstrap.roc)bootstrap.roc(grid))
        sds <- apply(tpr.bs,1,sd)
        q <- qnorm(1-alpha/2)
        gg.df <- mutate(gg.df,upper=tpr+q*sds,lower=tpr-q*sds)
        ggplot(gg.df,aes(x=fpr,y=tpr))+geom_line()+geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.5,fill='lightblue')+geom_abline(slope=1,col='lightblue',linetype=2)+ylim(0,1)
    }
}

## surv.positive <- Vectorize(surv.positive,'quantile')
## ## surv.positive(5,c(0,1),1,matrix(rep(.5,4),nrow=2),T=30,pi0=c(.5,.5))
## surv.positive(2,theta,sigma,rho,p,T,pi0=rep(1/2,2))

## renamed from "ppatient" 1/9/17--reflects ms notation
ppatient.inf <- function(quantile,xi,sigma,rho,status='both',lower.tail=T) {
    theta <- xi[1:2]; p <- xi[3:4]
    p <- matrix(c(p[1],1-p[2],1-p[1],p[2]),nrow=2)
    alpha0 <- -diff(theta^2)/(2*sigma^2)+diff(log(p[1,]))
    alpha1 <- log(prod(diag(p))/(p[1,2]*p[2,1]))
    alpha2 <- -rho*diff(theta)/sigma^2
    alpha3 <- diff(theta)/sigma^2

    ## eigen0 <- eigen(p)
    ## prior <- eigen0$vectors %*% diag(eigen0$values^(T-1)) %*% solve(eigen0$vectors)
    ## prior <- pi0%*%prior
    pi <- c(p[2,1],p[1,2])
    pi <- pi/sum(pi)

    tpr <- pnorm(((quantile-alpha0-alpha1*c(0,1))/alpha3 - theta[2])/sigma,lower.tail=diff(theta)<0) %*% (p[,2]*pi)/pi[2]
    fpr <- pnorm(((quantile-alpha0-alpha1*c(0,1))/alpha3 - theta[1])/sigma,lower.tail=diff(theta)<0) %*% (p[,1]*pi)/pi[1]

    if(lower.tail) {
        if(status=='both') {
            return(c(fnr=1-tpr,tnr=1-fpr))
        } else {
            if(status=='positive') return(c(fnr=1-tpr)) else return(c(tnr=1-fpr))
        }
    } else {
        if(status=='both') {
            return(c(tpr=tpr,fpr=fpr))
        } else {
            if(status=='positive') return(c(tpr=tpr)) else return(c(fpr=fpr))
        }
    }
}

ppatient.inf <- Vectorize(ppatient.inf,'quantile')

qpatient.inf <- function(p,xi,sigma,rho,status='positive',lower.tail=T,tol=1e-4) {
    upper <- 1; lower <- -1
    if(!lower.tail) p <- 1-p
    repeat {
        p.star <- ppatient.inf(upper,xi,sigma,rho,status=status)
        if(p.star<p) {
            lower <- upper
            upper <- upper*2
        } else break
    }
    repeat {
        p.star <- ppatient.inf(lower,xi,sigma,rho,status=status)
        if(p.star>p) {
            lower <- lower*2
        } else break
    }
    ## upper <- 20
    ## lower <- -20
    repeat {
        q.star <- (upper+lower)/2
        p.star <- ppatient.inf(q.star,xi,sigma,rho,status=status)
        if(abs(p.star-p)<tol) break
        if(p.star>p) upper <- q.star else lower <- q.star
    }
    return((upper+lower)/2)
}
qpatient.inf <- Vectorize(qpatient.inf,'p')


## conditioning on subject's history
## ppatient2 <- function(quantile,S.prev,theta,p,sigma,rho) {
##     ## p <- matrix(c(p[1],1-p[2],1-p[1],p[2]),nrow=2)
##     theta <- as.numeric(theta)
##     alpha0 <- -diff(theta^2)/(2*sigma^2)+diff(log(p[1,]))
##     alpha1 <- log(prod(diag(p))/(p[1,2]*p[2,1]))
##     alpha3 <- diff(theta)/sigma^2
##     gamma <- (quantile-alpha0-alpha1*S.prev)/alpha3
##     fpr <- pnorm((gamma-theta[1])/sigma,lower.tail=theta[1]>theta[2])
##     tpr <- pnorm((gamma-theta[2])/sigma,lower.tail=theta[1]>theta[2])
##     c(fpr=fpr,tpr=tpr)
## }
## ppatient2 <- Vectorize(ppatient2,'quantile')

## curve(ppatient(x,theta=c(0,1),sigma,rho,p,status='positive',lower.tail=F),-10,5)

## renamed from "dpatient" 1/9/17--reflects ms notation
dpatient.inf <- function(quantile,xi,sigma,rho,status=0) {
    theta <- xi[1:2]; p <- xi[3:4]
    p <- matrix(c(p[1],1-p[2],1-p[1],p[2]),nrow=2)
    alpha0 <- -diff(theta^2)/(2*sigma^2)+diff(log(p[1,]))
    alpha1 <- log(prod(diag(p))/(p[1,2]*p[2,1]))
    alpha2 <- -rho*diff(theta)/sigma^2
    alpha3 <- diff(theta)/sigma^2

    pi <- c(p[2,1],p[1,2])
    pi <- pi/sum(pi)

    positive.status <-  dnorm(((quantile-alpha0-alpha1*c(0,1))/alpha3 - theta[2])/sigma) %*% (p[,2]*pi)/pi[2] / (sigma*alpha3)
    negative.status <-  dnorm(((quantile-alpha0-alpha1*c(0,1))/alpha3 - theta[1])/sigma) %*% (p[,1]*pi)/pi[1] / (sigma*alpha3)

    return(sign(diff(theta))*c(negative.status=negative.status,positive.status=positive.status)[status+1])
}
dpatient.inf <- Vectorize(dpatient.inf,'quantile')
## theta <- c(1.18,1.19)
## curve(dpatient(x,theta,sigma,rho,p,0),-3,3)

biomarker <- function(dat,theta,p,rho,sigma,T) {
    sapply(1:length(dat),function(i) {
        p.i <- matrix(c(p[i,1],1-p[i,2],1-p[i,1],p[i,2]),nrow=2)
        theta.i <- theta[i,]
        alpha0 <- -diff(theta.i^2)/(2*sigma^2)+diff(log(p.i[1,]))
        alpha1 <- log(prod(diag(p.i))/(p.i[1,2]*p.i[2,1]))
        alpha2 <- -rho*diff(theta.i)/sigma^2
        alpha3 <- diff(theta.i)/sigma^2

        alpha0 + alpha1*dat[[i]]$S[T-1] + alpha2*dat[[i]]$X[T-1] + alpha3*dat[[i]]$X[T]
    })
}

## renamed from "auc" 1/9/17--reflects ms notation
## returns auc of 0 when theta1<theta0
auc.2 <- function(xi,sigma,rho,tolerance=5e-10  ) {
    xi <- matrix(xi,ncol=4)
    theta <- matrix(xi[,1:2],ncol=2); p <- matrix(xi[,3:4],ncol=2)
    ## if(dim(theta)==1) theta <- matrix(theta,nrow=1)
    ## if(dim(p)==1) p <- matrix(p,nrow=1)
    sapply(1:nrow(xi),function(i) {
        if(abs(diff(theta[i,]))<tolerance){return(NA)} ## tolerance for numerical integra
        ## if(diff(theta[i,])<tolerance)return(.5) ## tolerance for numerical integration
        integrate(Vectorize(function(x)(1-ppatient.inf(x,xi[i,],sigma,rho,status='positive'))*dpatient.inf(x,xi[i,],sigma,rho,0),'x'),-Inf,Inf)$value
    })
}

mc.subj.posterior <- function(x,...) {
    UseMethod('mc.subj.posterior',x)
    }
mc.subj.posterior.default <- function(n.0,n.1,n.00,n.11,mu.p.prior,sigma.p.prior) {
    post <- function(par){##,n.0,n.1,n.00,n.11) {
        p.00 <- par[1]; p.11 <- par[2]
        if(p.00>=1||p.11>=1||p.00<=0||p.11<=0)return(-.Machine$double.xmax)
        (n.0-n.00)*log(1-p.00)+n.00*log(p.00)+(n.1-n.11)*log(1-p.11)+n.11*log(p.11)-.5*t(qlogis(par)-mu.p.prior)%*%solve(sigma.p.prior)%*%(qlogis(par)-mu.p.prior)
    }
    ## cat(n.0,' ',n.1,' ',n.00,' ',n.11,'\n')
    optim(par=c(.5,.5),fn=post,control=list(fnscale=-1))$par
}
## mc.subj.posterior(10,10,5,8,sigma.p.prior,mu.p.prior)

mc.subj.posterior.PatientClass <- function(dat,mu.p.prior,sigma.p.prior) {
    sapply(dat,function(dat.i) {
        ## if (sum(dat.i$S)==0) return(c(1,0))
        ## if (sum(dat.i$S)==length(dat.i$S)) return(c(0,1))
        S.prev <- dat.i$S[-length(dat.i$S)]
        S.next <- dat.i$S[-1]
        ## print(S.prev);print(S.next)
        n.0 <- sum(S.prev==0)
        n.1 <- sum(S.prev==1)
        n.00 <- sum((1-S.prev)*(1-S.next))
        n.11 <- sum(S.prev*S.next)
        mc.subj.posterior(n.0,n.1,n.00,n.11,mu.p.prior,sigma.p.prior)
        })
}

ar.subj.posterior <- function(dat,rho,sigma,mu.theta.prior,sigma.theta.prior) {
    t <- length(dat[[1]]$X)
    P <- diag(t)
    P <- rho^(row(P)-col(P))
    P[col(P)>row(P)] <- 0
    R <- P%*%t(P)*sigma^2

    sapply(dat,function(dat.i) {
        ## X <- dat.i$X[1:t]; S <- dat.i$S[1:t]
        ## X <- dat.t[[1]]$X; S <- dat.t[[1]]$S
        X <- dat.i$X; S <- dat.i$S
        S <- matrix(as.numeric(c(S==0,S==1)),ncol=2)
        sig.inv <- solve(sigma.theta.prior)
        D <- t(S)%*%solve(P)%*%X + sig.inv%*%mu.theta.prior
        C.inv <- solve(t(S)%*%S+sig.inv)
        C.inv%*%D
    })
}

subj.fit <- function(dat,rho,sigma,mu.theta.prior,sigma.theta.prior,
                                     mu.p.prior,sigma.p.prior,Ts) {


    res <- sapply(Ts, function(t) {
        ## t <- 5
        ## param <- c(0,3,1,.3,1)

        ## dat.t <- structure(lapply(dat,function(dat.i)list(X=dat.i$X[1:t],S=dat.i$S[1:t])), class='PatientClass')
        dat.t <- dat[,1:t]

        ## ar.params <- tryCatch(ar.subj.posterior(dat.t,rho,sigma,mu.theta.prior,sigma.theta.prior),
        ##                       error=function(e){print(e);rep(NA,5)},finally={})
        ## mc.params <- tryCatch(mc.subj.posterior(dat.t,mu.p.prior,sigma.p.prior),
        ##                       error=function(e){print(e);rep(NA,5)},finally={})
        ar.params <- ar.subj.posterior(dat.t,rho,sigma,mu.theta.prior,sigma.theta.prior)
        mc.params <- mc.subj.posterior(dat.t,mu.p.prior,sigma.p.prior)

        rbind(ar.params,mc.params)

    })

    res <- lapply(1:length(Ts),function(t)matrix(res[,t],byrow=T,ncol=4))
    lapply(res,function(m){colnames(m) <- c('theta0','theta1','p00','p11');m})

    ## colnames(res) <- paste0('T=',Ts)
    ## rownames(res) <- rep(c('theta0','theta1','p00','p11'),length(dat))
    ## ## res <- lapply(1:length(dat),function(n)res[(4*(n-1)+1):(4*n),])
    ## res <- lapply(1:length(dat),function(n)t(res[(4*(n-1)+1):(4*n),]))
    ## return(res)


}

subj.fit <- function(dat,rho,sigma,eta,Ts) {

    mu.theta.prior <- eta['ar',1:2]; mu.p.prior <- eta['mc',1:2]
    sigma.theta.prior <- matrix(eta['ar',c(3,4,4,5)],nrow=2);
    sigma.p.prior <- matrix(eta['mc',c(3,4,4,5)],nrow=2)

    res <- sapply(Ts, function(t) {
        dat.t <- dat[,1:t]

        ar.params <- ar.subj.posterior(dat.t,rho,sigma,mu.theta.prior,sigma.theta.prior)
        mc.params <- mc.subj.posterior(dat.t,mu.p.prior,sigma.p.prior)

        rbind(ar.params,mc.params)

    })

    res <- lapply(1:length(Ts),function(t)matrix(res[,t],byrow=T,ncol=4))
    lapply(res,function(m){colnames(m) <- c('theta0','theta1','p00','p11');m})

}

## no longer vectorized in t/T
subj.fit <- function(dat,rho,sigma,eta,t) {

    mu.theta.prior <- eta['ar',1:2]; mu.p.prior <- eta['mc',1:2]
    sigma.theta.prior <- matrix(eta['ar',c(3,4,4,5)],nrow=2);
    sigma.p.prior <- matrix(eta['mc',c(3,4,4,5)],nrow=2)

    ## res <- sapply(Ts, function(t) {
        dat.t <- dat[,1:t]

        ar.params <- ar.subj.posterior(dat.t,rho,sigma,mu.theta.prior,sigma.theta.prior)
        mc.params <- mc.subj.posterior(dat.t,mu.p.prior,sigma.p.prior)

        params <- t(rbind(ar.params,mc.params))

    colnames(params) <- c('theta0','theta1','p00','p11')
    return(params)
    ## })

    ## res <- lapply(1:length(Ts),function(t)matrix(res[,t],byrow=T,ncol=4))
    ## lapply(res,function(m){colnames(m) <- ;m})
}

mc.pop.mle <- function(data,weights=rep(1,length(data))) {
    ## T <- attr(data,'T'); if(T<=0) return
    T <- sapply(data,function(dat.i)length(dat.i$X))
    N <- length(data); if(N<=0) return
    weights <- rep(weights,T-1)

    S.next <- lapply(data,function(data.i)data.i$S[-1])
    S <- lapply(data,function(data.i)data.i$S[-length(data.i$S)])
    data.glmer <- data.frame(S.next=unlist(S.next),S=unlist(S),patient=rep(1:N,T-1))
    tryCatch({
        glmer0 <- suppressWarnings(glmer(S.next ~ S + (S|patient),family=binomial,data=data.glmer,weights=weights))
        ## summary(glmer0)
        A <- matrix(c(-1,1,0,1),nrow=2)
        sigma.p.mle <- A%*%matrix(as.numeric(VarCorr(glmer0)$patient),nrow=2)%*%t(A)
        ## print(eigen(sigma.p.mle))
        mu.p.mle <- A%*%fixef(glmer0)
        c(theta0=mu.p.mle[1],theta1=mu.p.mle[2],sigma1=sigma.p.mle[1],sigma12=sigma.p.mle[2],sigma2=sigma.p.mle[4])
    },
    ## warning=function(w){rep(NA,5)},
    error=function(e){rep(NA,5)}
    )
}



## method of moments
mc.pop.mom <- function(data) {
    T <- attr(data,'T')
    ans <- sapply(data,function(data.i) {
        S.next <- data.i$S[-1]
        S.curr <- data.i$S[-T]
        p00.i <- sum(S.curr==0 & S.next==0)/sum(S.curr==0)
        p11.i <- sum(S.curr==1 & S.next==1)/sum(S.curr==1)
        c(p00.i,p11.i)
    })
    ans <- qlogis(ans)
    return(c(theta0=mean(ans[1,]),theta1=mean(ans[2,]),sigma1=var(ans[1,]),sigma12=cov(ans[1,],ans[2,]),sigma2=var(ans[2,])))
}

## ar population mle calling optimizer directly
ar.pop.mle <- function(dat,rho,sigma,init.param=c(0,1,1,.3,1)) {

    loglik <- function(param,dat,R) {
        mu.theta.prior <- matrix(param[1:2],ncol=1)
        sigma.theta.prior <- matrix(c(param[3],param[4],param[4],param[5]),nrow=2)
        sig.inv <- solve(sigma.theta.prior)
        if(det(sigma.theta.prior)<=0)return(.Machine$double.xmin)

        ans <- sapply(dat,function(dat.i) {
            X <- dat.i$X; S <- dat.i$S
            S <- matrix(as.numeric(c(S==0,S==1)),ncol=2)
            D <- t(S)%*%solve(P)%*%X + sig.inv%*%mu.theta.prior
            C.inv <- solve(t(S)%*%S+sig.inv)

            -log(det(sigma.theta.prior)) - t(mu.theta.prior)%*%sig.inv%*%mu.theta.prior + log(det(C.inv)) + t(D)%*%C.inv%*%D
        })
        return(sum(ans))
    }

    T <- attr(dat,'T'); if(T<=0) return
    P <- diag(T)
    P <- rho^(row(P)-col(P))
    P[col(P)>row(P)] <- 0
    R <- P%*%t(P)*sigma^2

    res <- optim(init.param,fn=loglik,dat=dat,R=R,lower=c(-Inf,-Inf,.1,-Inf,.1),method='L-BFGS-B',control=list(fnscale=-1))$par
    names(res) <- c('theta0','theta1','sigma1','sigma12','sigma2')
    return(res)
}

## ar population mle using lme
ar.pop.mle <- function(data,rho,weights=rep(1,length(data))) {
    ## T <- attr(data,'T');
    N <- length(data)
    Ts <- sapply(data,function(data.i)length(data.i$X))
    ## weights <- rep(weights,each=T-1)
    weights <- rep(weights,Ts-1)
    X <- sapply(data,function(data.i)data.i$X[-1]) %>% unlist
    S <- sapply(data,function(data.i)data.i$S[-1]) %>% unlist
    X.prev <- sapply(data,function(data.i)data.i$X[-length(data.i$X)]) %>% unlist
    data.lmer <- data.frame(Z=c(X-rho*X.prev),S=c(S),patient=rep(1:N,Ts-1))
    A <- matrix(c(1,1,0,1),nrow=2)
    tryCatch({
        lmer0 <- lmer(Z ~ S + (S|patient),data=data.lmer,weights=weights)
        sigma.theta.mle <- A%*%matrix(as.numeric(VarCorr(lmer0)$patient),nrow=2)%*%t(A)
        mu.theta.mle <- A%*%fixef(lmer0)
        c(theta0=mu.theta.mle[1],theta1=mu.theta.mle[2],sigma1=sigma.theta.mle[1],sigma12=sigma.theta.mle[2],sigma2=sigma.theta.mle[4])
    },
    warning=function(w){rep(NA,5)},
    error=function(e){rep(NA,5)}
    )
}

##parameter 'sigma' removed 1/9/17
pop.fit <- function(dat,rho,weights=rep(1,length(dat))) {
    rbind(ar=ar.pop.mle(dat,rho,weights),mc=mc.pop.mle(dat,weights))
}




## AUC <- sapply(1:length(Ts),function(n) {
##     positive <- biomarkers[S[,i],i]
##     negative <- biomarkers[!S[,i],i]
##     mean(outer(negative,positive,`<`))
## })

roc <- function(x,...) {
    UseMethod('roc',x)
}

roc.PatientClass <- function(data,xi,rho,sigma,t) {
    biomarkers <- biomarker(data,theta=xi[,1:2],p=xi[,3:4],rho,sigma,t)
    statuses <- sapply(data,function(d)d$S[t])
    roc(biomarkers,statuses,interpolator=T)
}

roc.default <- function(marker,status,thresholds=marker,interpolator=F) {
    total.pos <- sum(status)
    total.neg <- sum(1-status)
    res <- sapply(1:length(thresholds),function(i) {
        idx <- which(marker>thresholds[i])
        c(tpr=sum(status[idx])/total.pos,fpr=sum(1-status[idx])/total.neg)
    })
    res <- res[,order(-marker)]
    if(!interpolator) {
        return(res)
    } else {

        fpr <- function(new.fpr) {
            i <- which(new.fpr>=res['fpr',])

            if(length(i)==0) {

                tpr <- mean(c(0,res['tpr',1]),na.rm=T)
            } else {
                i.lower <- max(i)
                if(res['fpr',i.lower]==new.fpr) {
                    tpr <- res['tpr',i.lower]
                } else {
                    if(i.lower==ncol(res)) {
                        tpr <- mean(c(res['tpr',i.lower],1),na.rm=T)
                    } else {
                        tpr <- mean(res['tpr',c(i.lower,min(i.lower+1,ncol(res)))])
                    }
                }
            }
            return(tpr)
        }

        fpr <- Vectorize(fpr)
        return(structure(fpr,class='PatientROC'))
    }
}



## roc.est <- function(dat,rho,sigma,pop.fits=pop.fit(dat,rho=rho,sigma),Ts,weights=rep(1,attr(dat,'N')),interpolator=F,auc=F) {
##     subj.fits <- subj.fit(dat,rho,sigma,mu.theta.prior=pop.fits['ar',c('theta0','theta1')],
##                           sigma.theta.prior=matrix(pop.fits[1,c(3,4,4,5)],nrow=2),
##                           mu.p.prior=pop.fits['mc',c('theta0','theta1')],
##                           sigma.p.prior=matrix(pop.fits[2,c(3,4,4,5)],nrow=2),Ts-1)
##     ## theta <- do.call(rbind,subj.fits)[,c('theta0','theta1')]
##     ## p <- do.call(rbind,subj.fits)[,c('p00','p11')]
##     ## theta <- do.call(cbind,subj.fits)[,c('theta0','theta1')]
##     ## p <- do.call(cbind,subj.fits)[,c('p00','p11')]


##     ## biomarkers <- as.matrix(biomarker(dat,theta,p,rho,sigma,Ts))
##     biomarkers <- sapply(1:length(Ts),function(i)biomarker(dat,subj.fits[[i]][,1:2],subj.fits[[i]][,3:4],rho,sigma,Ts[i])) %>% matrix(ncol=length(Ts))
##     S <- sapply(dat,function(dat.patient)as.logical(dat.patient$S[Ts])) %>% matrix(ncol=length(Ts))
##     ## roc(biomarkers,S)

##     if(auc) {#!allow for weights
##         aucs <- sapply(1:length(Ts),function(i){
##             positive <- biomarkers[S[,i],i]
##             negative <- biomarkers[!S[,i],i]
##             mean(outer(negative,positive,`<`))
##         })
##         return(aucs)
##     } else {
##         rocs <- lapply(1:length(Ts),function(i)roc=roc(biomarkers[,i],S[,i],weights=weights,interpolator=interpolator))
##         return(rocs)
##     }
##     ## rocs <- do.call(rbind,rocs)
## }


## roc.bootstrap <- function(dat,rho,sigma,bootstrap.reps,alpha) {
##     dat.bs <- list()
##     ans <- replicate(bootstrap.reps, {
##         dat.bs <<- dat[sample(N,N,replace=T)]
##         ## roc.est(dat.bs,rho,sigma,pop.fits,Ts)
##         roc.est(dat=dat.bs,rho=rho,sigma=sigma,Ts=Ts)
##     })

##     rocs <- data.frame(tpr=unlist(ans['tpr',]),fpr=unlist(ans['fpr',]),sample=rep(1:bootstrap.reps,each=N))
##     rocs <- arrange(rocs,sample,desc(fpr))

##     roc.observed <- roc.est(dat=dat,rho=rho,sigma=sigma,Ts=Ts)
##     ## subj.fits <- subj.fit(dat,rho,sigma,mu.theta.prior=pop.fits['ar',c('theta0','theta1')],
##     ##                       sigma.theta.prior=matrix(c(pop.fits[1,3],pop.fits[1,4],pop.fits[1,4],pop.fits[1,5]),nrow=2),
##     ##                       mu.p.prior=pop.fits['mc',c('theta0','theta1')],
##     ##                       sigma.p.prior=matrix(c(pop.fits[2,3],pop.fits[2,4],pop.fits[2,4],pop.fits[2,5]),nrow=2),Ts-1)
##     ## theta <- do.call(rbind,subj.fits)[,c('theta0','theta1')]
##     ## p <- do.call(rbind,subj.fits)[,c('p00','p11')]

##     ## biomarkers <- as.matrix(biomarker(dat,theta,p,rho,sigma,Ts))
##     ## S <- as.matrix(sapply(dat,function(dat.patient)as.logical(dat.patient$S[Ts])))
##     ## roc.observed <- roc(biomarkers,S)%>%t%>%as.data.frame


##     ## roc.observed$tpr <- qlogis(roc.observed$tpr)
##     ## rocs$tpr <- qlogis(rocs$tpr)
##     grid <- roc.observed$fpr
##     ans <- sapply(grid,function(x) {
##         sapply(seq_len(bootstrap.reps),function(sample.no) {
##             idx <- min(which(rocs[rocs$sample==sample.no,'fpr']<=x))
##             rocs[rocs$sample==sample.no,'tpr'][idx]
##         }
##         )
##     })
##     sds <- sqrt(apply(ans,2,var))

##     q <- qnorm(alpha)
##     gg.df <- cbind(roc.observed,upper=roc.observed$tpr+q*sds,lower=roc.observed$tpr-q*sds,sds=sds)
##     for(fpr in unique(gg.df$fpr)) {
##         gg.df$upper[gg.df$fpr==fpr] <- max(gg.df$upper[gg.df$fpr==fpr])
##         gg.df$lower[gg.df$fpr==fpr] <- min(gg.df$lower[gg.df$fpr==fpr])
##     }
##     gg.df$upper <- pmin(gg.df$upper,1); gg.df$lower <- pmax(gg.df$lower,0)
##     return(gg.df)
## }



## ## test mc.pop.mle
## rho <- .5
## sigma <- 1#sigmas[2]
## mu.theta.prior <- c(0,10)
## sigma.theta.prior <- matrix(c(5,3,3,5),nrow=2)
## mu.p.prior <- c(0,0)
## sigma.p.prior <- matrix(c(1,-.6,-.6,1),nrow=2)
## N <- 3e1
## T <- 5e1

## ans <- replicate(3e1,{
## theta <- rmvnorm(N,mean=mu.theta.prior,sigma=sigma.theta.prior)
## p <- rmvnorm(N,mean=mu.p.prior,sigma=sigma.p.prior) %>% plogis
## data <- rpatient(rho=rho,sigma=sigma,theta=theta,p=p,T=T)
## mc.pop.mle(data)
## })

## gg.df <- gather(data.frame(t(ans)))
## gg.df <- mutate(gg.df,param=0)
## gg.df$param[gg.df$key=='theta0'] <- mu.p.prior[1]
## gg.df$param[gg.df$key=='theta1'] <- mu.p.prior[2]
## gg.df$param[gg.df$key=='sigma1'] <- sigma.p.prior[1,1]
## gg.df$param[gg.df$key=='sigma12'] <- sigma.p.prior[1,2]
## gg.df$param[gg.df$key=='sigma2'] <- sigma.p.prior[2,2]
## ggplot(gg.df,aes(x=value))+geom_histogram()+geom_vline(aes(xintercept=param),color='red')+facet_wrap(~key)


## ## test ar.pop.mle
## rho <- .5
## sigma <- 5
## mu.theta.prior <- c(0,10)
## sigma.theta.prior <- matrix(c(5,3,3,5),nrow=2)
## mu.p.prior <- c(0,0)
## sigma.p.prior <- matrix(c(1,-.6,-.6,1),nrow=2)
## N <- 3e1
## T <- 5e1

## ans <- replicate(3e1,{
## theta <- rmvnorm(N,mean=mu.theta.prior,sigma=sigma.theta.prior)
## p <- rmvnorm(N,mean=mu.p.prior,sigma=sigma.p.prior) %>% plogis
## data <- rpatient(rho=rho,sigma=sigma,theta=theta,p=p,T=T)
## ar.pop.mle(data,rho)
## })

## gg.df <- gather(data.frame(t(ans)))
## gg.df <- mutate(gg.df,param=0)
## gg.df$param[gg.df$key=='theta0'] <- mu.theta.prior[1]
## gg.df$param[gg.df$key=='theta1'] <- mu.theta.prior[2]
## gg.df$param[gg.df$key=='sigma1'] <- sigma.theta.prior[1,1]
## gg.df$param[gg.df$key=='sigma12'] <- sigma.theta.prior[1,2]
## gg.df$param[gg.df$key=='sigma2'] <- sigma.theta.prior[2,2]
## ggplot(gg.df,aes(x=value))+geom_histogram()+geom_vline(aes(xintercept=param),color='red')+facet_wrap(~key)


## ## test qpatient.inf
## status <- 'positive'
## xi <- c(0,1,.5,.5)
## sigma <- 1
## rho <- .6
## p <- .3
## grid <- seq(0,1,length.out=100)
## fprs <- qpatient.inf(grid,xi,sigma,rho,status='negative',lower.tail=F)
## tprs <- ppatient.inf(fprs,xi,sigma,rho,status='positive',lower.tail=F)
## plot(tprs ~ grid)

## grid <- seq(-10,10,length.out=1e2)
## res <- ppatient.inf(grid,xi,sigma,rho,lower.tail=F)
## lines(res['fpr',],res['tpr',],col='red')


roc.hat <- function(xi,eta,rho,sigma,t,reps) {
    xi <- matrix(rep(xi,reps),byrow=TRUE,ncol=4)
    data <- rpatient(rho,sigma,xi[,1:2],xi[,3:4],t)
    xi <- subj.fit(data,rho,sigma,eta,t-1)
    biomarkers <- biomarker(data,theta=xi[,1:2],p=xi[,3:4],rho,sigma,t)
    statuses <- sapply(data,function(d)d$S[t])

    roc(biomarkers,statuses,interpolator=T)
}

roc.1 <- function(eta,rho,sigma,ts,B,roc.hat.reps) {
    mu.theta.prior <- eta['ar',1:2]; mu.p.prior <- eta['mc',1:2]
    sigma.theta.prior <- matrix(eta['ar',c(3,4,4,5)],nrow=2);
    sigma.p.prior <- matrix(eta['mc',c(3,4,4,5)],nrow=2)
    theta <- rmvnorm(B,mean=mu.theta.prior,sigma=sigma.theta.prior)
    p <- rmvnorm(B,mean=mu.p.prior,sigma=sigma.p.prior) %>% plogis
    xi <- cbind(theta,p)
    roc1s <- lapply(ts,function(t) {
        rocs <- lapply(1:nrow(xi),function(i)roc.hat(xi=xi[i,],eta,rho=rho,sigma=sigma,t=t,reps=roc.hat.reps))
        roc <- function(fpr)
            sapply(rocs,function(roc.i)roc.i(fpr)) %>% apply(1,mean,na.rm=T)
        return(structure(roc,class='PatientROC'))
    })
    if(length(roc1s)==1)return(roc1s[[1]]) else return(roc1s)
}


roc.3 <- function(x,...) {
    UseMethod('roc.3',x)
}

roc.3.PatientClass <- function(data,rho,sigma,ts,weights=rep(1,length(data))) {
    ## N <- attr(data,'N')
    eta <- pop.fit(data,rho,weights)
    rocs <- lapply(ts,function(t) {
        xi <- subj.fit(data,rho,sigma,eta,t-1)
        roc(data,xi,rho,sigma,t)
    })
    if(length(rocs)==1)return(rocs[[1]]) else return(rocs)
}

## parametric version
roc.3.default <- function(eta,rho,sigma,t,B) {
    mu.theta.prior <- eta['ar',1:2]; mu.p.prior <- eta['mc',1:2]
    sigma.theta.prior <- matrix(eta['ar',c(3,4,4,5)],nrow=2);
    sigma.p.prior <- matrix(eta['mc',c(3,4,4,5)],nrow=2)
    theta <- rmvnorm(B,mean=mu.theta.prior,sigma=sigma.theta.prior)
    p <- rmvnorm(B,mean=mu.p.prior,sigma=sigma.p.prior) %>% plogis
    data <- rpatient(rho,sigma,theta,p,t)
    xi <- subj.fit(data,rho,sigma,eta,t-1)
    roc(data,xi=xi,rho=rho,sigma=sigma,t=t)
}

roc.2 <- function(xi,sigma,rho) {
    roc0 <- function(t) {
        fprs <- qpatient.inf(t,xi=xi,sigma=sigma,rho=rho,status='negative',lower.tail=F)
        ppatient.inf(fprs,xi=xi,sigma=sigma,rho=rho,status='positive',lower.tail=F)
    }
    return(structure(roc0,class='PatientROC'))
}
