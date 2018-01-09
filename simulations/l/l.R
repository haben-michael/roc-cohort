.libPaths(c('/afs/ir/users/h/a/habnice/Desktop/R/lib/R/library',"/usr/lib/R/library"))
source('../../ragon.R')

## 2. HIV data

## set.seed(1)
hiv.data <- get.hiv.data('../../final.csv')
hiv.params <- get.rho.sigma(hiv.data)

T <- attr(hiv.data,'T')
N <- attr(hiv.data,'N')
B <- 1e2


rho.sigma.bs <- replicate(B, {
    cat('.')
    hiv.data.star <- hiv.data[sample(1:N,replace=TRUE)]
    hiv.params.star <- get.rho.sigma(hiv.data.star)
})
rho.sigma.obs <- get.rho.sigma(hiv.data)
gg.df <- data.frame(value=as.numeric(t(rho.sigma.bs)),param=rep(c('rho','intercept','sigma'),each=B),observed=rep(unlist(rho.sigma.obs),each=B))
ggplot(gg.df,aes(x=value,y=1))+geom_point(alpha=.2)+facet_wrap(~param,scales='free_x')+geom_vline(aes(xintercept=observed))
bootstrap.reps <- 3e0


t <- 5
dfs <- list()
eta.bs <- replicate(B, {
    cat('.')
    hiv.data.star <- hiv.data[sample(1:N,replace=TRUE)]
    dfs <<- c(dfs,list(hiv.data.star))
    ## hiv.params.star <- get.eta(hiv.data.star)
    eta <- pop.fit(hiv.data.star,rho=rho.sigma.obs$rho,t=t)
})
eta.obs <- pop.fit(hiv.data,rho=rho.sigma.obs$rho,t=t)

labels <- apply(expand.grid(rownames(eta.obs),colnames(eta.obs)),1,paste,collapse='.')
gg.df <- data.frame(value=as.numeric(eta.bs),param=rep(labels,each=B),observed=rep(as.numeric(eta.obs),each=B))
ggplot(gg.df,aes(x=value,y=1))+geom_point(alpha=.2)+facet_wrap(~param,scales='free_x')+geom_vline(aes(xintercept=observed))

hist(eta.bs['mc','sigma1',])
idx <- which.max(eta.bs['mc','sigma1',])
pop.fit(dfs[[idx]],rho=rho.sigma.obs$rho,t=t)
