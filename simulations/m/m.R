.libPaths(c('/afs/ir/users/h/a/habnice/Desktop/R/lib/R/library',"/usr/lib/R/library"))
source('../../ragon.R')


## glm estimate of hiv data
require(grid)
require(gridExtra)
dat <- read.csv('../../final.csv')
dat$blip <- as.numeric(dat$VL>1e3)
dat <- subset(dat,select=c('ID','Blip_YN50','CD4','CD4P','Visits','blip','Sex','Transmission','trtm','AgeDiag','TimeDiag','AGE'))
dat <- subset(dat,select=c('ID','CD4','CD4P','Visits','blip','Blip_YN50'))
dat <- rename(dat,Visit=Visits)
dat$CD4[dat$CD4>6000] <- NA
dat$CD4P[dat$CD4P>=99] <- NA
dat$Visit <- unlist(sapply(rle(dat$ID)$lengths,function(x)seq_len(x)))
dat$blip <- factor(dat$blip)
dat$ID <- as.character(dat$ID)
dat <- droplevels(dat)

glmer.fit <- function(dat) {
glmer0 <- glmer(blip ~ CD4P|ID,family=binomial, data=dat)
fits <- fitted(glmer0)
df <- data.frame(marker=fits,status=as.numeric(dat$blip[as.numeric(names(fits))])-1,row=as.numeric(names(fits)))
## idx <- which(dat$Visit>=0)
## with(df[df$row%in%idx,],roc(marker,status,interpolator=T))
with(na.omit(df),roc(marker,status,interpolator=TRUE))
}

IDs <- unique(dat$ID)
bootstrap.reps <- 1e1
bootstrap.rocs <- replicate(bootstrap.reps, {
    dat.star <- dat[dat$ID %in% sample(IDs,replace=TRUE),]
    rownames(dat.star) <- 1:nrow(dat.star)
    tryCatch(glmer.fit(dat.star), error=function(e)NA)
},simplify=FALSE)

save.image(paste0('m-',round(runif(1)*1e8),'.RData'))

## from corn
source('../../ragon.R')
filelist <- dir()
filelist <- filelist[grep('m-[0-9]{8}.RData$',filelist)]
bootstrap.rocs.full <- list()
## k <- 1
for(file in filelist) {
    ## print(k)
    ## k <- k+1
    load(file)
    bootstrap.rocs.full <- c(bootstrap.rocs.full,bootstrap.rocs)
}
bootstrap.rocs <- bootstrap.rocs.full[!is.na(bootstrap.rocs.full)]
glmer.roc <- glmer.fit(dat)
## plot(glmer.roc)
plt.glmer <- plot(roc.observed=glmer.roc,bootstrap.rocs=bootstrap.rocs,alpha=.05,resolution=1e2)
plt.glmer
## save.image('m.RData')
load('m.RData')

