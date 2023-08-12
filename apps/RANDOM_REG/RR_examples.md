```{r}

library(BGLR)
library(data.table)
Y=fread('/mnt/research/quantgen/projects/G2F/data/PHENO.csv',data.table=FALSE)
load('/mnt/research/quantgen/projects/G2F/data/G.RData')
W=as.matrix(read.csv('/mnt/research/quantgen/projects/G2F/data/ECOV.csv',row.names=1))

Y=Y[Y$region=='North',]
G=G[rownames(G)%in%Y$genotype,rownames(G)%in%Y$genotype]
W=W[rownames(W)%in%Y$year_loc,]

# Creating the incidence matrix for genotypes
 EVD=eigen(G)
 PC=EVD$vectors[,EVD$values>1e-5]
 PC=sweep(x=PC,MARGIN=2,STAT=sqrt(EVD$values[EVD$values>1e-5]),FUN='*')
 rownames(PC)=rownames(G)
 colnames(PC)=paste0('PC',1:ncol(PC))
 X=PC[match(Y$genotype,rownames(G)),]
 rm(EVD)


# Creating the incidence matrix for the random regression
 Z=as.matrix(scale(W[,'yield_pMatHar'])) #this is the crop-model predicted yield
 colnames(Z)='ecov_yield'
 Z.EXPANDED=Z[match(Y$year_loc,rownames(Z)),,drop=FALSE]

# Setting the linear predictor
 LP=RR.set(X,Z.EXPANDED,model='BRR')

# Fitting the model
 setwd(tempdir())
 y=scale(Y$yield,center=FALSE,scale=TRUE)#this is not needed

# lmer (no G matrix)
 ID=Y$genotype
 library(lme4)
 fm0=lmer(y~Z.EXPANDED+(Z.EXPANDED|ID))

 fm=BGLR(y=y,ETA=LP ,verbose=TRUE)
 # fm=BLRXy(y=y,ETA=LP ,verbose=TRUE) # will be much faster in this data set

 RR.variances(fm)

 # compare with summary(fm0)

 B=RR.coef(fm=fm,X=PC)
 B0=coef(fm0)$ID

 plot(B$rand[,1]+B$fixed[1],B0[,1],main='Intercept',xlab='Bayesian (with G-matrix)',ylab='lmer (no SNPs)')
 abline(a=0,b=1,col=4,lwd=2)

 plot(B$rand[,2]+B$fixed[2],B0[,2],main='Slope',xlab='Bayesian (with G-matrix)',ylab='lmer (no SNPs')
 abline(a=0,b=1,col=4,lwd=2)


 ## Let's predict performance of lines at the .1,.25,.5,.75,.9 percentile of the Env. Cov

 Z=as.matrix(quantile(scale(W[,'yield_pMatHar']),prob=c(.1,.25,.5,.75,.9)))
 
 YHat=RR.predict(fm,PC,Z)
 head(YHat)

```
