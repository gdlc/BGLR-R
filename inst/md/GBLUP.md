


``` R

library(BGLR)
data(mice)

nIter=1200
burnIn=200
h2=.8
nQTL=20


X=scale(mice.X)
n=nrow(X)


QTL=seq(from=500,by=500,length=nQTL)
b=rnorm(n=nQTL,sd=sqrt(h2/nQTL))
signal=X[,QTL]%*%b
error=rnorm(sd=sqrt(1-h2),n)
y=signal+error

X=X/sqrt(ncol(X))


fm1=BGLR(	y=y,ETA=list(mrk=list(X=X,model='BRR')),
			nIter=nIter,burnIn=burnIn,saveAt='brr_'
		)


G=tcrossprod(X)

fm2=BGLR(	y=y,ETA=list(G=list(K=G,model='RKHS')),
			nIter=nIter,burnIn=burnIn,saveAt='eig_'
		)


EVD=eigen(G)

fm2b=BGLR(	y=y,ETA=list(G=list(V=EVD$vectors,d=EVD$values,model='RKHS')),
			nIter=nIter,burnIn=burnIn,
			burnIn=burnIn,saveAt='eigb_')

PC=EVD$vectors
for(i in 1:ncol(PC)){
	PC[,i]=PC[,i]*sqrt(EVD$values[i])
}

PC=PC[,EVD$values>1e-5]

fm3=BGLR(	y=y,ETA=list(pc=list(X=PC,model='BRR')),nIter=nIter,
			burnIn=burnIn,saveAt='pc_')

```
