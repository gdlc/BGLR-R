
### GxE using marker-by-environment interactions

The following examples illustrate how to implement marker-by-environments interaction models using BGLR, for further details about these models see [Lopez-Cruz et al. (2015)](https://doi.org/10.1534/g3.114.016097).


**(1) As a random regression on markers**

``` R
 library(BGLR)
 data(wheat)
 Y=wheat.Y # grain yield evaluated in 4 different environments
 round(cor(Y),2)

 X=scale(wheat.X)/sqrt(ncol(wheat.X))
 
 y2=Y[,2]
 y3=Y[,3]
 y=c(y2,y3)

 X0=matrix(nrow=nrow(X),ncol=ncol(X),0) # a matrix full of zeros

 X_main=rbind(X,X)
 X_1=rbind(X,X0)
 X_2=rbind(X0,X)

 fm=BGLR( y=y,ETA=list(             
                       main=list(X=X_main,model='BRR'),
                       int1=list(X=X_1,model='BRR'),
                       int2=list(X=X_2,model='BRR')
                      ),
	   nIter=6000,burnIn=1000,saveAt='GxE_',groups=rep(1:2,each=nrow(X))
 	 )
 varU_main=scan('GxE_ETA_main_varB.dat')[-c(1:200)]
 varU_int1=scan('GxE_ETA_int1_varB.dat')[-c(1:200)]
 varU_int2=scan('GxE_ETA_int2_varB.dat')[-c(1:200)]
 varE=read.table('GxE_varE.dat',header=FALSE)[-c(1:200),]
 varU1=varU_main+varU_int1
 varU2=varU_main+varU_int2
 h2_1=varU1/(varU1+varE[,1])
 h2_2=varU2/(varU2+varE[,2])
 COR=varU_main/sqrt(varU1*varU2)
 mean(h2_1)
 mean(h2_2)
 mean(COR)
```



**(2) Using genomic relationships**

A model equivalent to the one presented above can be implemented using G-matrices (or factorizations of it) with off-diagnoal blocks zeroed out for interactions, for further detials see [Lopez-Cruz et al., 2015](https://doi.org/10.1534/g3.114.016097).



```r
 library(BGLR)
 data(wheat)
 Y=wheat.Y # grain yield evaluated in 4 different environments
 round(cor(Y),2)
 y2=Y[,2]
 y3=Y[,3]
 y=c(y2,y3)

 library(BGData)
 
 G=getG(wheat.X,center=TRUE,scaleG=TRUE,scale=TRUE)
 EVD=eigen(G)
 PC=EVD$vectors[,EVD$values>1e-5]
 for(i in 1:ncol(PC)){ PC[,i]=EVD$vectors[,i]*sqrt(EVD$values[i]) }
 
 XMain=rbind(PC,PC)
 X0=matrix(nrow=nrow(X),ncol=ncol(PC),0) # a matrix full of zeros
 X1=rbind(PC,X0)
 X2=rbind(X0,PC)
 
 LP=list(main=list(X=XMain,model='BRR'), int1=list(X=X1,model='BRR'),int2=list(X=X2,model='BRR'))
 fmGRM=BGLR(y=y,ETA=LP,nIter=12000,burnIn=2000,saveAt='GRM_',groups=rep(1:2,each=nrow(X)))
 
 plot(fm$yHat,fmGRM$yHat)
 
 rbind(fm$varE,fmGRM$varE)
  
  
```


**(3) Incomplete, unbalanced designs**

In the example presented above, as in [Lopez-Cruz et al. (2015)](https://doi.org/10.1534/g3.114.016097), the data is assumed to be complete and balanced (i.e., all gentoypes on all environment), but this is not needed. The example below illustrates how to fit the interaction model with unbalanced data.


```R
library(BGLR)
data(wheat)
Y=wheat.Y
X=scale(wheat.X)/sqrt(ncol(wheat.X))
rownames(X)<-1:nrow(X)# use your IDs, don't need to be 1:n
rownames(Y)<-1:nrow(Y) # Y has data from 4 env (each col=1 env) 
n1=150
n2=305
y=c(sample(Y[,1],size=n1),sample(Y[,2],size=n2)) #unbalanced data from 2 env
env=c(rep(1,n1),rep(2,n2))

# Method 1: using SNPs explicitly

X0=X[names(y),] # Matrix for main effects
stopifnot(all(rownames(X0)==names(y)))

# now interactions
X1=X0
X2=X0
for(i in 1:nrow(X0)){
	X1[i,]<-(env[i]==1)*X0[i,]
	X2[i,]<-(env[i]==2)*X0[i,]	
}
fm=BGLR(y=y,ETA=list(list(X=X0,model='BRR'),
                     list(X=X1,model='BRR'),
                     list(X=X2,model='BRR'))
        ,groups=env)

 fm$varE
 fm$ETA[[1]]$varB
 fm$ETA[[2]]$varB
 fm$ETA[[3]]$varB
 
```

[Back to examples](https://github.com/gdlc/BGLR-R/blob/master/README.md)
