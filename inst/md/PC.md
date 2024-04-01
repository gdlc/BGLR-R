### Estimating the proportion of variance explained by PC

For further details see Janss et al. [Genetics, 2012](https://doi.org/10.1534/genetics.112.141143) and references therein.

```R
 
   library(BGLR)
   data(wheat)
   X=scale(wheat.X)/sqrt(ncol(X))
   y=wheat.Y[,1]
   G=tcrossprod(X)
   EVD=eigen(G)
   PC=EVD$vectors[,-599] 
   for(i in 1:ncol(PC)){ PC[,i]=PC[,i]*sqrt(EVD$values[i])}

   fm=BGLR(y=y,ETA=list( list(X=PC,model='BRR',saveEffects=T)),
	  nIter=6000,burnIn=1000)
   B=readBinMat('ETA_1_b.bin')

   VAR=matrix(nrow=nrow(B),ncol=ncol(B))
   for(i in 1:nrow(B)){
	for(j in 1:ncol(B)){ VAR[i,j]=var(PC[,j]*B[i,j]) }
        print(i)
   }

   mean(rowSums(VAR))+fm$varE
   varExplained=colMeans(VAR)
  plot(cumsum(varExplained)/sum(varExplained),type='o',col=4)
  lines(x=1:598,y=cumsum(EVD$values[1:598])/sum(EVD$values[1:598]),col=2,lty=2)

   plot(varExplained~EVD$values[1:598])
```

[Return to examples](https://github.com/gdlc/BGLR-R/blob/master/README.md)



