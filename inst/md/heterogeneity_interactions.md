** Estimating the proportion of variance explained by PC**

For further details see de los Campos et al. (JABES, 2015) and references therein.

```R
  ## Clustering using 5-PCs and kmeans
   library(BGLR)
   data(wheat)
   X=scale(wheat.X)/sqrt(ncol(X))
   y=wheat.Y[,1]
   G=tcrossprod(X)/ncol(X)
   EVD=eigen(G)
   group=kmeans(x=EVD$vectors[,1:5],centers=2,nstart=100)$cluster
   plot(EVD$vectors[,1:2],col=group)
```

```R
  X0=X # for main effects
  X1=X; X1[group==2,]=0 #interactions
  X2=X; X2[group==1,]=0 #interactions
  
  fm=BGLR(y=y,ETA=list(
                  main=list(X=X0,model='BRR'),
                  int1=list(X=X1,model='BRR'),
                  int2=list(X=X2,model='BRR')
		),
	  nIter=6000,burnIn=1000,groups=group)

  varU1=fm$ETA[[1]]$varB+fm$ETA[[2]]$varB
  varU2=fm$ETA[[1]]$varB+fm$ETA[[3]]$varB
  varE1=fm$varE[1] 
  varE2=fm$varE[2]

  h2_1=varU1/(varU1+varE1)
  h2_2=varU2/(varU2+varE2)
  
  fm$ETA[[1]]$varB/sqrt(varU1*varU2) #correlation of effects

```

**NOTE**: Similar models can be fitted using G-matrices or PCs.

