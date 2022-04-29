## Modeling Genetic Heterogenity using Random-effects Interactions

For further details see [de los Campos et al., (JABES, 2015)](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4666286), [Veturi et al. (Genetics, 2019)](https://www.genetics.org/content/211/4/1395.abstract), and [this example](https://github.com/gdlc/BGLR-R/blob/master/inst/md/GxE_usingInteractions.md). For a similar application involving GxE see [Lopez-Cruz et al. (G3, 2015)](https://www.g3journal.org/content/5/4/569.short), and for one involving multi-breed and GxE see [Yao et al. (J. D. Sci., 2017)](https://www.sciencedirect.com/science/article/pii/S0022030217300450).

**Clustering and displaying PCs**
```R
  ## Clustering using 5-PCs and kmeans
   library(BGLR)
   data(wheat)
   X=scale(wheat.X)/sqrt(ncol(wheat.X))
   y=wheat.Y[,1]
   PC=svd(X,nu=2,nv=0)$u
   group=kmeans(x=PC,centers=2,nstart=100)$cluster
   plot(PC,col=group*2,cex=.5)
```
**Preparing inputs**
```R
  X0=X # for main effects
  X1=X; X1[group==2,]=0 #interactions
  X2=X; X2[group==1,]=0 #interactions
  Z2=as.matrix(as.integer(group==2)) # a dummy variable for group 2
```

Next we will fit 3 models: (1) a model assuming that effects are homogeneous across clusters, (2) a statified analysis, and (3) a model with interactions. 

**(1) Homogeneous effects model**

```R
 fm0=BGLR( y=y,ETA=list( 
 		  int=list(X=Z2, model='FIXED'),
 		  list(X=X0,model='BRR') 
 	      ),
	  nIter=6000,burnIn=1000,saveAt='m0_')
```

**(2) Stratified Analysis**
```R
 fm1=BGLR( y=y[group==1],
 	   ETA=list( list(X=X0[group==1,],model='BRR') ),
	   nIter=6000,burnIn=1000, saveAt='m1_')
	   
 fm2=BGLR( y=y[group==2],
 	   ETA=list( list(X=X0[group==2,],model='BRR') ),
	   nIter=6000,burnIn=1000, saveAt='m2_')
```

**(3) Interaction model**
```R
  fm12=BGLR(y=y,ETA=list(
  		  int=list(X=Z2, model='FIXED'),
                  main=list(X=X0,model='BRR'),
                  int1=list(X=X1,model='BRR'),
                  int2=list(X=X2,model='BRR')
		),
	     nIter=6000,burnIn=1000,groups=group, saveAt='m12_')

  varU1=fm12$ETA[[2]]$varB+fm12$ETA[[3]]$varB
  varU2=fm12$ETA[[2]]$varB+fm12$ETA[[4]]$varB
  varE1=fm12$varE[1] 
  varE2=fm12$varE[2]

  h2_1=varU1/(varU1+varE1)
  h2_2=varU2/(varU2+varE2)
  
  fm12$ETA[[2]]$varB/sqrt(varU1*varU2) #correlation of effects
```

**NOTES**:

  - The above example uses a Gaussian prior, for information on other priors implemented in BGLR see [Bayesian Alphabet](https://github.com/gdlc/BGLR-R/blob/master/inst/md/BayesianAlphabet.md)
  - Within Gaussian contexts, the models described above can also be fitted using kinship matrices or marker-derived PCs, for this see [GBLUP](https://github.com/gdlc/BGLR-R/blob/master/inst/md/GBLUP.md)
  - Random effects interactions can also be used to deal with [GxE](https://github.com/gdlc/BGLR-R/blob/master/inst/md/GxE_usingInteractions.md)


[Back to examples](https://github.com/gdlc/BGLR-R/blob/master/README.md)
