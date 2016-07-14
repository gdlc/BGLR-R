## Modeling Genetic Heterogenity using Random-effects Interactions

For further details see [de los Campos et al., (JABES, 2015)](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4666286) and references therein.

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

**Evaluation in trainint-testing**

Note: since this involves 50 training-testing partitions, running the example will take time.

```R

  nRep=50
  COR=matrix(nrow=nRep,ncol=3); colnames(COR)=c('across','within','interaction')
  group1=which(groups==1)
  group2=which(groups==2)
  
  for(i in 1:nRep){
    yNA=y
    tst1=sample(group1,size=86)
    tst2=sample(group2,size=64)
    yNA[c(tst1,tst2)]=NA
    
    # Across
    fm0=BGLR( y=yNA,ETA=list( 
 		  int=list(X=Z2, model='FIXED'),
 		  list(X=X0,model='BRR') 
 	      ),
	  nIter=6000,burnIn=1000,saveAt='m0_',verbose=FALSE)
    # Within
    fm1=BGLR( y=yNA[group==1],
 	   ETA=list( list(X=X0[group==1,],model='BRR') ),
	   nIter=6000,burnIn=1000, saveAt='m1_',verbose=FALSE)
	   
    fm2=BGLR( y=yNA[group==2],
 	   ETA=list( list(X=X0[group==2,],model='BRR') ),
	   nIter=6000,burnIn=1000, saveAt='m2_',verbose=FALSE)

    # Interaction	  
    fm12=BGLR(y=yNA,ETA=list(
  		  int=list(X=Z2, model='FIXED'),
                  main=list(X=X0,model='BRR'),
                  int1=list(X=X1,model='BRR'),
                  int2=list(X=X2,model='BRR')
		),
	     nIter=6000,burnIn=1000,groups=group, saveAt='m12_',verbose=FALSE)
     COR1[i,1]=cor(y[tst1],fm0$yHat[tst1])
     COR1[i,2]=cor(y[tst1],fm1$yHat[which(group1%in%tst1)])
     COR1[i,3]=cor(y[tst1],fm12$yHat[tst1])

     COR2[i,1]=cor(y[tst2],fm0$yHat[tst2])
     COR2[i,2]=cor(y[tst2],fm2$yHat[which(group2%in%tst2)])
     COR2[i,3]=cor(y[tst2],fm12$yHat[tst2])
     cat('Done with parititon ',i,'\n')   
    }

```
**NOTES**:

  - The above example uses a Gaussian prior, for information on other priors implemented in BGLR see [Bayesian Alphabet](https://github.com/gdlc/BGLR-R/blob/master/inst/md/BayesianAlphabet.md)
  - Within Gaussian contexts, the models described above can also be fitted using kinship matrices or marker-derived PCs, for this see [GBLUP](https://github.com/gdlc/BGLR-R/blob/master/inst/md/GBLUP.md)
  - Random effects interactions can also be used to deal with [GxE](https://github.com/gdlc/BGLR-R/blob/master/inst/md/GxE_usingInteractions.md)


[Back to examples](https://github.com/gdlc/BGLR-R/blob/master/README.md)
