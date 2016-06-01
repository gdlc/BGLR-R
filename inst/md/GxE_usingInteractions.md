
### GxE using marker-by-environment interactions

The following examples illustrate how to implement marker-by-environments interaction models using BGLR, for further details about these models see [Lopez-Cruz et al., 2015](http://www.g3journal.org/content/5/4/569.full?sid=81d404b6-7d0f-4ace-8556-936393eb829d).


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

**Note**: similar models can be fitted using G-matrices (or factorizations of it) with off-diagnoal blocks zeroed out for interactions, for further detials see [Lopez-Cruz et al., 2015](http://www.g3journal.org/content/5/4/569.full?sid=81d404b6-7d0f-4ace-8556-936393eb829d).


[Back to examples](https://github.com/gdlc/BGLR-R/blob/master/README.md)
