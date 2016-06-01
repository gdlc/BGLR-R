
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
	   nIter=12999,burnIn=2000,saveAt='GxE_'
 	 )

```

**(2) Providing the G-matrix**



[Back to examples](https://github.com/gdlc/BGLR-R/blob/master/README.md)
