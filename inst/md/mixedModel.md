#### Using BGLR to fitmodels with multiple sets of effects

```R

  library(BGLR)
  data(mice)

  X=scale(mice.X) # markers
  A=mice.A        # pedigree-relationships
  pheno=mice.pheno

  fm=BGLR(y=pheno$Obesity.BMI,
          ETA=list(
			fixed=list(~factor(GENDER)+factor(Litter),data=pheno,model='FIXED'),
                    	cage=list(~factor(cage),data=pheno,model='BRR'),
                    	ped=list(K=A,model='RKHS'),
                    	mrk=list(X=X,model='BayesB')
		)
	)

    fm$ETA$fixed$b
    fm$ETA$cage$b
    fm$ETA$ped$u
    fm$ETA$mrk$b
```
[Back to examples](https://github.com/gdlc/BGLR-R/blob/master/README.md)
