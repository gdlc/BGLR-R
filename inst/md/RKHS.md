#### Semi-parametric regression

Examples (modified) from [de los Campos et al., Genetics Research, 2010](https://doi.org/10.1017/s0016672310000285).

**(1) RKHS with a Gaussian Kernel and fixed bandwidth parameter**

```R
 library(BGLR)
 data(wheat); X=wheat.X; Y=wheat.Y

### DISTANCE MATRIX #############################
  D<-as.matrix(dist(X,method="euclidean"))^2
  D<-D/mean(D)
  h<-1

### KERNEL ######################################
  K<-exp(-h*D)
  
### GENERATES TESTING SET #######################
 set.seed(12345)
 tst<-sample(1:599,size=100,replace=FALSE)
 y<-Y[,4]
 yNA<-y
 yNA[tst]<-NA
  
### MODEL FITTING #################################

 fm=BGLR(y=yNA,ETA=list(list(K=K,model='RKHS')),nIter=6000,burnIn=1000)
 fm$varE
 fm$ETA[[1]]$varU
 cor(y[tst],fm$yHat[tst])

```

**(2) Fitting the model over a grid of values of the bandwidth parameter**

```R
 h=c(.01,.1,.4,.8,1.5,3,5)
 PMSE<-numeric(); VARE<-numeric(); VARU<-numeric();
 pD<-numeric(); DIC<-numeric()
 fmList<-list()

 for(i in 1:length(h)){
 	print(paste('Working with h=',h[i],sep=''))
    # COMPUTES THE KERNEL
    K<-exp(-h[i]*D)
    # FITS THE MODEL
    ETA<-list(list(K=K,model='RKHS'))
    prefix<- paste(h[i], "_",sep="")
    fm<-BGLR(y=yNA,ETA=ETA,
           nIter=5000,burnIn=1000,df0=5,S0=2,saveAt=prefix)
    fmList[[i]]<-fm
    PMSE[i]<-mean((y[tst]-fm$yHat[tst])^2)
    VARE[i]<-fm$varE
    VARU[i]<-fm$ETA[[1]]$varU
    DIC[i]<-fm$fit$DIC
    pD[i]<-fm$fit$pD
 }
  
 R2<-1-PMSE/mean((y[tst]-mean(y[-tst]))^2)

 ### PLOTS ############################### 
 plot(VARE~h,xlab="Bandwidth", 
      ylab="Residual Variance",type="o",col=4)
  
 plot(I(VARE/VARU)~h,xlab="Bandwidth",
      ylab="variance ratio (noise/signal)",type="o",col=4)

 plot(pD~h,xlab="Bandwidth", ylab="pD",type="o",col=2)

 plot(DIC~h,xlab="Bandwidth", ylab="DIC",type="o",col=2)

 plot(R2~h,xlab="Bandwidth", ylab="R-squared",type="o",col=2)

```


**(3) Kernel Averaging**

```R
 PMSE<-numeric()
 VARE<-numeric()
 KList<-list()
 for(i in 1:length(h)){
    KList[[i]]<-list(K=exp(-h[i]*D),model='RKHS')
 }

 fmKA<-BGLR(y=yNA,ETA=KList,thin=10,
            nIter=103000,burnIn=3000 ,saveAt="KA_")

 VARG<-numeric()
 for(i in 1:length(KList)){  VARG[i]<-fmKA$ETA[[i]]$varU }
 weights<-round(VARG/sum(VARG),5)

 PMSE<-mean((y[tst]-fmKA$yHat[tst])^2)
 R2_KA<-1-PMSE/mean((y[tst]-mean(y[-tst]))^2)
```

[Return to examples](https://github.com/gdlc/BGLR-R/blob/master/README.md)
