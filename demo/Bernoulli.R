#Finney, D. J. (1947). The estimation from Individual Records of the Relationship Between Dose and Quantal Response. 
#Biometrika, 34, 320-334
#Albert, J., Chib, S. (1993). Bayesian Analysis of Binary and Polychotomus Response Data.
#JASA, 88, 669-679.

rm(list=ls())
setwd(tempdir())

y=c(1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,1,1,1,0,1,0,0,0,0,1,0,1,0,1,0,1,0,0,1,1,1,0,0,1)
v=c(3.7,3.5,1.25,.75,.8,.7,.6,1.1,.9,.9,.8,.55,.6,1.4,.75,2.3,3.2,.85,1.7,1.8,
    0.4,.95,1.35,1.5,1.6,.6,1.8,.95,1.9,1.6,2.7,2.35,1.1,1.1,1.2,.8,.95,.75,1.3)
r=c(.825,1.09,2.5,1.5,3.2,3.5,.75,1.7,.75,.45,.57,2.75,3,2.33,3.75,1.64,1.6,1.415,
    1.06,1.8,2,1.36,1.35,1.36,1.78,1.5,1.5,1.9,.95,.4,.75,.03,1.83,2.2,2,3.33,1.9,1.9,1.625)
X=cbind(v,r)

nIter=5000;
burnIn=2500;
thin=10;
saveAt='';
ETA=list(list(X=X,model='FIXED'))

fit_Bernoulli_BGLR=BGLR(y=y,response_type='ordinal',ETA=ETA,nIter=nIter,burnIn=burnIn,
                   thin=thin,saveAt=saveAt)
fit_Bernoulli_BGLR$ETA[[1]]$b
fit_Bernoulli_BGLR$ETA[[1]]$SD.b

#In our parameterization, mu is set to 0 and it is estimated as a threshold
-fit_Bernoulli_BGLR$threshold
fit_Bernoulli_BGLR$SD.threshold 

fit_mle=glm(y~v+r,family=binomial(link="probit")) 
summary(fit_mle)

readline("Press <return> to continue with next example: ") 

#Example of prediction for missing values
rm(list=ls())
setwd(tempdir())

#data
data(wheat)

#libraries
library(pROC)
 
# extracts phenotypes
#continous

y=wheat.Y[,1]  

#binary                
yBin=ifelse(y>0,1,0)

# generates testing dataset
tst=sample(1:599,size=100,replace=FALSE)
yNA=yBin 
yNA[tst]=NA
  
nIter=5000;
burnIn=2500;
thin=10;
saveAt='';
ETA=list(list(X=wheat.X,model='BRR'))

fit_Bernoulli_BGLR=BGLR(y=yNA,response_type='ordinal',ETA=ETA,nIter=nIter,burnIn=burnIn,
                   thin=thin,saveAt=saveAt)

mean((yBin[tst]-pnorm(fit_Bernoulli_BGLR$yHat[tst]))^2) # mean-sq. error
auc(response=yBin[tst],predictor=fit_Bernoulli_BGLR$yHat[tst])


