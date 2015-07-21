#A Cheese testing experiment

#The following data kindly provided by Dr. Graeme Newell were obtained from an 
#experiment concerning the effect on taste of various cheese additives.  
#The so-called hedonic scale has nine response categories, ranging from 
#'strong dislike' (1) to 'excellent taste' (9). In this instance, four additives 
#labeled A, B, C and D were tested. 

#Here the effects are so great that the qualitative ordering (D, A, C, B) can easily 
#be deduced from visual inspection. Nevertheless it is of some interest to check 
#whether the models described earlier are capable of describing these differences 
#and evaluating the statistical significance of the differences observed

#References
#Albert, J., Chib, S. (1993). Bayesian Analysis of Binary and Polychotomus 
#Response Data. JASA, 88, 669-679.

rm(list=ls())
setwd(tempdir())

#polr function
library(MASS)

#Function to expand a data.frame using column  Freq
#x is a data.frame, it should inclde a column named Freq
#the function returns another data.frame
expand.dft=function(x, na.strings = "NA", as.is = FALSE, dec = ".")
{
  DF=sapply(1:nrow(x), function(i) x[rep(i, each = x$Freq[i]), ],
               simplify = FALSE)
  DF=subset(do.call("rbind", DF), select = -Freq)
  for (i in 1:ncol(DF))
  {
    DF[[i]]=type.convert(as.character(DF[[i]]),
                         na.strings = na.strings,
                         as.is = as.is, dec = dec)                                           
  }
  DF
}  


#Data
Freq=c(0,0,1,7,8,8,19,8,1,6,9,12,11,7,6,1,0,0,1,1,6,8,23,7,5,1,0,0,0,0,1,3,7,14,16,11)
response=gl(9,1,36,labels=c("I","II","III","IV","V","VI","VII","VIII","IX"))
additive=gl(4,9,labels=c("A","B","C","D"))
cheese=data.frame(Freq,response,additive)
cheese


#a)Bayesian model

#Design matrix without intercept and expand it using frequencies
data=expand.dft(as.data.frame(cbind(response,model.matrix(~additive)[,-1],Freq)))

#The response should be ordered
data=data[order(data[,1]),]

#Response
y=as.vector(data[,1])

#Design matrix
X=as.matrix(data[,c(2:4)])

nIter=5000;
burnIn=2500;
thin=10;
saveAt='';
ETA=list(list(X=X,model='FIXED'))

fit_ordinal_BGLR=BGLR(y=y,response_type='ordinal',ETA=ETA,nIter=nIter,burnIn=burnIn,
                   thin=thin,saveAt=saveAt)

fit_ordinal_BGLR$ETA[[1]]$b
fit_ordinal_BGLR$threshold
fit_ordinal_BGLR$SD.threshold

#b)Frequentist

#polr function
fitted_mle=polr(response ~ additive, weights = Freq, data = cheese,method="probit")
fitted_mle
