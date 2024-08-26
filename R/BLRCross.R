#Bayes A, Mewissen et al. (2001).
#Prediction of Total Genetic Value Using Genome-Wide Dense Marker Maps
#Genetics 157: 1819-1829, Modified so that the Scale parameter is estimated from data (a gamma prior is assigned)

setLT.BayesA.Cross=function(prior,vy,j,p,idColumns,sumVarX,R2,nLT,verbose,
                      saveAt,rmExistingFiles,thin,nIter,burnIn)
{	
	#Just a copy of values provided by user
	LT=list()
	LT$Name=prior$Name
	LT$R2=prior$R2
	LT$df0=prior$df0
	LT$S0=prior$S0
	LT$shape0=prior$shape0
	LT$rate0=prior$rate0
	LT$p=p
	LT$idColumns=idColumns
	LT$saveEffects=prior$saveEffects
	LT$storageMode=prior$storageMode
		
	LT$MSx=sumVarX
	
	if(is.null(LT$R2))
    {
    	LT$R2=R2/nLT
    	if(verbose)
    	{
      		message("R2 in LP ",j, " was missing and was set to ",LT$R2)
    	}
  	}
  	
  	#Default value for the degrees of freedom associated with the distribution 
  	#assigned to the variance of betas
  	if(is.null(LT$df0))
  	{
    	LT$df0= 5
    	if(verbose)
    	{
    		message("DF in LP ",j," was missing and was set to ",LT$df0)
    	}
  	}
  		
  	#Default value for the scale parameter associated with the distribution 
  	#assigned to the variance of betas
  	if(is.null(LT$S0))
  	{
  		if(LT$df0<=0) stop("df0>0 in ",LT$model," in order to set S0");
     	LT$S0=vy*LT$R2/(LT$MSx)*(LT$df0+2)
     	if(verbose)
     	{
     		message("Scale parameter in LP ",j," was missing and was set to ",LT$S0)
     	}
  	}
  	
  	# Improvement: Treat Scale as random, assign a gamma density 
  	message("Scale parameter is treated as a random, we have assigned a gamma density")
  	if(is.null(LT$shape0))
  	{
  		LT$shape0=1.1
  		message("shape parameter for the Scale in LP ",j, " was missing and was set to ",LT$shape0)
  	}

  	if(is.null(LT$rate0))
  	{
  		LT$rate0=(LT$shape0-1)/LT$S0		
  		message("rate parameter for the Scale in LP ",j, " was missing and was set to ",LT$rate0)
  	}
  	
  	LT$S=LT$S0
  	LT$b=rep(0,LT$p)
  	LT$varB=rep(LT$S0/(LT$df0+2),LT$p)
  	
  	# Add one file when S0 is treated as random.
  	fname=paste(saveAt,LT$Name,"_ScaleBayesA.dat",sep="") 
  	
  	if(rmExistingFiles)
  	{ 
    	unlink(fname) 
  	}
  	
  	LT$fileOut=file(description=fname,open="w")
  	LT$NamefileOut=fname;
  	
  	#Objects for saving posterior means for MCMC
    LT$post_b=rep(0,LT$p)
    LT$post_b2=rep(0,LT$p)
    LT$post_varB=0                 
    LT$post_varB2=0
    LT$post_S=0
    LT$post_S2=0
    
    #*#
    if(is.null(LT$saveEffects)){LT$saveEffects=FALSE}
    if(LT$saveEffects){
        if(is.null(LT$storageMode)){LT$storageMode="double"}
        if(!LT$storageMode%in%c("single","double")) {
            stop("storageMode of LP ",j," can either be 'single' or 'double' (default)")
        }
    	if(is.null(LT$thin)){ LT$thin=thin }
    	fname=paste(saveAt,LT$Name,"_b.bin",sep="")
    	if(rmExistingFiles){ unlink(fname) }
    	LT$fileEffects=file(fname,open='wb')
    	nRow=floor((nIter-burnIn)/LT$thin)
    	writeBin(object=c(nRow,LT$p),con=LT$fileEffects,size=ifelse(LT$storageMode=="single",4,8))
    }#*#
    
  	return(LT)
}

##########################################################################################
#Set linear term for BayesB
##########################################################################################

setLT.BayesB.Cross=function(prior,vy,j,p,idColumns,sumVarX,R2,nLT,verbose,
					  saveAt,rmExistingFiles,thin,nIter,burnIn)
{	
	#Just a copy of values provided by user
	LT=list()
	LT$Name=prior$Name
	LT$R2=prior$R2
	LT$df0=prior$df0
	LT$rate0=prior$rate0
	LT$shape0=prior$shape0
	LT$probIn=prior$probIn
	LT$counts=prior$counts
	LT$p=p
	LT$idColumns=idColumns
	LT$saveEffects=prior$saveEffects
	LT$storageMode=prior$storageMode
		
	LT$MSx=sumVarX
	
	if(is.null(LT$R2))
    {
    	LT$R2=R2/nLT
    	if(verbose)
    	{
      		message("R2 in LP ",j, " was missing and was set to ",LT$R2)
    	}
  	}
  	
  	#Default value for the degrees of freedom associated with the distribution 
  	#assigned to the variance of betas
  	if(is.null(LT$df0))
  	{
    	LT$df0= 5
    	if(verbose)
    	{
    		message("DF in LP ",j," was missing and was set to ",LT$df0)
    	}
  	}
  	
  	#Default value for a predictor being "in" the model
  	if(is.null(LT$probIn))
  	{
    	LT$probIn=0.5
    	if(verbose)
    	{	
       		message("probIn in LP ",j," was missing and was set to ",LT$probIn)
    	}
  	}
  	
  	#Default value for prior counts
  	if(is.null(LT$counts))
  	{
    	LT$counts=10
    	if(verbose)
    	{
       		message("Counts in LP ",j," was missing and was set to ",LT$counts)
    	}
  	}
  	
  	LT$countsIn=LT$counts * LT$probIn
  	LT$countsOut=LT$counts - LT$countsIn
  	
  	#Default value for the scale parameter associated with the distribution assigned to the variance of 
    #marker effects
    if(is.null(LT$S0))
    {
    	if(LT$df0<=0) stop("df0>0 in ",LT$model," in order to set S0")
    	LT$S0=vy*LT$R2/(LT$MSx)*(LT$df0+2)/LT$probIn
    	if(verbose)
    	{
    		message("Scale parameter in LP ",j," was missing and was set to ",LT$S0)
    	}
    }
    
    if(is.null(LT$shape0))
    {
    	LT$shape0=1.1
    	message("shape0 in LP ",j," was missing and was set to ",LT$shape0)
    }
    
    if(is.null(LT$rate0))
    {
    	LT$rate0=(LT$shape0-1)/LT$S0
    	message("rate0 in LP ",j," was missing and was set to ",LT$rate0)
    	
    }
    
    LT$S=LT$S0
  	
  	LT$a=rep(0, LT$p)	
  	LT$d=rbinom(n = LT$p, size = 1, prob = LT$probIn)
  	LT$b=LT$a*LT$d  #b=a*d, for compatibility with BGLR we use b instead of beta in linear terms
  	
  	LT$varB = rep(LT$S0/(LT$df0+2),LT$p)
  	
  	fname=paste(saveAt,LT$Name,"_parBayesB.dat",sep="")
  	
  	if(rmExistingFiles)
  	{
		unlink(fname) 
  	}
  	
  	LT$fileOut=file(description=fname,open="w")
  	LT$NamefileOut=fname
  	
  	tmp=c('probIn','scale')
   	write(tmp, ncolumns = 2, file = LT$fileOut, append = TRUE)
   	
  	
  	#Objects for saving posterior means for MCMC
    LT$post_b=rep(0,LT$p)
    LT$post_b2=rep(0,LT$p)
    LT$post_varB=0                 
    LT$post_varB2=0
    
    LT$post_d=0
    LT$post_probIn=0
    LT$post_probIn2=0
    
    LT$post_S=0
    LT$post_S2=0
     
    #*#
    if(is.null(LT$saveEffects)){LT$saveEffects=FALSE}
    if(LT$saveEffects){
        if(is.null(LT$storageMode)){LT$storageMode="double"}
        if(!LT$storageMode%in%c("single","double")) {
            stop("storageMode of LP ",j," can either be 'single' or 'double' (default)")
        }
    	if(is.null(LT$thin)){ LT$thin=thin }
    	fname=paste(saveAt,LT$Name,"_b.bin",sep="")
    	if(rmExistingFiles){ unlink(fname) }
    	LT$fileEffects=file(fname,open='wb')
    	nRow=floor((nIter-burnIn)/LT$thin)
    	writeBin(object=c(nRow,LT$p),con=LT$fileEffects,size=ifelse(LT$storageMode=="single",4,8))
    }#*#
      	
  	return(LT)
}

##########################################################################################
#Set linear term for BayesC
##########################################################################################

setLT.BayesC.Cross=function(prior,vy,j,p,idColumns,sumVarX,R2,nLT,verbose,
                      saveAt,rmExistingFiles,thin,nIter,burnIn)
{	
	#Just a copy of values provided by user
	LT=list()
	LT$Name=prior$Name
	LT$R2=prior$R2
	LT$df0=prior$df0
	LT$S0=prior$S0
	LT$probIn=prior$probIn
	LT$counts=prior$counts
	LT$p=p
	LT$idColumns=idColumns
	LT$saveEffects=prior$saveEffects
	LT$storageMode=prior$storageMode
		
	LT$MSx=sumVarX
		
	if(is.null(LT$R2))
    {
    	LT$R2=R2/nLT
    	if(verbose)
    	{
      		message("R2 in LP ",j, " was missing and was set to ",LT$R2)
    	}
  	}
  	
  	#Default value for the degrees of freedom associated with the distribution 
  	#assigned to the variance of betas
  	if(is.null(LT$df0))
  	{
    	LT$df0= 5
    	if(verbose)
    	{
    		message("DF in LP ",j," was missing and was set to ",LT$df0)
    	}
  	}
  	
  	#Default value for a predictor being "in" the model
  	if(is.null(LT$probIn))
  	{
    	LT$probIn=0.5
    	if(verbose)
    	{	
       		message("probIn in LP ",j," was missing and was set to ",LT$probIn)
    	}
  	}
  	
  	#Default value for prior counts
  	if(is.null(LT$counts))
  	{
    	LT$counts=10
    	if(verbose)
    	{
       		message("Counts in LP ",j," was missing and was set to ",LT$counts)
    	}
  	}
  	
  	LT$countsIn=LT$counts * LT$probIn
  	LT$countsOut=LT$counts - LT$countsIn
  	
  	#Default value for the scale parameter associated with the distribution 
  	#assigned to the variance of betas
  	if(is.null(LT$S0))
  	{
  		if(LT$df0<=0) stop("df0>0 in ",LT$model," in order to set S0");
     	LT$S0=vy*LT$R2/(LT$MSx)*(LT$df0+2)/LT$probIn
     	if(verbose)
     	{
     		message("Scale parameter in LP ",j," was missing and was set to ",LT$S0)
     	}
  	}
  	
  	LT$a=rep(0, LT$p)	
  	LT$d=rbinom(n = LT$p, size = 1, prob = LT$probIn)
  	LT$b=LT$a*LT$d  #b=a*d, for compatibility with BGLR we use b instead of beta in linear terms
  	
  	LT$varB = LT$S0
  	
  	fname=paste(saveAt,LT$Name,"_parBayesC.dat",sep="")
  	
  	if(rmExistingFiles)
  	{
		unlink(fname) 
  	}
  	
  	LT$fileOut=file(description=fname,open="w")
  	LT$NamefileOut=fname
  	
  	tmp=c('probIn','varB')
   	write(tmp, ncolumns = 2, file = LT$fileOut, append = TRUE)
  	
  	#Objects for saving posterior means for MCMC
    LT$post_b=rep(0,LT$p)
    LT$post_b2=rep(0,LT$p)
    LT$post_varB=0                 
    LT$post_varB2=0
    
    LT$post_d=0
    LT$post_probIn=0
    LT$post_probIn2=0
    
    #*#
    if(is.null(LT$saveEffects)){LT$saveEffects=FALSE}
    if(LT$saveEffects){
        if(is.null(LT$storageMode)){LT$storageMode="double"}
        if(!LT$storageMode%in%c("single","double")) {
            stop("storageMode of LP ",j," can either be 'single' or 'double' (default)")
        }
    	if(is.null(LT$thin)){ LT$thin=thin }
    	fname=paste(saveAt,LT$Name,"_b.bin",sep="")
    	if(rmExistingFiles){ unlink(fname) }
    	LT$fileEffects=file(fname,open='wb')
    	nRow=floor((nIter-burnIn)/LT$thin)
    	writeBin(object=c(nRow,LT$p),con=LT$fileEffects,size=ifelse(LT$storageMode=="single",4,8))
    }#*#
      	
  	return(LT)
}

##########################################################################################
#Set linear term for SSVS
#George, E. I. and McCulloch, R. E. 1993. Variable selection via Gibbs Sampling, 
#Journal of the American Statistical Association, 88(423): 881-889

##########################################################################################

#Evaluates the logarithm of p(c|else)
LogCondc=function(c,varB,b,d,shape1,shape2)
{
	p=length(b)
	bs=b[d!=1]
	-1/(2*varB*c^2)*sum(bs^2)+(shape1-1)*log(c)+(shape2-1)*log(1-c)-(p-sum(d))*log(c)
}

#Metropolis sampler for c|else
metropc=function(c,varB,b,d,shape1,shape2)
{

	flag=TRUE
	while(flag)
	{
    	c_new=c+rnorm(1,0,sd=0.05)
    	if(c_new>0 & c_new<1)
    	{
    		flag=FALSE
    	}
    }
    
    logP_old = LogCondc(c,varB,b,d,shape1,shape2)
    logP_new = LogCondc(c_new,varB,b,d,shape1,shape2) 
    
    accept = (logP_new - logP_old) > log(runif(1))
    if (accept) {
        return(c_new)
    }else{
    	return(c)
    }
}



setLT.SSVS.Cross=function(prior,vy,j,p,idColumns,sumVarX,R2,nLT,verbose,
                    saveAt,rmExistingFiles,thin,nIter,burnIn)
{	
	#Just a copy of values provided by user
	LT=list()
	
	LT$Name=prior$Name
	
	LT$cprobIn=prior$cprobIn
	LT$ccounts=prior$ccounts
	
	LT$R2=prior$R2
	LT$df0=prior$df0
	LT$S0=prior$S0
	LT$probIn=prior$probIn
	LT$counts=prior$counts
	LT$p=p
	LT$idColumns=idColumns
	LT$saveEffects=prior$saveEffects
	LT$storageMode=prior$storageMode
		
	LT$MSx=sumVarX
	
	
	if(is.null(LT$cprobIn))
	{
		LT$cprobIn=1/100
		if(verbose)
		{	
			message("cprobIn in LP ",j," was missing and was set to ",LT$cprobIn)
		}
	}
		
	if(is.null(LT$ccounts))
	{
		LT$ccounts=1E3
		if(verbose)
		{
    		message("ccounts in LP ",j," was missing and was set to ",LT$ccounts)
		}
	}
		
	LT$ccountsIn=LT$cprobIn*LT$ccounts        #Shape1
	LT$ccountsOut=(1-LT$cprobIn)*LT$ccounts   #Shape2
	LT$c=rbeta(n=1,LT$ccountsIn,LT$ccountsOut)
	
	
	if(is.null(LT$R2))
    {
    	LT$R2=R2/nLT
    	if(verbose)
    	{
      		message("R2 in LP ",j, " was missing and was set to ",LT$R2)
    	}
  	}
  	
  	#Default value for the degrees of freedom associated with the distribution 
  	#assigned to the variance of betas
  	if(is.null(LT$df0))
  	{
    	LT$df0= 5
    	if(verbose)
    	{
    		message("DF in LP ",j," was missing and was set to ",LT$df0)
    	}
  	}
  	
  	#Default value for a predictor being "in" the model
  	if(is.null(LT$probIn))
  	{
    	LT$probIn=0.5
    	if(verbose)
    	{	
       		message("probIn in LP ",j," was missing and was set to ",LT$probIn)
    	}
  	}
  	
  	#Default value for prior counts
  	if(is.null(LT$counts))
  	{
    	LT$counts=10
    	if(verbose)
    	{
       		message("Counts in LP ",j," was missing and was set to ",LT$counts)
    	}
  	}
  	
  	LT$countsIn=LT$counts * LT$probIn
  	LT$countsOut=LT$counts - LT$countsIn
  	
  	#Default value for the scale parameter associated with the distribution 
  	#assigned to the variance of betas
  	if(is.null(LT$S0))
  	{
  		if(LT$df0<=0) stop("df0>0 in ",LT$model," in order to set S0");
     	LT$S0=vy*LT$R2/(LT$MSx)*(LT$df0+2)/LT$probIn
     	if(verbose)
     	{
     		message("Scale parameter in LP ",j," was missing and was set to ",LT$S0)
     	}
  	}
  	
  	LT$b=rep(0, LT$p)	
  	LT$d=rbinom(n = LT$p, size = 1, prob = LT$probIn)
  	LT$a=rep(1, LT$p)  #a=1 if d==1 and a=c if d==0, here the values are all set to 1 initially
  	LT$varB = LT$S0
  	
  	
  	fname=paste(saveAt,LT$Name,"_parSSVS.dat",sep="")
  	
  	if(rmExistingFiles)
  	{
		unlink(fname) 
  	}
  	
  	LT$fileOut=file(description=fname,open="w")
  	LT$NamefileOut=fname
  	
  	tmp=c('probIn','varB','c')
   	write(tmp, ncolumns = 3, file = LT$fileOut, append = TRUE)
   	
  	#Objects for saving posterior means for MCMC
    LT$post_b=rep(0,LT$p)
    LT$post_b2=rep(0,LT$p)
    LT$post_varB=0                 
    LT$post_varB2=0
    
    LT$post_d=0
    LT$post_probIn=0
    LT$post_probIn2=0
    
    LT$post_c=0
    LT$post_c2=0
    
    #*#
    if(is.null(LT$saveEffects)){LT$saveEffects=FALSE}
    if(LT$saveEffects){
        if(is.null(LT$storageMode)){LT$storageMode="double"}
        if(!LT$storageMode%in%c("single","double")) {
            stop("storageMode of LP ",j," can either be 'single' or 'double' (default)")
        }
    	if(is.null(LT$thin)){ LT$thin=thin }
    	fname=paste(saveAt,LT$Name,"_b.bin",sep="")
    	if(rmExistingFiles){ unlink(fname) }
    	LT$fileEffects=file(fname,open='wb')
    	nRow=floor((nIter-burnIn)/LT$thin)
    	writeBin(object=c(nRow,LT$p),con=LT$fileEffects,size=ifelse(LT$storageMode=="single",4,8))
    }#*#
      	
  	return(LT)
}

##########################################################################################
#Set linear term for Bayesian Ridge Regression
##########################################################################################

setLT.BRR.Cross=function(prior,vy,j,p,idColumns,sumVarX,R2,nLT,verbose,
                   saveAt,rmExistingFiles,thin,nIter,burnIn)
{	
	#Just a copy of values provided by user
	LT=list()
	LT$Name=prior$Name
	LT$R2=prior$R2
	LT$df0=prior$df0
	LT$S0=prior$S0
	LT$p=p
	LT$idColumns=idColumns
	LT$saveEffects=prior$saveEffects
	LT$storageMode=prior$storageMode
		
	LT$MSx=sumVarX
	
	if(is.null(LT$R2))
    {
    	LT$R2=R2/nLT
    	if(verbose)
    	{
      		message("R2 in LP ",j, " was missing and was set to ",LT$R2)
    	}
  	}
  	
  	#Default value for the degrees of freedom associated with the distribution 
  	#assigned to the variance of betas
  	if(is.null(LT$df0))
  	{
    	LT$df0=5
    	if(verbose)
    	{
    		message("DF in LP ",j," was missing and was set to ",LT$df0)
    	}
  	}
  	
  	#Default value for the scale parameter associated with the distribution 
  	#assigned to the variance of betas
  	if(is.null(LT$S0))
  	{
  		if(LT$df0<=0) stop("df0>0 in ",LT$model," in order to set S0");
     	LT$S0=vy*LT$R2/(LT$MSx)*(LT$df0+2)
     	if(verbose)
     	{
     		message("Scale parameter in LP ",j," was missing and was set to ",LT$S0)
     	}
  	}
  	
  	LT$b=rep(0, LT$p)
  	LT$varB = LT$S0/(LT$df0+2)
  	
  	fname=paste(saveAt,LT$Name,"_varB.dat",sep=""); 
    
    if(rmExistingFiles)
    { 
       unlink(fname) 
    }
    
    LT$NamefileOut=fname
    LT$fileOut=file(description=fname,open="w")
    
  	#Objects for saving posterior means for MCMC
    LT$post_b=rep(0,LT$p)
    LT$post_b2=rep(0,LT$p)
    LT$post_varB=0                 
    LT$post_varB2=0
    
    #*#
    if(is.null(LT$saveEffects)){LT$saveEffects=FALSE}
    if(LT$saveEffects){
        if(is.null(LT$storageMode)){LT$storageMode="double"}
        if(!LT$storageMode%in%c("single","double")) {
            stop("storageMode of LP ",j," can either be 'single' or 'double' (default)")
        }
    	if(is.null(LT$thin)){ LT$thin=thin }
    	fname=paste(saveAt,LT$Name,"_b.bin",sep="")
    	if(rmExistingFiles){ unlink(fname) }
    	LT$fileEffects=file(fname,open='wb')
    	nRow=floor((nIter-burnIn)/LT$thin)
    	writeBin(object=c(nRow,LT$p),con=LT$fileEffects,size=ifelse(LT$storageMode=="single",4,8))
    }#*#
    
  	return(LT)
}

##########################################################################################
#FIXED Effects (Ridge Regression with huge variance  for regression coefficients, 
#which efective leads to a flat prior)
##########################################################################################

setLT.Fixed.Cross=function(p,idColumns,Name,saveAt,rmExistingFiles)
{	
	#Just a copy of values provided by user
	LT=list()
	
	LT$Name=Name
	LT$p=p
	LT$idColumns=idColumns
  	
  	LT$b=rep(0, LT$p)
  	LT$varB = 1e10
  	
  	
  	if(LT$Name=="___INTERCEPT___")
  	{
  		#Dirty hack to save samples from the intercept as "mu.dat"
        #when the intercept is as fixed effect included automatically in routine BLRXy
  		fname=paste(saveAt,"mu.dat",sep="")
  	}else{
  		fname=paste(saveAt,LT$Name,"_b.dat",sep="")
  	}

    LT$NamefileOut=fname

    if(rmExistingFiles)
    { 
       unlink(fname) 
    }

    LT$fileOut=file(description=fname,open="w")
  	
  	#Objects for saving posterior means for MCMC
    LT$post_b=rep(0,LT$p)
    LT$post_b2=rep(0,LT$p)
        
  	return(LT)
}

#RKHS
setLT.RKHS.Cross=function(LT,j)
{
	if(is.null(LT$EVD) & is.null(LT$K))
	{
		text<-"Either variance co-variance matrix K or its eigen-value decomposition\n"
		text<-paste(text,"must be provided for linear term ",j,"\n")
		text<-paste(text,"To specify the variance covariance matrix K use:\n")
		text<-paste(text,"list(K=?,model='RKHS'), where ? is the user defined (between subjects) co-variance matrix\n")
		text<-paste(text,"To specify the eigen-value decomposition for K use:\n")
		text<-paste(text,"list(EVD=?,model='RKHS'), where ? is the output from eigen function for a user defined (between subjects) co-variance matrix\n")
		stop(text)
	}
	
	if((!is.null(LT$K)) & (!is.null(LT$EVD)))
	{
		message("Variance covariance matrix K and its eigen-value decomposition for linear term ",j, " was provided")
		message("ONLY EVD will be used")
		LT$K<-NULL
	}
	
	if((!is.null(LT$K)) & is.null(LT$EVD))
	{
		message("Checking variance co-variance matrix K  for linear term ",j)
		if(nrow(LT$K)!=ncol(LT$K)) stop("variance covariance matrix must be square")
		LT$EVD <- eigen(LT$K,symmetric=TRUE)
		message("Ok")
	}
	
	if(is.null(LT$K) & (!is.null(LT$EVD)))
	{
		message("Checking EVD provided for linear term ",j)
		if(!is.matrix(LT$EVD$vectors)) stop("eigen-vectors must be a matrix\n")
		if(nrow(LT$EVD$vectors)!=ncol(LT$EVD$vectors)) stop("eigen-vectors must be a square matrix\n")
		if(!is.numeric(LT$EVD$values)) stop("eigen-values must be a numeric vector\n")
		message("Ok")
	}
	
	keep <- LT$EVD$values>1e-10
	LT$EVD$vectors <- LT$EVD$vectors[,keep]
	LT$EVD$values <- LT$EVD$values[keep]
	
	#X=Gamma*Lambda^{1/2}
	LT$X<-sweep(x=LT$EVD$vectors,MARGIN=2,STATS=sqrt(LT$EVD$values),FUN="*")
	
	#Return the goodies
	return(LT)
	
}


BLRCross=function(y=NULL,my=NULL, vy=NULL, n=NULL,
                 XX,Xy,nIter=1500,burnIn=500,
                 thin=5,R2=0.5,
                 S0=NULL,df0=5,
                 priors=NULL,
                 idPriors=NULL,
                 verbose=TRUE,
                 saveAt="",
                 rmExistingFiles = TRUE)
{
	if(verbose)
    {
		welcome()
    }
	
	#Check burning and thin
	if(burnIn>=nIter)
	{
		burnIn=as.integer(nIter/2)
		message("burnIn was set to ",burnIn, " because burnIn can not be bigger than nIter")
	}
	
	if(is.null(y))
	{
		message("Since y was not provided, you must provide my, vy and n")
		message("my: sample mean of y")
		message("vy: sample variance of y")
		message("n: sample size")
		
		if(is.null(my)) stop("you must provide my")
		if(is.null(vy)) stop("you must provide vy")
		if(is.null(n)) stop("you must provide n")
		
		#Sum of squares of y
		ssy=(n-1)*vy + n*my^2
	}else{
	
		message("you provided y, so my, vy and n will be ignored (if provided)")
	
		ssy=sum(y^2)
		my=mean(y)
		vy=var(y)
		n=length(y)
	}

	#Assuming all efects are zero, including the intercept, if it exists...
    RSS=ssy
    
    p=ncol(XX)
    
    #Number of predictors in each group
    nCols=table(idPriors)
    
    if(p!=sum(nCols)) stop("The number of columns in X'X is different to the number of elements in idPriors\n")
        
    varE=vy*(1-R2)
    
    if(is.null(S0))
    {
    	S0=varE*(df0+2)
    	message("S0 was missing and was set to ",S0)
    }
    
    if(is.null(priors)) stop("priors can not be NULL\n")
    
    if(!is.list(priors)) stop("priors should be a list\n")
    
    nLT = length(priors)
    
    if(!(nLT>0)) stop("priors should have at least one component\n")
    
    if(is.null(names(priors)))
    {
    	names(priors)=rep("",nLT)
    }
    
    
    #Setting the linear terms
    
    #Create an empty list
    ETA=list()
    	
    #Loop over the components for the linear terms and set up 
    #hyperparameters

    	
    for(j in 1:nLT)
    {
    		diagonal=XX[1 + 0L:(p - 1L) * (p + 1)]
    		sumVarX=sum(diagonal[idPriors==j])/n
    		idColumns=which(idPriors==j)
    		
    		if(!(priors[[j]]$model %in% c("FIXED", "BRR", "BayesA", "BayesB","BayesC","SSVS","RKHS"))) 
            {
                stop("Error in priors[[", j, "]]", " model ", priors[[j]]$model, " not implemented (note: evaluation is case sensitive)")
            }
            
            if(names(priors)[j]=="")
            {
            	priors[[j]]$Name=paste("ETA_",j,sep="")
            }else{
            	if(names(priors)[j]!="___INTERCEPT___")
            	{
            		#Dirty hack to avoid changing Name, see setLT.Fixed.Cross
            		priors[[j]]$Name=paste("ETA_",names(priors)[j],sep="")
            	}
            }
    		
    		ETA[[j]]=switch(priors[[j]]$model,
    						BayesA=setLT.BayesA.Cross(prior=priors[[j]],vy=vy,j=j,p=nCols[j],
    											idColumns=idColumns,sumVarX=sumVarX,
    											R2=R2,nLT=nLT,verbose=verbose,
    											saveAt=saveAt,
    											rmExistingFiles=rmExistingFiles,
    											thin=thin,
    											nIter=nIter,
    											burnIn=burnIn),
    						BayesB=setLT.BayesB.Cross(prior=priors[[j]],vy=vy,j=j,p=nCols[j],
    											idColumns=idColumns,sumVarX=sumVarX,
    											R2=R2,nLT=nLT,verbose=verbose,
    											saveAt=saveAt,
    											rmExistingFiles=rmExistingFiles,
    											thin=thin,
    											nIter=nIter,
    											burnIn=burnIn),
    						BayesC=setLT.BayesC.Cross(prior=priors[[j]],vy=vy,j=j,p=nCols[j],
    											idColumns=idColumns,sumVarX=sumVarX,
    											R2=R2,nLT=nLT,verbose=verbose,
    											saveAt=saveAt,
    											rmExistingFiles=rmExistingFiles,
    											thin=thin,
    											nIter=nIter,
    											burnIn=burnIn),
    						BRR=setLT.BRR.Cross(prior=priors[[j]],vy=vy,j=j,p=nCols[j],
    											idColumns=idColumns,sumVarX=sumVarX,
    											R2=R2,nLT=nLT,verbose=verbose,
    											saveAt=saveAt,
    											rmExistingFiles=rmExistingFiles,
    											thin=thin,
    											nIter=nIter,
    											burnIn=burnIn),
    						RKHS=setLT.BRR.Cross(prior=priors[[j]],vy=vy,j=j,p=nCols[j],
    											idColumns=idColumns,sumVarX=sumVarX,
    											R2=R2,nLT=nLT,verbose=verbose,
    											saveAt=saveAt,
    											rmExistingFiles=rmExistingFiles,
    											thin=thin,
    											nIter=nIter,
    											burnIn=burnIn),
    						FIXED=setLT.Fixed.Cross(p=nCols[j],idColumns=idColumns,
    						                  Name=priors[[j]]$Name,
    						                  saveAt=saveAt,
    						                  rmExistingFiles=rmExistingFiles),
    						SSVS=setLT.SSVS.Cross(prior=priors[[j]],vy=vy,j=j,p=nCols[j],
    											idColumns=idColumns,sumVarX=sumVarX,
    											R2=R2,nLT=nLT,verbose=verbose,
    											saveAt=saveAt,
    											rmExistingFiles=rmExistingFiles,
    											thin=thin,
    											nIter=nIter,
    											burnIn=burnIn)
    		               )
    }
    
    #File for storing sample for varE 

    fname = paste(saveAt, "varE.dat", sep = "")

    if(rmExistingFiles) {
        unlink(fname)
    }

    fileOutVarE = file(description = fname, open = "w")
    
    post_varE = 0
    post_varE2 = 0
    
        
    #Gibbs sampler
    #Loop over iterations
    
    nSums=0		#For running means
    
    #VERY IMPORTANT DO NOT CHANGE THIS LINE IN THE INITIALIZATION
    beta=rep(0,p)
    dots=rep(0,p) 
        
    for(i in 1:nIter)
    {	
    	start=proc.time()[3]
    
    	#Loop over linear predictors
    	for(j in 1:nLT)
    	{
    		
    		if(priors[[j]]$model=="BayesA")
    		{
    			ans=.Call("sampler_others",p, XX, Xy, ETA[[j]]$idColumns, ETA[[j]]$p, beta,
    			           rep(ETA[[j]]$varB,nCols[j]), varE,RSS)
    			
    			beta=ans[[1]]
    			ETA[[j]]$b=beta[ETA[[j]]$idColumns]
    			RSS=ans[[2]]
    			
    			#Update variances
    			SS = ETA[[j]]$S + ETA[[j]]$b^2
                DF = ETA[[j]]$df0 + 1
                ETA[[j]]$varB = SS/rchisq(n = ETA[[j]]$p, df = DF)
                
                tmpShape=ETA[[j]]$p*ETA[[j]]$df0/2+ETA[[j]]$shape0
                tmpRate=sum(1/ETA[[j]]$varB)/2+ETA[[j]]$rate0
                ETA[[j]]$S=rgamma(shape=tmpShape,rate=tmpRate,n=1)
                
    		}#End of BayesA
    		
    		if(priors[[j]]$model=="BayesB")
    		{
    			#ans=.Call("sampler_DiracSS",p, XX, Xy, ETA[[j]]$idColumns, ETA[[j]]$p, ETA[[j]]$a, beta,
    			#           ETA[[j]]$d, ETA[[j]]$varB, varE, ETA[[j]]$probIn,RSS)
                       
 			ans=.Call("sampler_DiracSS_v2",p, XX, Xy, ETA[[j]]$idColumns, ETA[[j]]$p, ETA[[j]]$a, beta,
    			           ETA[[j]]$d, ETA[[j]]$varB, varE, ETA[[j]]$probIn,RSS,dots)

    			
    			ETA[[j]]$a=ans[[1]]
    			ETA[[j]]$d=ans[[2]]
    			beta=ans[[3]]
    			ETA[[j]]$b=beta[ETA[[j]]$idColumns]
    			RSS=ans[[4]]
    			
    			#Sampling hyper-parameters   
    			SS=(ETA[[j]]$a)^2 + ETA[[j]]$S
    			DF=ETA[[j]]$df0+1
    			ETA[[j]]$varB=SS/rchisq(n=ETA[[j]]$p,df=DF)
    			
    			#Update the scale
    			tmpShape=ETA[[j]]$p*ETA[[j]]$df0/2+ETA[[j]]$shape0
                	tmpRate=sum(1/ETA[[j]]$varB)/2+ETA[[j]]$rate0
                	ETA[[j]]$S=rgamma(shape=tmpShape,rate=tmpRate,n=1)
                                
                	#Update inclusion probabilities            
                	mrkIn = sum(ETA[[j]]$d)
                	ETA[[j]]$probIn=rbeta(shape1 = (mrkIn + ETA[[j]]$countsIn + 1),
                                                shape2 = (ETA[[j]]$p - mrkIn + ETA[[j]]$countsOut + 1), n = 1)
    		                
    		}#End of BayesB
    		
    		#BayesC case
    		if(priors[[j]]$model=="BayesC")
    		{
    			
    			#ans=.Call("sampler_DiracSS",p, XX, Xy, ETA[[j]]$idColumns, ETA[[j]]$p, ETA[[j]]$a, beta,
    			#           ETA[[j]]$d, rep(ETA[[j]]$varB,nCols[j]), varE, ETA[[j]]$probIn,RSS)
	
    			ans=.Call("sampler_DiracSS_v2",p, XX, Xy, ETA[[j]]$idColumns, ETA[[j]]$p, ETA[[j]]$a, beta,
    			           ETA[[j]]$d, rep(ETA[[j]]$varB,nCols[j]), varE, ETA[[j]]$probIn,RSS,dots)
    		
    			ETA[[j]]$a=ans[[1]]
    			ETA[[j]]$d=ans[[2]]
    			beta=ans[[3]]
    			ETA[[j]]$b=beta[ETA[[j]]$idColumns]
    			RSS=ans[[4]]
    			
    			#Sampling hyper-parameters   
    			S=sum((ETA[[j]]$a)^2)+ETA[[j]]$S0
    			ETA[[j]]$varB=S/rchisq(n=1,df=ETA[[j]]$df0+nCols[j])
			                
                mrkIn = sum(ETA[[j]]$d)
                ETA[[j]]$probIn=rbeta(shape1 = (mrkIn + ETA[[j]]$countsIn + 1),
                                                shape2 = (ETA[[j]]$p - mrkIn + ETA[[j]]$countsOut + 1), n = 1)
    		} #End of BayesC
    		
    		if(priors[[j]]$model=="SSVS")
    		{	    			
    			ans=.Call("sampler_ACSS", p, XX, Xy, ETA[[j]]$idColumns, ETA[[j]]$p, ETA[[j]]$a, beta,
    			                 ETA[[j]]$d, rep(ETA[[j]]$varB, nCols[j]), varE, ETA[[j]]$probIn, 
    			                 RSS, ETA[[j]]$c)
    			
    			ETA[[j]]$a=ans[[1]]
    			ETA[[j]]$d=ans[[2]]
    			beta=ans[[3]]
    			ETA[[j]]$b=beta[ETA[[j]]$idColumns]
    			RSS=ans[[4]]
    			
    			#Sampling hyper-parameters
    			S=sum((ETA[[j]]$b/ETA[[j]]$a)^2)+ETA[[j]]$S0
    			ETA[[j]]$varB=S/rchisq(n=1,df=ETA[[j]]$df0+nCols[j])
			                
                mrkIn = sum(ETA[[j]]$d)
                ETA[[j]]$probIn=rbeta(shape1 = (mrkIn + ETA[[j]]$countsIn + 1),
                                                shape2 = (ETA[[j]]$p - mrkIn + ETA[[j]]$countsOut + 1), n = 1)
                                                
            	#Update c
            	ETA[[j]]$c=metropc(ETA[[j]]$c,ETA[[j]]$varB,ETA[[j]]$b,ETA[[j]]$d,ETA[[j]]$ccountsIn,ETA[[j]]$ccountsOut)
            	#message("c=",ETA[[j]]$c)
    			
    		} #End of SSVS
    		
    		#BRR case
    		if(priors[[j]]$model%in%c("BRR","RKHS"))
    		{
    			
    			
    			ans=.Call("sampler_others",p, XX, Xy, ETA[[j]]$idColumns, ETA[[j]]$p, beta,
    			           rep(ETA[[j]]$varB,nCols[j]), varE,RSS)
    			
    			beta=ans[[1]]
    			ETA[[j]]$b=beta[ETA[[j]]$idColumns]
    			RSS=ans[[2]]
    			
    			#Sampling hyper-parameters
    			DF = ETA[[j]]$df0 + ETA[[j]]$p
                SS = sum(ETA[[j]]$b^2) + ETA[[j]]$S0
                ETA[[j]]$varB = SS/rchisq(df = DF, n = 1)	
    		} #End of BRR
    		
    		if(priors[[j]]$model=="FIXED")
    		{
    			ans=.Call("sampler_others",p, XX, Xy, ETA[[j]]$idColumns, ETA[[j]]$p, beta,
    			           rep(ETA[[j]]$varB,nCols[j]), varE,RSS)
    			
    			beta=ans[[1]]
    			#print(beta)
    			ETA[[j]]$b=beta[ETA[[j]]$idColumns]
    			RSS=ans[[2]]
    			#message("RSS=",RSS)
    		} #End of FIXED
    	}
    	
    	#Update residual
    	varE=(RSS+S0)/rchisq(n=1,df=n+df0)
    	
    	#Saving samples and computing running means
    	if(i%%thin==0)
    	{
    		
    		for(j in 1:nLT)
    		{
    			if(priors[[j]]$model == "BayesA") 
    			{
                    tmp=ETA[[j]]$S
                    write(tmp, ncolumns = 1, file = ETA[[j]]$fileOut, append = TRUE)
                }
                
                if(priors[[j]]$model == "BayesB")
                {
                        tmp=c(ETA[[j]]$probIn,ETA[[j]]$S)
                        write(tmp, ncolumns = 2, file = ETA[[j]]$fileOut, append = TRUE)
                }
                
                if(priors[[j]]$model == "BayesC") 
                {
                    tmp = c(ETA[[j]]$probIn, ETA[[j]]$varB)
                    write(tmp, ncolumns = 2, file = ETA[[j]]$fileOut, append = TRUE)
                }
                
                if(priors[[j]]$model == "SSVS")
                {
                	tmp = c(ETA[[j]]$probIn, ETA[[j]]$varB, ETA[[j]]$c)
                	write(tmp, ncolumns = 3, file = ETA[[j]]$fileOut, append = TRUE)
                	
                }
                
                if (priors[[j]]$model %in% c("BRR","RKHS")) 
                {
                    write(ETA[[j]]$varB, file = ETA[[j]]$fileOut, append = TRUE)
                }
                
                if (priors[[j]]$model == "FIXED") {
                    write(ETA[[j]]$b,ncolumns=ETA[[j]]$p, file = ETA[[j]]$fileOut, append = TRUE)
                }
                
    		}
    		
    		#Write output file for varE
    		write(x = varE, ncolumns=1,file = fileOutVarE, append = TRUE)
    		
    	
    		if(i>burnIn)
    		{
    			nSums = nSums + 1
                k = (nSums - 1)/(nSums)
    		
    			for(j in 1:nLT)
    			{
    				#BayesA
    				if(priors[[j]]$model=="BayesA")
    				{
    					ETA[[j]]$post_b = ETA[[j]]$post_b * k + ETA[[j]]$b/nSums
                      	ETA[[j]]$post_b2 = ETA[[j]]$post_b2 * k + (ETA[[j]]$b^2)/nSums
                      	ETA[[j]]$post_varB = ETA[[j]]$post_varB * k + (ETA[[j]]$varB)/nSums
                      	ETA[[j]]$post_varB2 = ETA[[j]]$post_varB2 * k + (ETA[[j]]$varB^2)/nSums
                      	ETA[[j]]$post_S = ETA[[j]]$post_S * k + (ETA[[j]]$S)/nSums
		      		  	ETA[[j]]$post_S2 = ETA[[j]]$post_S2 * k + (ETA[[j]]$S^2)/nSums
		      		  	
		      		  	if(ETA[[j]]$saveEffects&&(i%%ETA[[j]]$thin)==0){
                            writeBin(object=ETA[[j]]$b,con=ETA[[j]]$fileEffects,size=ifelse(ETA[[j]]$storageMode=="single",4,8))
                        }#*#
    				}
    				
    				#BayesB case
    				if(priors[[j]]$model=="BayesB")
    				{
    					ETA[[j]]$post_b = ETA[[j]]$post_b * k + ETA[[j]]$b/nSums
                      	ETA[[j]]$post_b2 = ETA[[j]]$post_b2 * k + (ETA[[j]]$b^2)/nSums
                      	ETA[[j]]$post_varB = ETA[[j]]$post_varB * k + (ETA[[j]]$varB)/nSums
                      	ETA[[j]]$post_varB2 = ETA[[j]]$post_varB2 * k + (ETA[[j]]$varB^2)/nSums
                      	
                      	ETA[[j]]$post_d = ETA[[j]]$post_d * k + (ETA[[j]]$d)/nSums
                        ETA[[j]]$post_probIn = ETA[[j]]$post_probIn * k + (ETA[[j]]$probIn)/nSums
                        ETA[[j]]$post_probIn2 = ETA[[j]]$post_probIn2 * k + (ETA[[j]]$probIn^2)/nSums  	
                        
                        ETA[[j]]$post_S = ETA[[j]]$post_S * k + (ETA[[j]]$S)/nSums
						ETA[[j]]$post_S2 = ETA[[j]]$post_S2 * k + (ETA[[j]]$S^2)/nSums
						
						if(ETA[[j]]$saveEffects&&(i%%ETA[[j]]$thin)==0){
                            writeBin(object=ETA[[j]]$b,con=ETA[[j]]$fileEffects,size=ifelse(ETA[[j]]$storageMode=="single",4,8))
                        }#*#
    				}
    				
    				#BayesC case
    				if(priors[[j]]$model=="BayesC")
    				{
    					ETA[[j]]$post_b = ETA[[j]]$post_b * k + ETA[[j]]$b/nSums
                      	ETA[[j]]$post_b2 = ETA[[j]]$post_b2 * k + (ETA[[j]]$b^2)/nSums
                      	ETA[[j]]$post_varB = ETA[[j]]$post_varB * k + (ETA[[j]]$varB)/nSums
                      	ETA[[j]]$post_varB2 = ETA[[j]]$post_varB2 * k + (ETA[[j]]$varB^2)/nSums
                      	
                      	ETA[[j]]$post_d = ETA[[j]]$post_d * k + (ETA[[j]]$d)/nSums
                        ETA[[j]]$post_probIn = ETA[[j]]$post_probIn * k + (ETA[[j]]$probIn)/nSums
                        ETA[[j]]$post_probIn2 = ETA[[j]]$post_probIn2 * k + (ETA[[j]]$probIn^2)/nSums  
                        
                        if(ETA[[j]]$saveEffects&&(i%%ETA[[j]]$thin)==0){
                            writeBin(object=ETA[[j]]$b,con=ETA[[j]]$fileEffects,size=ifelse(ETA[[j]]$storageMode=="single",4,8))
                        }#*#
    				}
    				
    				#SSVS case
    				if(priors[[j]]$model=="SSVS")
    				{
    					ETA[[j]]$post_b = ETA[[j]]$post_b * k + ETA[[j]]$b/nSums
                      	ETA[[j]]$post_b2 = ETA[[j]]$post_b2 * k + (ETA[[j]]$b^2)/nSums
                      	ETA[[j]]$post_varB = ETA[[j]]$post_varB * k + (ETA[[j]]$varB)/nSums
                      	ETA[[j]]$post_varB2 = ETA[[j]]$post_varB2 * k + (ETA[[j]]$varB^2)/nSums
                      	
                      	ETA[[j]]$post_d = ETA[[j]]$post_d * k + (ETA[[j]]$d)/nSums
                        ETA[[j]]$post_probIn = ETA[[j]]$post_probIn * k + (ETA[[j]]$probIn)/nSums
                        ETA[[j]]$post_probIn2 = ETA[[j]]$post_probIn2 * k + (ETA[[j]]$probIn^2)/nSums
                        
                        ETA[[j]]$post_c = ETA[[j]]$post_c * k + ETA[[j]]$c/nSums
                      	ETA[[j]]$post_c2 = ETA[[j]]$post_c2 * k + (ETA[[j]]$c^2)/nSums
                        
                        if(ETA[[j]]$saveEffects&&(i%%ETA[[j]]$thin)==0){
                            writeBin(object=ETA[[j]]$b,con=ETA[[j]]$fileEffects,size=ifelse(ETA[[j]]$storageMode=="single",4,8))
                        }#*#
    				}
    			
    				#BRR case
    				if(priors[[j]]$model%in%c("BRR","RKHS"))
    				{
						ETA[[j]]$post_b = ETA[[j]]$post_b * k + ETA[[j]]$b/nSums
                      	ETA[[j]]$post_b2 = ETA[[j]]$post_b2 * k + (ETA[[j]]$b^2)/nSums
                      	ETA[[j]]$post_varB = ETA[[j]]$post_varB * k + (ETA[[j]]$varB)/nSums
                      	ETA[[j]]$post_varB2 = ETA[[j]]$post_varB2 * k + (ETA[[j]]$varB^2)/nSums  
                      	
                      	if(ETA[[j]]$saveEffects&&(i%%ETA[[j]]$thin)==0){
                            writeBin(object=ETA[[j]]$b,con=ETA[[j]]$fileEffects,size=ifelse(ETA[[j]]$storageMode=="single",4,8))
                        }#*#
    				}
    				
    				#FIXED case
    				if(priors[[j]]$model=="FIXED")
    				{
    					ETA[[j]]$post_b = ETA[[j]]$post_b * k + ETA[[j]]$b/nSums
                      	ETA[[j]]$post_b2 = ETA[[j]]$post_b2 * k + (ETA[[j]]$b^2)/nSums
    				}
    			}
    			
    			post_varE = post_varE * k + varE/nSums
                post_varE2 = post_varE2 * k + (varE^2)/nSums
    		}
    		
    	}#End if loop for checking that i is a multiple of thin
    	
    	end=proc.time()[3]
    	
    	if(verbose)
    	{
    		message("Iter=",i," Time/Iter=",round(end-start,3))
    		message("varE=",round(varE,4))
    	}
    
    }#End of loop for Gibbs sampler
    
    
    #close output files
    close(fileOutVarE)
    
    for (j in 1:nLT) 
    {
    	if (!is.null(ETA[[j]]$fileOut)) 
    	{
            flush(ETA[[j]]$fileOut)                
            close(ETA[[j]]$fileOut)
            ETA[[j]]$fileOut = NULL
        }
        
        if(!is.null(ETA[[j]]$fileEffects))
        {
            flush(ETA[[j]]$fileEffects)
            close(ETA[[j]]$fileEffects)
            ETA[[j]]$fileEffects = NULL
        }  
    }
    
    #Return the goodies
    
    out=list()
    
    #Renaming/removing objects in ETA ...
    
    for(j in 1:nLT)
    {
    	ETA[[j]]$b = ETA[[j]]$post_b
        ETA[[j]]$SD.b = sqrt(ETA[[j]]$post_b2 - ETA[[j]]$post_b^2)
        ETA[[j]]$varB = ETA[[j]]$post_varB
        ETA[[j]]$SD.varB = sqrt(ETA[[j]]$post_varB2 - (ETA[[j]]$post_varB^2))	
        tmp = which(names(ETA[[j]]) %in% c("post_b", "post_b2","post_varB", "post_varB2"))
        ETA[[j]] = ETA[[j]][-tmp]
        
        if(priors[[j]]$model%in%c("BayesB","BayesC","SSVS"))
        {
        		ETA[[j]]$d=ETA[[j]]$post_d
                ETA[[j]]$probIn=ETA[[j]]$post_probIn
				ETA[[j]]$SD.probIn=sqrt(ETA[[j]]$post_probIn2 - (ETA[[j]]$post_probIn^2))
                tmp = which(names(ETA[[j]]) %in% c("a","post_d", "post_probIn","post_probIn2"))
                ETA[[j]] = ETA[[j]][-tmp]
        }
        
        if(priors[[j]]$model=="SSVS")
        {
        		ETA[[j]]$c=ETA[[j]]$post_c
        		ETA[[j]]$SD.c=sqrt(ETA[[j]]$post_c2 - (ETA[[j]]$post_c^2))
        		tmp = which(names(ETA[[j]]) %in% c("post_c", "post_c2"))
        		ETA[[j]]=ETA[[j]][-tmp]
        }
        
    }
    
    out$ETA=ETA
    
    out$varE = post_varE
    out$SD.varE = sqrt(post_varE2 - post_varE^2)
    
    class(out)="BLRCross"
    
    return(out)
}


#No missing values allowed
#For RKHS we assume that X=Gamma*Lambda^{1/2} is provided...

BLRXy_no_miss<-function(y, intercept=TRUE, ETA, 
               nIter = 1500, burnIn = 500, thin = 5, 
               S0 = NULL, df0 = 5, R2 = 0.5, 
               verbose = TRUE, saveAt="",rmExistingFiles = TRUE) 
{
		
		nLT <- ifelse(is.null(ETA), 0, length(ETA))
		
		if(!(nLT>0)) stop("ETA should have at least one component\n")
		
		if(is.null(names(ETA)))
    	{
    		names(ETA)<-rep("",nLT)
    	}
		
		#Check NA's
		if(any(is.na(y))) stop("NA's not allowed in this implementation")
		
		n<-length(y)
	
		ps<-rep(NA,nLT)
		ns<-rep(NA,nLT)
		
		#Check supported models, number of columns in each element of the list
		for(j in 1:nLT)
		{
			if(!(ETA[[j]]$model %in% c("FIXED", "BRR", "BayesA", "BayesB","BayesC","SSVS","RKHS"))) 
			{
				stop("Error in ETA[[", j, "]]", " model ", ETA[[j]]$model, " not implemented (note: evaluation is case sensitive)")
			}
			
			if(!is.null(ETA[[j]]$X))
			{
				if(is.matrix(ETA[[j]]$X))
				{
					ps[j]<-ncol(ETA[[j]]$X)
					ns[j]<-nrow(ETA[[j]]$X)
				}else{
					stop("Error in ETA[[", j, "]], X should be a matrix")
				}
			}else{
				stop("Error in ETA[[", j, "]], X is NULL")
			}
						
		}#End of checking inputs
		
		if(any(is.na(ps))) stop("Check that every X in ETA is a matrix at least one column and ",n," rows")
		if(any(is.na(ns))) stop("Check that every X in ETA is a matrix at least one column and ",n," rows")
		if(any(n!=ns)) stop("Every X in ETA must have the same number of rows")
		
		low<-cumsum(ps)-ps+1
		up<-cumsum(ps)

		X<-matrix(NA,nrow=n,ncol=sum(ps))
		idPriors<-rep(NA,sum(ps))
		priors<-list()
		
		for(j in 1:nLT)			
		{
			X[1:n,c(low[j]:up[j])]<-ETA[[j]]$X
			idPriors[c(low[j]:up[j])]<-j
			priors[[j]]<-ETA[[j]]
			priors[[j]]$X<-NULL	
		}
		
		tmp<-names(ETA)
		
		#Dirty hack
		if(intercept)
		{
			X<-cbind(X,1)
			idPriors<-c(idPriors,nLT+1)
			priors[[nLT+1]]<-list(model="FIXED")
			priors[[nLT+1]]$Name="___INTERCEPT___"
			tmp<-c(tmp,"___INTERCEPT___")
		}
		
		names(priors)<-tmp
		
		#Crossprod
		XX<-crossprod(X)
		Xy<-as.vector(crossprod(X,y))
		
		out<-BLRCross(y=y,XX=XX,Xy=Xy,nIter=nIter,burnIn=burnIn,
                      thin=thin,R2=R2,
                      S0=S0,df0=df0,
                      priors=priors,
                      idPriors=idPriors,
                      verbose=verbose,
                      saveAt=saveAt,
                      rmExistingFiles = rmExistingFiles)
                      
        if(intercept)
        {
        	out$mu<-out$ETA[[nLT+1]]$b
        	out$SD.mu<-out$ETA[[nLT+1]]$SD.b
        	out$ETA[[nLT+1]]<-NULL
        }
        
        #Compute yHat
        yHat<-rep(0,n)
        
        for(j in 1:nLT)
        {
        	yHat<-yHat+as.vector(ETA[[j]]$X%*%out$ETA[[j]]$b)
        }
        
        if(intercept)
        {
        	yHat<-yHat+out$mu
        }
        
        out$yHat<-yHat
                 
		return(out)
}

#Missing values allowed
BLRXy<-function(y, intercept=TRUE, ETA, 
               nIter = 1500, burnIn = 500, thin = 5, 
               S0 = NULL, df0 = 5, R2 = 0.5, 
               verbose = TRUE, saveAt="",rmExistingFiles = TRUE)

{
	nLT <- ifelse(is.null(ETA), 0, length(ETA))
	
	n<-length(y)
	
	if(!(nLT>0)) stop("ETA should have at least one component\n")
	
	if(is.null(names(ETA)))
    {
    		names(ETA)<-rep("",nLT)
    }
    
	missings<-any(is.na(y))
	
	if(missings)
	{
		whichNA<-is.na(y)
	}
	
	for(j in 1:nLT)
	{		
			if(!(ETA[[j]]$model %in% c("FIXED", "BRR", "BayesA", "BayesB","BayesC","SSVS","RKHS"))) 
			{
				stop("Error in ETA[[", j, "]]", " model ", ETA[[j]]$model, " not implemented (note: evaluation is case sensitive)")
			}
			
			#We set linear term to fit RKHS as BRR
			if(ETA[[j]]$model=="RKHS")
			{
				ETA[[j]]<-setLT.RKHS.Cross(ETA[[j]],j)
			}
			
			if(!is.null(ETA[[j]]$X))
			{
				if(is.matrix(ETA[[j]]$X))
				{
					if(nrow(ETA[[j]]$X)!=n) stop("Error in ETA[[", j, "]], X should have ",n," rows")
				}else{
					stop("Error in ETA[[", j, "]], X should be a matrix")
				}
				
				#Keep only the rows without missing values
				if(missings)
				{
					ETA[[j]]$X_NAs<-ETA[[j]]$X[whichNA,,drop=FALSE]
					ETA[[j]]$X<-ETA[[j]]$X[!whichNA,,drop=FALSE]
				}
			}else{
				stop("Error in ETA[[", j, "]], X is NULL")
			}
	}
		
	if(missings)
	{
		out<-BLRXy_no_miss(y=y[!whichNA],intercept=intercept,ETA=ETA,
		           		   nIter=nIter, burnIn=burnIn,thin=thin,
		                   S0 = S0, df0=df0, R2=R2,
		                   verbose=verbose,saveAt=saveAt,
		                   rmExistingFiles=rmExistingFiles)
		yHat<-rep(0,n)                   
		yHat[!whichNA]<-out$yHat
		
		for(j in 1:nLT)
        {
        	yHat[whichNA]<-yHat[whichNA]+as.vector(ETA[[j]]$X_NAs%*%out$ETA[[j]]$b)
        }
        
        if(intercept)
        {
        	yHat[whichNA]<-yHat[whichNA]+out$mu
        }
		
		out$yHat<-yHat
	}else{
	
		out<-BLRXy_no_miss(y=y,intercept=intercept,ETA=ETA,
		           		   nIter=nIter, burnIn=burnIn,thin=thin,
		                   S0 = S0, df0=df0, R2=R2,
		                   verbose=verbose,saveAt=saveAt,
		                   rmExistingFiles=rmExistingFiles)
		                   
	}
	
	return(out)
	
}

