#This function creates an incidence matrix that will be included in the 
#linear term of the model
#Arguments: LT, Linear term, an object of the class "formula" that also includes
#optionally a data.frame to obtain the information
#It returns the incidence matrix
set.X=function(LT)
{

	flag=TRUE
	n_elements=length(LT)
	i=0
	while(i<=n_elements & flag)
	{
	   i=i+1;
	   if(is(LT[[i]],"formula"))
	   {
	   		flag=FALSE
			rhs=LT[[i]]
			if(is.null(LT$data))
			{
				mf = model.frame(formula=rhs)
			}else{
				mf = model.frame(formula=rhs,data=LT$data)
			}
    		X = model.matrix(attr(mf, "terms"), data=mf)
    		Xint = match("(Intercept)", colnames(X), nomatch=0L)
    		if(Xint > 0L) X = X[, -Xint, drop=FALSE]
	   }
	}
	if(flag) stop("Unable to build incidence matrix, wrong formula or data")
	return(X)
}

## Fixed Effects ##################################################################
#Function for initializing regression coefficients for Fixed effects.
#All the arguments are defined in the function BGLR
setLT.Fixed=function(LT,n,j,y,weights,nLT,saveAt,rmExistingFiles,groups,nGroups)
{

    if(is.null(LT$X)) LT$X=set.X(LT)

    LT$X=as.matrix(LT$X)
    LT$p=ncol(LT$X)
    LT$colNames=colnames(LT$X)
	
    if(any(is.na(LT$X)))
    { 
	stop("LP ",j," has NAs in X")
    }

    if(nrow(LT$X)!=n)
    {
        stop("Number of rows of LP ",j,"  not equal to the number of phenotypes")
    }

    #weight inputs if necessary
        
    LT$X=sweep(LT$X,1L,weights,FUN="*")        #weights

    if(!is.null(groups))
    {
        x2=matrix(NA,nrow=nGroups,ncol=ncol(LT$X))
        for(g in 1:nGroups)
        {
                x2[g,]=apply(LT$X[groups==g,,drop=FALSE],2L,function(x) sum(x^2)) #the sum of the square of each of the columns for each group
        }
        LT$x2=x2;
    }else{
        LT$x2=apply(LT$X,2L,function(x) sum(x^2))  #the sum of the square of each of the columns
    }	
    
    #Objects for saving posterior means from MCMC
    LT$b=rep(0,LT$p)
    LT$post_b=rep(0,LT$p)
    LT$post_b2=rep(0,LT$p)   

    fname=paste(saveAt,LT$Name,"_b.dat",sep="")

    LT$NamefileOut=fname; 

    if(rmExistingFiles)
    { 
       unlink(fname) 
    }

    LT$fileOut=file(description=fname,open="w")
    tmp=LT$colNames
    write(tmp, ncolumns = LT$p, file = LT$fileOut, append = TRUE)
    LT$X=as.vector(LT$X)
    LT$x2=as.vector(LT$x2)
    LT$varB=1e10
    return(LT)
}

## Gaussian Regression ############################################################
#Function for initializing regression coefficients for Ridge Regression.
#All the arguments are defined in the function BGLR
setLT.BRR=function(LT,y,n,j,weights,nLT,R2,saveAt,rmExistingFiles,groups,nGroups,verbose,thin,nIter,burnIn,lower_tri){ #*#

    #Check inputs

    if(is.null(LT$X)) LT$X=set.X(LT)
       
    LT$X=as.matrix(LT$X)
    LT$p=ncol(LT$X)
    LT$colNames=colnames(LT$X)
	
    if(any(is.na(LT$X)))
    { 
      stop("LP ",j," has NAs in X")
    }
    
    if(nrow(LT$X)!=n)
    {
      stop("Number of rows of LP ",j,"  not equal to the number of phenotypes")
    }   

    #Weight inputs if necessary
    LT$X=sweep(LT$X,1L,weights,FUN="*")  #weights

    if(!is.null(groups))
    {
    
	x2=matrix(NA,nrow=nGroups,ncol=ncol(LT$X))
	for(g in 1:nGroups)
	{
		x2[g,]=apply(LT$X[groups==g,,drop=FALSE],2L,function(x) sum(x^2)) #the sum of the square of each of the columns for each group
	}
        LT$x2=x2;
    }else{
	LT$x2=apply(LT$X,2L,function(x) sum(x^2))  #the sum of the square of each of the columns
    }  

    sumMeanXSq = sum((apply(LT$X,2L,mean))^2)

    #Default df for the prior assigned to the variance of marker effects
    if(is.null(LT$df0))
    {
	LT$df0=5

	if(verbose)
	{
		message("Degree of freedom of LP ",j,"  set to default value (",LT$df0,")")
	}
    }

    if(is.null(LT$R2))
    { 
        LT$R2=R2/nLT
    }

 
    #Default scale parameter for the prior assigned to the variance of marker effects
    if(is.null(LT$S0))
    {
        if(LT$df0<=0) stop("df0>0 in BRR in order to set S0")

	LT$MSx=sum(LT$x2)/n-sumMeanXSq       
	LT$S0=((var(y,na.rm=TRUE)*LT$R2)/(LT$MSx))*(LT$df0+2)  
	
	if(verbose)
	{
		message("Scale parameter of LP ",j,"  set to default value (",LT$S0,")")
	}
    }

    
    #Objects for saving posterior means from MCMC
    LT$b=rep(0,LT$p)
    LT$post_b=rep(0,LT$p)
    LT$post_b2=rep(0,LT$p)
    LT$varB=LT$S0/(LT$df0+2)
    LT$post_varB=0                 
    LT$post_varB2=0
    
    fname=paste(saveAt,LT$Name,"_varB.dat",sep=""); 
    
    if(rmExistingFiles)
    { 
       unlink(fname) 
    }

    LT$NamefileOut=fname
    LT$fileOut=file(description=fname,open="w")

    if(is.null(LT$lower_tri)) LT$lower_tri=FALSE;

    if(LT$lower_tri)
    {
	message("You have provided a lower triangular matrix for LP ", j)
	message("Checking dimmensions...")
	if(ncol(LT$X)==nrow(LT$X))
	{
		message("Ok.")
		LT$X=LT$X[lower.tri(LT$X,diag=TRUE)]
        }
    }else{
   	 LT$X=as.vector(LT$X)
    }

    LT$x2=as.vector(LT$x2)
    
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

#Ridge regression using sets of markers
#This is just a Ridge Regression set-specific variances, 
#LT has an extra attribute: sets

setLT.BRR_sets=function(LT,y,n,j,weights,nLT,R2,saveAt,rmExistingFiles,verbose,thin,nIter,burnIn)
{

	
    #Check the inputs
    if(is.null(LT$X)) LT$X=set.X(LT)
   
    LT$X=as.matrix(LT$X)
    LT$p=ncol(LT$X)
    LT$colNames=colnames(LT$X)

    if(is.null(LT$sets)) stop("Argument sets (a vector grouping effects into sets) is required in BRR_sets");
  	if(length(LT$sets)!=LT$p){ stop("The length of sets must be equal to the number of predictors") }
	
	LT$sets<-as.integer(factor(LT$sets,ordered=TRUE,levels=unique(LT$sets)))
		
	LT$n_sets=length(unique(LT$sets))	
  	
  	if(LT$n_sets>=LT$p){ stop("The number of sets is greater or equal than the number of effects!") }
    
    if(any(is.na(LT$X))){ 
      stop("LP ",j," has NAs in X")
    }
    
    if(nrow(LT$X)!=n){
      stop("Number of rows of LP ",j,"  not equal to the number of phenotypes")
    }   

    #Weight inputs if necessary
    LT$X=sweep(LT$X,1L,weights,FUN="*")  #weights
    LT$x2=apply(LT$X,2L,function(x) sum(x^2))  #the sum of the square of each of the columns
    sumMeanXSq = sum((apply(LT$X,2L,mean))^2)
    

    if(is.null(LT$df0)){
		LT$df0=5
		if(verbose){
			message("Degree of freedom of LP ",j,"  set to default value (",LT$df0,")")
		}
    }

    if(is.null(LT$R2)){ 
        LT$R2=R2/nLT
    }

    if(is.null(LT$S0))    {
        if(LT$df0<=0) stop("df0 must be greater than 0")

		LT$MSx=sum(LT$x2)/n-sumMeanXSq       
		LT$S0=((var(y,na.rm=TRUE)*LT$R2)/(LT$MSx))*(LT$df0+2) 
		if(verbose){ 
			message("Scale parameter of LP ",j,"  set to default value (",LT$S0,")")
		}
    }
    
    LT$DF1=table(LT$sets)+LT$df0
    
    LT$b=rep(0,LT$p)
    LT$post_b=rep(0,LT$p)
    LT$post_b2=rep(0,LT$p)
    LT$varB=rep(LT$S0/(LT$df0+2),LT$p)
    LT$varSets=rep(0,LT$n_sets)
    LT$post_varSets=rep(0,LT$n_sets)
    LT$post_varSets2<-rep(0,LT$n_sets)
    LT$post_varB=rep(0 ,LT$p)                
    LT$post_varB2=rep(0,LT$p)

    fname=paste(saveAt,LT$Name,"_varB.dat",sep=""); 
    
    if(rmExistingFiles){ 
       unlink(fname) 
    }

    LT$NamefileOut=fname
    LT$fileOut=file(description=fname,open="w")
    LT$X=as.vector(LT$X)

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
    }
    return(LT)
}

## Bayesian LASSO ############################################################
## The well known Bayesian LASSO (Park and Casella, 2008) and 
## de los Campos et al (2009). 
#  This functions simply sets hyper-parameters for quantities involved in BL regression

setLT.BL=function(LT,y,n,j,weights,nLT,R2,saveAt,rmExistingFiles,verbose,thin,nIter,burnIn)
{
    #Check the inputs
    if(is.null(LT$minAbsBeta)) LT$minAbsBeta=1e-9

    if(is.null(LT$X)) LT$X=set.X(LT)

    LT$X=as.matrix(LT$X)
    LT$p=ncol(LT$X)
    LT$colNames=colnames(LT$X)
		
    if(any(is.na(LT$X)))
    {
       stop("LP ",j," has NAs in X")
    }

    if(nrow(LT$X)!=n)
    {
        stop("Number of rows of LP ",j,"  not equal to the number of phenotypes")
    } 

    #Wheight inputs if necessary
    LT$X=sweep(LT$X,1L,weights,FUN="*")  #weights
    LT$x2=apply(LT$X,2L,function(x) sum(x^2))  #the sum of the square of each of the columns
    sumMeanXSq = sum((apply(LT$X,2L,mean))^2)	

    LT$MSx=sum(LT$x2)/n-sumMeanXSq

    # Prior
    if(is.null(LT$R2))
    { 
       LT$R2=R2/nLT
    }
	
    # Setting default value of lambda
    if(!is.null(LT$lambda))
    { 
    	if(LT$lambda<0)
        {
		stop("lambda should be positive")
	}  
    }
    if(is.null(LT$lambda))
    {
		LT$lambda2=2*(1-R2)/(LT$R2)*LT$MSx
		LT$lambda=sqrt(LT$lambda2)
		if(verbose)
		{
			message("Initial value of lambda in LP ",j," was set to default value (",LT$lambda,")")
		}
    }else{
	if(LT$lambda<0) stop("lambda should be positive");
        LT$lambda2=LT$lambda^2
    }
	
    # Checking lambda-type
    if(is.null(LT$type))
    {
		LT$type="gamma"
		if(verbose)
		{
			message("By default, the prior density of lambda^2 in the LP ",j,"  was set to gamma")
		}
    }else{
		if(!LT$type%in%c("gamma","beta","FIXED")) stop("The prior for lambda^2 should be gamma, beta or a point of mass (i.e., fixed lambda)")
    }
    if(LT$type=="gamma")
    {
		if(is.null(LT$shape))
                {
			 LT$shape=1.1
			 if(verbose)
			 {
			 	message("shape parameter in LP ",j," was missing and was set to ",LT$shape)
			 }
		}
		
		if(is.null(LT$rate))
                {
			 LT$rate=(LT$shape-1)/LT$lambda2
			 if(verbose)
			 {
			 	message("rate parameter in LP ",j," was missing and was set to ",LT$rate)
			 }
		}	
    }
    
    if(LT$type=="beta")
    {
                if(is.null(LT$probIn))
  		{
    			LT$probIn=0.5
			if(verbose)
			{
    				message("probIn in LP ",j," was missing and was set to ",LT$probIn)
			}
  		}

  		if(is.null(LT$counts))
  		{
    			LT$counts=2
			if(verbose)
			{
    				message("Counts in LP ",j," was missing and was set to ",LT$counts)
			}
  		} 

                LT$shape1=LT$probIn*LT$counts;
                LT$shape2=(1-LT$probIn)*LT$counts;

		if(is.null(LT$max))
		{
		    LT$max=10*LT$lambda
		    if(verbose)
		    {
		    	message("max parameter in LP ",j," was missing and was set to ",LT$max)
		    }
		}
    }

    #Objects to storing information for MCMC iterations

    LT$b=rep(0,LT$p)
    LT$post_b=rep(0,LT$p)
    LT$post_b2=rep(0,LT$p)
    
    tmp=((var(y,na.rm=TRUE)*R2/nLT)/(LT$MSx))
    LT$tau2=rep(tmp,LT$p)
    LT$post_tau2=0  
    LT$post_lambda=0
    
    fname=paste(saveAt,LT$Name,"_lambda.dat",sep="");
    
    if(rmExistingFiles)
    { 
       unlink(fname) 
    }
    
    LT$NamefileOut=fname
    LT$fileOut=file(description=fname,open="w")

    LT$X=as.vector(LT$X)
    
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


#Reproducing kernel Hilbert spaces
#This function simply sets hyperparameters and prepares inputs 
#for Reproducing Kernel Hilbert Spaces. The algorithm used here is 
#Fully described in de los Campos et al (2010).

setLT.RKHS=function(LT,y,n,j,weights,saveAt,R2,nLT,rmExistingFiles,verbose)  
{

    #Checking inputs
    if(is.null(LT$V))
    {
        if(is.null(LT$K)) stop("Kernel for linear term ",j, " was not provided, specify it with list(K=?,model='RKHS'), where ? is the kernel matrix")

		if(!is.matrix(LT$K)) stop("Kernel for linear term ",j, " should be a matrix, the kernel provided is of class ", class(LT$K))
		
		LT$K = as.matrix(LT$K)

        if(nrow(LT$K)!=ncol(LT$K)) stop("Kernel for linear term ",j, " is not a square matrix")

		#This code was rewritten to speed up computations
    	#T = diag(weights)   
    	#LT$K = T %*% LT$K %*% T 
        
    	#Weight kernels
		#for(i in 1:nrow(LT$K))
    	#{
		#	#j can not be used as subindex because its value is overwritten
		#	for(m in i:ncol(LT$K))
    	#    {    
		#			LT$K[i,m]=LT$K[i,m]*weights[i]*weights[m];
    	#            LT$K[m,i]=LT$K[i,m]
		#	}
		#}
		
		#Added January 10/2020
		#This is faster than the for loop
		
		LT$K=sweep(sweep(LT$K,1L,weights,"*"),2L,weights,"*")
    
    	tmp =eigen(LT$K,symmetric=TRUE)
    	LT$V =tmp$vectors
    	LT$d =tmp$values
		rm(tmp)
	
    }else{
		if(any(weights!=1))
        { 
			warning("Eigen decomposition for LT",j," was provided and the model involves weights. Note: You should have weighted the kernel before computing eigen(K)") 
        }
    }
    
    #Defaul value for tolD
    #Only those eigenvectors whose eigenvalues> tolD are kept.
    if (is.null(LT$tolD)) 
    {
       LT$tolD = 1e-10
       if(verbose)
       {
       		message("Default value of minimum eigenvalue in LP ",j," was set to ",LT$tolD)
       }
    }
    
    #Removing elements whose eigenvalues < tolD
    tmp= LT$d > LT$tolD
    LT$levelsU = sum(tmp)
    LT$d = LT$d[tmp]
    LT$V = LT$V[, tmp]
    
    #Default degrees of freedom and scale parameter associated with the variance component for marker effect
    if (is.null(LT$df0)) 
    {
      LT$df0 = 5
      if(verbose)
      {
      	message("default value of df0 in LP ",j," was missing and was set to ",LT$df0)
      }
    }
   
    if(is.null(LT$R2))
    { 
           LT$R2=R2/nLT
    }

    if (is.null(LT$S0)) 
    {
          if(LT$df0<=0) stop("df0>0 in RKHS in order to set S0");

	  	  LT$S0=((var(y,na.rm=TRUE)*LT$R2)/(mean(LT$d)))*(LT$df0+2)

	  if(verbose)
	  {
             message("default value of S0 in LP ",j," was missing and was set to ",LT$S0)
	  }
    }
    
    LT$u=rep(0,nrow(LT$V))
    
    LT$varU=LT$S0/(LT$df0+2)
       
    LT$uStar=rep(0, LT$levelsU)
    
    #Output files
    fname=paste(saveAt,LT$Name,"_varU.dat",sep="")
    LT$NamefileOut=fname; 

    if(rmExistingFiles)
    { 
       unlink(fname) 
    }

    #Objects for storing information for MCMC iterations
    LT$fileOut=file(description=fname,open="w")
    LT$post_varU=0
    LT$post_varU2=0
    LT$post_uStar = rep(0, LT$levelsU)
    LT$post_u = rep(0, nrow(LT$V))
    LT$post_u2 = rep(0,nrow(LT$V))
    
    #return object
    return(LT)
}

###Bayes B and C########################################################################################################################################                 

setLT.BayesBandC=function(LT,y,n,j,weights,saveAt,R2,nLT,rmExistingFiles, groups, nGroups,verbose,thin,nIter,burnIn)
{

  model=LT$model
 
  if(is.null(LT$X)) LT$X=set.X(LT)

  #Be sure that your X is a matrix
  LT$X=as.matrix(LT$X)  
  LT$p=ncol(LT$X)
  LT$colNames=colnames(LT$X)

  #Weight inputs if necessary
  LT$X=sweep(LT$X,1L,weights,FUN="*")  #weights

  if(!is.null(groups))
  {

        x2=matrix(NA,nrow=nGroups,ncol=ncol(LT$X))
        for(g in 1:nGroups)
        {
                x2[g,]=apply(LT$X[groups==g,,drop=FALSE],2L,function(x) sum(x^2)) #the sum of the square of each of the columns for each group
        }
        LT$x2=x2;
    }else{
        LT$x2=apply(LT$X,2L,function(x) sum(x^2))  #the sum of the square of each of the columns
  }

  sumMeanXSq = sum((apply(LT$X,2L,mean))^2)
  LT$MSx=sum(LT$x2)/n-sumMeanXSq

  
  if(any(is.na(LT$X))){ stop("LP ",j," has NAs in X") }
  if(nrow(LT$X)!=n){ stop("Number of rows of LP ",j,"  not equal to the number of phenotypes") }
  
  if(is.null(LT$R2))
  {
    LT$R2=R2/nLT
    if(verbose)
    {
      message("R2 in LP ",j," was missing and was set to ",LT$R2)
    }
  }

  #Default value for the degrees of freedom associated with the distribution assigned to the variance
  #of marker effects
  if(is.null(LT$df0))
  {
    LT$df0= 5
    if(verbose)
    {
    	message("DF in LP ",j," was missing and was set to ",LT$df0)
    }
  }


  #Default value for a marker being "in" the model
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
     if(LT$df0<=0) stop("df0>0 in ",model," in order to set S0");
     LT$S0=var(y, na.rm = TRUE)*LT$R2/(LT$MSx)*(LT$df0+2)/LT$probIn
     if(verbose)
     {
     	message("Scale parameter in LP ",j," was missing and was set to ",LT$S0)
     }
  }
 
  LT$b=rep(0, LT$p)
  LT$d=rbinom(n = LT$p, size = 1, prob = LT$probIn)

  if(model=="BayesB")
  {
        if(is.null(LT$shape0))
  	{
        	LT$shape0=1.1
  	}
  	if(is.null(LT$rate0)){
    		LT$rate0=(LT$shape0-1)/LT$S0
  	}
        LT$S=LT$S0
  	LT$varB=rep(LT$S0/(LT$df0+2),LT$p)
        fname=paste(saveAt,LT$Name,"_parBayesB.dat",sep="")
  }else{
	LT$varB = LT$S0
	fname=paste(saveAt,LT$Name,"_parBayesC.dat",sep="")
        
  }

  LT$X=as.vector(LT$X)
  LT$x2=as.vector(LT$x2)

  if(rmExistingFiles)
  { 
      unlink(fname) 
  }

  LT$fileOut=file(description=fname,open="w")
  LT$NamefileOut=fname;
  
  if(model=="BayesB")
  {
	tmp=c('probIn','scale')
   	write(tmp, ncolumns = LT$p, file = LT$fileOut, append = TRUE)
  }

  #Objects for storing MCMC information 
  LT$post_varB=0
  LT$post_varB2=0
  LT$post_d=0
  LT$post_probIn=0
  LT$post_probIn2=0
  LT$post_b=rep(0,LT$p)
  LT$post_b2=rep(0,LT$p)

  if(model=="BayesB")
  {
     LT$post_S=0
     LT$post_S2=0
  }
  
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
  
  #return object
  return(LT) 
}

#Bayes A, Mewissen et al. (2001).
#Prediction of Total Genetic Value Using Genome-Wide Dense Marker Maps
#Genetics 157: 1819-1829, Modified so that the Scale parameter is estimated from data (a gamma prior is assigned)

setLT.BayesA=function(LT,y,n,j,weights,saveAt,R2,nLT,rmExistingFiles,verbose,thin,nIter,burnIn)
{

  #Ckecking inputs 
  if(is.null(LT$X)) LT$X=set.X(LT)
 
  LT$X=as.matrix(LT$X)
  LT$p=ncol(LT$X)
  LT$colNames=colnames(LT$X)

  #Weight inputs if necessary
  LT$X=sweep(LT$X,1L,weights,FUN="*")  #weights
  LT$x2=apply(LT$X,2L,function(x) sum(x^2))  #the sum of the square of each of the columns
  sumMeanXSq = sum((apply(LT$X,2L,mean))^2)
  LT$MSx=sum(LT$x2)/n-sumMeanXSq

  #Default degrees of freedom for the prior assigned to the variance of markers
  if(is.null(LT$df0))
  {
     LT$df0 = 5 
     if(verbose)
     { 
     	message("DF in LP ",j," was missing and was set to ",LT$df0)
     }
  }
  if(is.null(LT$R2))
  {
    LT$R2=R2/nLT
    if(verbose)
    {
    	message("R2 in LP ",j," was missing and was set to ",LT$R2)
    }
  }

  #Defuault scale parameter for the prior assigned to the variance of markers
  if(is.null(LT$S0))
  {
     if(LT$df0<=0) stop("df0>0 in BayesA in order to set S0")
     LT$S0 = var(y, na.rm = TRUE)*LT$R2/(LT$MSx)*(LT$df0+2)
     if(verbose)
     {
     	message("Scale parameter in LP ",j," was missing and was set to ",LT$S0)
     }
  }

  # Improvement: Treat Scale as random, assign a gamma density 
  if(is.null(LT$shape0))
  {
     LT$shape0=1.1
  }

  if(is.null(LT$rate0))
  {    
     LT$rate0=(LT$shape0-1)/LT$S0
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
    
  LT$X=as.vector(LT$X)
  
  #Objects for storing information generated during MCMC iterations
  LT$post_varB=0
  LT$post_varB2=0
  
  LT$post_b=rep(0,LT$p)
  LT$post_b2=rep(0,LT$p)
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
  
  #return object
  return(LT)
}

##################################################################################################
#Just the welcome function that will appear every time that your run the program
welcome=function()
{
  message("\n");
  message("#--------------------------------------------------------------------#");
  message("#        _\\\\|//_                                                     #");
  message("#       (` o-o ')      BGLR v1.1.2                                   #");
  message("#------ooO-(_)-Ooo---------------------------------------------------#");
  message("#                      Bayesian Generalized Linear Regression        #");
  message("#                      Gustavo de los Campos, gdeloscampos@gmail.com #");
  message("#    .oooO     Oooo.   Paulino Perez-Rodriguez, perpdgo@gmail.com    #");
  message("#    (   )     (   )   April, 2024                                   #");
  message("#_____\\ (_______) /_________________________________________________ #");
  message("#      \\_)     (_/                                                   #");
  message("#                                                                    #");
  message("#------------------------------------------------------------------- #");
  message("\n");
}
##################################################################################################

##################################################################################################
#The density of a scaled inverted chi-squered distribution
#df: degrees of freedom, S: Scale parameter
dScaledInvChisq=function (x, df, S)
{
    tmp = dchisq(S/x, df = df)/(x^2)
    return(tmp)
}


##################################################################################################
#The density function for lambda
#Density function for Regularization parameter in Bayesian LASSO
#Rate: rate parameter, shape: the value for the shape parameter
dLambda=function (rate, shape, lambda) 
{
    tmp = dgamma(x = I(lambda^2), rate = rate, shape = shape) * 2 * lambda
    return(tmp)
}

##################################################################################################
#Metropolis sampler for lambda in the Bayesian LASSO
metropLambda=function (tau2, lambda, shape1 = 1.2, shape2 = 1.2, max = 200, ncp = 0)
{
    lambda2 = lambda^2
    l2_new = rgamma(rate = sum(tau2)/2, shape = length(tau2),
        n = 1)
    l_new = sqrt(l2_new)
    logP_old = sum(dexp(x = tau2, log = TRUE, rate = (lambda2/2))) +
        dbeta(x = lambda/max, log = TRUE, shape1 = shape1, shape2 = shape2) -
        dgamma(shape = sum(tau2)/2, rate = length(tau2), x = (2/lambda2),
            log = TRUE)
    logP_new = sum(dexp(x = tau2, log = TRUE, rate = (l2_new/2))) +
        dbeta(x = l_new/max, log = TRUE, shape1 = shape1, shape2 = shape2) -
        dgamma(shape = sum(tau2)/2, rate = length(tau2), x = (2/l2_new),
            log = TRUE)
    accept = (logP_new - logP_old) > log(runif(1))
    if (accept) {
        lambda = l_new
    }
    return(lambda)
}

##################################################################################################
#Startup function
#this function is executed once the library is loaded
.onAttach = function(library, pkg)
{
  if(interactive())
  {
    packageStartupMessage("# Gustavo de los Campos & Paulino Perez-Rodriguez")
    packageStartupMessage("# Support provided by the U.S., National Institutes of Health (NIH)")
    packageStartupMessage("# (Grant: R01GM101219, NIGMS)")
    packageStartupMessage("# and by the International Maize and Wheat Improvement Center (CIMMyT).")                        
    packageStartupMessage("# Type 'help(BGLR)' for summary information")
  }
  invisible()
}


##################################################################################################
#rtrun draws from a truncated univariate normal distribution using the inverse CDF algorithm
#Arguments:
#mu: mean
#sigma: standard deviation
#a: lower bound
#b: upper bound
#NOTES: 1) This routine was taken from bayesm package, December 18, 2012
#       2) The inputs are not checked, 
#It is assumed that are ok.
#rtrun=function (mu, sigma, a, b) 
#{
#    FA = pnorm(((a - mu)/sigma))
#    FB = pnorm(((b - mu)/sigma))
#    return(mu + sigma * qnorm(runif(length(mu)) * (FB - FA) + FA))
#}

#Using the rtruncnorm function in the truncnorm package
rtrun=function(mu,sigma,a,b)
{
        n=max(c(length(mu),length(sigma),length(a),length(b)))
	rtruncnorm(n,a,b,mu,sigma)
}


#Extract the values of z such that y[i]=j
#z,y vectors, j integer
#extract=function(z,y,j) subset(as.data.frame(z,y),subset=(y==j))
extract=function(z,y,j) z[y==j]

#This routine was adapted from rinvGauss function from S-Plus
# Random variates from inverse Gaussian distribution
# Reference:
#      Chhikara and Folks, The Inverse Gaussian Distribution,
#      Marcel Dekker, 1989, page 53.
# GKS  15 Jan 98
#n: Number of samples
#nu: nu parameter
#lambda: lambda parameter

rinvGauss=function(n, nu, lambda)
{
	if(any(nu<=0)) stop("nu must be positive")
	if(any(lambda<=0)) stop("lambda must be positive")
	if(length(n)>1) n = length(n)
	if(length(nu)>1 && length(nu)!=n) nu = rep(nu,length=n)
	if(length(lambda)>1 && length(lambda)!=n) lambda = rep(lambda,length=n)
        tmp = rnorm(n)
	y2 = tmp*tmp
	u = runif(n)
	r1 = nu/(2*lambda) * (2*lambda + nu*y2 - sqrt(4*lambda*nu*y2 + nu*nu*y2*y2))
	r2 = nu*nu/r1
	ifelse(u < nu/(nu+r1), r1, r2)
}


#log-likelihood for ordinal data
#y: response vector
#predicted response vector, yHat=X%*%beta
#threshold
loglik_ordinal=function(y,yHat,threshold)
{
	sum=0
        n=length(y)
	for(i in 1:n)
        {
           sum=sum + log(pnorm(threshold[y[i] + 1]-yHat[i])-pnorm(threshold[y[i]]-yHat[i]))
        }
        return(sum)
}

##################################################################################################

#Arguments:
#y: data vector, NAs allowed
#response_type: can be "gaussian", "ordinal",
#ETA: The linear predictor
#nIter: Number of MCMC iterations
#burnIn: burnIn
#thin: thin
#saveAt: string, where to save the information
#Se: Scale parameter for the prior for varE
#dfe: Degrees of freedom for the prior for varE
#weights: 
#R2
#Note: The function was designed to work with gaussian responses, some changes were made to deal binary and ordinal responses


#To add new method: 
#(a) create setLT, 
#(b) add it to the switch statement,
#(c) add code to update parameters in the Gibbs sampler, 
#(d) add code to save samples
#(e) add code to compute posterior means
#(f) Test:
#(f1) Test simple example without hyper-paramaeters, evaluate how  
#        default values were set
	#(f2)  Check posterior means and files
	#(f3)  Test simple example with some hyper-parameters give and 
	#         some set by default
#(f4) Run an example with a few missing values, compare to BLR 
#       example, check: (i) residual variance, (ii) plot of effects, (iii) plot 
#        of predictions in trn, (iv) plot of prediction in tst.


BGLR=function (y, response_type = "gaussian", a = NULL, b = NULL, 
    ETA = NULL, nIter = 1500, burnIn = 500, thin = 5, saveAt = "", 
    S0 = NULL, df0 = 5, R2 = 0.5, weights = NULL, 
    verbose = TRUE, rmExistingFiles = TRUE, groups=NULL) 
{
   
    if(verbose)
    {
	welcome()
    }

    IDs=names(y)
    if (!(response_type %in% c("gaussian", "ordinal")))  stop("Only gaussian and ordinal responses are allowed")

    if (saveAt == "") {
        saveAt = paste(getwd(), "/", sep = "")
    }

    y=as.vector(y)
    y0=y
    a = as.vector(a)
    b = as.vector(b)
    n = length(y)

    nGroups=1
    if(!is.null(groups))
    {
		groups<-as.character(groups)  #Groups as character and then as factor to avoid dummy levels
		groups<-as.factor(groups)
		#Number of records by group
		countGroups=table(groups)
		nGroups=length(countGroups)
		groupLabels=names(countGroups)
                groups=as.integer(groups)
                ggg=as.integer(groups-1);  #In C we begin to count in 0
		if(sum(countGroups)!=n) stop("length of groups and y differs, NA's not allowed in groups");
    }

    if(response_type=="ordinal")
    {

    	y=factor(y,ordered=TRUE)
        lev=levels(y)
        nclass=length(lev)
        if(nclass==n) stop("The number of classes in y must be smaller than the number of observations");

        y=as.integer(y)
        z=y  

		fname = paste(saveAt, "thresholds.dat", sep = "")
        fileOutThresholds = file(description = fname, open = "w")
    }

    
    if (is.null(weights)) 
    {
        weights = rep(1, n)
    }

    if(!is.null(groups))
    {
      sumW2=tapply(weights^2,groups,"sum")
    }else{
      sumW2 = sum(weights^2)
    }

    nSums = 0

    whichNa = which(is.na(y))
    nNa = length(whichNa)

    Censored = FALSE
     
    if (response_type == "gaussian") 
    {
        if ((!is.null(a)) | (!is.null(b))) 
        {
            Censored = TRUE
            if ((length(a) != n) | (length(b) != n)) stop("y, a and b must have the same dimension")
            if (any(weights != 1)) stop("Weights are only implemented for Gausian uncensored responses")
        }
        mu = weighted.mean(x = y, w = weights, na.rm = TRUE)
    }
    post_mu = 0
    post_mu2 = 0

    fname = paste(saveAt, "mu.dat", sep = "")
    if (rmExistingFiles) 
    {
        unlink(fname)
    }
    else {
        message("Note: samples will be appended to existing files")
    }

    fileOutMu = file(description = fname, open = "w")

    if (response_type == "ordinal") {
        if(verbose){message("Prior for residual is not necessary, if you provided it, it will be ignored")}
        if (any(weights != 1)) stop("Weights are not supported")
       
        countsZ=table(z)

        if (nclass <= 1) stop("Data vector y has only ", nclass, " differente values, it should have at least 2 different values")
        threshold=qnorm(p=c(0,cumsum(as.vector(countsZ)/n)))
          
        y = rtrun(mu =0, sigma = 1, a = threshold[z], b = threshold[ (z + 1)])

        mu=0
        #posterior for thresholds
        post_threshold = 0
        post_threshold2 = 0
        
	post_prob=matrix(nrow=n,ncol=nclass,0)
        post_prob2=post_prob
    }

    post_logLik = 0

    # yStar & yHat
    yStar = y * weights
    yHat = mu * weights
    
    if (nNa > 0) {
        yStar[whichNa] = yHat[whichNa]
    }

    post_yHat = rep(0, n)
    post_yHat2 = rep(0, n)

    # residual and residual variance
    e = (yStar - yHat)

    varE = var(e, na.rm = TRUE) * (1 - R2)

    if (is.null(S0)) {
        S0 = varE * (df0 + 2)
    }

    if(!is.null(groups))
    {
        varE=rep(varE/nGroups,nGroups)
        names(varE)=groupLabels
    }

    sdE = sqrt(varE)


    post_varE = 0
    post_varE2 = 0

    #File for storing sample for varE 

    fname = paste(saveAt, "varE.dat", sep = "")

    if (rmExistingFiles) {
        unlink(fname)
    }

    fileOutVarE = file(description = fname, open = "w")

    nLT = ifelse(is.null(ETA), 0, length(ETA))
    

    #Setting the linear terms
    if (nLT > 0) {
	
	if(is.null(names(ETA)))
    	{ 
             names(ETA)<-rep("",nLT)
    	}

        for (i in 1:nLT) {  

	    if(names(ETA)[i]=="")
	    {
	       	ETA[[i]]$Name=paste("ETA_",i,sep="")
	    }else{
               ETA[[i]]$Name=paste("ETA_",names(ETA)[i],sep="")
	    }

            if (!(ETA[[i]]$model %in% c("FIXED", "BRR", "BL", "BayesA", "BayesB","BayesC", "RKHS","BRR_sets"))) 
            {
                stop("Error in ETA[[", i, "]]", " model ", ETA[[i]]$model, " not implemented (note: evaluation is case sensitive)")
                
            }

            if(!is.null(groups))
            {
		if(!(ETA[[i]]$model %in%  c("BRR","FIXED","BayesB","BayesC"))) stop("Error in ETA[[", i, "]]", " model ", ETA[[i]]$model, " not implemented for groups")
            }

            ETA[[i]] = switch(ETA[[i]]$model, 
			      FIXED = setLT.Fixed(LT = ETA[[i]],  n = n, j = i, weights = weights, y = y, nLT = nLT, saveAt = saveAt, rmExistingFiles = rmExistingFiles,groups=groups,nGroups=nGroups), 
                              BRR = setLT.BRR(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles,groups=groups,nGroups=nGroups,verbose=verbose,thin=thin,nIter=nIter,burnIn=burnIn,lower_tri=ETA[[i]]$lower_tri),#*#
                              BL = setLT.BL(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles,verbose=verbose,thin=thin,nIter=nIter,burnIn=burnIn), 
                              RKHS = setLT.RKHS(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles,verbose=verbose), 
                              BayesC = setLT.BayesBandC(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles,groups=groups,nGroups=nGroups,verbose=verbose,thin=thin,nIter=nIter,burnIn=burnIn),
                              BayesA = setLT.BayesA(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles,verbose=verbose,thin=thin,nIter=nIter,burnIn=burnIn),
                              BayesB = setLT.BayesBandC(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles,groups=groups,nGroups=nGroups,verbose=verbose,thin=thin,nIter=nIter,burnIn=burnIn),
                              BRR_sets = setLT.BRR_sets(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles,verbose=verbose,thin=thin,nIter=nIter,burnIn=burnIn)
                              )
        }
    }

    # Gibbs sampler

    time = proc.time()[3]

    for (i in 1:nIter) {
        # intercept
	if(!is.null(groups))
	{
		e = e + weights * mu 
                varEexpanded=varE[groups]
		#rhs = sum(tapply(e*weights,groups,"sum")/varE)
                rhs = as.numeric(crossprod(e/varEexpanded,weights));
                C = sum(sumW2/varE)
                sol = rhs/C
                mu = rnorm(n = 1, sd = sqrt(1/C)) + sol;
	}else{
        	e = e + weights * mu
        	rhs = sum(weights * e)/varE
        	C = sumW2/varE
        	sol = rhs/C
        	mu = rnorm(n = 1, sd = sqrt(1/C)) + sol
	}
        if (response_type == "ordinal") {
            mu=0 
        }

        e = e - weights * mu
        
        #deltaSS and deltadf for updating varE
        deltaSS = 0
        deltadf = 0

        if (nLT > 0) {
            for (j in 1:nLT) {
                ## Fixed effects ####################################################################
                if (ETA[[j]]$model == "FIXED") {
                  #cat("varB=",ETA[[j]]$varB,"\n");
                  varBj = rep(ETA[[j]]$varB, ETA[[j]]$p)
                  if(!is.null(groups)){
                        ans = .Call("sample_beta_groups", n, ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b,
                                             e, varBj, varE, 1e-9,ggg,nGroups)
		  }else{
                  	ans = .Call("sample_beta", n, ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, 
                                	     e, varBj, varE, 1e-9)
		  }
                  ETA[[j]]$b = ans[[1]]
                  e = ans[[2]]
                }#End of fixed effects

                ## Ridge Regression ##################################################################
                if (ETA[[j]]$model == "BRR") {
                  varBj = rep(ETA[[j]]$varB, ETA[[j]]$p)

                  if(!is.null(groups))
		  {
                        ans = .Call("sample_beta_groups",n, ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, 
                                    e, varBj, varE, 1e-9,ggg,nGroups)
	          }else{
			if(!(ETA[[j]]$lower_tri))
			{
                  		ans = .Call("sample_beta", n, ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, 
                                	             e, varBj, varE, 1e-9)
			}else{
				ans = .Call("sample_beta_lower_tri", n, ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b,
                                                     e, ETA[[j]]$varB, varE, 1e-9)
			}
		  }
                  ETA[[j]]$b = ans[[1]]
                  e = ans[[2]]

                  DF = ETA[[j]]$df0 + ETA[[j]]$p
                  SS = sum(ETA[[j]]$b^2) + ETA[[j]]$S0
                  ETA[[j]]$varB = SS/rchisq(df = DF, n = 1)
                }# END BRR
                
                if(ETA[[j]]$model=="BRR_sets"){
                   	ans = .Call("sample_beta", n, ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b,
                                             e, ETA[[j]]$varB, varE, 1e-9)
		  		   ETA[[j]]$b = ans[[1]]
                   	e = ans[[2]]
			SS=tapply(X=ETA[[j]]$b^2,INDEX=ETA[[j]]$sets,FUN=sum)+ETA[[j]]$S0
				   
                   	ETA[[j]]$varSets=SS/rchisq(df=ETA[[j]]$DF1,n=ETA[[j]]$n_sets)
                   	ETA[[j]]$varB=ETA[[j]]$varSets[ETA[[j]]$sets]
                }


                ## Bayesian LASSO ####################################################################
                if (ETA[[j]]$model == "BL") {
                  
                   varBj = ETA[[j]]$tau2 * varE
                   ans = .Call("sample_beta", n, ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, 
                                e, varBj, varE, ETA[[j]]$minAbsBeta)

                  ETA[[j]]$b = ans[[1]]
                  e = ans[[2]]

                  nu = sqrt(varE) * ETA[[j]]$lambda/abs(ETA[[j]]$b)
                  tmp = NULL
                  try(tmp <- rinvGauss(n = ETA[[j]]$p, nu = nu, lambda = ETA[[j]]$lambda2))
                  if (!is.null(tmp) && !any(tmp<0)) {
                    if (!any(is.na(sqrt(tmp)))) {
                      ETA[[j]]$tau2 = 1/tmp
                    }
                    else {
                      warning(paste("tau2 was not updated in iteration",i, "due to numeric problems with beta\n",sep=" "),immediate. = TRUE)
                    }
                  }
                  else {
                    warning(paste("tau2 was not updated  in iteration",i,"due to numeric problems with beta\n",sep=" "),immediate. = TRUE)
                  }

                  #Update lambda 
                  if (ETA[[j]]$type == "gamma") {
                    rate = sum(ETA[[j]]$tau2)/2 + ETA[[j]]$rate
                    shape = ETA[[j]]$p + ETA[[j]]$shape
                    ETA[[j]]$lambda2 = rgamma(rate = rate, shape = shape, n = 1)
                    if (!is.na(ETA[[j]]$lambda2)) {
                      ETA[[j]]$lambda = sqrt(ETA[[j]]$lambda2)
                    }
                    else {
                      warning(paste("lambda was not updated in iteration",i, "due to numeric problems with beta\n",sep=" "),immediate. = TRUE)
                    }
                  }

                  if (ETA[[j]]$type == "beta") {
                    ETA[[j]]$lambda = metropLambda(tau2 = ETA[[j]]$tau2, 
                                                   lambda = ETA[[j]]$lambda, shape1 = ETA[[j]]$shape1, shape2 = ETA[[j]]$shape2, 
                                                   max = ETA[[j]]$max)
                    ETA[[j]]$lambda2 = ETA[[j]]$lambda^2
                  }

                  deltaSS = deltaSS + sum((ETA[[j]]$b/sqrt(ETA[[j]]$tau2))^2)
                  deltadf = deltadf + ETA[[j]]$p
                }#END BL

                ## RKHS ####################################################################
                if (ETA[[j]]$model == "RKHS") {
                  #error
                  e = e + ETA[[j]]$u
                  rhs = crossprod(ETA[[j]]$V, e)/varE
                  varU = ETA[[j]]$varU * ETA[[j]]$d
                  C = as.numeric(1/varU + 1/varE)
                  SD = 1/sqrt(C)
                  sol = rhs/C
                  tmp = rnorm(n = ETA[[j]]$levelsU, mean = sol, sd = SD)
                  ETA[[j]]$uStar = tmp
                  ETA[[j]]$u = as.vector(ETA[[j]]$V %*% tmp)
		  
                  #update error
                  e = e - ETA[[j]]$u
                   
                  #update the variance
                  tmp = ETA[[j]]$uStar/sqrt(ETA[[j]]$d)
                  SS = as.numeric(crossprod(tmp)) + ETA[[j]]$S0
                  DF = ETA[[j]]$levelsU + ETA[[j]]$df0
                  ETA[[j]]$varU = SS/rchisq(n = 1, df = DF)
                }#END RKHS

                ## BayesA ##############################################################################
                if (ETA[[j]]$model == "BayesA") {
                  varBj = ETA[[j]]$varB
                  ans = .Call("sample_beta", n, ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, 
                                             e, varBj, varE, 1e-9)
                  ETA[[j]]$b = ans[[1]]
                  e = ans[[2]]
                   
                  #Update variances
                  SS = ETA[[j]]$S + ETA[[j]]$b^2
                  DF = ETA[[j]]$df0 + 1
                  ETA[[j]]$varB = SS/rchisq(n = ETA[[j]]$p, df = DF)

                  tmpShape=ETA[[j]]$p*ETA[[j]]$df0/2+ETA[[j]]$shape0
                  tmpRate=sum(1/ETA[[j]]$varB)/2+ETA[[j]]$rate0
                  ETA[[j]]$S=rgamma(shape=tmpShape,rate=tmpRate,n=1)

                }#End BayesA

		#BayesB and BayesC
		if(ETA[[j]]$model %in% c("BayesB","BayesC"))
		{
			#Update marker effects
                      	mrkIn=ETA[[j]]$d==1
                      	pIn=sum(mrkIn)
		        
                        if(ETA[[j]]$model=="BayesB")
                        {
                          if(!is.null(groups))
                          {
                             ans=.Call("sample_beta_BB_BCp_groups",n,ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, ETA[[j]]$d, e, ETA[[j]]$varB, varE, 1e-9, ETA[[j]]$probIn,ggg,nGroups);
                          }else{
                             ans=.Call("sample_beta_BB_BCp",n,ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, ETA[[j]]$d, e, ETA[[j]]$varB, varE, 1e-9, ETA[[j]]$probIn);
                          }
                        }else{
                          if(!is.null(groups))
                          {
                             ans=.Call("sample_beta_BB_BCp_groups",n,ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, ETA[[j]]$d, e, rep(ETA[[j]]$varB,ETA[[j]]$p), varE, 1e-9, ETA[[j]]$probIn,ggg,nGroups);
                          }else{   
                             ans=.Call("sample_beta_BB_BCp",n,ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, ETA[[j]]$d, e, rep(ETA[[j]]$varB,ETA[[j]]$p), varE, 1e-9, ETA[[j]]$probIn);
                          }
                        }

                        ETA[[j]]$d=ans[[1]]
                        e=ans[[2]]
                        ETA[[j]]$b=ans[[3]] 

			#Update the variance component associated with the markers
			if(ETA[[j]]$model=="BayesB")
			{
				SS = ETA[[j]]$b^2 + ETA[[j]]$S
                      		DF = ETA[[j]]$df0+1
                      		ETA[[j]]$varB = SS/rchisq(df=DF, n = ETA[[j]]$p)

                                # Update scale
                                tmpShape=ETA[[j]]$p*ETA[[j]]$df0/2+ETA[[j]]$shape0
                                tmpRate=sum(1/ETA[[j]]$varB)/2+ETA[[j]]$rate0
                                ETA[[j]]$S=rgamma(shape=tmpShape,rate=tmpRate,n=1)

			}else{
				SS = sum(ETA[[j]]$b^2) + ETA[[j]]$S0
				DF = ETA[[j]]$df0 + ETA[[j]]$p
                  		ETA[[j]]$varB = SS/rchisq(df = DF, n = 1)
			}
                      	mrkIn = sum(ETA[[j]]$d)
                      	ETA[[j]]$probIn = rbeta(shape1 = (mrkIn + ETA[[j]]$countsIn + 1),
                                                shape2 = (ETA[[j]]$p - mrkIn + ETA[[j]]$countsOut + 1), n = 1)
		}
            }#Loop for
        }#nLT
        
        # yHat
        yHat = yStar - e
        
       #4#
         # residual variance # missing values
        if (response_type == "gaussian") {
	    
            if(!is.null(groups))
	    {	
        	for(g in 1:nGroups)
        	{
			SS=sum(e[groups==g]^2)+ S0 + deltaSS
                        DF=countGroups[g]+df0+deltadf
			varE[g]=SS/rchisq(n=1,df=DF)
		}
	     }else{
            		SS = sum(e * e) + S0 + deltaSS
            		DF = n + df0 + deltadf
            		varE = SS/rchisq(n = 1, df = DF)
	    }
            sdE = sqrt(varE)
        
          if (nNa > 0) {
            if (Censored) {
                if(!is.null(groups))
                {
                   #FIXME: Double check this, I was testing it and is ok
                   sdEexpanded=sdE[groups]
                   yStar[whichNa] = rtrun(mu = yHat[whichNa], a = a[whichNa], b = b[whichNa], sigma = sdEexpanded)

                }else{
                  yStar[whichNa] = rtrun(mu = yHat[whichNa], a = a[whichNa], b = b[whichNa], sigma = sdE)
                }
            }
            else{
                 if(!is.null(groups))
                 {
                    #FIXME: Double check this, I was testing it and is ok
                    sdEexpanded=sdE[groups]
                    yStar[whichNa] = yHat[whichNa] + rnorm(n = nNa, sd = sdEexpanded)
                 }else{
                    yStar[whichNa] = yHat[whichNa] + rnorm(n = nNa, sd = sdE)
                 }
            }
            e[whichNa] = yStar[whichNa] - yHat[whichNa]
          }
        }else{  #ordinal
            varE = 1
            sdE = 1
            
            #Update yStar, this is the latent variable
            if(nNa==0){
               yStar=rtrun(mu = yHat, sigma = 1, a = threshold[z], b = threshold[(z + 1)])
            }else{
               yStar[-whichNa]=rtrun(mu = yHat[-whichNa], sigma = 1, a = threshold[z[-whichNa]], b = threshold[(z[-whichNa] + 1)])
               yStar[whichNa]=yHat[whichNa] + rnorm(n = nNa, sd = sdE)           
            }

            #Update thresholds           
            if(nNa==0){ 
              for (m in 2:nclass) {
            
                lo = max(max(extract(yStar, z, m - 1)), threshold[m - 1])
                hi = min(min(extract(yStar, z, m)), threshold[m + 1])
                threshold[m] = runif(1, lo, hi)
              }
            }else{

              for (m in 2:nclass) {
                tmpY=yStar[-whichNa]
                tmpZ=z[-whichNa]
                lo = max(max(extract(tmpY, tmpZ, m - 1)), threshold[m - 1])
                hi = min(min(extract(tmpY, tmpZ, m)), threshold[m + 1])
                threshold[m] = runif(1, lo, hi)
              }
            }
            
            #Update error
            e = yStar - yHat
        }

        # Saving samples and computing running means
        if ((i%%thin == 0)) {
            if (nLT > 0) {
                for (j in 1:nLT) {

                  if (ETA[[j]]$model == "FIXED") {
                    write(ETA[[j]]$b,ncolumns=ETA[[j]]$p, file = ETA[[j]]$fileOut, append = TRUE)
                  }

                  if (ETA[[j]]$model == "BRR") {
                    write(ETA[[j]]$varB, file = ETA[[j]]$fileOut, append = TRUE)
                  }
                  if (ETA[[j]]$model == "BRR_sets") {
                    write(ETA[[j]]$varSets, ncolumns=ETA[[j]]$n_sets,file = ETA[[j]]$fileOut, append = TRUE)
                  }
                  if (ETA[[j]]$model == "BL") {
                    write(ETA[[j]]$lambda, file = ETA[[j]]$fileOut, append = TRUE)
                  }

                  if (ETA[[j]]$model == "RKHS") {
                    write(ETA[[j]]$varU, file = ETA[[j]]$fileOut, append = TRUE)
                  }

                  if (ETA[[j]]$model == "BayesC") {
                    tmp = c(ETA[[j]]$probIn, ETA[[j]]$varB)
                    write(tmp, ncolumns = 2, file = ETA[[j]]$fileOut, append = TRUE)
                  }

                  if (ETA[[j]]$model == "BayesA") {
                    tmp=ETA[[j]]$S
                    write(tmp, ncolumns = 1, file = ETA[[j]]$fileOut, append = TRUE)
                  }
                  
                  if(ETA[[j]]$model=="BayesB")
                  {
                        tmp=c(ETA[[j]]$probIn,ETA[[j]]$S)
                        write(tmp, ncolumns = 2, file = ETA[[j]]$fileOut, append = TRUE)
                  }
                }
            }

            #Output files
            write(x = mu, file = fileOutMu, append = TRUE)
            write(x = varE, ncolumns=nGroups,file = fileOutVarE, append = TRUE)

	    if (response_type == "ordinal") {
		  write(x=threshold[2:nclass],ncolumns=nclass-1,file=fileOutThresholds,append=TRUE)
	    }


            if (i > burnIn) {
                nSums = nSums + 1
                k = (nSums - 1)/(nSums)
                if (nLT > 0) {
                  for (j in 1:nLT) {
                    if (ETA[[j]]$model == "FIXED") {
                      ETA[[j]]$post_b = ETA[[j]]$post_b * k + ETA[[j]]$b/nSums
                      ETA[[j]]$post_b2 = ETA[[j]]$post_b2 * k + (ETA[[j]]$b^2)/nSums
                    }

                    if (ETA[[j]]$model == "BRR") {
                      ETA[[j]]$post_b = ETA[[j]]$post_b * k + ETA[[j]]$b/nSums
                      ETA[[j]]$post_b2 = ETA[[j]]$post_b2 * k + (ETA[[j]]$b^2)/nSums
                      ETA[[j]]$post_varB = ETA[[j]]$post_varB * k + (ETA[[j]]$varB)/nSums
                      ETA[[j]]$post_varB2 = ETA[[j]]$post_varB2 * k + (ETA[[j]]$varB^2)/nSums
                      if(ETA[[j]]$saveEffects&&(i%%ETA[[j]]$thin)==0){
                          writeBin(object=ETA[[j]]$b,con=ETA[[j]]$fileEffects,size=ifelse(ETA[[j]]$storageMode=="single",4,8))
                      }#*#
                    }

                    if (ETA[[j]]$model == "BRR_sets") {
                      ETA[[j]]$post_b = ETA[[j]]$post_b * k + ETA[[j]]$b/nSums
                      ETA[[j]]$post_b2 = ETA[[j]]$post_b2 * k + (ETA[[j]]$b^2)/nSums
                      ETA[[j]]$post_varB = ETA[[j]]$post_varB * k + (ETA[[j]]$varB)/nSums
                      ETA[[j]]$post_varB2 = ETA[[j]]$post_varB2 * k + (ETA[[j]]$varB^2)/nSums
                      ETA[[j]]$post_varSets<-ETA[[j]]$post_varSets*k+ETA[[j]]$varSets/nSums
                      ETA[[j]]$post_varSets2<-ETA[[j]]$post_varSets2*k+(ETA[[j]]$varSets^2)/nSums
                      if(ETA[[j]]$saveEffects&&(i%%ETA[[j]]$thin)==0){
                          writeBin(object=ETA[[j]]$b,con=ETA[[j]]$fileEffects,size=ifelse(ETA[[j]]$storageMode=="single",4,8))
                      }#*#
                    }
                    
                    if (ETA[[j]]$model == "BL") {
                      ETA[[j]]$post_b = ETA[[j]]$post_b * k + ETA[[j]]$b/nSums
                      ETA[[j]]$post_b2 = ETA[[j]]$post_b2 * k + (ETA[[j]]$b^2)/nSums
                      ETA[[j]]$post_tau2 = ETA[[j]]$post_tau2 * k + (ETA[[j]]$tau2)/nSums
                      ETA[[j]]$post_lambda = ETA[[j]]$post_lambda * k + (ETA[[j]]$lambda)/nSums
                      if(ETA[[j]]$saveEffects&&(i%%ETA[[j]]$thin)==0){
                          writeBin(object=ETA[[j]]$b,con=ETA[[j]]$fileEffects,size=ifelse(ETA[[j]]$storageMode=="single",4,8))
                      }#*#
                    }

                    if (ETA[[j]]$model == "RKHS") {
                      ETA[[j]]$post_varU = ETA[[j]]$post_varU * k + ETA[[j]]$varU/nSums
                      ETA[[j]]$post_varU2 = ETA[[j]]$post_varU2 * k + (ETA[[j]]$varU^2)/nSums
                      ETA[[j]]$post_uStar = ETA[[j]]$post_uStar * k + ETA[[j]]$uStar/nSums
                      ETA[[j]]$post_u = ETA[[j]]$post_u * k + ETA[[j]]$u/nSums
                      ETA[[j]]$post_u2 = ETA[[j]]$post_u2 * k + (ETA[[j]]$u^2)/nSums
                    }

                    if (ETA[[j]]$model == "BayesC") {
                      ETA[[j]]$post_b = ETA[[j]]$post_b * k + ETA[[j]]$b/nSums
                      ETA[[j]]$post_b2 = ETA[[j]]$post_b2 * k + (ETA[[j]]$b^2)/nSums
                      ETA[[j]]$post_varB = ETA[[j]]$post_varB * k + (ETA[[j]]$varB)/nSums
                      ETA[[j]]$post_varB2 = ETA[[j]]$post_varB2 * k + (ETA[[j]]$varB^2)/nSums
                      ETA[[j]]$post_d = ETA[[j]]$post_d * k + (ETA[[j]]$d)/nSums
                      ETA[[j]]$post_probIn = ETA[[j]]$post_probIn * k + (ETA[[j]]$probIn)/nSums
                      ETA[[j]]$post_probIn2 = ETA[[j]]$post_probIn2 * k + (ETA[[j]]$probIn^2)/nSums
                      if(ETA[[j]]$saveEffects&&(i%%ETA[[j]]$thin)==0){
                          writeBin(object=ETA[[j]]$b*ETA[[j]]$d,con=ETA[[j]]$fileEffects,size=ifelse(ETA[[j]]$storageMode=="single",4,8))
                      }#*#
                    }

                    if (ETA[[j]]$model == "BayesA") {
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

                    if(ETA[[j]]$model=="BayesB")
                    {
                        ETA[[j]]$post_b=ETA[[j]]$post_b*k+ETA[[j]]$b/nSums
                        ETA[[j]]$post_b2=ETA[[j]]$post_b2*k+(ETA[[j]]$b^2)/nSums
                        ETA[[j]]$post_varB=ETA[[j]]$post_varB*k+(ETA[[j]]$varB)/nSums
                        ETA[[j]]$post_varB2=ETA[[j]]$post_varB2*k+(ETA[[j]]$varB^2)/nSums
                        ETA[[j]]$post_d = ETA[[j]]$post_d * k + (ETA[[j]]$d)/nSums
                        ETA[[j]]$post_probIn = ETA[[j]]$post_probIn * k + (ETA[[j]]$probIn)/nSums
                        ETA[[j]]$post_probIn2 = ETA[[j]]$post_probIn2 * k + (ETA[[j]]$probIn^2)/nSums
                        ETA[[j]]$post_S = ETA[[j]]$post_S * k + (ETA[[j]]$S)/nSums
						ETA[[j]]$post_S2 = ETA[[j]]$post_S2 * k + (ETA[[j]]$S^2)/nSums
						if(ETA[[j]]$saveEffects&&(i%%ETA[[j]]$thin)==0){
                            writeBin(object=ETA[[j]]$b*ETA[[j]]$d,con=ETA[[j]]$fileEffects,size=ifelse(ETA[[j]]$storageMode=="single",4,8))
                        }#*#
                    }
                  }
                }

                post_mu = post_mu * k + mu/nSums
                post_mu2 = post_mu2 * k + (mu^2)/nSums

                post_yHat = post_yHat * k + yHat/nSums
                post_yHat2 = post_yHat2 * k + (yHat^2)/nSums

                post_varE = post_varE * k + varE/nSums
                post_varE2 = post_varE2 * k + (varE^2)/nSums

                if (response_type == "ordinal") {
                  post_threshold = post_threshold * k + threshold/nSums
                  post_threshold2 = post_threshold2 * k + (threshold^2)/nSums

                  TMP=matrix(nrow=n,ncol=nclass,0)

                  TMP[,1]=pnorm(threshold[2]-yHat)

                  if(nclass>2){
                     for(m in 2:(nclass-1)){
                       TMP[,m]=pnorm(threshold[(m+1)]-yHat)-rowSums(as.matrix(TMP[,1:(m-1)]))
                     }
                  }
                  TMP[,nclass]=1-rowSums(TMP)

                  post_prob=post_prob*k+TMP/nSums
                  post_prob2=post_prob2*k+(TMP^2)/nSums
                 
                  if(nNa==0){
                    logLik=loglik_ordinal(z,yHat,threshold)
		          }else{
		            logLik=loglik_ordinal(z[-whichNa],yHat[-whichNa],threshold)
		          }
                }

                if(response_type == "gaussian") {
                  
                  tmpE = e/weights
                  if(!is.null(groups)){
                    tmpSD=rep(NA,n)
                    for(g in 1:nGroups)
                    {
                       index=(groups==g)
                       tmpSD[index]=sqrt(varE[g])/weights[index]
                    }
                  }else{
                    tmpSD = sqrt(varE)/weights
		  		  }

                  if (nNa > 0) {
                    tmpE = tmpE[-whichNa]
                    tmpSD = tmpSD[-whichNa]
                  }
                  
	          logLik = sum(dnorm(tmpE, sd = tmpSD, log = TRUE))

                }#end gaussian

                post_logLik = post_logLik * k + logLik/nSums
            }
        }#end of saving samples and computing running means

        if (verbose) {
            message("---------------------------------------")
            tmp = proc.time()[3]
            message("Iter=",i," Time/Iter=",round(tmp-time,3))
            #message("VarE=",round(varE,3))
	    #In the case of variance by groups
	    message("varE=",paste(round(varE,3),collapse=", "))
            time = tmp
        }
    }#end of Gibbs sampler

    #Closing files
    close(fileOutVarE)
    close(fileOutMu)
    if(response_type == "ordinal") close(fileOutThresholds)

    if (nLT > 0) {
        for (i in 1:nLT) {
            if (!is.null(ETA[[i]]$fileOut)) {
                flush(ETA[[i]]$fileOut)                
                close(ETA[[i]]$fileOut)
                ETA[[i]]$fileOut = NULL
            }
            if(!is.null(ETA[[i]]$fileEffects)){
                flush(ETA[[i]]$fileEffects)
                close(ETA[[i]]$fileEffects)
                if(!is.null(ETA[[i]]$compressEffects)&&ETA[[i]]$compressEffects==TRUE){
                    compressFile(paste0(saveAt,ETA[[i]]$Name,"_b.bin"))
                }
                ETA[[i]]$fileEffects = NULL
            }
            
        }
    }
    
    #return goodies

    out = list(y = y0, a=a,b=b,whichNa = whichNa, saveAt = saveAt, nIter = nIter, 
               burnIn = burnIn, thin = thin, 
               weights = weights, verbose = verbose, 
               response_type = response_type, df0 = df0, S0 = S0)

    out$yHat = post_yHat

    names(out$yHat)=IDs
    names(out$y)=IDs

    out$SD.yHat = sqrt(post_yHat2 - (post_yHat^2))
    out$mu = post_mu
    out$SD.mu = sqrt(post_mu2 - post_mu^2)
    out$varE = post_varE
    out$SD.varE = sqrt(post_varE2 - post_varE^2)
    
    #goodness of fit 
    out$fit = list()
    
    if(response_type=="gaussian")
    {
    	tmpE = (yStar - post_yHat)/weights

        if(!is.null(groups))
        {
                    tmpSD=rep(NA,n)
                    for(g in 1:nGroups)
                    {
                       index=(groups==g)
                       tmpSD[index]=sqrt(varE[g])/weights[index]
                    }
         }else{
    		tmpSD = sqrt(post_varE)/weights
	 }
    
    	if (nNa > 0) {
        	tmpE = tmpE[-whichNa]
        	tmpSD = tmpSD[-whichNa]
    	}
    	out$fit$logLikAtPostMean = sum(dnorm(tmpE, sd = tmpSD, log = TRUE))

	if (Censored) {
            cdfA = pnorm(q = a[whichNa], sd = sqrt(post_varE), mean = post_yHat[whichNa])
            cdfB = pnorm(q = b[whichNa], sd = sqrt(post_varE), mean = post_yHat[whichNa])
            out$fit$logLikAtPostMean = out$fit$logLikAtPostMean + sum(log(cdfB - cdfA))
        }
    }

    if(response_type=="ordinal")
    {
         out$probs=post_prob
         out$SD.probs=sqrt(post_prob2-post_prob^2)
         colnames(out$probs)=lev
         colnames(out$SD.probs)=lev
         out$threshold = post_threshold[-c(1, nclass + 1)]
         out$SD.threshold = sqrt(post_threshold2 - post_threshold^2)[-c(1, nclass + 1)] 
         
         #out$fit$logLikAtPostMean = loglik_ordinal(y,post_yHat,post_threshold)#*#
        tmp=0
         for(i in 1:nclass){
            tmp=tmp+sum(ifelse(y0==lev[i],log(out$probs[,i]),0))
         }
         
         out$fit$logLikAtPostMean=tmp
         out$levels=lev
         out$nlevels=nclass
    }

    out$fit$postMeanLogLik = post_logLik
    out$fit$pD = -2 * (post_logLik - out$fit$logLikAtPostMean)
    out$fit$DIC = out$fit$pD - 2 * post_logLik

    # Renaming/removing objects in ETA and appending names
    if (nLT > 0) {
        for (i in 1:nLT) {

            if (ETA[[i]]$model != "RKHS") {
                ETA[[i]]$b = ETA[[i]]$post_b
                ETA[[i]]$SD.b = sqrt(ETA[[i]]$post_b2 - ETA[[i]]$post_b^2)
                names(ETA[[i]]$b)=ETA[[i]]$colNames
                names(ETA[[i]]$SD.b)=ETA[[i]]$colNames
                tmp = which(names(ETA[[i]]) %in% c("post_b", "post_b2","X","x2"))
                ETA[[i]] = ETA[[i]][-tmp]
            }
            
            if(ETA[[i]]$model=="RKHS")
            {
               ETA[[i]]$SD.u=sqrt(ETA[[i]]$post_u2 - ETA[[i]]$post_u^2)
               ETA[[i]]$u=ETA[[i]]$post_u
               ETA[[i]]$uStar=ETA[[i]]$post_uStar
               ETA[[i]]$varU=ETA[[i]]$post_varU
               ETA[[i]]$SD.varU=sqrt(ETA[[i]]$post_varU2 - ETA[[i]]$post_varU^2)
               tmp=which(names(ETA[[i]])%in%c("post_varU","post_varU2","post_uStar","post_u","post_u2"))
               ETA[[i]]=ETA[[i]][-tmp]
            }

            if (ETA[[i]]$model %in% c("BRR","BRR_sets", "BayesA", "BayesC","BayesB")) {
                ETA[[i]]$varB = ETA[[i]]$post_varB
                ETA[[i]]$SD.varB = sqrt(ETA[[i]]$post_varB2 - (ETA[[i]]$post_varB^2))
                tmp = which(names(ETA[[i]]) %in% c("post_varB", "post_varB2"))
                ETA[[i]] = ETA[[i]][-tmp]
            }
           	if(ETA[[i]]$model=="BRR_sets"){
				ETA[[i]]$varSets=ETA[[i]]$post_varSets
				ETA[[i]]$SD.varSets=sqrt(ETA[[i]]$post_varSets2-(ETA[[i]]$post_varSets^2))
				tmp<-which(names(ETA[[i]])%in%c("post_varSets","post_varSets2"))
				ETA[[i]]=ETA[[i]][-tmp]
			}

            if(ETA[[i]]$model %in% c("BayesB","BayesC"))
            {
	        	ETA[[i]]$d=ETA[[i]]$post_d
                ETA[[i]]$probIn=ETA[[i]]$post_probIn
				ETA[[i]]$SD.probIn=sqrt(ETA[[i]]$post_probIn2 - (ETA[[i]]$post_probIn^2))
                tmp = which(names(ETA[[i]]) %in% c("post_d", "post_probIn","post_probIn2"))
                ETA[[i]] = ETA[[i]][-tmp]
            }
            
            if(ETA[[i]]$model %in% c("BayesA","BayesB"))
            {
                ETA[[i]]$S=ETA[[i]]$post_S
                ETA[[i]]$SD.S=sqrt( ETA[[i]]$post_S2 - (ETA[[i]]$post_S^2))
                tmp=which(names(ETA[[i]])%in%c("post_S","post_S2"))
                ETA[[i]]=ETA[[i]][-tmp]
            }

            if(ETA[[i]]$model=="BL")
            {
                ETA[[i]]$tau2=ETA[[i]]$post_tau2
                ETA[[i]]$lambda=ETA[[i]]$post_lambda
                tmp = which(names(ETA[[i]]) %in% c("post_tau2", "post_lambda","lambda2"))
                ETA[[i]] = ETA[[i]][-tmp]
            }
        }
        out$ETA = ETA
    }
    class(out) = "BGLR"
    return(out)
}

#This function will be a wrapper for BGLR
#the idea is to maintain the compatibility with the function BLR in 
#the package BLR that was released in 2010, updated in 2011 and 2012

#NOTE: thin2 parameter is missing in BGLR, so it will be removed

BLR=function (y, XF = NULL, XR = NULL, XL = NULL, GF = list(ID = NULL, 
    A = NULL), prior = NULL, nIter = 1100, burnIn = 100, thin = 10, 
    thin2 = 1e+10, saveAt = "", minAbsBeta = 1e-09, weights = NULL) 
{

    ETA = NULL
    ETA = list()
    nLT = 0

    message("This implementation is a simplified interface for the more general")
    message("function BGLR, we keep it for backward compatibility with our package BLR")

    warning("thin2 parameter is not used any more and will be deleted in next releases\n",immediate. = TRUE);
   
    message("Setting parameters for BGLR...")
    if (is.null(prior)) {
        message("===============================================================")
        message("No prior was provided, BGLR will be running with improper priors.")
        message("===============================================================")
        prior = list(varE = list(S = NULL, df = 1), varBR = list(S = 0, 
            df = 0), varU = list(S = 0, df = 0), lambda = list(shape = 0, 
            rate = 0, type = "random", value = 50))
    }
    if (!is.null(XF)) {
        nLT = nLT + 1
        ETA[[nLT]] = list(X = XF, model = "FIXED")
    }
    if (!is.null(XR)) {
        nLT = nLT + 1
        ETA[[nLT]] = list(X = XR, model = "BRR", df0 = prior$varBR$df, 
            S0 = prior$varBR$S)
    }
    if (!is.null(XL)) {
        nLT = nLT + 1
        if (prior$lambda$type == "random") {
            if (is.null(prior$lambda$rate)) {
                message("Setting prior for lambda^2 to beta")
                prior$lambda$type = "beta"
                prior$lambda$shape = NULL
                prior$lambda$rate = NULL
            }
            else {
                message("Setting prior for lambda^2 to gamma")
                prior$lambda$type = "gamma"
                prior$lambda$max = NULL
                prior$lambda$shape1 = NULL
                prior$lambda$shape2 = NULL
            }
        }
        ETA[[nLT]] = list(X = XL, model = "BL", type = prior$lambda$type, 
                          rate = prior$lambda$rate, shape = prior$lambda$shape, 
                          max = prior$lambda$max, shape1 = prior$lambda$shape1, 
                          shape2 = prior$lambda$shape2, lambda = prior$lambda$value,
                          minAbsBeta=minAbsBeta)
    }

    #NOTE: In original BLR IDS are used to buid A matrix Z,
    #and then run the model y=Zu+e, u~MN(0,varU*A), and then using the Cholesky factorization
    #it was possible to fit the model. The algorithm used here is different (Orthogonal variables)
    #and it may be that the IDS are not longer necessary

    if (!is.null(GF[[1]])) {
        nLT = nLT + 1
        ETA[[nLT]] = list(K = GF$A, model = "RKHS", df0 = prior$varU$df, 
            S0 = prior$varU$S)
        warning("IDs are not used any more and will be deleted in next releases...\n",immediate. = TRUE) 
    }

    message("Finish setting parameters for BGLR")
    message("Fitting model using BGLR...")
    out = BGLR(y = y, ETA = ETA, df0 = prior$varE$df, S0 = prior$varE$S, 
               nIter = nIter, burnIn = burnIn, thin = thin, saveAt = saveAt, 
               weights = weights)

    #Backward compatibility with BLR
    if (nLT > 0) {
        for (j in 1:nLT) {
            if (ETA[[j]]$model == "FIXED") {
                out$bF = out$ETA[[j]]$b
                out$SD.bF = out$ETA[[j]]$SD.b
            }
            if (ETA[[j]]$model == "BL") {
                out$bL = out$ETA[[j]]$b
                out$SD.bL = out$ETA[[j]]$SD.b
                out$lambda = out$ETA[[j]]$lambda
            }
            if (ETA[[j]]$model == "BRR") {
                out$bR = out$ETA[[j]]$b
                out$SD.bR = out$ETA[[j]]$SD.b
                out$varBR = out$ETA[[j]]$varB
                out$SD.bR = out$ETA[[j]]$SD.varB
            }
            if (ETA[[j]]$model == "RKHS") {
                out$u = out$ETA[[j]]$u
                out$SD.u =out$ETA[[j]]$SD.u
                out$varU = out$ETA[[j]]$varU
            }
        }
    }
    out$ETA = NULL
    class(out) = "BLR"
    return(out)
}

