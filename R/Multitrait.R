
#Auxiliary functions
#Converts a logical vector to decimal representation
#eg x=c(TRUE,FALSE), means that we have a binary number "10" which in decimal is 2

logicalToDec<-function(x)
{
	sum(x * 2^(rev(seq_along(x)) - 1))
}

#generating pseudo-random sequences from Wishart distribution
#v Degrees of freedom (scalar).
#S Inverse scale matrix pxp.
#Function included in the MCMCpack, version 1.4-5

rwish<-function(v, S) 
{
    if (!is.matrix(S))
      S <- matrix(S)
    if (nrow(S) != ncol(S)) 
    {
      stop(message="S not square in rwish().\n")
    }
    if (v < nrow(S)) 
    {
      stop(message="v is less than the dimension of S in rwish().\n")
    }
    p <- nrow(S)
    CC <- chol(S)
    Z <- matrix(0, p, p)
    diag(Z) <- sqrt(rchisq(p, v:(v-p+1)))
    if(p > 1) 
    {
      pseq <- 1:(p-1)
      Z[rep(p*pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p*(p-1)/2)
    }
    return(crossprod(Z %*% CC))
}

#generating pseudo-random sequences from inverse Wishart distribution
#v Degrees of freedom (scalar).
#S Inverse scale matrix p x p
#Function included in the MCMCpack, version 1.4-5
riwish<-function(v, S) 
{
	return(solve(rwish(v,solve(S))))
}

#Extract Lower Triangular Elements from a Symmetric Matrix.
#The elements are stored in column major order, it does not check for symmetry. 
#Function included in the MCMCpack, version 1.4-5
vech<-function (x) 
{
    x <- as.matrix(x)
    
    if (dim(x)[1] != dim(x)[2]) 
    {
      stop("Non-square matrix passed to vech().\n")
    }
    
    output <- x[lower.tri(x, diag = TRUE)]
    dim(output) <- NULL
    return(output)
}

#Takes vector x and returns an nrow times nrow symmetric matrix,
#this will recycle the elements of x as needed to fill the matrix.
#The number of rows is computed autmatically is not given.
#Function included in the MCMCpack, version 1.4-5
xpnd <-function (x, nrow = NULL) 
{
    dim(x) <- NULL
    if(is.null(nrow)) nrow <- (-1 + sqrt(1 + 8 * length(x))) / 2
    output <- matrix(0, nrow, nrow)
    output[lower.tri(output, diag = TRUE)] <- x
    hold <- output
    hold[upper.tri(hold, diag=TRUE)] <- 0
    output <- output + t(hold)
    return(output)
}
  
#UNstructured covariance matrix
#Cov, list
#traits, integer
#mo, Mode
setCov.UN<-function(Cov,traits,j,mo,saveAt)
{	
	
	message("UNstructured covariance matrix")
	
	if(is.null(Cov$df0))
	{
		Cov$df0<-traits+1
		message("df0 was set to ",Cov$df0)
	}
	
	if(is.null(Cov$S0))
	{
		#Cov$S0<-diag(traits)
		#message("S0 set to an identity matrix")
		Cov$S0<-mo*(Cov$df0+traits+1)
		message("S0 set to ")
		print(Cov$S0)	
	}
		
	#Omega=cov(b)
	Cov$Omega<-riwish(v=Cov$df0,S=Cov$S0)
	Cov$Omegainv<-solve(Cov$Omega)
	
	#Objects for saving posterior means for MCMC
	Cov$post_Omega<-matrix(0,nrow=traits,ncol=traits)
	Cov$post_Omega2<-matrix(0,nrow=traits,ncol=traits)
	
	#Output files
	Cov$fName_Omega<-paste(saveAt,"Omega_",j,".dat",sep="")
	Cov$f_Omega<-file(description=Cov$fName_Omega,open="w")
	
	return(Cov)

}

#DIAGonal covariance matrix
#Cov, list
#traits, integer
#mo, Mode
setCov.DIAG<-function(Cov,traits,j,mo,saveAt)
{
	message("DIAGonal covariance matrix")
	
	if(is.null(Cov$df0))
	{
		Cov$df0<-rep(traits+1,traits)
		message("df0 set to  ",traits+1," for all the traits")
	}
	
	if(is.null(Cov$S0))
	{
		#Cov$S0 <- rep(1,traits)
		#message("S0 was set to 1 for all the traits")
		Cov$S0<-mo*(Cov$df0+2)
		message("S0 was set to: ")
		print(Cov$S0)	
	}
	
	Cov$Omega<-matrix(0,nrow=traits,ncol=traits)
		
	for(k in 1:traits)
	{
		Cov$Omega[k,k]<-Cov$S0[k]/rchisq(n=1,df=Cov$df0[k])
	}
		
	Cov$Omegainv<-solve(Cov$Omega)
	
	#Objects for saving posterior means for MCMC
	Cov$post_Omega<-matrix(0,nrow=traits,ncol=traits)
	Cov$post_Omega2<-matrix(0,nrow=traits,ncol=traits)
	
	#Output files
	Cov$fName_Omega<-paste(saveAt,"Omega_",j,".dat",sep="")
	Cov$f_Omega<-file(description=Cov$fName_Omega,open="w")
	
	return(Cov)	
	
}

#RECursive covariance matrix
setCov.REC<-function(Cov,traits,j,mo,saveAt)
{
	message("RECursive covariance matrix")
	
	if(is.null(Cov$df0))
	{
		Cov$df0<-rep(traits+1,traits)
		message("df0 set to ", traits+1, " for all the traits")
	}
	
	if(is.null(Cov$S0))
	{
		#Cov$S0 <- rep(1,traits)
		#message("S0 was set to 1 for all the traits")
		Cov$S0<-mo*(Cov$df0+2)
		message("S0 was set to: ")
		print(Cov$S0)	
	}
			
	if(is.null(Cov$var))
	{
		Cov$var<-100
		message("var was set to 100")
	}

	if(is.null(Cov$M)) stop("M can not be null")
		
	if(!is.logical(Cov$M)) stop("M must be logical matrix (with entries being TRUE/FALSE)")
		
	if(!is.matrix(Cov$M)) stop("M must be a matrix")
		
	if(nrow(Cov$M)!=ncol(Cov$M)) stop("M must be a square matrix") 
		
	if(nrow(Cov$M)!=traits) stop("M must have ", traits, " rows and columns")
		
	if(any(diag(Cov$M)==TRUE)) stop("All diagonal entries of M must be set to FALSE")
		
	Cov$M[upper.tri(Cov$M)]<-FALSE
		

	Cov$W<-matrix(0,nrow=traits,ncol=traits)
	Cov$PSI<-rep(NA,traits)
		
	for(k in 1:traits)
	{
		Cov$PSI[k]<-Cov$S0[k]/rchisq(n=1,df=Cov$df0[k])
	}
		
	#Omega=cov(b)
	Cov$Omega<-riwish(v=traits,S=diag(Cov$S0))
	Cov$Omegainv<-solve(Cov$Omega)
		
	#Objects for saving posterior means for MCMC
	Cov$post_Omega<-matrix(0,nrow=traits,ncol=traits)
	Cov$post_Omega2<-matrix(0,nrow=traits,ncol=traits)
	Cov$post_W<-matrix(0,nrow=traits,ncol=traits)
	Cov$post_W2<-matrix(0,nrow=traits,ncol=traits)
	Cov$post_PSI<-rep(0,traits)
	Cov$post_PSI2<-rep(0,traits)
	
	#Output files
	Cov$fName_W<-paste(saveAt,"W_",j,".dat",sep="")
	Cov$fName_PSI<-paste(saveAt,"PSI_",j,".dat",sep="")
	
	Cov$f_W<-file(description=Cov$fName_W,open="w")
	Cov$f_PSI<-file(description=Cov$fName_PSI,open="w")
	
	return(Cov)
}


#FA covariance matrix
setCov.FA<-function(Cov,traits,nD,j,mo,saveAt)
{
	message("FA covariance matrix")
	
	if(is.null(Cov$df0))
	{
		Cov$df0<-rep(traits+1,traits)
		message("df0 set to ",traits+1," for all the traits")
	}

	
	if(is.null(Cov$S0))
	{
		#Cov$S0 <- rep(1/100,traits)
		#message("S0 was set to 1/100 for all the traits")
		Cov$S0<-mo*(Cov$df0+2)
		message("S0 was set to: ")
		print(Cov$S0)
	}
			
	if(is.null(Cov$var))
	{
		Cov$var<-100
		message("var was set to 100")
	}
			
	
	if(is.null(Cov$M)) stop("M can not be null")
		
	if(!is.logical(Cov$M)) stop("M must be logical matrix (with entries being TRUE/FALSE)")
		
	if(!is.matrix(Cov$M)) stop("M must be a matrix")
		
	if(nrow(Cov$M)!=traits) stop("M must have ", traits, " rows")
		
	if(ncol(Cov$M)>traits) stop("Number of columns of M must be smaller than ", traits)
		
	Cov$nF<-ncol(Cov$M)
	
	Cov$nD<-nD
		
	#Omega=cov(b)
	Cov$Omega<-riwish(v=traits,S=diag(Cov$S0))
	
	sdU <- sqrt(diag(Cov$Omega))
    FA <- factanal(covmat = Cov$Omega, factors = Cov$nF)
    Cov$W <- matrix(nrow = traits, ncol = Cov$nF, 0)
    Cov$W[Cov$M] <- (diag(sdU) %*% FA$loadings)[Cov$M]
    Cov$PSI <- (sdU^2) * FA$uniquenesses + 1e-04
    Cov$Omega <- tcrossprod(Cov$W) + diag(Cov$PSI)
    Cov$Omegainv<-solve(Cov$Omega)
    
    Cov$F <- matrix(nrow = nD, ncol = Cov$nF, 0)
        
    #Objects for saving posterior means for MCMC
	Cov$post_Omega<-matrix(0,nrow=traits,ncol=traits)
	Cov$post_Omega2<-matrix(0,nrow=traits,ncol=traits)
	Cov$post_W<-matrix(0,nrow=traits,ncol=Cov$nF)
	Cov$post_W2<-matrix(0,nrow=traits,ncol=Cov$nF)
	Cov$post_PSI<-rep(0,traits)
	Cov$post_PSI2<-rep(0,traits)
	
	#Output files
	Cov$fName_W<-paste(saveAt,"W_",j,".dat",sep="")
	Cov$fName_PSI<-paste(saveAt,"PSI_",j,".dat",sep="")
	Cov$f_W<-file(description=Cov$fName_W,open="w")
	Cov$f_PSI<-file(description=Cov$fName_PSI,open="w")
		
	return(Cov)
}


#Set linear term for DiracSS_mt
setLT.DiracSS_mt<-function(LT,traits,j,Sy,nLT,R2,saveAt,nRow)
{
	
	message("Setting linear term ",j)
	
	#Inclusion probabilities
	
	if(is.null(LT$inclusionProb))
	{
		LT$inclusionProb<-list(probIn=rep(0.5,traits),
		                       counts=rep(2,traits))
		message("probIn set to 0.5 for all the traits")
		message("counts set to 2 for all the traits")
	
	}else{
		
		if(is.null(LT$inclusionProb$probIn))
		{
			LT$inclusionProb$probIn<-rep(0.5,traits)
			message("probIn set to 0.5 for all the traits")
		}
		
		if(is.null(LT$inclusionProb$counts))
		{
			LT$inclusionProb$probIn$counts<-rep(2,traits)
			message("counts set to 2 for all the traits")
		}
	}

	#Compute countsIn and countsOut
	LT$inclusionProb$countsIn <- LT$inclusionProb$counts * LT$inclusionProb$probIn
	LT$inclusionProb$countsOut <- LT$inclusionProb$counts - LT$inclusionProb$countsIn
	
	#X
	if(is.null(LT$X)) stop("X can not be NULL\n")
	if(!is.matrix(LT$X))  stop("X must be a matrix\n")
	if(any(is.na(LT$X))) stop("X has NAs\n")
	
	LT$x2<-as.vector(colSums(LT$X^2))
	
	sumMeanXSq<-sum((apply(LT$X,2L,mean))^2)
	MSx<-sum(LT$x2)/nrow(LT$X)-sumMeanXSq
	message("MSx=",MSx)
	
	#Initialize b, d, beta, beta=b#d
	LT$p<-ncol(LT$X)
	
	LT$b<-matrix(0,nrow=LT$p,ncol=traits)
	LT$d<-matrix(1,nrow=LT$p,ncol=traits)

	LT$beta<-LT$b*LT$d
	
	#Distribution for Omega, which is related to (b_1j,b_2j,...,b_tj)
	#j=1,...,p, where p is the number of columns of X
	#t the number of traits
	
	if(is.null(LT$Cov))
	{
		#Cov is null
		LT$Cov<-list()
		LT$Cov$type<-"UN"
		
	}else{
		#Cov is not null
		if(is.null(LT$Cov$type))
		{
		
			LT$Cov$type<-"UN"
			
		}else{
		
			if(!(LT$Cov$type %in% c("UN","REC","FA","DIAG")))
			{
				stop("Error '", LT$Cov$type, "' not implemented (note: evaluation is case sensitive)")
			}
			 
		}
	}
	
	#Beyond this point Cov is not NULL and we already know the covariance structure

	#Select appropriate covariance structure
	
	LT$Cov<-switch(LT$Cov$type,
					UN=setCov.UN(Cov=LT$Cov,traits=traits,j=j,mo=(R2/nLT)*Sy/MSx,saveAt=saveAt),
					DIAG=setCov.DIAG(Cov=LT$Cov,traits=traits,j=j,mo=(R2/nLT)*diag(Sy)/MSx,saveAt),
					REC=setCov.REC(Cov=LT$Cov,traits=traits,j=j,mo=(R2/nLT)*diag(Sy)/MSx,saveAt=saveAt),
					FA=setCov.FA(Cov=LT$Cov,traits=traits,nD=LT$p,j=j,mo=(R2/nLT)*diag(Sy)/MSx,saveAt=saveAt)
	               )
	
	
	#It is not working very well when probIn is small
	#Dm<-diag(1/LT$inclusionProb$probIn)
	
	#LT$Cov<-switch(LT$Cov$type,
	#				UN=setCov.UN(Cov=LT$Cov,traits=traits,mo=(R2/nLT)*Dm%*%Sy%*%Dm/MSx),
	#				DIAG=setCov.DIAG(Cov=LT$Cov,traits=traits,mo=(R2/nLT)*Dm%*%diag(Sy)%*%Dm/MSx),
	#				REC=setCov.REC(Cov=LT$Cov,traits=traits,j=j,mo=(R2/nLT)*Dm%*%diag(Sy)%*%Dm/MSx,saveAt=saveAt),
	#				FA=setCov.FA(Cov=LT$Cov,traits=traits,nD=LT$p,j=j,mo=(R2/nLT)*Dm%*%diag(Sy)%*%Dm/MSx,saveAt=saveAt)
	#               )
	
	
	
	#Add a new object to compute covariance between entries of
	#beta=b*d with MCMC output, Sigma=Cov(beta,beta'), beta is 
	#a vector with marker effects for one locus, dimmension 1*traits
	LT$Cov$Sigma=matrix(0,nrow=traits,ncol=traits)   
	LT$Cov$post_Sigma=matrix(0,nrow=traits,ncol=traits)
	LT$Cov$post_Sigma2=matrix(0,nrow=traits,ncol=traits)            
	
	#Objects for saving posterior means for MCMC
	LT$post_b<-matrix(0,nrow=LT$p,ncol=traits)
	LT$post_b2<-matrix(0,nrow=LT$p,ncol=traits)
	
	LT$post_d<-matrix(0,nrow=LT$p,ncol=traits)
	LT$post_d2<-matrix(0,nrow=LT$p,ncol=traits)
	
	LT$post_beta<-matrix(0,nrow=LT$p,ncol=traits)
	LT$post_beta2<-matrix(0,nrow=LT$p,ncol=traits)
	
	#Files to save binary files with betas
	if(is.null(LT$saveEffects))
	{
		LT$saveEffects<-FALSE
	}
	
    if(LT$saveEffects)
    {
        if(is.null(LT$storageMode))
        {
        	LT$storageMode<-"double"
        }
        
        if(!LT$storageMode%in%c("single","double")) 
        {
            stop("storageMode of LP ",j," can either be 'single' or 'double' (default)")
        }
         	
    	fname<-paste(saveAt,LT$Name,"_beta.bin",sep="")
    	LT$fileEffects<-file(fname,open='wb')
    	writeBin(object=c(nRow,traits,LT$p),con=LT$fileEffects,size=ifelse(LT$storageMode=="single",4,8))
    	
    }#*#
    
    #NEW, 
    #perhaps we need to remove, this is only useful for computing genomic relationship
    #matrix when using Cheng (2018) method.
    #Files to save indicator variables
    #binary file are saved in single mode (4 bytes for each number)
    if(is.null(LT$saveIndicators))
    {
    	LT$saveIndicators<-FALSE
    }
    
    if(LT$saveIndicators)
    {

    	fname2<-paste(saveAt,LT$Name,"_d.bin",sep="")
    	LT$fileIndicators<-file(fname2,open='wb')
    	#nrow, traits and p stored as single
    	writeBin(object=c(nRow,traits,LT$p),con=LT$fileIndicators,size=4)
    }
	
	return(LT)
}

#Set linear term for Ridge Regression
setLT.BRR_mt<-function(LT,traits,j,Sy,nLT,R2,saveAt,nRow)
{		

	message("Setting linear term ",j)
		
	#X
	if(is.null(LT$X)) stop("X can not be NULL\n")
	if(!is.matrix(LT$X))  stop("X must be a matrix\n")
	if(any(is.na(LT$X))) stop("X has NAs\n")
	
	LT$x2<-as.vector(colSums(LT$X^2))
	
	sumMeanXSq<-sum((apply(LT$X,2L,mean))^2)
	MSx<-sum(LT$x2)/nrow(LT$X)-sumMeanXSq
	message("MSx=",MSx)
	
	#Initialize beta
	LT$p<-ncol(LT$X)
	
	LT$beta<-matrix(0,nrow=LT$p,ncol=traits)
	
	#Distribution for Omega, which is related to (b_1j,b_2j,...,b_tj)
	#j=1,...,p, where p is the number of columns of X
	#t the number of traits
	
	if(is.null(LT$Cov))
	{
		#Cov is null
		LT$Cov<-list()
		LT$Cov$type<-"UN"
		
	}else{
		#Cov is not null
		if(is.null(LT$Cov$type))
		{
		
			LT$Cov$type<-"UN"
			
		}else{
		
			if(!(LT$Cov$type %in% c("UN","REC","FA","DIAG")))
			{
				stop("Error '", LT$Cov$type, "' not implemented (note: evaluation is case sensitive)")
			}
			 
		}
	}
	
	#Beyond this point Cov is not NULL and we already know the covariance structure

	#Select appropriate covariance structure
	
	LT$Cov<-switch(LT$Cov$type,
					UN=setCov.UN(Cov=LT$Cov,traits=traits,j=j,mo=(R2/nLT)*Sy/MSx,saveAt=saveAt),
					DIAG=setCov.DIAG(Cov=LT$Cov,traits=traits,j=j,mo=(R2/nLT)*diag(Sy)/MSx,saveAt=saveAt),
					REC=setCov.REC(Cov=LT$Cov,traits=traits,j=j,mo=(R2/nLT)*diag(Sy)/MSx,saveAt=saveAt),
					FA=setCov.FA(Cov=LT$Cov,traits=traits,nD=LT$p,j=j,mo=(R2/nLT)*diag(Sy)/MSx,saveAt=saveAt)
	               )
	
	#Objects for saving posterior means for MCMC
	LT$post_beta<-matrix(0,nrow=LT$p,ncol=traits)
	LT$post_beta2<-matrix(0,nrow=LT$p,ncol=traits)
		
	#Files to save binary files with betas
	if(is.null(LT$saveEffects))
	{
		LT$saveEffects<-FALSE
	}
	
    if(LT$saveEffects)
    {
        if(is.null(LT$storageMode))
        {
        	LT$storageMode<-"double"
        }
        
        if(!LT$storageMode%in%c("single","double")) 
        {
            stop("storageMode of LP ",j," can either be 'single' or 'double' (default)")
        }
         	
    	fname<-paste(saveAt,LT$Name,"_beta.bin",sep="")
    	LT$fileEffects<-file(fname,open='wb')
    	writeBin(object=c(nRow,traits,LT$p),con=LT$fileEffects,size=ifelse(LT$storageMode=="single",4,8))
    	
    }#*#
    
	return(LT)
}


#Set linear term for u_t ~ N(0, \sigma^2_r K)
#Note also  (u_1',...,u_t') ~ N(0, G_0 x K), where x represents the Kronecker product
#Internally we represent this using the eigen-value decomposition, 
#using as incidence matrix X=Gamma*Lambda^{1/2}, where K=Gamma*Lambda*Gamma' 
setLT.RKHS_mt<-function(LT,traits,j,Sy,nLT,R2,saveAt)
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
	
	
	#We do not save effects in RKHS
	LT$saveEffects<-FALSE
	LT<-setLT.BRR_mt(LT=LT,traits=traits,j=j,Sy=Sy,nLT=nLT,R2=R2,saveAt=saveAt,nRow=0)
	
	return(LT)
	
}

#Set linear term for Fixed effects
#Modified by Gustavo to support saving fixed effects, April 13, 2022
setLT.FIXED_mt<-function(LT,traits,j,saveAt,nRow)
{

        message("Setting linear term ",j)

        if(is.null(LT$common))
        {
                LT$common<-TRUE
                message("matrix for fixed effects X is the same for all the traits,")
                message("so the same effects are assumed for all the traits")
        }else{
                if(LT$common)
                {
                        message("matrix for fixed effects X is the same for all the traits,")
                        message("so the same effects are assumed for all the traits")
                }else{
                        message("each trait has its own set of predictors, we assumme")
                        message("X=[X_1,...,X_t], k=1,...,t (traits)")
                        
                    	if(is.null(LT$idColumns))
                    	{
                    		stop("You need to provide the vector idColumns \n which indicates which columns in X affect each trait")
                    	}
                }
        }

        #X
        if(is.null(LT$X)) stop("X can not be NULL\n")
        if(!is.matrix(LT$X))  stop("X must be a matrix\n")
        if(any(is.na(LT$X))) stop("X has NAs\n")

        if(LT$common)
        {
       		 	#Omega
        		LT$Cov<-list()
        		LT$Cov$Omega<-diag(rep(1E6,traits))
        		LT$Cov$Omegainv<-solve(LT$Cov$Omega)
        		
                #check rank
                if(qr(LT$X)$rank<ncol(LT$X)) stop("X is rank deficient")

                LT$x2<-as.vector(colSums(LT$X^2))


                #Initialize beta
                LT$p<-ncol(LT$X)

                LT$beta<-matrix(0,nrow=LT$p,ncol=traits)

                #Objects for saving posterior means for MCMC
                LT$post_beta<-matrix(0,nrow=LT$p,ncol=traits)
                LT$post_beta2<-matrix(0,nrow=LT$p,ncol=traits)

        }else{

				if(ncol(LT$X)!=length(LT$idColumns))
				{
					stop("Number of columns in X does not match with length of idColumns")
				}
				
				#Check the rank of each matrix
				for(k in 1:traits)
				{
						Ck<-which(LT$idColumns==k)
						if(length(Ck)>0)
						{
							if(qr(LT$X[,Ck])$rank<length(Ck)) stop("X_",k, " is rank deficient")	
						}else{
							message("covariates for fixed effects for trait ",k, " were not provided")
						}
				}
				
		        #Do not move this code out here!!!
                #It appears to be repeated but is not the case

                LT$x2<-as.vector(colSums(LT$X^2))

				#Stack the betas for each trait in one single vector, it is not possible
				#to use a matrix because each matrix in [X_1,...,X_t], k=1,...,t (traits)
				#can have different number of columns 
				
                LT$beta<-rep(0,ncol(LT$X))

                #Objects for saving posterior means for MCMC
                LT$post_beta<-rep(0,ncol(LT$X))
                LT$post_beta2<-rep(0,ncol(LT$X))
        }

     	#*#
     	
     	#Files to save binary files with betas     	
		if(is.null(LT$saveEffects))
		{
			LT$saveEffects<-FALSE
		} 
     	       
        if(LT$saveEffects)
        {
        	if(is.null(LT$storageMode))
        	{
        		LT$storageMode<-"double"
        	}

        	if(!LT$storageMode%in%c("single","double"))
        	{
            	stop("storageMode of LP ",j," can either be 'single' or 'double' (default)")
        	}

        	fname<-paste(saveAt,LT$Name,"_beta.bin",sep="")
        	LT$fileEffects<-file(fname,open='wb')
        	
        	if(LT$common)
        	{
        		writeBin(object=c(nRow,traits,LT$p),con=LT$fileEffects,size=ifelse(LT$storageMode=="single",4,8))
        	}else{
        		writeBin(object=c(nRow,ncol(LT$X)),con=LT$fileEffects,size=ifelse(LT$storageMode=="single",4,8))
        	}

    	}#*#

		return(LT)	
}

#Initialize residual covariance structure
setResCov<-function(resCov,traits,error,Sy,R2,saveAt)
{

	message("Initializing resCov")
	
	resCov$R<-var(error)/2
	resCov$Rinv<-solve(resCov$R)
	
	if(is.null(resCov$type))
	{
		resCov$type<-"UN"
		message("Modelling R as UNstructured")
		
	}else{
	
		if(!(resCov$type %in% c("UN","DIAG","REC","FA")))
		{
			stop("Error '", resCov$type, "' not implemented (note: evaluation is case sensitive)")
		} 
	}
	
	if(resCov$type=="UN")
	{
		message("Setting hyperparameters for UNstructured R")
		
		if(is.null(resCov$df0))
		{
			resCov$df0<-5
			message("df0 was set to ", resCov$df0)
		}
		
		if(is.null(resCov$S0))
		{
			#resCov$S0 <- diag(traits)
			resCov$S0<-(1-R2)*Sy*(resCov$df0+traits+1)
			message("S0 was set to ")
			print(resCov$S0)	
		}	
	}
	
	if(resCov$type=="DIAG")
	{
		message("Setting hyperparameters for DIAG R")
		
		if(is.null(resCov$df0))
		{
			resCov$df0<-rep(5,traits)
			message("df0 set to  5 for all the traits")
		}
		
		if(is.null(resCov$S0))
		{
			#resCov$S0 <- rep(1,traits)
			#message("S0 was set to 1 for all the traits")
			resCov$S0<-(1-R2)*diag(Sy)*(resCov$df0+2)
			message("S0 was set to ")
			print(resCov$S0)	
		}
		

	}
	
	if(resCov$type=="REC")
	{
		message("Setting hyperparameters for REC R")
		
		if(is.null(resCov$M)) stop("M can not be null")
		
		if(!is.logical(resCov$M)) stop("M must be logical matrix (with entries being TRUE/FALSE)")
		
		if(!is.matrix(resCov$M)) stop("M must be a matrix")
		
		if(nrow(resCov$M)!=ncol(resCov$M)) stop("M must be a square matrix") 
		
		if(nrow(resCov$M)!=traits) stop("M must have ", traits, " rows and columns")
		
		if(any(diag(resCov$M)==TRUE)) stop("All diagonal entries of M must be set to FALSE")
		
		resCov$M[upper.tri(resCov$M)]<-FALSE
		
		
		if(is.null(resCov$df0))
		{
			resCov$df0<-rep(5,traits)
			message("df0 set to  5 for all the traits")
		}
		
		if(is.null(resCov$S0))
		{
			#resCov$S0 <- rep(1,traits)
			#message("S0 was set to 1 for all the traits")
			resCov$S0<-(1-R2)*diag(Sy)*(resCov$df0+2)
			message("S0 was set to ")
			print(resCov$S0)	
		}
			
		if(is.null(resCov$var))
		{
			resCov$var<-100
			message("var was set to 100")
		}
			
		resCov$W<-matrix(0,nrow=traits,ncol=traits)
		resCov$PSI<-rep(NA,traits)
		
		for(k in 1:traits)
		{
			resCov$PSI[k]<-resCov$S0[k]/rchisq(n=1,df=resCov$df0[k])
		}
			
		#Objects for saving posterior means for MCMC
		resCov$post_W<-matrix(0,nrow=traits,ncol=traits)
		resCov$post_W2<-matrix(0,nrow=traits,ncol=traits)
		resCov$post_PSI<-rep(0,traits)
		resCov$post_PSI2<-rep(0,traits)

		#Output files
		resCov$fName_W<-paste(saveAt,"W_R.dat",sep="")
		resCov$fName_PSI<-paste(saveAt,"PSI_R.dat",sep="")
		resCov$f_W<-file(description=resCov$fName_W,open="w")
		resCov$f_PSI<-file(description=resCov$fName_PSI,open="w")
		
	}
	
	if(resCov$type=="FA")
	{
		message("Setting hyperparameters for FA R")
		
		if(is.null(resCov$M)) stop("M can not be null")
		
		if(!is.logical(resCov$M)) stop("M must be logical matrix (with entries being TRUE/FALSE)")
		
		if(!is.matrix(resCov$M)) stop("M must be a matrix")
		
		if(nrow(resCov$M)!=traits) stop("M must have ", traits, " rows")
		
		if(ncol(resCov$M)>traits) stop("Number of columns of M must be smaller than ", traits)
		
		resCov$nF<-ncol(resCov$M)
		
		
		if(is.null(resCov$df0))
		{
			resCov$df0<-rep(5,traits)
			message("df0 set to  5 for all the traits")
		}
		
		if(is.null(resCov$S0))
		{
			#resCov$S0 <- rep(1/100,traits)
			#message("S0 was set to 1/100 for all the traits")
			resCov$S0<-(1-R2)*diag(Sy)*(resCov$df0+2)
			message("S0 was set to ")
			print(resCov$S0)	
		}
		

		if(is.null(resCov$var))
		{
			resCov$var<-100
			message("var was set to 100")
		}
				
		sdU <- sqrt(diag(resCov$R))
        FA <- factanal(covmat = resCov$R, factors = resCov$nF)
        resCov$W <- matrix(nrow = traits, ncol = resCov$nF, 0)
        resCov$W[resCov$M] <- (diag(sdU) %*% FA$loadings)[resCov$M]
        resCov$PSI <- (sdU^2) * FA$uniquenesses + 1e-04
        resCov$R <- tcrossprod(resCov$W) + diag(resCov$PSI)
        resCov$Rinv<-solve(resCov$R)
        resCov$F <- matrix(nrow = nrow(error), ncol = resCov$nF, 0)
        
        #Objects for saving posterior means for MCMC
		resCov$post_W<-matrix(0,nrow=traits,ncol=resCov$nF)
		resCov$post_W2<-matrix(0,nrow=traits,ncol=resCov$nF)
		resCov$post_PSI<-rep(0,traits)
		resCov$post_PSI2<-rep(0,traits)
		
		#Output files
		resCov$fName_W<-paste(saveAt,"W_R.dat",sep="")
		resCov$fName_PSI<-paste(saveAt,"PSI_R.dat",sep="")
		resCov$f_W<-file(description=resCov$fName_W,open="w")
		resCov$f_PSI<-file(description=resCov$fName_PSI,open="w")
	}
		
	#Objects for saving posterior means for MCMC
	resCov$post_R<-matrix(0,nrow=traits,ncol=traits)
	resCov$post_R2<-matrix(0,nrow=traits,ncol=traits)
	
	resCov$fName_R<-paste(saveAt,"R.dat",sep="")
	resCov$f_R<-file(description=resCov$fName_R,open="w")
	
	message("Done")
	
	return(resCov)
	
}

#Evaluates
#partial Loglikelihood for complete and partially observed records
#-0.5 * n * log(det(R)) - 0.5 * sum (error_i' R^{-1} error_i) 
#error a matrix with errors, R residual variance covariance matrix

partialLogLik<-function(error,R)
{
	if(is.matrix(error))
	{
		n<-nrow(error)
		Linv<-solve(chol(R))
		ans<- -0.5*n*log(det(R)) - 0.5 * sum(crossprod(t(error), Linv)^2)
	}else{
		stop("error must be a matrix")
	}
	return(ans)
}


#Computes the  Deviance Information Criterion and Effective Number of Parameters 
#parameters: 
#y.back: matrix of dimension n x traits, NA's for missing values
#ETAHat the posterior mean of the conditional expectation 
#meanLogLik the posterior mean of logLik
#cte = -(n_complete_observed + n_partially_observed)*traits/2*log(2*pi)

getDIC<-function(y.back, ETAHat, meanLogLik, cte, complete_records, R,
                 missings=FALSE,
                 missing_records=NULL,
				 Dpatterns=NULL,
				 Upatterns=NULL,
				 dUpatterns=NULL,
				 dAllMissings=NULL)
{
		
		error <- y.back-ETAHat
		
		logLikAtPostMean <- cte
		
		LogLikAtPostMean <- logLikAtPostMean + partialLogLik(error[complete_records,,drop=FALSE],R)
		
		if(missings)
		{
		
			for (q in 1:length(dUpatterns))
			{
				#Some traits observed
				if(dUpatterns[q]!=dAllMissings)
				{
					#1=missing, 2=observed
					S22<-R[!Upatterns[q,],!Upatterns[q,],drop=FALSE]
					
					index<-missing_records[Dpatterns==dUpatterns[q]]
					
					#logLik
					logLikAtPostMean <- logLikAtPostMean + partialLogLik(error[index,!Upatterns[q,],drop=FALSE],S22)
				}
			
			}
		}
		
		pD <- -2 * (meanLogLik - logLikAtPostMean)
    	DIC <- pD - 2 * meanLogLik
    	
    	return(list(pD=pD,DIC=DIC,logLikAtPostMean=logLikAtPostMean))		
}


#Recursive structures for variance covariance matrix
#Model: U=UB+D (ONLY RECURSIVE ALLOWED!)  Current sample of random effects ('data') 
#M a traits x traits matrix with TRUE/FALSE indicating position of non-null
#recursive effects (FALSE in diagonal!)  PSI px1 the variance of the orthogonal
#shocks ...

#Function taken from MTM package
    
sample_G0_REC <- function(U, M, PSI, traits, priorVar = 100, 
                          df0 = rep(0, traits),S0 = rep(0, traits)) 
{

    B <- matrix(nrow = traits, ncol = traits, 0)
    
    for (i in 1:traits)
    {

        dimX <- sum(M[i, ])

        if (dimX > 0) {
            tmpX <- U[, M[i, ]]
            tmpY <- U[, i]

            C <- crossprod(tmpX)/PSI[i] + 1/priorVar
            CInv <- chol2inv(chol(C))
            rhs <- crossprod(tmpX, tmpY)/PSI[i]

            sol <- crossprod(CInv, rhs)
            L <- chol(CInv)
            shock <- crossprod(L, rnorm(dimX))
            tmpB <- as.numeric(sol + shock)
            B[i, M[i, ]] <- tmpB
            uStar <- tmpY - matrix(tmpX, ncol = dimX) %*% (tmpB)
            SS <- as.numeric(crossprod(uStar)) + S0[i]
            df <- nrow(U) + df0[i]
            PSI[i] <- SS/rchisq(n = 1, df = df)
        }else{
            SS <- as.numeric(crossprod(U[, i])) + S0[i]
            df <- nrow(U) + df0
            PSI[i] <- SS/rchisq(n = 1, df = df)
        }
    }

    tmp <- solve(diag(traits) - B)
    G <- tmp %*% diag(PSI) %*% t(tmp)

    out <- list(B = B, PSI = PSI, G = G)

    return(out)
}


#Sampler for FA
#Model U=BF+D
#G_0 = B B' + PSI
#B is a matrix of loadings (regressions of the original random effects into common factors)
#F is a matrix of common factors
#M a logical matrix of the same dimensions that F, with TRUE for the loadings that the 
#user wants to estimate
#PSI is a diagonal matrix whose non-null entries give the variances of factors that 
#are trait-specific
#nF number of common factors, is equal to the number of columns of B
#nD: number of rows of U
#df0: degrees of freedom associated to the prior of PSI
#S0: scale parameter associated to the prior for PSI
#priorVar: prior variance if the Gaussian prior assigned to the unknown loadings

#Function taken from MTM package


sample_G0_FA <- function(U, F, M, B, PSI, traits, nF, nD, 
                       df0 = rep(1, traits), S0 = rep(1/100,traits), priorVar = 100) 
{
    ## sampling common factors LOOP OVER FACTORS
    for (i in 1:nF) 
    {
        tmpY <- U - F[, -i] %*% matrix((B[, -i]), ncol = traits)
        rhs <- tmpY %*% matrix(B[, i]/PSI, ncol = 1)
        CInv <- 1/(sum((B[, i]^2)/PSI) + 1)
        sol <- CInv * rhs
        SD <- sqrt(CInv)
        F[, i] <- rnorm(n = nD, sd = SD, mean = sol)
    }

    # sampling loadings LOOP OVER TRAITS LOOP OVER FACTORS
    for (i in 1:traits) 
    {
    
        for (j in 1:nF) 
        {
            if (M[i, j]) 
            {
                tmpY <- U[, i] - F[, -j] %*% matrix(B[i, -j], ncol = 1)
                CInv <- 1/as.numeric(crossprod(F[, j])/PSI[i] + 1/priorVar)
                rhs <- as.numeric(crossprod(F[, j], tmpY)/PSI[i])
                sol <- CInv * rhs
                SD <- sqrt(CInv)
                B[i, j] <- rnorm(n = 1, mean = sol, sd = SD)
            }
        }
        
        D <- U[, i] - F %*% B[i, ]
        df <- df0[i] + nD
        SS <- S0[i] + crossprod(D)
        PSI[i] <- SS/rchisq(df = df, n = 1)
    }
    
    if (nF > 1) 
    {
        B <- varimax(B)$loadings[]
    }
    
    G <- tcrossprod(B) + diag(PSI)
    
    out <- list(F = F, PSI = PSI, B = B, G=G)
    
    return(out)
}

#Gibbs sampler for DIAG covariance matrices
#Function taken from MTM package

sample_G0_DIAG <- function(U, traits = ncol(U), n = nrow(U), 
                           df0 = rep(0, traits),
                           S0 = diag(0, traits)) 
{

    G <- matrix(nrow = traits, ncol = traits, 0)
    
    ## sampling common factors LOOP OVER FACTORS
    for (i in 1:traits) 
    {
        tmp_SS <- sum(U[, i]^2) + S0[i]
        tmp_df <- n + df0[i]
        G[i, i] <- tmp_SS/rchisq(df = tmp_df, n = 1)
    }

    return(G)
}

#Sample the vector mu (intercept)
#The prior for mu is non informative
#For trait k, y_k = mu_k * 1 + eta_j + e_j
#ystar_k = y_k - eta_j = mu_k * 1 + e_j 
#mu | else ~ MN(m,R/n), where m is the vector of sample means obtained obtained from ystar
#Arguments: 
#ystar: matrix, individuals in rows, traits in columns
#R: residual covariance matrix
#n: number of rows of y
#traits: number of traits

sample_mu <- function(ystar, R, n, traits) 
{
    sol <- colMeans(ystar)
    L <- chol(R)/sqrt(n)
    mu <- as.vector(crossprod(L, rnorm(n=traits,mean=0,sd=1)) + sol)
    return(mu)
}

# Function to read effects saved by Multitrait when ETA[[j]]$saveEffects=TRUE
# It returns a 3D array, with dim=c(nRow,p,traits)
# nRow number of MCMC samples saved,
# p number of predictors
# traits number of traits

readBinMatMultitrait<-function(filename,storageMode="double")
{

    if(!storageMode%in%c("single","double")){
        stop("storageMode can either be 'single' or 'double' (default)")
    }
    
  	fileIn<-gzfile(filename,open='rb')
 	nRow<-readBin(fileIn,n=1,what=numeric(),size=ifelse(storageMode=="single",4,8))
 	traits<-readBin(fileIn,n=1,what=numeric(),size=ifelse(storageMode=="single",4,8))
 	p<-readBin(fileIn,n=1,what=numeric(),size=ifelse(storageMode=="single",4,8))

 	Beta<-array(data=NA,dim=c(nRow,p,traits))
 	
 	for(j in 1:nRow)
 	{
 		for(k in 1:traits)
 		{
 			Beta[j,,k]<-readBin(fileIn,n=p,what=numeric(),size=ifelse(storageMode=="single",4,8))
 		}
 	}
 	 	
 	close(fileIn)
 	
 	return(Beta)
}

#Get genetic co-variance matrix
#Internal function, used by getGCovar

getG0i<-function(Z,Bi)
{
    U<-Z%*%Bi
    G0i<-cov(U)
    return(G0i[row(G0i)>=col(G0i)])
}

#Genetic co-variance matrix using MCMC samples
#Lehermeier et al., 2017. 
#Genomic Variance Estimates: With or without Disequilibrium Covariances?
#J Anim Breed Genet, 134(3):232-241.
#Arguments:
#X: matrix of covariates
#B: samples for regression coefficients, 3D array, with dim=c(nRow,p,traits)
# nRow number of MCMC samples saved,
# p number of predictors
# traits number of traits

getGCovar<-function(X,B)
{
        q<-dim(B)[3]
        G<-t(apply(FUN=getG0i,X=B,Z=X,MARGIN=1))
        return(G)
}

#Function to fit multi-trait model
#Arguments:
#y: A matrix of dimension n * t, where t is the number of traits, NAs allowed.
#ETA: Linear predictor A two level list to specify the linear predictors.
#intecept: Logical, TRUE an intercept is included and FALSE an intercept is not included.
#resCov: List to specify the prior for residual covariance matrix (R).
#nIter, burnIn, thin: Number of iterations, burnIn and thin.
#verbose: logical, TRUE/FALSE to print iteration history

Multitrait<-function(y,
					 ETA,
					 intercept=TRUE,
                     resCov = list(df0=5,S0=NULL,type="UN"),
                     R2=0.5,
                     nIter=1000,burnIn=500,thin=10,
                     saveAt="",verbose=TRUE)
{
	
	#Check inputs
	if(!is.matrix(y)) stop("y must be a matrix\n")
	
	traits<-ncol(y)
	
	if(traits<2) stop("y must hava at least 2 columns\n")
	
	n<-nrow(y)
	
	#Compute sample variance covariance matrix
	Sy<-cov(y,use="pairwise.complete.obs")
	Sy[is.na(Sy)]=0
	
	#Deep copy of y, DO NOT REPLACE with y.back<-y, it does not work, both objects 
	#share the same memory address
	y.back<-matrix(NA,nrow=n,ncol=traits)
	y.back<-y[]
	
	
	#Now check if there are missing values
	complete_records <- which(complete.cases(y))
	n_complete_observed <- length(complete_records)
	n_partially_observed <- 0                        #Just initial value
	n_complete_missing <- 0							 #Just initial value
	
	#Note: 
	#n=n_complete_observed + n_partially_observed + n_complete_missing
	
	missings<-any(is.na(y))

	if(missings)
	{
		
		missing_records<-which(!complete.cases(y))
	
		#patterns of missing records, should be a matrix
		patterns<-is.na(y[missing_records,])

		#Case when there is only one missing record
		if(!is.matrix(patterns))
		{
			patterns<-matrix(patterns,ncol=ncol(y))	
		}
	
		Dpatterns<-apply(patterns,1,logicalToDec)     #patterns in decimal
		Upatterns<-unique(patterns,drop=FALSE)	      #Unique patterns
		dUpatterns<-apply(Upatterns,1,logicalToDec)	  #decimal of unique patterns
		dAllMissings<-logicalToDec(rep(TRUE,traits))  #decimal representation of a record
		                                              #with all missing values
		
		m<-colMeans(y,na.rm=TRUE)
		for(k in 1:traits)
		{
			tmp<-is.na(y[,k])
			y[tmp,k]<-m[k]
			rm(tmp)
		}
		
		#Partially observed records
		n_partially_observed <- sum(Dpatterns!=dAllMissings)
		
		#Complete missing
		n_complete_missing <-  sum(Dpatterns==dAllMissings)
		
	}else{
		
		missing_records <- NULL
		Dpatterns <- NULL
		Upatterns <- NULL
		dUpatterns <- NULL
		dAllMissings <- NULL
		
	}
	
	#For likelihood
	cte<- -(n_complete_observed + n_partially_observed)*traits/2*log(2*pi)
	
	if(intercept)
	{
		mu<-colMeans(y)
		post_mu<-rep(0,traits)
		post_mu2<-rep(0,traits)
		f_mu<-file(description=paste(saveAt,"mu.dat",sep=""),open="w")
	}
	
	#Setting the linear terms
	
	nLT <- ifelse(is.null(ETA), 0, length(ETA))
	
	if(nLT<1) stop("Provide at least a linear predictor in ETA\n")
	
	#Names of linear terms
	if(is.null(names(ETA)))
    { 
		names(ETA)<-rep("",nLT)
    }
    
    nRow<-nIter/thin
    if(nRow<1) stop("Check nIter, thin\n")
	
	for(j in 1:nLT)
	{
		
		if(names(ETA)[j]=="")
	    {
	       	ETA[[j]]$Name=paste("ETA_",j,sep="")
	    }else{
            ETA[[j]]$Name=paste("ETA_",names(ETA)[j],sep="")
	    }
	
        if(!(ETA[[j]]$model %in% c("SpikeSlab","BRR","RKHS","FIXED"))) 
        {
        	stop("Error in ETA[[", j, "]]", " model ", ETA[[j]]$model, " not implemented (note: evaluation is case sensitive)")    
        }
        
		ETA[[j]]<-switch(ETA[[j]]$model,
						SpikeSlab=setLT.DiracSS_mt(LT=ETA[[j]],traits=traits,j=j,Sy=Sy,nLT=nLT,R2=R2,saveAt=saveAt,nRow=nRow),
						BRR=setLT.BRR_mt(LT=ETA[[j]],traits=traits,j=j,Sy=Sy,nLT=nLT,R2=R2,saveAt=saveAt,nRow=nRow),
						RKHS=setLT.RKHS_mt(LT=ETA[[j]],traits=traits,j=j,Sy=Sy,nLT=nLT,R2=R2,saveAt=saveAt),
						FIXED=setLT.FIXED_mt(LT=ETA[[j]],traits=traits,j=j,saveAt=saveAt,nRow=nRow))
	
	} 	#End of cycle for setting linear terms
	
	
	#error initialization assuming beta=0 in all the entries
	#Deep copy of y, DO NOT REPLACE with error<-y, it does not work, both objects 
	#share the same memory address
	
	error<-matrix(NA,nrow=n,ncol=traits)
	error<-y[]
	
	if(intercept)
	{
		error<-sweep(x=error,MARGIN=2,STATS=mu,FUN="-")
	}

	#Initialization of variance covariance matrix for errors (R)
	resCov<-setResCov(resCov=resCov,traits=traits,error=error,Sy=Sy,R2=R2,saveAt=saveAt)
	
	ETAHat<-matrix(0,nrow=n,ncol=traits)
	ETAHat2<-matrix(0,nrow=n,ncol=traits)
	
	post_logLik<-0
	
	#Objects for running means
	nSums<-0	
	
	#Iterations
	for(iter in 1:nIter)
	{
		start<-proc.time()[3]
		
		logLik <- cte 
		
		if(intercept)
		{
			error<-sweep(x=error,MARGIN=2,STATS=mu,FUN="+")
			mu<-sample_mu(ystar=error, R=resCov$R, n=n, traits=traits)
			error<-sweep(x=error,MARGIN=2,STATS=mu,FUN="-")
			#print(mu) 
		}
		
		for(j in 1:nLT)
		{
		
			#SpikeSlab
			if(ETA[[j]]$model=="SpikeSlab")
			{
					# for(k in 1:traits)
# 					{
# 						#cat("k=",k,"\n")	
# 		
# 						S11 <- ETA[[j]]$Cov$Omega[k,k,drop=FALSE]
# 						S22 <- ETA[[j]]$Cov$Omega[-k,-k,drop=FALSE]
# 						S12 <- ETA[[j]]$Cov$Omega[k,-k,drop=FALSE]
# 						tmp12 <- as.vector(S12%*%solve(S22))
# 						tmp11 <- tmp12%*%t(S12)
# 						sigma2 <- as.numeric(S11-tmp11)
# 		
# 						logPriorOdds<-log(ETA[[j]]$inclusionProb$probIn[k]/(1-ETA[[j]]$inclusionProb$probIn[k]))
# 		
# 						#b, d, beta and error are overwritten with this .Call
# 		
# 						.Call("sampler_DiracSS_mt", k, logPriorOdds, n, ETA[[j]]$p,
# 							  traits, resCov$Rinv, ETA[[j]]$X, error,
# 							  ETA[[j]]$beta,
# 							  ETA[[j]]$b,
# 							  ETA[[j]]$d,
# 							  ETA[[j]]$x2,
# 							  tmp12,
# 							  sigma2,
# 							  ETA[[j]]$Cov$Omegainv[k,-k],
# 							  ETA[[j]]$Cov$Omegainv[k,k])
# 		
# 						#Sampling inclusion probabilities | else
# 						mrkIn <- sum(ETA[[j]]$d[,k])
# 						shape1 <- mrkIn + ETA[[j]]$inclusionProb$countsIn[k]
# 						shape2 <- ETA[[j]]$p - mrkIn + ETA[[j]]$inclusionProb$countsOut[k]
# 						ETA[[j]]$inclusionProb$probIn[k]=rbeta(shape1 = shape1,
# 													           shape2 = shape2, 
# 													           n = 1)
# 												
# 					}#End of loop for traits
					
					#### BEGIN NEW code 

					logPriorOdds<-log(ETA[[j]]$inclusionProb$probIn/(1-ETA[[j]]$inclusionProb$probIn))
				
					.Call("sampler_DiracSS_mt_v2",logPriorOdds,n, ETA[[j]]$p,
						  traits,resCov$Rinv, ETA[[j]]$X, error,
							  ETA[[j]]$beta,
							  ETA[[j]]$b,
							  ETA[[j]]$d,
							  ETA[[j]]$x2,
							  ETA[[j]]$Cov$Omega,
							  ETA[[j]]$Cov$Omegainv)
					
					for(k in 1:traits)
					{
						#Sampling inclusion probabilities | else
						mrkIn <- sum(ETA[[j]]$d[,k])
						shape1 <- mrkIn + ETA[[j]]$inclusionProb$countsIn[k]
						shape2 <- ETA[[j]]$p - mrkIn + ETA[[j]]$inclusionProb$countsOut[k]
						ETA[[j]]$inclusionProb$probIn[k]=rbeta(shape1 = shape1,
													           shape2 = shape2, 
													           n = 1)
						
					}
					
					### END NEW code

					
					#Sampling from Omega | else						
					if(ETA[[j]]$Cov$type=="UN")
					{
						S4<-crossprod(ETA[[j]]$b)		
						ETA[[j]]$Cov$Omega<-riwish(v=ETA[[j]]$Cov$df0+traits+ETA[[j]]$p,
					                               S=S4+ETA[[j]]$Cov$S0)
					                               
						ETA[[j]]$Cov$Omegainv<-solve(ETA[[j]]$Cov$Omega)
					}
					
					if(ETA[[j]]$Cov$type=="DIAG")
					{
						ETA[[j]]$Cov$Omega<-sample_G0_DIAG(U=ETA[[j]]$b, traits=traits, 
						                                   n=nrow(ETA[[j]]$b),
						                                   df0=ETA[[j]]$Cov$df0,
                             			                   S0=ETA[[j]]$Cov$S0)
                             			                   
                        ETA[[j]]$Cov$Omegainv<-solve(ETA[[j]]$Cov$Omega)
					}
					
					if(ETA[[j]]$Cov$type=="REC")
					{
						tmp<-sample_G0_REC(U=ETA[[j]]$b, M=ETA[[j]]$Cov$M, 
							               PSI=ETA[[j]]$Cov$PSI, 
							               traits=traits, 
							               priorVar = ETA[[j]]$Cov$var, 
                             			   df0 = ETA[[j]]$Cov$df0,
                             			   S0 = ETA[[j]]$Cov$S0)
                            
                        ETA[[j]]$Cov$Omega<-tmp$G
                        ETA[[j]]$Cov$Omegainv<-solve(ETA[[j]]$Cov$Omega)
                        ETA[[j]]$Cov$W<-tmp$B
                        ETA[[j]]$Cov$PSI<-tmp$PSI
                            
                        rm(tmp)
					}
					
					if(ETA[[j]]$Cov$type=="FA")
					{
						tmp<-sample_G0_FA(U=ETA[[j]]$b, F=ETA[[j]]$Cov$F, M=ETA[[j]]$Cov$M, 
						    	          B=ETA[[j]]$Cov$W, PSI=ETA[[j]]$Cov$PSI, 
						        	      traits=traits, nF=ETA[[j]]$Cov$nF, 
						        	      nD=ETA[[j]]$Cov$nD, df0 = ETA[[j]]$Cov$df0, 
                       					  S0 = ETA[[j]]$Cov$S0, 
                       					  priorVar = ETA[[j]]$Cov$var) 

						ETA[[j]]$Cov$F<-tmp$F
						ETA[[j]]$Cov$PSI<-tmp$PSI
						ETA[[j]]$Cov$W<-tmp$B
						ETA[[j]]$Cov$Omega<-tmp$G
                    	ETA[[j]]$Cov$Omegainv<-solve(ETA[[j]]$Cov$Omega)
						
						rm(tmp)
						
					}
					
					#Compute the variance-covariance matrix for beta=b*d
					ETA[[j]]$Cov$Sigma<-cov(ETA[[j]]$beta)					
					
			} #End of SpikeSlab
			
			if(ETA[[j]]$model%in%c("BRR","RKHS"))
			{		

# 					for(k in 1:traits)
# 					{
# 						#cat("k=",k,"\n")	
# 						
# 						#beta and error are overwritten with this .Call
# 						.Call("sampler_BRR_mt", k, n, ETA[[j]]$p,
# 							  traits, resCov$Rinv, ETA[[j]]$X, error,
# 							  ETA[[j]]$beta,
# 							  ETA[[j]]$x2,
# 							  ETA[[j]]$Cov$Omegainv[k,-k],
# 							  ETA[[j]]$Cov$Omegainv[k,k])
# 												
# 					}#End of loop for traits

					### BEGIN NEW code
					
					.Call("sampler_BRR_mt_v2", n, ETA[[j]]$p,
						  traits, resCov$Rinv, ETA[[j]]$X, error,
						  ETA[[j]]$beta,
						  ETA[[j]]$x2,
						  ETA[[j]]$Cov$Omegainv)
						  
					### END NEW code
					
					#Sampling from Omega | else
					
					if(ETA[[j]]$Cov$type=="UN")
					{
						S4<-crossprod(ETA[[j]]$beta)		
						ETA[[j]]$Cov$Omega<-riwish(v=ETA[[j]]$Cov$df0+traits+ETA[[j]]$p,
					    	                   S=S4+ETA[[j]]$Cov$S0)
						ETA[[j]]$Cov$Omegainv<-solve(ETA[[j]]$Cov$Omega)
						
					}#End of UN
					
					if(ETA[[j]]$Cov$type=="DIAG")
					{
												
						ETA[[j]]$Cov$Omega<-sample_G0_DIAG(U=ETA[[j]]$beta, traits=traits, 
						                                   n=nrow(ETA[[j]]$beta),
						                                   df0=ETA[[j]]$Cov$df0,
                             			                   S0=ETA[[j]]$Cov$S0)
                             			                   
                        ETA[[j]]$Cov$Omegainv<-solve(ETA[[j]]$Cov$Omega)
                        
					}#End of DIAG
					
					if(ETA[[j]]$Cov$type=="REC")
					{
						tmp<-sample_G0_REC(U=ETA[[j]]$beta, M=ETA[[j]]$Cov$M, 
							               PSI=ETA[[j]]$Cov$PSI, 
							               traits=traits, 
							               priorVar = ETA[[j]]$Cov$var, 
                             			   df0 = ETA[[j]]$Cov$df0,
                             			   S0 = ETA[[j]]$Cov$S0)
                            
                        ETA[[j]]$Cov$Omega<-tmp$G
                        ETA[[j]]$Cov$Omegainv<-solve(ETA[[j]]$Cov$Omega)
                        ETA[[j]]$Cov$W<-tmp$B
                        ETA[[j]]$Cov$PSI<-tmp$PSI
                            
                        rm(tmp)
					}#End of REC
					
					if(ETA[[j]]$Cov$type=="FA")
					{
						tmp<-sample_G0_FA(U=ETA[[j]]$beta, F=ETA[[j]]$Cov$F, M=ETA[[j]]$Cov$M, 
						    	          B=ETA[[j]]$Cov$W, PSI=ETA[[j]]$Cov$PSI, 
						        	      traits=traits, nF=ETA[[j]]$Cov$nF, 
						        	      nD=ETA[[j]]$Cov$nD, df0 = ETA[[j]]$Cov$df0, 
                       					  S0 = ETA[[j]]$Cov$S0, 
                       					  priorVar = ETA[[j]]$Cov$var) 

						ETA[[j]]$Cov$F<-tmp$F
						ETA[[j]]$Cov$PSI<-tmp$PSI
						ETA[[j]]$Cov$W<-tmp$B
						ETA[[j]]$Cov$Omega<-tmp$G
                    	ETA[[j]]$Cov$Omegainv<-solve(ETA[[j]]$Cov$Omega)
						
						rm(tmp)
						
					}#End of FA
					
					
			}#End of BRR and RKHS
			
			if(ETA[[j]]$model=="FIXED")
			{
					if(ETA[[j]]$common)
					{
						for(k in 1:traits)
						{							
							#beta and error are overwritten with this .Call
							.Call("sampler_BRR_mt", k, n, ETA[[j]]$p,
								  traits, resCov$Rinv, ETA[[j]]$X, error,
								  ETA[[j]]$beta,
								  ETA[[j]]$x2,
								  ETA[[j]]$Cov$Omegainv[k,-k],
								  ETA[[j]]$Cov$Omegainv[k,k])
												
						}#End of loop for traits
						
					}else{
					
						for(k in 1:traits)
						{
							#X=[X_1,...,X_t], k=1,...,t (traits), we pass as
							
							Ck<-which(ETA[[j]]$idColumns==k)
							
							if(length(Ck)>0)
							{
								Xk<-ETA[[j]]$X[,Ck,drop=FALSE]
								x2k<-ETA[[j]]$x2[Ck]
								#beta and error are overwritten with this .Call
								.Call("sampler_BRR_mt_fixed", k, n, Ck,length(Ck),
								  	   traits, resCov$Rinv, Xk, error,
								  		ETA[[j]]$beta, x2k)
							}
													
						}#End of loop for traits
						
					}
			} #End of FIXED			
		} #End of loop linear terms
	
		#Sampling from R
		
		if(resCov$type=="UN")
		{
			CP<-crossprod(error)
			resCov$R<-riwish(v=resCov$df0+n,S=CP+resCov$S0)
		}
		
		if(resCov$type=="DIAG")
		{
		
			resCov$R<-sample_G0_DIAG(U=error, traits=traits, n=nrow(error),
			 			             df0=resCov$df0, S0=resCov$S0)
		}
		
		if(resCov$type=="REC")
		{
			tmp<-sample_G0_REC(U=error,M=resCov$M, PSI=resCov$PSI, 
							   traits=traits, priorVar = resCov$var, 
                               df0 = resCov$df0,
                               S0 = resCov$S0)
                            
            resCov$R<-tmp$G
            resCov$W<-tmp$B
            resCov$PSI<-tmp$PSI
            rm(tmp)
            
		}
		
		if(resCov$type=="FA")
		{
				tmp<-sample_G0_FA(U=error, F=resCov$F, M=resCov$M, 
						    	  B=resCov$W, PSI=resCov$PSI, 
						          traits=traits, nF=resCov$nF, nD=nrow(error),
                       			  df0 = resCov$df0, 
                       			  S0 = resCov$S0, 
                       			  priorVar = resCov$var)

                resCov$F<-tmp$F
				resCov$PSI<-tmp$PSI
				resCov$W<-tmp$B
				resCov$R<-tmp$G
		}
		
		resCov$Rinv<-solve(resCov$R)		
		
		#Impute missing values
				
		#Linear predictor
		lp <- y-error
		
		#Log likelihood for complete cases
		logLik <- logLik + partialLogLik(error[complete_records,,drop=FALSE],resCov$R)
		
		if(missings)
		{
			
			for (q in 1:length(dUpatterns))
			{
				#Some traits observed
				if(dUpatterns[q]!=dAllMissings)
				{
					#1=missing, 2=observed
					S11<-resCov$R[Upatterns[q,],Upatterns[q,],drop=FALSE]
					S12<-resCov$R[Upatterns[q,],!Upatterns[q,],drop=FALSE]
					S21<-resCov$R[!Upatterns[q,],Upatterns[q,],drop=FALSE]
					S22<-resCov$R[!Upatterns[q,],!Upatterns[q,],drop=FALSE]
			
					tmp12<-S12%*%solve(S22)
			
					index<-missing_records[Dpatterns==dUpatterns[q]]
					
					#logLik
					logLik <- logLik + partialLogLik(error[index,!Upatterns[q,],drop=FALSE],S22)
					
					mu1<-lp[index,Upatterns[q,],drop=FALSE]
					mu2<-lp[index,!Upatterns[q,],drop=FALSE]
			
					Sigma<-S11-tmp12%*%S21
			
					tmp3<-tmp12%*%t(y[index,!Upatterns[q,],drop=FALSE]-mu2)
					tmp3<-t(tmp3)
			
					mean<-mu1+tmp3
			
					#Impute y and overwrite the value
					y[index,Upatterns[q,]]<-mean+mvrnorm(n=length(index),mu=rep(0,nrow(Sigma)),Sigma=Sigma)
					
				
				}else{
				
					#Observations for all traits are missing for some records
					index<-missing_records[Dpatterns==dUpatterns[q]]
									
					#predict y and overwrite the value
					y[index,]<-lp[index,]+mvrnorm(n=length(index),mu=rep(0,nrow(resCov$R)),Sigma=resCov$R)
					
				}
				
			}
		
			#Update residuals
			for(k in 1:traits)
			{			
				index<-missing_records[patterns[,k]]
				error[index,k] <- y[index,k]-lp[index,k]
			}
		
		}
		
		#Saving files
		if(iter%%thin==0)
		{
			#mu
			if(intercept)
			{
				write(mu,ncolumns=length(mu),file=f_mu,append=TRUE,sep=" ")
			}
			
			#resCov
			tmp <- vech(resCov$R)
			write(tmp, ncolumns = length(tmp), file = resCov$f_R, append = TRUE, sep = " ")
			rm(tmp)
			
			if(resCov$type%in%c("REC","FA"))
			{
				if (sum(resCov$M) > 0) 
				{
					tmp <- resCov$W[resCov$M]
					write(tmp, ncolumns = length(tmp), file = resCov$f_W, append = TRUE, 
					      sep = " ")
					rm(tmp)
				}
								
				write(resCov$PSI, ncolumns = length(resCov$PSI), file = resCov$f_PSI, 
				      append = TRUE, sep = " ")
			}
			
			for(j in 1:nLT)
			{
				if(ETA[[j]]$model%in%c("SpikeSlab","BRR","RKHS"))
				{
				
					if(ETA[[j]]$Cov$type%in%c("REC","FA"))
					{
						if (sum(ETA[[j]]$Cov$M) > 0) 
						{
							tmp <- ETA[[j]]$Cov$W[ETA[[j]]$Cov$M]
							write(tmp, ncolumns = length(tmp), file = ETA[[j]]$Cov$f_W, 
							      append = TRUE, sep = " ")
							rm(tmp)
						}
								
						write(ETA[[j]]$Cov$PSI, ncolumns = length(ETA[[j]]$Cov$PSI), 
						      file = ETA[[j]]$Cov$f_PSI, append = TRUE, sep = " ")			
					}#End REC and FA
					
					if(ETA[[j]]$Cov$type%in%c("UN","DIAG"))
					{
						tmp <- vech(ETA[[j]]$Cov$Omega)
						write(tmp, ncolumns = length(tmp), file = ETA[[j]]$Cov$f_Omega, append = TRUE, sep = " ")
						rm(tmp)
					}
					
				}#End of SpikeSlab, BRR and RKHS
				
			}#End for
			
			#Saving beta effects and indicator variables
			for(j in 1:nLT)
			{
				if(ETA[[j]]$model%in%c("SpikeSlab","BRR","FIXED"))
				{
					if(ETA[[j]]$saveEffects)
					{
                        writeBin(object=as.vector(ETA[[j]]$beta),
                                 con=ETA[[j]]$fileEffects,
                                 size=ifelse(ETA[[j]]$storageMode=="single",4,8))
                    }
                    
                    if(ETA[[j]]$model=="SpikeSlab")                          
                    {
                        if(ETA[[j]]$saveIndicators)
                        {
                                #entries are saved in "single" mode
                                writeBin(object=as.vector(ETA[[j]]$d),
                                         con=ETA[[j]]$fileIndicators,
                                         size=4)
                        }
                    }
                }
			}
			
		}#End of saving files
		
		
		#Running means
		if((iter>burnIn) & (iter%%thin==0))
		{
			nSums<-nSums + 1
			fraction <- (nSums - 1)/nSums
			
			#mu
			if(intercept)
			{
				post_mu<-post_mu * fraction + mu/nSums
				post_mu2<-post_mu2 * fraction + mu^2/nSums
			}
			
			#Predictions
			ETAHat<-ETAHat * fraction + lp/nSums
			ETAHat2<-ETAHat2 * fraction + lp^2/nSums
			
			#Residual
			resCov$post_R <- resCov$post_R * fraction + resCov$R/nSums
			resCov$post_R2<- resCov$post_R2 * fraction + resCov$R^2/nSums
			
			if(resCov$type%in%c("REC","FA"))
			{
				resCov$post_W <- resCov$post_W * fraction + resCov$W/nSums
    			resCov$post_W2 <- resCov$post_W2 * fraction + resCov$W^2/nSums
    					
    			resCov$post_PSI <- resCov$post_PSI * fraction + resCov$PSI/nSums
    			resCov$post_PSI2 <- resCov$post_PSI2 * fraction + resCov$PSI^2/nSums
    			
			}		
			
			#Likelihood
			post_logLik<- post_logLik * fraction + logLik/nSums
			
			for(j in 1:nLT)
			{
				#All the models, post_beta and post_beta2
				ETA[[j]]$post_beta <- ETA[[j]]$post_beta * fraction + ETA[[j]]$beta/nSums
				ETA[[j]]$post_beta2 <- ETA[[j]]$post_beta2 * fraction + ETA[[j]]$beta^2/nSums
				
				
				#post_b, post_b2, post_d, post_d2
				if(ETA[[j]]$model=="SpikeSlab")
				{
					ETA[[j]]$post_b <- ETA[[j]]$post_b * fraction + ETA[[j]]$b/nSums
					ETA[[j]]$post_b2 <- ETA[[j]]$post_b2 * fraction + ETA[[j]]$b^2/nSums
					
					ETA[[j]]$post_d <- ETA[[j]]$post_d * fraction + ETA[[j]]$d/nSums
					ETA[[j]]$post_d2 <- ETA[[j]]$post_d2 * fraction + ETA[[j]]$d^2/nSums				
				}

				if(ETA[[j]]$model%in%c("SpikeSlab","BRR","RKHS"))
				{
					#post_Omega, post_Omega2
					ETA[[j]]$Cov$post_Omega <- ETA[[j]]$Cov$post_Omega * fraction + ETA[[j]]$Cov$Omega/nSums
    				ETA[[j]]$Cov$post_Omega2 <- ETA[[j]]$Cov$post_Omega2 * fraction + ETA[[j]]$Cov$Omega^2/nSums
    				
    				if(ETA[[j]]$model%in%"SpikeSlab")
    				{
    					ETA[[j]]$Cov$post_Sigma <- ETA[[j]]$Cov$post_Sigma * fraction + ETA[[j]]$Cov$Sigma/nSums
    					ETA[[j]]$Cov$post_Sigma2 <- ETA[[j]]$Cov$post_Sigma2 * fraction + ETA[[j]]$Cov$Sigma^2/nSums
    				}
    				
    				if(ETA[[j]]$Cov$type%in%c("REC","FA"))
    				{
    					ETA[[j]]$Cov$post_W <- ETA[[j]]$Cov$post_W * fraction + ETA[[j]]$Cov$W/nSums
    					ETA[[j]]$Cov$post_W2 <- ETA[[j]]$Cov$post_W2 * fraction + ETA[[j]]$Cov$W^2/nSums
    					
    					ETA[[j]]$Cov$post_PSI <- ETA[[j]]$Cov$post_PSI * fraction + ETA[[j]]$Cov$PSI/nSums
    					ETA[[j]]$Cov$post_PSI2 <- ETA[[j]]$Cov$post_PSI2 * fraction + ETA[[j]]$Cov$PSI^2/nSums
    					
    				}			
				}
			}#End for linear terms
		}
		
		end<-proc.time()[3]
		
		if(verbose)
		{
			cat(paste("Iter: ", iter, "time: ", round(end - start, 4)," s \n"))
		}
				
	}#End of loop for iterations
	
	if(intercept)
	{
		mu<-post_mu
		SD.mu<-sqrt(post_mu2-post_mu^2)
	}
		
	#Renaming/removing objects in ETA
	for(j in 1:nLT)
	{
			#beta and SD.beta
			if(ETA[[j]]$model%in%c("SpikeSlab","BRR","FIXED"))
			{
				ETA[[j]]$beta<-ETA[[j]]$post_beta
				ETA[[j]]$SD.beta<-sqrt(ETA[[j]]$post_beta2-ETA[[j]]$post_beta^2)
			}
			
			#u, FIXME: SD.u missing 
			if(ETA[[j]]$model=="RKHS")
			{
				ETA[[j]]$u<-ETA[[j]]$X%*%ETA[[j]]$post_beta
			}
			
			#post_b, post_b2, post_d, post_d2
			if(ETA[[j]]$model=="SpikeSlab")
			{
				ETA[[j]]$b<-ETA[[j]]$post_b
				ETA[[j]]$SD.b<-sqrt(ETA[[j]]$post_b2-ETA[[j]]$post_b^2)
		
				ETA[[j]]$d<-ETA[[j]]$post_d
				ETA[[j]]$SD.d<-sqrt(ETA[[j]]$post_d2-ETA[[j]]$post_d^2)
			}

			if(ETA[[j]]$model%in%c("SpikeSlab","BRR","RKHS"))
			{
				#Omega
				ETA[[j]]$Cov$Omega<-ETA[[j]]$Cov$post_Omega
				ETA[[j]]$Cov$SD.Omega<-sqrt(ETA[[j]]$Cov$post_Omega2-ETA[[j]]$Cov$post_Omega^2)
				
				if(ETA[[j]]$model=="SpikeSlab")
				{
					ETA[[j]]$Cov$Sigma<-ETA[[j]]$Cov$post_Sigma
					ETA[[j]]$Cov$SD.Sigma<-sqrt(ETA[[j]]$Cov$post_Sigma2-ETA[[j]]$Cov$post_Sigma^2)
				}
				
				if(ETA[[j]]$Cov$type%in%c("REC","FA"))
				{
					ETA[[j]]$Cov$W<-ETA[[j]]$Cov$post_W
					ETA[[j]]$Cov$SD.W<-sqrt(ETA[[j]]$Cov$post_W2-ETA[[j]]$Cov$post_W^2)
				
					ETA[[j]]$Cov$PSI<-ETA[[j]]$Cov$post_PSI
					ETA[[j]]$Cov$SD.PSI<-sqrt(ETA[[j]]$Cov$post_PSI2-ETA[[j]]$Cov$post_PSI^2)	
				}
				
				tmp<-which(names(ETA[[j]]$Cov)%in%c("post_Omega","post_Omega2","Omegainv",
				                                    "post_W","post_W2","post_PSI","post_PSI2",
				                                    "F",
				                                    "post_Sigma","post_Sigma2"))
				ETA[[j]]$Cov<-ETA[[j]]$Cov[-tmp] 
				
				rm(tmp)                                   
								
			}
			
			#Deep cleaning!!
			tmp<-which(names(ETA[[j]])%in%c("X","x2","post_beta","post_beta2",
			                                "post_b","post_b2","post_d","post_d2"))
			ETA[[j]]<-ETA[[j]][-tmp]
			
			rm(tmp)
			
	}#End for linear terms
	
	#resCov
	resCov$R<-resCov$post_R
	resCov$SD.R<-sqrt(resCov$post_R2-resCov$post_R^2)
	
	if(resCov$type%in%c("REC","FA"))
	{
		resCov$W<-resCov$post_W
		resCov$SD.W<-sqrt(resCov$post_W2-resCov$post_W^2)
				
		resCov$PSI<-resCov$post_PSI
		resCov$SD.PSI<-sqrt(resCov$post_PSI2-resCov$post_PSI^2)
		
		close(resCov$f_W)
		resCov$f_W<-NULL
		close(resCov$f_PSI)
		resCov$f_PSI<-NULL
		
	}
	
	#Deep cleaning!
	tmp<-which(names(resCov)%in%c("post_R","post_R2","Rinv","post_W","post_W2",
			                      "post_PSI","post_PSI2","F"))
	resCov<-resCov[-tmp]
			    
	rm(tmp)
	
	#Closing files
	
	if(intercept)
	{
		close(f_mu)
		f_mu<-NULL
	}
	
	close(resCov$f_R)
	resCov$f_R<-NULL
	
	#Covariance matrices
	for(j in 1:nLT)
	{
		if(ETA[[j]]$model%in%c("SpikeSlab","BRR","RKHS"))
		{
			if(ETA[[j]]$Cov$type%in%c("REC","FA"))
			{
				close(ETA[[j]]$Cov$f_W)
				ETA[[j]]$Cov$f_W<-NULL
				close(ETA[[j]]$Cov$f_PSI)
				ETA[[j]]$Cov$f_PSI<-NULL
				
			}#End of if REC, FA
			
			if(ETA[[j]]$Cov$type%in%c("UN","DIAG"))
			{
				close(ETA[[j]]$Cov$f_Omega)
				ETA[[j]]$Cov$f_Omega<-NULL
			}
			
		}#End of if SpikeSlab, BRR, RKHS
	}#End of for
	
	#Effect files & indicators files
	for(j in 1:nLT)
	{
		if(!is.null(ETA[[j]]$fileEffects))
		{
        	flush(ETA[[j]]$fileEffects)
            close(ETA[[j]]$fileEffects)
            ETA[[j]]$fileEffects<-NULL
            
            if(!is.null(ETA[[j]]$compressEffects) && ETA[[j]]$compressEffects)
            {
            	message("Compressing binary file for effects for term ", j)
            	compressFile(paste(saveAt,ETA[[j]]$Name,"_beta.bin",sep=""))
            	message("Done")
            }        
        }
        
        if(!is.null(ETA[[j]]$fileIndicators))
        {
        	flush(ETA[[j]]$fileIndicators)
            close(ETA[[j]]$fileIndicators)
            ETA[[j]]$fileIndicators<-NULL
            
            #Compress file by default to save a lot of space
            message("Compressing binary file for indicators for term ", j)
            compressFile(paste(saveAt,ETA[[j]]$Name,"_d.bin",sep=""))
            message("Done")
        }
    }
	
	SD.ETAHat<-sqrt(ETAHat2-ETAHat^2)
	
	#Fit
	fit<-getDIC(y.back, ETAHat, post_logLik, cte, complete_records, resCov$R, 
				missings, missing_records, Dpatterns, Upatterns, 
				dUpatterns,dAllMissings)
	
	fit$postMeanLogLik<-post_logLik
	
	#Return the goodies
	out<-list(ETA=ETA,resCov=resCov,ETAHat=ETAHat,SD.ETAHat=SD.ETAHat,
	          fit=fit)
	          
	if(intercept)
	{
		out$mu<-mu
		out$SD.mu<-SD.mu
	}
	
	if(missings)
	{
		out$missing_records<-missing_records
		out$patterns<-patterns
	}
	
	return(out)
}


