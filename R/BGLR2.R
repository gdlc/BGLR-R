
BGLR2=function (y, response_type = "gaussian", a = NULL, b = NULL, 
    ETA = NULL, nIter = 1500, burnIn = 500, thin = 5, saveAt = "", 
    S0 = NULL, df0 = 5, R2 = 0.5, weights = NULL, 
    verbose = TRUE, rmExistingFiles = TRUE, groups=NULL,
    saveEnv=FALSE,BGLR_ENV=NULL,newChain=TRUE) {
   
  if(verbose){welcome()}
    
  if(is.null(BGLR_ENV)){  #*#
    GIBBS_start=1 #*#
    
    IDs=names(y)
    if (!(response_type %in% c("gaussian", "ordinal")))  stop(" Only gaussian and ordinal responses are allowed\n")

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
		if(sum(countGroups)!=n) stop("length of groups and y differs, NA's not allowed in groups\n");	
    }

    if(response_type=="ordinal")
    {

    	y=factor(y,ordered=TRUE)
        lev=levels(y)
        nclass=length(lev)
        if(nclass==n) stop("The number of classes in y must be smaller than the number of observations\n");

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
            if ((length(a) != n) | (length(b) != n)) stop(" y, a and b must have the same dimension\n")
            if (any(weights != 1)) stop(" Weights are only implemented for Gausian uncensored responses\n")
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
        if(verbose) {cat(" Note: samples will be appended to existing files. \n") }
    }

    fileOutMu = file(description = fname, open = "w")

    if (response_type == "ordinal") {
        if(verbose){ cat(" Prior for residual is not necessary, if you provided it, it will be ignored\n")}
        if (any(weights != 1)) stop(" Weights are not supported \n")
       
        countsZ=table(z)

        if (nclass <= 1) stop(paste(" Data vector y has only ", nclass, " differente values, it should have at least 2 different values\n"))
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
                stop(paste(" Error in ETA[[", i, "]]", " model ", ETA[[i]]$model, " not implemented (note: evaluation is case sensitive).", sep = ""))
                
            }

            if(!is.null(groups))
            {
		if(!(ETA[[i]]$model %in%  c("BRR","FIXED","BayesB","BayesC"))) stop(paste(" Error in ETA[[", i, "]]", " model ", ETA[[i]]$model, " not implemented for groups\n", sep = ""))
            }

            if(!is.null(ETA[[i]]$sizeEffects)){
                if(!(ETA[[i]]$sizeEffects %in% c(4L,8L))) stop("Error in ETA[[", i, "]]", " sizeEffects can be either '4' (single) or '8' (double, default)")
            }else{
                ETA[[i]]$sizeEffects=8L
            }

            ETA[[i]] = switch(ETA[[i]]$model, 
			      FIXED = setLT.Fixed(LT = ETA[[i]],  n = n, j = i, weights = weights, y = y, nLT = nLT, saveAt = saveAt, rmExistingFiles = rmExistingFiles,groups=groups,nGroups=nGroups), 
                              BRR = setLT.BRR(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles,groups=groups,nGroups=nGroups,verbose=verbose,thin=thin,nIter=nIter,burnIn=burnIn),#*# 
                              BL = setLT.BL(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles,verbose=verbose,thin=thin,nIter=nIter,burnIn=burnIn), 
                              RKHS = setLT.RKHS(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles,verbose=verbose), 
                              BayesC = setLT.BayesBandC(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles,groups=groups,nGroups=nGroups,verbose=verbose,thin=thin,nIter=nIter,burnIn=burnIn),
                              BayesA = setLT.BayesA(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles,verbose=verbose,thin=thin,nIter=nIter,burnIn=burnIn),
                              BayesB = setLT.BayesBandC(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles,groups=groups,nGroups=nGroups,verbose=verbose,thin=thin,nIter=nIter,burnIn=burnIn),
                              BRR_sets = setLT.BRR_sets(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles,verbose=verbose,thin=thin,nIter=nIter,burnIn=burnIn)
                              )
        }
    }

    }else{ #*# Block of code for the case when the environment is re-loaded
        callParameters=list(nIter=nIter, burnIn=burnIn, thin=thin, saveAt=saveAt, 
        			verbose=verbose,rmExistingFiles=rmExistingFiles,newChain=newChain)

    	load(BGLR_ENV)
    	
    	# restoring some call parameters
    	newChain=callParameters$newChain
    	if(newChain){
    		nIter=callParameters$nIter
  		GIBBS_start=1
  		saveAt=callParameters$saveAt
    	}else{
    		GIBBS_start=nIter+1
    		nIter=GIBBS_start+callParameters$nIter-1
    	}

    	# Restore call parameters
    	 burnIn=callParameters$burnIn
    	 thin=callParameters$thin
    	 verbose=callParameters$verbose
    	 rmExisitingFiles=callParameters$rmExistingFiles
    	 rm(callParameters)
  
    	# Reseting Running means and connections
    	if(newChain){

 		nSums=0   		
    		# Running means
	    	tmp=ls(pattern='post_')
		for(i in 1:length(tmp)){
			eval(parse(text=paste(tmp[i],'[]<-0')))
		}
		
		for(i in 1:length(ETA)){
			tmp=names(ETA[[i]])[grep(names(ETA[[i]]),pattern='post_')]
			for(j in 1:length(tmp)){
				eval(parse(text=paste0('ETA[[i]]$',tmp[j],"[]<-0")))
			}
		}
    	
    		# Resets Connections
    		fname = paste(saveAt, "mu.dat", sep = "")

    		if (rmExistingFiles) {
        		unlink(fname)
    		}
		fileOutMu = file(description = fname, open = "w")
    		
    		fname = paste(saveAt, "varE.dat", sep = "")
    		if (rmExistingFiles) {
        		unlink(fname)
    		}
		fileOutVarE = file(description = fname, open = "w")    		

    		if(response_type=="ordinal"){
			fname = paste(saveAt, "thresholds.dat", sep = "")
    			if (rmExistingFiles) {
        			unlink(fname)
    			}			
        		fileOutThresholds = file(description = fname, open = "w")
    		}
		for(i in 1:length(ETA)){        	
          		fname=paste0("\'",saveAt, basename(normalizePath(ETA[[i]]$NamefileOut)),"\'")
          		ETA[[i]]$fileOut=fname
			if(rmExistingFiles){ 
       				unlink(fname) 
    			}
    			ETA[[i]]$fileOut=file(description=fname,open="w")
    			
    			if(ETA[[i]]$model=='FIXED'){
    				tmp=ETA[[i]]$colNames
    				write(tmp, ncolumns = ETA[[i]]$p, file = ETA[[i]]$fileOut, append = TRUE)
          		}else{
          			if(ETA[[i]]$saveEffects){
    					fname=paste(saveAt,ETA[[i]]$Name,"_b.bin",sep="")
    					if(rmExistingFiles){ unlink(fname) }
    					ETA[[i]]$fileEffects=file(fname,open='wb')
    					nRow=floor((nIter-burnIn)/thin)
    					writeBin(object=c(nRow,ETA[[i]]$p),con=ETA[[i]]$fileEffects,size=ETA[[i]]$sizeEffects)
    				}
    			}
         	}
    	}else{
    		# if we just continue the chain, we re-open connections in mode 'append'
    		fname = paste(saveAt, "mu.dat", sep = "")
	    	fileOutMu = file(description = fname, open = "a")

    		fname = paste(saveAt, "varE.dat", sep = "")
	    	fileOutVarE = file(description = fname, open = "a")    		
    
    		if(response_type=="ordinal"){
			fname = paste(saveAt, "thresholds.dat", sep = "")
        		fileOutThresholds = file(description = fname, open = "a")
    		}
    		
    		for(i in 1:length(ETA)){
    			fname=switch(ETA[[i]]$model, 
    				FIXED=paste(saveAt,ETA[[i]]$Name,"_b.dat",sep=""),
    				BRR=paste(saveAt,ETA[[i]]$Name,"_varB.dat",sep=""),
    				BRR_sets=paste(saveAt,ETA[[i]]$Name,"_varB.dat",sep=""),
    				BL=paste(saveAt,ETA[[i]]$Name,"_varB.dat",sep=""),
    				RKHS=paste(saveAt,ETA[[i]]$Name,"_varU.dat",sep=""),
    				BayesA=paste(saveAt,ETA[[i]]$Name,"_ScaleBayesA.dat",sep=""),
    				BayesB=paste(saveAt,ETA[[i]]$Name,"_parBayesB.dat",sep=""),    				
    				BayesC=paste(saveAt,ETA[[i]]$Name,"_parBayesC.dat",sep=""),   
    			      )
    			ETA[[i]]$NamefileOut=fname
			if(rmExistingFiles){  unlink(fname)  }
			ETA[[i]]$fileOut=file(description=fname,open="a")
			
			if(ETA[[i]]$saveEffects){
				fname=paste(saveAt,ETA[[i]]$Name,"_b.bin",sep="")
    				ETA[[i]]$fileEffects=file(fname,open='ab')
			}
		}
	}
    }
    # Gibbs sampler

    time = proc.time()[3]

       	# Restore seed
    if(!newChain){ .Random.seed=seed }
    
    for (i in GIBBS_start:nIter) {
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
                  	ans = .Call("sample_beta", n, ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, 
                                             e, varBj, varE, 1e-9)
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
				   
                   tmp=SS/rchisq(df=ETA[[j]]$DF1,n=ETA[[j]]$n_sets)
              
                   ETA[[j]]$varB=tmp[ETA[[j]]$sets]
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
                          writeBin(object=ETA[[j]]$b,con=ETA[[j]]$fileEffects,size=ETA[[j]]$sizeEffects)
                      }#*#
                    }

                    if (ETA[[j]]$model == "BRR_sets") {
                      ETA[[j]]$post_b = ETA[[j]]$post_b * k + ETA[[j]]$b/nSums
                      ETA[[j]]$post_b2 = ETA[[j]]$post_b2 * k + (ETA[[j]]$b^2)/nSums
                      ETA[[j]]$post_varB = ETA[[j]]$post_varB * k + (ETA[[j]]$varB)/nSums
                      ETA[[j]]$post_varB2 = ETA[[j]]$post_varB2 * k + (ETA[[j]]$varB^2)/nSums
                      ETA[[j]]$post_varSets<-ETA[[j]]$post_varSets*k+tmp/nSums
                      ETA[[j]]$post_varSets2<-ETA[[j]]$post_varSets2*k+(tmp^2)/nSums
                      if(ETA[[j]]$saveEffects&&(i%%ETA[[j]]$thin)==0){
                          writeBin(object=ETA[[j]]$b,con=ETA[[j]]$fileEffects,size=ETA[[j]]$sizeEffects)
                      }#*#
                    }
                    
                    if (ETA[[j]]$model == "BL") {
                      ETA[[j]]$post_b = ETA[[j]]$post_b * k + ETA[[j]]$b/nSums
                      ETA[[j]]$post_b2 = ETA[[j]]$post_b2 * k + (ETA[[j]]$b^2)/nSums
                      ETA[[j]]$post_tau2 = ETA[[j]]$post_tau2 * k + (ETA[[j]]$tau2)/nSums
                      ETA[[j]]$post_lambda = ETA[[j]]$post_lambda * k + (ETA[[j]]$lambda)/nSums
                      if(ETA[[j]]$saveEffects&&(i%%ETA[[j]]$thin)==0){
                          writeBin(object=ETA[[j]]$b,con=ETA[[j]]$fileEffects,size=ETA[[j]]$sizeEffects)
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
                          writeBin(object=ETA[[j]]$b*ETA[[j]]$d,con=ETA[[j]]$fileEffects,size=ETA[[j]]$sizeEffects)
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
                          writeBin(object=ETA[[j]]$b,con=ETA[[j]]$fileEffects,size=ETA[[j]]$sizeEffects)
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
                            writeBin(object=ETA[[j]]$b*ETA[[j]]$d,con=ETA[[j]]$fileEffects,size=ETA[[j]]$sizeEffects)
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
            cat("---------------------------------------\n")
            tmp = proc.time()[3]
            cat(c(paste(c("  Iter=", "Time/Iter="), round(c(i, c(tmp - time)), 3), sep = "")), "\n")
            cat("  VarE=",round(varE,3),"\n")
            time = tmp
        }
    }#end of Gibbs sampler
    seed=.Random.seed
    
    #Closing files
    close(fileOutVarE)
    close(fileOutMu)
    if(response_type == "ordinal") close(fileOutThresholds)

    if (nLT > 0) {
        for (i in 1:nLT) {
            if (!is.null(ETA[[i]]$fileOut)) {
                close(ETA[[i]]$fileOut)
                ETA[[i]]$fileOut = NULL
            }
            if(!is.null(ETA[[i]]$fileEffects)){
            		close(ETA[[i]]$fileEffects)
            		ETA[[i]]$fileEffects = NULL
            }
            
        }
    }
    
    if(saveEnv){
    	save(list=ls(),file=paste0(saveAt,'BGLR_ENV.RData'))
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

