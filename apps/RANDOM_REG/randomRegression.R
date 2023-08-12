
RR.set=function(X,Z,...){
    # X could be a matrix with SNPs or a factorization of a GRM such that G=XX'
    # Z can be a matrix, data.frame, or a formula for the covariates for the random regression.
    
    # Making Z a matrix
    if(!is.matrix(Z)){
        if(is.data.frame(Z)){
            Z=as.matrix(Z)
        }else{
            if(class(Z)=='formula'){
              Z=as.matrix(model.matrix(Z,...))
            }else{ 
              stop('Z must be either a matrix, a data.frame, or  a formula')
            }
        }
    }
    
    tmp=apply(FUN=var,X=Z,MARGIN=2)
    hasInt=any(tmp==0)

    if(hasInt){
        colInt=which(tmp==0)
    }else{
        Z=cbind('Int'=1,Z)
        colInt=1
    }

    ETA=list()
    ETA[[1]]=list(X=Z[,-colInt,drop=FALSE],model='FIXED')
    for(i in 1:ncol(Z)){
        W=sweep(X,FUN='*',STAT=Z[,i],MARGIN=1)
        ETA[[i+1]]=list(X=W,...)
    }
    names(ETA)=c('FIXED',colnames(Z))
    return(ETA)
}

RR.coef=function(fm,X){
    nLevels=nrow(X)
    nCov=length(fm$ETA$FIXED$b)
    B=list( fixed=c(fm$mu,fm$ETA$FIXED$b),
            random=matrix(nrow=nLevels,ncol=nCov+1))
   

    for(i in 1:(nCov+1)){
        B$random[,i]=X%*%fm$ETA[[i+1]]$b
        
    }
    colnames(B$random)=names(fm$ETA[-1])
    rownames(B$random)=rownames(X)
    return(B)
}


RR.variances=function(fm){
    n=length(fm$ETA)
    OUT=data.frame(estimate=rep(NA,n),sd=rep(NA,n))
    rownames(OUT)=c(names(fm$ETA)[-1],'Error')

    for(i in 2:length(fm$ETA)){
        OUT$estimate[i-1]=fm$ETA[[i]]$varB
        OUT$sd[i-1]=fm$ETA[[i]]$SD.varB
    }
    
    OUT$estimate[n]=fm$varE
    OUT$sd[n]=fm$SD.varE
    return(OUT)
}

RR.predict=function(fm,X,Z,returnUnique=TRUE){

    # Making Z a matrix
    if(!is.matrix(Z)){
        if(is.data.frame(Z)){
            Z=as.matrix(Z)
        }else{
            if(class(Z)=='formula'){
              Z=as.matrix(model.matrix(Z,...))
            }else{ 
              stop('Z must be either a matrix, a data.frame, or  a formula')
            }
        }
    }
    
    tmp=apply(FUN=var,X=Z,MARGIN=2)
    hasInt=any(tmp==0)

    if(hasInt){
        colInt=which(tmp==0)
        X=cbind('Int'=1,Z[,-colInt])
    }else{
        Z=cbind('Int'=1,Z)
        colInt=1
    }

    tmp=RR.coef(fm,X)

    B=tmp$random
    for (i in 1:ncol(B)){
        B[,i]=B[,i]+tmp$fixed[i]
    }
    YHat=t(tcrossprod(Z,B))

    return(YHat)

}

