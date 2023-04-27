RKHS.Groups=function(y,K,group,...){
 ## Group-specific intercepts
  Z=as.matrix(model.matrix(~factor(group)))[,-1,drop=FALSE]
  ETA=list(int=list(X=Z,model='FIXED'))
 
 ## Main effect
  EVD0=eigen(K,symmetric=TRUE)
  PC0=sweep(x=EVD0$vectors[,EVD0$values>1e-8],STATS=sqrt(EVD0$values[EVD0$values>1e-8]),FUN='*',MARGIN=2)
 
  ETA[['main']]=list(X=PC0,model='BRR',saveEffects=TRUE)
  groups=unique(group)
  nGroups=length(groups)
 
  ## Interactions
  n=nrow(K)
  for(i in 1:nGroups){
 	 tmp=which(group==groups[i])
 	 EVD=eigen(K[tmp,tmp,drop=FALSE],symmetric=TRUE)
 	 PC=sweep(x=EVD$vectors[,EVD$values>1e-8],STATS=sqrt(EVD$values[EVD$values>1e-8]),FUN='*',MARGIN=2L)
 	
     X=matrix(nrow=n,ncol=ncol(PC),0)
 	 X[group==groups[i],]=PC
 	 ETA[[paste0('group_',groups[i])]]=list(X=X,model='BRR',saveEffects=TRUE)
 }
 
 tmp=paste0(tempfile(),'_')
 fm=BGLR(y=y,ETA=ETA,group=group,saveAt=tmp)# group=group accommodates group-specific error variances
 
 ## Recovering breeding values from PC and effects
 for(i in 1:nGroups){
 	tmp=paste0('group_',groups[i])
 	fm$ETA[[tmp]]$u=ETA[[tmp]]$X%*%fm$ETA[[tmp]]$b
 }
 
 return(fm)
}
