getInt=function(X,eID,gID=rep(1:nrow(X)),sparse=FALSE,centerX=FALSE,...){
	stopifnot(all(gID%in%rownames(X)))
	
	X=X[gID,]
	if(centerX){
		X=scale(X,center=TRUE,scale=FALSE)
	}
	
	ETA=list()
	# Means
	 W=as.matrix(model.matrix(~factor(eID))[,-1])
         ETA[[1]]=list(X=W,model='FIXED')
	
	# Main effects
	 ETA[[2]]=list(X=X,...)

	# Interacctions
	levels=unique(eID)
	nLevels=length(levels)
	for(i in 1:nLevels){
		Z=X
		Z[eID!=levels[i],]=0
		ETA[[i+2]]=list(X=Z,...)
		if(sparse){
			ETA[[i+2]]$X=as(ETA[[i+2]]$X,"CsparseMatrix") 
			ETA[[i+2]]$model=paste0(ETA[[i+2]]$model,'_sparse')
		}
	}
	names(ETA)=c('meansDif','main',paste0('int_',levels))

	return(ETA)
}
