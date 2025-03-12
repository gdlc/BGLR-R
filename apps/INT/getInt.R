getInt=function(X,eID,gID=rep(1:nrow(X)),sparse=TRUE,centerX=FALSE,...){
	stopifnot(all(GID%in%rownames(X)))
	ETA=list()
	X=X[gID,]
	if(centerX){
		X=scale(X,center=TRUE,scale=FALSE)
	}
	# Main effects
	 ETA[[1]]=list(X=X,...)

	# Interacctions
	levels=unique(eID)
	nLevels=length(levels)
	for(i in 1:nLevels){
		Z=X
		Z[eID!=levels[i],]=0
		ETA[[i+1]]=list(X=Z,...)
		if(sparse){
			ETA[[i+1]]$X=as(ETA[[i+1]]$X,"CsparseMatrix") 
			ETA[[i+1]]$model=paste0(ETA[[i+1]]$model,'_sparse')
		}
	}
	names(ETA)=c('main',paste0('int_',levels))
	W=as.matrix(model.matrix(~factor(eID)))[,-1]
        ETA=c(meanDif=list(X=W,model='FIXED'),ETA)
	return(ETA)
}
