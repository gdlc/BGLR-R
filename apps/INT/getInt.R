
getInt=function(X,z,sparse=TRUE,centerX=FALSE,...){
	ETA=list()
	if(centerX){
		X=scale(X,center=TRUE,scale=FALSE)
	}

	# Main effects
	 ETA[[1]]=list(X=X,...)

	# Interacctions
	levels=unique(z)
	nLevels=length(levels)
	for(i in 1:nLevels){
		Z=X
		Z[z!=levels[i],]=0

		ETA[[i+1]]=list(X=Z,...)

		if(sparse){
			ETA[[i+1]]$X=as(ETA[[i+1]]$X,"CsparseMatrix") 
			ETA[[i+1]]$model=paste0(ETA[[i+1]]$model,'_sparse')
		}
	}
	names(ETA)=c('main',paste0('int_',levels))
	return(ETA)
}
