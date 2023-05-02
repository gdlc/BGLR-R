

RKHS.Hybrid_prediction<-function(y, location, id1, id2, idH,
								 G1, G2, H, 
                                 nIter=10000,burnIn=5000,thin=10,
                                 verbose=FALSE)
{
	location<-as.character(location)
	location<-factor(location)
	
	ZE<-model.matrix(~location-1)
	
	id1<-as.character(id1)
	id1<-factor(id1,levels=rownames(G1))
	Z1<-model.matrix(~id1-1)
	
	id2<-as.character(id2)
	id2<-factor(id2,levels=rownames(G2))
	Z2<-model.matrix(~id2-1)
	
	idH<-as.character(idH)
	idH<-factor(idH,levels=rownames(H))
	Z3<-model.matrix(~idH-1)

	G1star<-Z1%*%G1%*%t(Z1)
	G2star<-Z2%*%G2%*%t(Z2)
	Hstar<-Z3%*%H%*%t(Z3)

	ETA<-list(list(X=ZE,model="FIXED"),
	     list(K=G2star,model="RKHS"),
         list(K=G2star,model="RKHS"),
         list(K=Hstar,model="RKHS"))

	fm<-BGLR(y=y,ETA=ETA,nIter=nIter,burnIn = burnIn,thin=thin,
        verbose=verbose)

}

