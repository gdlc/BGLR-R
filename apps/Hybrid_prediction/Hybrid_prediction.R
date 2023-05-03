

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

	eigen_G1<-eigen(G1,symmetric=TRUE)
	index<-eigen_G1$values>1e-10
	eigen_G1$vectors<-eigen_G1$vectors[,index]
	eigen_G1$values<-eigen_G1$values[index]


	eigen_G2<-eigen(G2,symmetric=TRUE)
	index<-eigen_G2$values>1e-10
	eigen_G2$vectors<-eigen_G2$vectors[,index]
	eigen_G2$values<-eigen_G2$values[index]

	eigen_H<-eigen(H,symmetric=TRUE)
	index<-eigen_H$values>1e-10
	eigen_H$vectors<-eigen_H$vectors[,index]
	eigen_H$values<-eigen_H$values[index]


	Z1star<-Z1%*%eigen_G1$vectors%*%sqrt(diag(eigen_G1$values))
	Z2star<-Z2%*%eigen_G2$vectors%*%sqrt(diag(eigen_G2$values))
	Z3star<-Z3%*%eigen_H$vectors%*%sqrt(diag(eigen_H$values))

	#Linear predictor
	ETA<-list(list(X=ZE,model="BRR"),
              list(X=Z1star,model="BRR"),
       	      list(X=Z2star,model="BRR"),
              list(X=Z3star,model="BRR"))
           
	fm<-BGLR(y=y,ETA=ETA,nIter=nIter,burnIn = burnIn,thin=thin,
        verbose=verbose)

	#Getting the mean of random effects and add the results to the output	
	u1<-as.vector(eigen_G1$vectors%*%sqrt(diag(eigen_G1$values))%*%fm$ETA[[2]]$b)
	names(u1)<-rownames(G1)
	fm$ETA[[2]]$u<-u1

	u2<-as.vector(eigen_G2$vectors%*%sqrt(diag(eigen_G2$values))%*%fm$ETA[[3]]$b)
	names(u2)<-rownames(G2)
	fm$ETA[[3]]$u<-u2

	h<-as.vector(eigen_H$vectors%*%sqrt(diag(eigen_H$values))%*%fm$ETA[[4]]$b)
	names(h)<-rownames(H)
	fm$ETA[[4]]$u<-h

	return(fm)
	
}




