rm(list=ls())

#Fits a GBLUP model
#y: response variable, NAs allowed
#FEffects: formula for fixed effects
#GID: Genotype identifier
#K: Relationship matrix, the names or rows columns should correspond to the GIDs
#data: A data.frame
 
RKHS.GBLUP<-function(y,FEffects=NULL,GID,K,data=NULL,
                     nIter=10000,burnIn=5000,thin=10,verbose=FALSE)
{
	   
	   if(!is.vector(y)) stop("y must be a vector")
	   
	   nRecords<-length(y)
	   
	   
	   XF<-NULL
	   
	   if(is(FEffects,"formula"))
	   {
	   		if(is.null(data))
	   		{
	   			mf <- model.frame(formula=FEffects)
	   			
	   		}else{
	   			
	   			mf <- model.frame(formula=FEffects,data=data)	
	   		}
	   		
	   		XF <- model.matrix(attr(mf, "terms"), data=mf)
    		Xint <- match("(Intercept)", colnames(XF), nomatch=0L)
    		if(Xint > 0L) XF <- XF[, -Xint, drop=FALSE]
	   }
	   
	   if(is.matrix(FEffects))
	   {
	   		warning("BGLR includes by default one intercept, if your matrix as one intercept remove it\n")
	   		XF<-FEffects
	   }
	   
	   if(!is.null(XF))
	   {
	   		if(nrow(XF)!=nRecords) stop("Mismatch between number of rows of matrix of fixed effects and number of phenotypes\n")
	   }
	   
	   if(!is.matrix(K))
	   {
	   		stop("K must be a matrix\n")
	   }
	   
	   if(nrow(K)!=ncol(K))
	   {
	   		stop("K must be a square matrix\n")
	   }
	   
	   if(any(rownames(K)!=colnames(K))) stop("Row/Columns names mismatch in matrix K\n")
	   
	   if(any(duplicated(rownames(K)))) stop("Row names of K with duplicates\n")
	   
	   if(any(duplicated(colnames(K)))) stop("Column names of K with duplicates\n")
	   
	   common<-intersect(unique(as.character(GID)),rownames(K))
	   
	   nGenotypes<-length(common)
	   
	   if(sum(as.character(GID)%in%common)!=nRecords) stop("Some GIDs in phenotypes are not in K\n")
	   
	   if(sum(rownames(K)%in%common)!=nGenotypes) stop("Some GIDs in K not present in the phenotypes\n")
	   
	   #Design matrix for genotypes
	   GID<-as.character(GID)
	   GID<-factor(GID,levels=rownames(K))
	   Zg<-model.matrix(~GID-1)
	   
	   
	   eigen_K<-eigen(K,symmetric=TRUE)
	   index<-eigen_K$values>1e-10
	   eigen_K$vectors<-eigen_K$vectors[,index]
	   eigen_K$values<-eigen_K$values[index]
	   
	   #Zg* = Zg*Gamma*Lambda^0.5
	   Zgstar<-Zg%*%sweep(x=eigen_K$vectors,MARGIN=2,STATS=sqrt(eigen_K$values),FUN="*")
	   
	   if(!is.null(XF))
	   {
	   		ETA<-list(FEffects=list(X=XF,model="FIXED"),
	   		          REffects=list(X=Zgstar,model="BRR"))
	   }else{
	   		ETA<-list(REffects=list(X=Zgstar,model="BRR"))
	   }
	   
	   fm<-BGLR(y=y,ETA=ETA,
	            nIter=nIter, burnIn = burnIn, thin=thin, verbose=verbose)
	   
	   #Getting the mean of random effects and add the results to the output	
	   u<-as.vector(sweep(x=eigen_K$vectors,MARGIN=2,STATS=sqrt(eigen_K$values),FUN="*")%*%fm$ETA$REffects$b)
	   names(u)<-rownames(K)
	   fm$ETA$REffects$u<-u
	   
	   #Adding information about the variance associated to random effect
	   fm$ETA$REffects$varU<-fm$ETA$REffects$varB
	   fm$ETA$REffects$SD.varU<-fm$ETA$REffects$SD.varB
	   
	   return(fm)
}
           
