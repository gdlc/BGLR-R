#Function to fit Finlay-Wilkinson model in 2 steps using BGLR
#The model is:
# y[ijk] = mu + g[i] + h[j] + b[i] * h[j] + e[ijk]
#i indexed genotypes, j indexes environments and k indexes replicate.
# We assume g ~ N(0, sigma^2_g G), b ~ N(0, sigma^2_ G) 
#
#Function arguments
#pheno: A data.frame with 3 columns
# $VAR: variety
# $ENV: environment
# $y: response variable
# X: matrix of markers
# G: Genomic relationship matrix


FW.BGLR<-function(pheno,X=NULL,G=NULL,model1='BRR',model2='BRR',saveAt='',...){

    varieties<-unique(as.character(pheno$VAR))
    environments<-unique(as.character(pheno$ENV))
    common<-intersect(rownames(G),varieties)
    pheno<-pheno[pheno$VAR%in%common,]
    
    if(!is.null(G)){
    	tmp<-rownames(G)%in%common
    	G<-G[tmp,tmp]
    }
    
    if(!is.null(X)){
   	tmp<-rownames(X)%in%common
    	X<-X[tmp,]
    	X<-scale(X,center=TRUE,scale=FALSE)
		G<-tcrossprod(X)/ncol(X)
    }
    
    #Incidence matrix for environments
    pheno$ENV<-as.factor(pheno$ENV)
    ZE<-model.matrix(~pheno$ENV-1)
    ZE<-scale(ZE,center=T,scale=F) 

    #Incidenve matrix for main effect of Genotypes
    pheno$VAR<-factor(x=pheno$VAR,levels=rownames(G),ordered=TRUE)
    ZVAR<-model.matrix(~pheno$VAR-1)

    if(!is.null(G))
    {
	#Eigen-value decomposition for G
	evdG<-eigen(G)

	PC<-sweep(evdG$vector,STATS=sqrt(evdG$values),MARGIN=2,FUN="*")
	 XVAR<-ZVAR%*%PC
	 colnames(XVAR)=colnames(ZVAR)
    }

    
    #Finlay-Wilkinson
    #Step 1 Compute environmental means
    
    eta1<-list(list(X=ZE,model="BRR"),#*# centering ZE should help convergence
               list(X=XVAR,model=model1))
           
           
    fm1<-BGLR(y=pheno$y,ETA=eta1,saveAt=paste0(saveAt,'step1_'),...)
    hHat1<-fm1$ETA[[1]]$b #*# Why don't we add intercept? #because these are interpreted as deviations from general mean
                                                          #due to environmental effects?
    names(hHat1)<-levels(pheno$ENV)

    #Step 2
    hHatExpanded<-as.vector(ZE%*%hHat1)
    yc<-pheno$y-hHatExpanded # why not centering yc?      #Once that we estimated the environmental effect we fit the model
                                                          #y[ijk]-hat(h[j]) = mu + g[i] + b[i]*hat(h[j]) + e[ijk] 

    #Weight each row in ZVAR
    XInteraction<-sweep(ZVAR,1L,hHatExpanded,"*")
    XInteraction<-XInteraction%*%PC

    eta2<-list( int=list(X=XVAR,model=model2,saveEffects=TRUE),
	       		slope=list( X=XInteraction,model=model2,saveEffects=TRUE))

	saveAt=paste0(saveAt,'step2_')
    fm2<-BGLR(y=yc,ETA=eta2,saveAt=saveAt,...)
	
	BInt=readBinMat(paste0(saveAt,'ETA_int_b.bin'))
	BSlope=readBinMat(paste0(saveAt,'ETA_slope_b.bin'))
	
	
	INT=tcrossprod(PC,BInt)
	SLOPE=tcrossprod(PC,BSlope)
	
    INT<-INT + mean(pheno$y,na.rm=TRUE)
 	SLOPE<-SLOPE+1   
    rownames(INT)<-levels(pheno$VAR)
    rownames(SLOPE)<-levels(pheno$VAR)

	Int=rowMeans(INT)
	Slope=rowMeans(SLOPE)
    yHat=fm2$mu+ZE%*%hHat1+ZVAR%*%(Int-mean(pheno$y,na.rm=TRUE))+XInteraction%*%Slope
    
    #return the goodies
    out<-list(yHat=yHat,VAR=data.frame( ID=rownames(INT),int=Int,intSD=apply(FUN=sd,X=INT,MARGIN=1),
    									slope=Slope,slopeSD=apply(FUN=sd,X=SLOPE,MARGIN=1),stringsAsFactors=FALSE),ENV=hHat1)
    
    return(out)
}


plot.FW <- function(out,pheno)
{
 	## set par settings
 	par(mar=c(5,5,2,2))

  	xlim <- range(out$ENV)
 	ylim <- range(out$VAR$int)
 	ylim[1]<-ylim[1]-0.35*ylim[1]
 	ylim[2]<-1.35*ylim[2]

  	## draw an empty plot
 	plot(as.numeric()~1, ylab = "Genotype performance", 
	     xlab = "Environment effect", cex.axis=1.2, cex.lab=1.5, xlim=xlim, ylim=ylim)

  	## draw regression lines for all the genotypes
 	for(i in 1:nrow(out$VAR))
 	{
 		#we can add intercept to match scale of original data if needed
 		abline(a=out$VAR$int[i],b=out$VAR$slope[i],col="lightgrey")
 	}

  	## add vertical line @ hHat=0 to see G effect
 	abline(v=0, lty=2)

  	## add dashed line for b=1
 	abline(mean(pheno$y, na.rm=TRUE), 1, lty=2, lwd=1)

  	## select genotypes to plot
 	selgeno <- out$VAR$slope
 	names(selgeno) <- rownames(out$VAR)
 	selgeno <- selgeno[order(selgeno)]
 
 	## plot genotypes ith min/close to 1/max slope
 	selgeno <- c(selgeno[1], selgeno[floor(length(selgeno)/2)], selgeno[length(selgeno)])

  	mycols <- c("blue", "violetred", "orange")
 	count <- 1
 	for (i in names(selgeno))
	{
 		points(out$ENV[pheno$ENV[pheno$VAR == i]], 
		       pheno$y[pheno$VAR == i], 
		       pch=19, col=mycols[count])
 		abline(a=out$VAR$int[rownames(out$VAR) == i],
		       b=out$VAR$slope[ rownames(out$VAR) == i], 
		       col=mycols[count], lwd=2)
 		count <- count+1
 	}

  	## add legend
 	z <- match(names(selgeno), rownames(out$VAR))
 	mylegend <- paste0(names(selgeno), ": ", "g=", round(out$VAR$int[z], 2), "; b=", round(out$VAR$slope[z], 2))
 	legend("topleft", bty="n", cex=1.2, legend=mylegend, col=mycols, pch=rep(19, length(selgeno)))

  	## just to get nice and neat box
 	box()

  	## reset par settings
 	par(mar=c(5,4,4,2))
}
