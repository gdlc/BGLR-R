###

 # A function to find segments with elevated PIP
 CS.SEGMENTS=function(B,setBFDR,chr,pos,maxD,minProb=0.05){
   if(!is.logical(B[1,1])){
   	B=B!=0
   }
   PIP=colMeans(B)
   LFDR=1-PIP
   
   DS=segments(LFDR,chr=chr,bp=pos,threshold=1-minProb,gap=maxD)
   DS$setPIP=NA
   for(i in 1:nrow(DS)){
		DS$setPIP[i]=mean(apply(FUN=any,X=B[,DS$start[i]:DS$end[i]]!=0,MARGIN=1))
   }
   DS$BFDR=cumsum(1-DS$setPIP)/(1:nrow(DS))
   DS=DS[DS$BFDR<=setBFDR,]
   return(DS)
}

###
updateSamples=function(B0,CS){
     if(length(CS)>0){
        p=ncol(B0)
        B=B0
        tmp=apply(FUN=any,X=B[,CS,drop=FALSE],MARGIN=1)
        CS_FALSE=!tmp
        for(i in (1:p)){
                B[,i]=B[,i]&(CS_FALSE)
        }
        return(B)
    }else{
      return(B0)
    }
}



## This function is similar to findCS_OLD() but it restrictes the search to SNPs within
## a maximum distance to the leading variant (i.e., the one with highest posterior probability of inclussion)
## The function is meant to be called by findCS(), see function below

nextCS=function(B,minProbIn=.05,maxSize=min(ncol(B),10),
                 maxSetProb=.98,maxD=100,bp=1:ncol(B), verbose=FALSE){
  if(is.null(colnames(B))){
    colnames(B)=1:nrow(B)	  
  }
  colNames=colnames(B)
  CS=0

  PROB=0

  nActive=0
  p=ncol(B)

  prob_in=colMeans(B)

  ready=FALSE
  B0=B

  cycle=1

  while(!ready){

     B=updateSamples(B0,CS)

     if(cycle==1){
             prob_in=colMeans(B)
     }else{
        if(cycle==2){
            searchSet=which((abs(bp[CS[2]]-bp)<=maxD))
                prob_in=rep(0,ncol(B))
        }

        prob_in[searchSet]=colMeans(B[,searchSet,drop=F])

     }

     if(any(prob_in>=minProbIn)){
        tmp=which.max(prob_in)
        if(length(tmp)>1){
          tmp=tmp[1] # could be random
        }
        CS=c(CS,tmp)
        PROB=c(PROB,max(PROB)+prob_in[tmp])

        nActive=length(CS)-1
     }
     cycle=cycle+1
     ready=(all(prob_in<minProbIn)|(nActive>=maxSize)|(max(PROB)>=maxSetProb))
  }

  return(cbind('SNPs'=colNames[CS[-1]],'cumProb'=PROB[-1]))
}

## This function finds CS within chromosome by recursively calling nextCS

 findCS=function(B,minProbIn,maxSize,minSetProb,maxSetProb,maxD,chr,bp,verbose=FALSE){
 
   if(!is.logical(B)){
      B=(B!=0)
   }
 
 
    SETS=list()
    
    nChr=length(unique(chr))
    count_chr=1
    for(i in unique(chr)){
        ready=FALSE
	    count_sets=1
	    
		SETS[[count_chr]]=list()
		
		
		
		    tmp<-(chr==i)
		    B_CHR=B[,tmp]
		while(!ready){
   			TMP=nextCS(B=B_CHR,minProbIn=minProbIn,maxSize=maxSize,maxSetProb=maxSetProb,maxD=maxD,bp=bp[tmp],verbose=verbose)
   			ready=ifelse(nrow(TMP)==0,TRUE,max(TMP[,2,drop=FALSE])<minSetProb)

   			if(!ready){
        		B_CHR[,TMP[,1]]=FALSE
        		SETS[[count_chr]][[count_sets]]=cbind('set'=count_sets,TMP)
        		count_sets=count_sets+1
   			}
  
 		}
 		count_chr=count_chr+1
 		if(verbose){ 
 			 message('==> chr ',i, ' done','') 
 		}
 	
 	}
 	names(SETS)=unique(chr)
 	return(SETS)
}




####################################
## Old code
####################################
if(FALSE){
 ### now replaced with nextCS()
findCS_OLD=function(B,lfdr=.01,maxSize=min(ncol(B),10),maxProb=1-lfdr){

  CS=0

  PROB=0

  nActive=0
  p=ncol(B)


  prob_in=colMeans(B)

  ready=FALSE
  B0=B

  while(!ready){

     B=updateSamples(B0,CS)
     prob_in=colMeans(B)

     if(any(prob_in>=lfdr)){
        tmp=which.max(prob_in)
        if(length(tmp)>1){
                tmp=tmp[1] # could be random
        }
        CS=c(CS,tmp)
        PROB=c(PROB,max(PROB)+prob_in[tmp])

        nActive=length(CS)-1

        message('========' ,length(CS)-1 ,'==========')
        print(cbind(CS[-1],PROB[-1]))
     }
     ready=(all(prob_in<lfdr)|(nActive>=maxSize)|(max(PROB)>=maxProb))
  }
  message('========' ,'Done!','==========')

  return(cbind(CS[-1],PROB[-1]))
}

    
}


## Backward elimination


dropOne0=function(B,j){
	mean(apply(X=B[,-j,drop=FALSE],MARGIN=1,FUN=any))
}
dropOne=function(B,set=1:ncol(B)){
	probs=unlist(lapply(FUN=dropOne0,X=set,B=B))
	return(which.min(probs))
}


BKW=function(B,minProb=.8){
    p=ncol(B)
	DS=1:p
	RS=integer()
		
	colnames(B)=1:p
	
	p=ncol(B)
	
	probs=mean(apply(X=B,MARGIN=1,FUN=any))
	ready=probs < minProb
	counter=0
	
	while(!ready){
		counter=counter+1		
		tmp=dropOne(B)
		DS=DS[DS!=as.integer(colnames(B)[tmp])]
		RS=c(RS,as.integer(colnames(B)[tmp]))
		B=B[,-tmp,drop=FALSE]
	    	probs=c(probs,mean(apply(X=B,MARGIN=1,FUN=any)))
		ready= (min(probs) < minProb)| ncol(B)==1
		if(ready & (length(RS)>0)){
		    tmp=length(RS)
		    DS=c(DS,RS[tmp])
		    RS=RS[-tmp]
		    probs=probs[-length(probs)]
		}
	}
	return(list(DS=DS,RS=RS,prob=probs[-1]))
}


#set.seed(122345)
#B=matrix(nrow=100,ncol=10,runif(1000)>.5)
#BKW(B)

