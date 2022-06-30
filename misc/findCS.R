

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

  return(cbind('SNPs'=CS[-1],'cumProb'=PROB[-1]))
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

