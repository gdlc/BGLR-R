

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


###

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


## This function is similar to findCS() but it restrictes the search to SNPs within
## a maximum distance to the leading variant (i.e., the one with highest posterior probability of inclussion)

nextCS=function(B,minProbIn=.05,maxSize=min(ncol(B),10),
                 maxSetProb=.98,maxD=100,bp=1:ncol(B)){

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

        message('========' ,length(CS)-1 ,'==========')
        print(cbind(CS[-1],PROB[-1]))
     }

     cycle=cycle+1

     ready=(all(prob_in<minProbIn)|(nActive>=maxSize)|(max(PROB)>=maxSetProb))
  }
  message('========' ,'Done!','==========')

  return(cbind(CS[-1],PROB[-1]))
}

## This function finds CS recursively


