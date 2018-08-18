
 readBinMat=function(filename,byrow=TRUE){
 	## Function to read effects saved by BGLR when ETA[[j]]$saveEffects=TRUE
  	fileIn=file(filename,open='rb')
 	n=readBin(fileIn,n=1,what=numeric())
 	p=readBin(fileIn,n=1,what=numeric())
 	tmp=readBin(fileIn,n=(n*p),what=numeric())
 	X=matrix(data=tmp,nrow=n,ncol=p,byrow=byrow)
 	close(fileIn)
 	return(X)
 }
 
 writeBinMat=function(x,file,byrow=T){
    n=nrow(x)
    p=ncol(x)
    x=as.vector(x)
    fileOut<-file(file,open='rb')
    writeBin(object=n,con=fileOut)
    writeBin(object=p,con=fileOut)
    writeBin(object=x,con=fileOut)
    close(fileOut)
 }
 
getVariances<-function(X,B,sets,verbose=TRUE)
{

	nSets=length(unique(sets))
	n=nrow(X)
	setLabels=unique(sets)
	nIter=nrow(B)
	VAR=matrix(nrow=nIter,ncol=(nSets+1),NA)
	colnames(VAR)<-c(setLabels,'total')
	XList=list()

	for(i in 1:nSets)
	{
		index<-sets==setLabels[i]
		XList[[i]]<-list(index=index,X=X[,index,drop=F])
	}

	rm(X)

	for(i in 1:nIter)
	{
		
		yHat<-rep(0,n)

		for(j in 1:nSets)
		{
			uHat<-XList[[j]]$X%*%as.vector(B[i,XList[[j]]$index])
			VAR[i,j]<-var(uHat)
			yHat=yHat+uHat
		}

		if(verbose){ cat(' Working iteration',i,'(out of',nIter,')\n')}
		VAR[i,(nSets+1)]<-var(yHat)
	}
	return(VAR)
}
 


#
#Rotines for Plink support
#http://pngu.mgh.harvard.edu/~purcell/plink/
#

#The documentation for the C functions can be found in src/util_plink.c

#This function will read a bed file (binary file for genotypes in plink format)
#Arguments: 
#bed_file: Name for bed file
#bim_file: Name for bim file
#fam_file: Name for fam file
#na_strings: Missing value indicator
#verbose: Logical, it TRUE it will print information generated when reading the bed file.
#it returns a list with 3 components, n: number of individuals, p: number of markers, 
#x, integer vector of dimensions n*p with genotypic information.
#see demo/read_bed.R for an example
read_bed=function(bed_file,bim_file,fam_file,na.strings=c("0","-9"),verbose=FALSE)
{
	#Extended map file (this gives the number of snps)
	bim=read.table(bim_file, comment.char="", as.is=TRUE, na.strings=na.strings)
	snps=as.character(bim[,2])

	#First 6 columns of ped file (this gives the number of individuals)
	fam=read.table(fam_file, comment.char="", as.is=TRUE, na.strings=na.strings)

	n=nrow(fam)
	p=length(snps)

	out=rep(0,n*p)

        if(verbose)
        {
           verbose=1
        }else 
        {   
           verbose=0
        }

	out=.C("read_bed_",as.character(bed_file),as.integer(n),as.integer(p),as.integer(out),as.integer(verbose))[[4]]
        return(list(n=n,p=p,x=out))
}

#This function will read a ped file
#Note: It assumes that the missing value is 0
#It returns a list with 3 components, n: number of individuals, p: number of markers, 
#x, integer vector of dommensions n*p with genotypic information.
#see demo/read_ped.R for an example
read_ped=function(ped_file)
{
	out=.Call("read_ped_",ped_file)
	return(out)
}

#This function will write a bed file (binary file for genotypes in plink format)
#x: vector with genotypic information
#n: number of individuals
#p: number of markers
#bed_file: Output file in bed format
#See demo/write_bed.R for an example
write_bed=function(x,n,p,bed_file)
{
   	#Check inputs

   	if(!is.numeric(x)) stop("x should be an integer and numeric\n");
   	if(min(x)<0) stop("Supported codes are 0,1,2,3\n");
   	if(max(x)>3) stop("Supported codes are 0,1,2,3\n");
   	if(n<=0) stop("n should be bigger than 0");
   	if(p<=0) stop("p should be bigger than 0");
   	if(length(x)!=n*p) stop("length of x is not equal to n*p");
          
	#Function call
	.C("write_bed_",as.character(bed_file), as.integer(n), as.integer(p), as.integer(x)) 
}
