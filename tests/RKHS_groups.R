library(BGLR)

# Gaussian conditional distributions and kernel-factor equivalence are
# derived in inst/md/RKHS_groups_derivation.md.
local({
    outputDir=tempfile("RKHS_groups_")
    dir.create(outputDir)
    on.exit(unlink(outputDir,recursive=TRUE))

    equal=function(x,y,tolerance=1e-8)
    {
        stopifnot(isTRUE(all.equal(unname(x),unname(y),
                                  tolerance=tolerance,check.attributes=FALSE)))
    }

    fit=function(y,ETA,groups=NULL,weights=NULL,seed=107,nIter=200,burnIn=50)
    {
        set.seed(seed)
        BGLR(y=y,ETA=ETA,groups=groups,weights=weights,
             nIter=nIter,burnIn=burnIn,thin=5,df0=5,S0=1,
             saveAt=tempfile(tmpdir=outputDir),verbose=FALSE)
    }

    # Compare one complete coordinate sweep with a dense precision-matrix
    # calculation. The reference uses Q directly, independently of x2.
    set.seed(918)
    n=30L
    p=5L
    V=qr.Q(qr(matrix(rnorm(n*p),n,p)))
    d=seq(0.4,2,length.out=p)
    groups=rep(1:3,length.out=n)
    varE=c(0.2,1,3)
    varU=0.7
    y=rnorm(n)
    a=rnorm(p)
    residual=y-as.vector(V%*%a)
    Q=crossprod(V,V/varE[groups])+diag(1/(varU*d))
    h=as.vector(crossprod(V,y/varE[groups]))
    reference=a+0
    set.seed(53)
    for(j in seq_len(p))
    {
        m=(h[j]-sum(Q[j,-j]*reference[-j]))/Q[j,j]
        reference[j]=rnorm(1,mean=m,sd=sqrt(1/Q[j,j]))
    }
    x2=matrix(0,3,p)
    for(g in 1:3) x2[g,]=colSums(V[groups==g,,drop=FALSE]^2)
    set.seed(53)
    sampled=.Call("sample_beta_groups",n,p,V,as.vector(x2),a,residual,
                  varU*d,varE,0,as.integer(groups-1L),3L,PACKAGE="BGLR")
    equal(sampled[[1]],reference,1e-12)
    equal(sampled[[2]],y-as.vector(V%*%reference),1e-12)

    # Multiple kernels, fixed environment effects, unequal weights, unused
    # factor levels, and missing responses in different variance groups.
    set.seed(404)
    n=36
    markers=matrix(rnorm(n*5),n,5)
    K=tcrossprod(markers)/5
    environment=factor(rep(c("dry","wet","normal"),length.out=n),
                       levels=c("wet","dry","unused","normal"))
    kernels=list(G=K,GE=K*outer(as.character(environment),
                              as.character(environment),"=="))
    weights=seq(0.7,1.4,length.out=n)
    y=rnorm(n)
    y[c(8,10,19)]=NA
    X_E=model.matrix(~factor(as.character(environment)))[,-1,drop=FALSE]
    ETA_K=list(E=list(X=X_E,model="FIXED"))
    ETA_V=ETA_X=ETA_K
    factors=list()
    for(term in names(kernels))
    {
        eig=eigen(sweep(sweep(kernels[[term]],1L,weights,"*"),
                        2L,weights,"*"),symmetric=TRUE)
        keep=eig$values>1e-10
        V=eig$vectors[,keep,drop=FALSE]
        d=eig$values[keep]
        factors[[term]]=sweep(sweep(V,2L,sqrt(d),"*"),1L,weights,"/")
        ETA_K[[term]]=list(K=kernels[[term]],model="RKHS",df0=5,S0=2)
        ETA_V[[term]]=list(V=V,d=d,model="RKHS",df0=5,S0=2)
        ETA_X[[term]]=list(X=factors[[term]],model="BRR",df0=5,S0=2)
    }
    kernelFit=fit(y,ETA_K,environment,weights)
    # Supplied eigenvectors already describe W K W; the existing warning is expected.
    eigenFit=suppressWarnings(fit(y,ETA_V,environment,weights))
    factorFit=fit(y,ETA_X,environment,weights)
    equal(kernelFit$yHat,factorFit$yHat)
    equal(kernelFit$SD.yHat,factorFit$SD.yHat)
    equal(kernelFit$varE,factorFit$varE)
    equal(kernelFit$SD.varE,factorFit$SD.varE)
    equal(kernelFit$yHat,eigenFit$yHat)
    equal(kernelFit$varE,eigenFit$varE)
    prediction=kernelFit$mu+as.vector(X_E%*%kernelFit$ETA$E$b)
    for(term in names(kernels))
    {
        equal(kernelFit$ETA[[term]]$u,
              as.vector(factors[[term]]%*%factorFit$ETA[[term]]$b))
        equal(kernelFit$ETA[[term]]$varU,factorFit$ETA[[term]]$varB)
        equal(kernelFit$ETA[[term]]$SD.varU,factorFit$ETA[[term]]$SD.varB)
        equal(kernelFit$ETA[[term]]$SD.u,eigenFit$ETA[[term]]$SD.u)
        prediction=prediction+kernelFit$ETA[[term]]$u
        stopifnot(is.null(kernelFit$ETA[[term]]$x2))
        samples=scan(paste0(kernelFit$saveAt,"ETA_",term,"_varU.dat"),quiet=TRUE)
        equal(mean(samples[11:40]),kernelFit$ETA[[term]]$varU,1e-6)
    }
    equal(kernelFit$yHat/weights,prediction)
    stopifnot(all(is.finite(kernelFit$yHat)),
              identical(names(kernelFit$varE),c("dry","normal","wet")))
    residualSamples=as.matrix(read.table(paste0(kernelFit$saveAt,"varE.dat")))
    equal(colMeans(residualSamples[11:40,,drop=FALSE]),kernelFit$varE,1e-6)
    logLik=sum(dnorm(y-kernelFit$yHat/weights,
                    sd=sqrt(kernelFit$varE[as.character(environment)])/weights,
                    log=TRUE),na.rm=TRUE)
    equal(kernelFit$fit$logLikAtPostMean,logLik)

    # One residual group must reduce to the homogeneous RKHS conditional.
    y[is.na(y)]=0
    homogeneous=fit(y,ETA_K,weights=weights)
    oneGroup=fit(y,ETA_K,groups=rep("all",n),weights=weights)
    equal(homogeneous$yHat,oneGroup$yHat)
    equal(homogeneous$varE,oneGroup$varE)
    for(term in names(kernels)) equal(homogeneous$ETA[[term]]$u,oneGroup$ETA[[term]]$u)

    # Row permutations preserve the model when supplied eigenvectors and
    # group labels are permuted with the observations.
    set.seed(1072)
    perm=sample.int(n)
    permETA=ETA_V
    permETA$E$X=permETA$E$X[perm,,drop=FALSE]
    for(term in names(kernels)) permETA[[term]]$V=permETA[[term]]$V[perm,,drop=FALSE]
    complete=suppressWarnings(fit(y,ETA_V,environment,weights))
    permuted=suppressWarnings(fit(y[perm],permETA,environment[perm],weights[perm]))
    equal(complete$yHat[perm],permuted$yHat)
    equal(complete$varE,permuted$varE)

    # Rank-one and full-rank kernels, including initialization with default priors.
    rankOne=tcrossprod(seq(-1,1,length.out=n))
    rankFit=fit(y,list(list(K=rankOne,model="RKHS")),environment)
    stopifnot(rankFit$ETA[[1]]$levelsU==1,
              is.matrix(rankFit$ETA[[1]]$V),all(is.finite(rankFit$yHat)))
    fullFit=fit(y,list(list(K=K+diag(0.2,n),model="RKHS")),environment)
    stopifnot(fullFit$ETA[[1]]$levelsU==n,all(is.finite(fullFit$yHat)))
    id=rep(seq_len(n/3),3)
    repeatedFit=fit(y,list(list(K=K[id,id],model="RKHS")),environment)
    stopifnot(all(is.finite(repeatedFit$yHat)))

    # An independent intercept-only Gibbs chain checks missing-response
    # augmentation and the group variance posterior, rather than comparing
    # two ETA parameterizations that share the same augmentation code.
    y=c(0.2,-1,NA,0.4,1.2,NA,-0.8,0.6,1.5,NA,0.1,-0.2)
    groups=rep(1:3,each=4)
    weights=seq(0.8,1.3,length.out=length(y))
    missing=which(is.na(y))
    nIter=100
    burnIn=25
    thin=5
    set.seed(107)
    mu=weighted.mean(y,weights,na.rm=TRUE)
    augmented=y*weights
    augmented[missing]=mu*weights[missing]
    varE=rep(var(augmented-mu*weights)*0.5/3,3)
    muSamples=numeric(nIter)
    varianceSamples=matrix(0,nIter,3)
    for(iter in seq_len(nIter))
    {
        precision=sum(weights^2/varE[groups])
        meanMu=sum(augmented*weights/varE[groups])/precision
        mu=rnorm(1,sd=sqrt(1/precision))+meanMu
        residual=augmented-mu*weights
        for(g in 1:3)
        {
            varE[g]=(sum(residual[groups==g]^2)+1)/rchisq(1,df=4+5)
        }
        augmented[missing]=mu*weights[missing]+rnorm(length(missing),
                                                    sd=sqrt(varE[groups[missing]]))
        muSamples[iter]=mu
        varianceSamples[iter,]=varE
    }
    grouped=fit(y,NULL,groups,weights,nIter=nIter,burnIn=burnIn)
    retained=seq(burnIn+thin,nIter,by=thin)
    equal(grouped$mu,mean(muSamples[retained]))
    equal(grouped$varE,colMeans(varianceSamples[retained,,drop=FALSE]))
    equal(grouped$yHat,mean(muSamples[retained])*weights)

    # Invalid kernel dimensions must fail before entering the compiled sampler.
    initializer=getFromNamespace("setLT.RKHS","BGLR")
    invalid=function(LT,pattern)
    {
        result=tryCatch(initializer(LT=LT,y=y,n=length(y),j=1,
                         weights=rep(1,length(y)),saveAt=tempfile(tmpdir=outputDir),
                         R2=0.5,nLT=1,rmExistingFiles=TRUE,verbose=FALSE,
                         groups=groups,nGroups=3),error=conditionMessage)
        stopifnot(is.character(result),grepl(pattern,result,fixed=TRUE))
    }
    invalid(list(K=diag(2)),"number of phenotypes")
    invalid(list(V=matrix(1,length(y),2),d=1),"number of eigenvectors")
    invalid(list(V=matrix(NA_real_,length(y),1),d=1),"must be finite")
    invalid(list(K=matrix(0,length(y),length(y))),"No eigenvalues")

    invalidCall=function(extra,pattern)
    {
        arguments=list(y=rep(c(0,1),6),groups=rep(1:3,4),
                       ETA=list(list(K=diag(12),model="RKHS")),
                       nIter=10,burnIn=5,thin=1,
                       saveAt=tempfile(tmpdir=outputDir),verbose=FALSE)
        arguments[names(extra)]=extra
        result=tryCatch(do.call(BGLR,arguments),error=conditionMessage)
        stopifnot(is.character(result),grepl(pattern,result,fixed=TRUE))
    }
    invalidCall(list(response_type="ordinal"),"uncensored Gaussian")
    invalidCall(list(a=rep(-Inf,12),b=rep(Inf,12)),"uncensored Gaussian")
    namedModel=list(list(K=diag(12),model=c(term="RKHS")))
    invalidCall(list(ETA=namedModel,response_type="ordinal"),"uncensored Gaussian")
    invalidCall(list(ETA=namedModel,a=rep(-Inf,12),b=rep(Inf,12)),"uncensored Gaussian")
    invalidCall(list(groups=c(rep(1,12),NA)),"NA's not allowed in groups")
})
