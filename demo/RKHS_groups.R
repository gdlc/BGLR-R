# Single-response genomic main effects and G x E interactions with
# environment-specific residual variances. See inst/md/RKHS_groups.md.

library(BGLR)
data(wheat)

# A subset keeps the demonstration short. The same construction applies to
# the full set of genotypes. The bundled matrices have matching row orders.
set.seed(2401)
selected=sort(sample(seq_len(nrow(wheat.Y)),120))
Y=wheat.Y[selected,,drop=FALSE]
M=wheat.X[selected,,drop=FALSE]
M=M[,apply(M,2,sd)>0,drop=FALSE]
M=scale(M,center=TRUE,scale=TRUE)
G=tcrossprod(M)/ncol(M)
G=G/mean(diag(G))

# One response vector, ordered by genotype within environment.
genotype=rep(seq_len(nrow(Y)),ncol(Y))
environment=factor(rep(colnames(Y),each=nrow(Y)))
y=as.vector(Y)
K_G=G[genotype,genotype]
K_GE=K_G*outer(as.character(environment),as.character(environment),"==")
X_E=model.matrix(~environment)[,-1,drop=FALSE]

# Withhold genotype-environment cells. Every row remains in both kernels.
set.seed(2402)
testing=sample(seq_along(y),floor(length(y)*0.1))
yNA=y
yNA[testing]=NA

ETA=list(environment=list(X=X_E,model="FIXED"),
         G=list(K=K_G,model="RKHS",df0=5,S0=1),
         GE=list(K=K_GE,model="RKHS",df0=5,S0=1))
outputDir=tempfile("MET_RKHS_")
dir.create(outputDir)
set.seed(2403)
fm=BGLR(y=yNA,ETA=ETA,groups=environment,
        nIter=6000,burnIn=1000,thin=5,df0=5,S0=1,
        saveAt=file.path(outputDir,"MET_"),verbose=FALSE)

fm$varE
c(G=fm$ETA$G$varU,GE=fm$ETA$GE$varU)
tapply((y[testing]-fm$yHat[testing])^2,environment[testing],
       function(x) sqrt(mean(x)))
outputDir

# These are demonstration settings. Examine variance traces and run
# independent chains before using posterior summaries for scientific inference.
