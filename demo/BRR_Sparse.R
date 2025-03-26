rm(list=ls())

setwd(tempdir())

#Libraries, load BGLR first... and then Matrix
library(Matrix)

#Example with sparse matrix

data(wheat)
nIter<-12000
burnIn<-2000
X<-scale(wheat.X)/sqrt(ncol(wheat.X))
y<-wheat.Y[,1]
G<-tcrossprod(X)
diag(G)<-diag(G)+1/1e4
L<-t(chol(G))
Ls<-as(L,"dgCMatrix")

object.size(L)
object.size(Ls)

#Non zero values
Ls@x

#Cumulative number of non zero elements
Ls@p

#row index of each element
Ls@i

ETAs<-list(list(X=Ls,model="BRR_sparse"))

set.seed(123)
fms<-BGLR(y=y,ETA=ETAs,nIter=10000,burnIn=5000)
unlink("*.dat")

plot(y,fms$yHat)

ETA<-list(list(X=as.matrix(Ls),model="BRR"))

set.seed(123)
fm<-BGLR(y=y,ETA=ETA,nIter=10000,burnIn=5000)
unlink("*.dat")

plot(y,fm$yHat)

plot(fms$yHat,fm$yHat)
cor(fms$yHat,fm$yHat)

plot(fms$ETA[[1]]$b,fm$ETA[[1]]$b)
cor(fms$ETA[[1]]$b,fm$ETA[[1]]$b)


