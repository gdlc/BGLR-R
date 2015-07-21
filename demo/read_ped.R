
ped_file = system.file("extdata/sample.ped", package="BGLR")

out=read_ped(ped_file)

p=out$p
n=out$n
out=out$x

#Recode snp to 0,1,2 format using allele 1
# 0 --> 0
# 1 --> 1
# 2 --> NA
# 3 --> 2

out[out==2]=NA
out[out==3]=2

X=matrix(out,nrow=p,ncol=n,byrow=TRUE)
X=t(X)


