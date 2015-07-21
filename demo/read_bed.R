
bed_file = system.file("extdata/sample.bed", package="BGLR")

#Extended map file (this gives the number of snps)
bim_file = system.file("extdata/sample.bim", package="BGLR")

#First 6 columns of ped file (this gives the number of individuals)
fam_file = system.file("extdata/sample.fam", package="BGLR")

out=read_bed(bed_file=bed_file,bim_file=bim_file,fam_file=fam_file,verbose=TRUE)

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
