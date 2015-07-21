bed_file = system.file("extdata/sample.bed", package="BGLR")

#Extended map file (this gives the number of snps)
bim_file = system.file("extdata/sample.bim", package="BGLR")

#First 6 columns of ped file (this gives the number of individuals)
fam_file = system.file("extdata/sample.fam", package="BGLR")

out=read_bed(bed_file=bed_file,bim_file=bim_file,fam_file=fam_file,verbose=TRUE)

#Now write the bed file using the internal routine
#Using the routine xxd compare both files, i.e. extdata/sample.bed and test.bed

new_bed_file="test.bed"
write_bed(x=out$x,n=out$n,p=out$p,bed_file=new_bed_file)


