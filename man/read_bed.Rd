\name{read_bed}
\alias{read_bed}
\title{read_bed}
\description{

  This function reads genotype information stored in  binary PED (BED) files used in plink.
  These files save space and time. The pedigree/phenotype information is stored in a separate file (*.fam) 
  and the map information is stored in an extededed MAP file (*.bim) that contains information about the allele names, 
  which would otherwise be lost in the BED file. More details \url{http://zzz.bwh.harvard.edu/plink/binary.shtml}.
}
\usage{
  read_bed(bed_file,bim_file,fam_file,na.strings,verbose)
}
\arguments{
   \item{bed_file}{binary file with genotype information.}
   \item{bim_file}{text file with pedigree/phenotype information.}
   \item{fam_file}{text file with extended map information.}
   \item{na.strings}{missing value indicators, default=c("0","-9").}
   \item{verbose}{logical, if true print hex dump of bed file.}
}


\value{

The routine will return a vector of dimension n*p (n=number of individuals, p=number of snps), with the snps(individuals) stacked, 
depending whether the BED file is in SNP-major or individual-major mode.

The vector contains integer codes:

\tabular{ll}{
   Integer code \tab Genotype \cr
   0 \tab 00 Homozygote "1"/"1" \cr
   1 \tab 01 Heterozygote \cr
   2 \tab 10 Missing genotype \cr
   3 \tab 11 Homozygote "2"/"2" \cr
 }

}

\author{
Gustavo de los Campos, Paulino Perez Rodriguez,
}

\examples{

\dontrun{

library(BGLR)
demo(read_bed)

}

}
