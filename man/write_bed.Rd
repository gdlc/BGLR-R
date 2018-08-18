\name{write_bed}
\alias{write_bed}
\title{write_bed}
\description{

  This function writes genotype information into  a binary PED (BED) filed used in plink.
  For more details about this format see  \url{http://zzz.bwh.harvard.edu/plink/binary.shtml}.
}
\usage{
  write_bed(x,n,p,bed_file)
}
\arguments{
   \item{n}{integer, number of individuals.}
   \item{p}{integer, number of SNPs.}
   \item{x}{integer vector that contains the genotypic information coded as 0,1,2 and 3 (see details below). 
            The information must be in snp major order. The vector should be of dimension n*p with the snps stacked.}
   \item{bed_file}{output binary file with genotype information.}
}

\details{

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
demo(write_bed)

}

}
