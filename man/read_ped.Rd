\name{read_ped}
\alias{read_ped}
\title{read_ped}
\description{
   This function reads genotype information stored in  PED format used in plink.
}
\usage{
  read_ped(ped_file)
}
\arguments{
   \item{ped_file}{ASCII file with genotype information.}
}

\details{

The PED file is a white-space (space or tab) delimited file: the first six columns are mandatory:

     Family ID
     Individual ID
     Paternal ID
     Maternal ID
     Sex (1=male; 2=female; other=unknown)
     Phenotype

The IDs are alphanumeric: the combination of family and individual ID should uniquely identify a person. 
A PED file must have 1 and only 1 phenotype in the sixth column. The phenotype can be either a quantitative 
trait or an affection status column.

}

\value{

The routine will return a vector of dimension n*p (n=number of individuals, p=number of snps), with the snps stacked.

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
demo(read_ped)

}

}
