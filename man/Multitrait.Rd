\name{Multitrait}
\alias{Multitrait}
\title{Multi trait models}

\description{
  The Multitrait  function fits Bayesian multitrait models with arbitrary number of random effects ussing a Gibbs 
  sampler. The data equation is as follows:
  
  \ifelse{html}{\out{y<sub>j</sub> = 1&mu;<sub>j</sub> + X<sub>Fj</sub> &beta;<sub>Fj</sub>+X<sub>1</sub>&beta;<sub>j1</sub> + ... +X<sub>k</sub>&beta;<sub>jk</sub> + u<sub>j1</sub>+ ... + u<sub>jq</sub> + e<sub>j</sub>,}}{\deqn{y_j= 1 \mu_j + X_{Fj} \beta_{Fj} + X_1 \beta_{j1} + ... + X_k \beta_{jk} + u_{j1}+...+u_{jq}+e_j,}}
 
  where:
  
  \itemize{
  	
  	\item{\ifelse{html}{\out{y<sub>j</sub> is a n-dimensional response vector of phenotypes (NAs allowed) with y<sub>ij</sub> representing the phenotypic record of the i-th subject for the j-th trait.}}{\eqn{y_j}{y_j} is a n-dimensional response vector of phenotypes (NAs allowed) with \eqn{y_{ij}}{y_ij} representing the phenotypic record of the i-th subject for the j-th trait.}}
  
   	\item{\ifelse{html}{\out{&mu;<sub>j</sub> is an intercept.}}{\eqn{\mu_j} is an intercept.}}
   	
   	\item{\ifelse{html}{\out{X<sub>Fj</sub> is a matrix of fixed effects.}}{\eqn{X_{Fj}} is a matrix of fixed effects.}}
   	
   	\item{\ifelse{html}{\out{&beta;<sub>Fj</sub> is a vector of fixed effects.}}{\eqn{\beta_{Fj}} is a vector of fixed effects.}}
   	
   	\item{\ifelse{html}{\out{X<sub>s</sub> is an incidence matrix for predictors that are common for all individuals (e.g. markers), s=1,...,k.}}{\eqn{X_s} is an incidence matrix for predictors that are common for all the individuals (e.g. markers), s=1,...,k.}}
   	
   	\item{\ifelse{html}{\out{&beta;<sub>js</sub>  is a vector of regression coefficients, s=1,..,k. Different priors can be assigned to these regression coefficients (spike-slab and Gaussian) and regression coefficients are correlated across traits.}}{\eqn{\beta_{js}} is a vector of regression coefficients, s=1,...,k. Different priors can be assigned to these regression coefficients (spike-slab and Gaussian) and regression coefficients are correlated across traits.}}
   	
   	\item{\ifelse{html}{\out{u<sub>jr</sub> is a n-dimensional vector of random effects.}}{\eqn{u_{jr}} is a n-dimensional vector of random effects.}}
   	
   	\item{\ifelse{html}{\out{e<sub>j</sub> is a n-dimensional vector of residuals.}}{\eqn{e_j} is a n-dimensional vector of residuals.}}
   	}
}

\usage{
 Multitrait(y, ETA, intercept=TRUE,resCov = list(df0=5,S0=NULL,type="UN"), 
   R2=0.5, nIter=1000,burnIn=500,thin=10, saveAt="",verbose=TRUE)
}

\arguments{
   \item{y}{a matrix of dimension \eqn{n \times t}{n x t}, where \eqn{n}{n} is the number of individials in each trait, 
   	           \eqn{t}{t} is the number of traits, NAs are allowed.}
   \item{ETA}{ (list) This is a two-level list used to specify the regression function (or linear predictor). 
              Regression on covariates and other types of random effects are specified in this two-level list. For instance:

\preformatted{ETA=list(list(X=W, model="FIXED"), 
      list(X=Z,model="BRR"), 
      list(K=G,model="RKHS"))} 

	          specifies that the linear predictor should include: a linear regression on W with regression 
              coefficients treated as fixed effects (i.e., flat prior), plus regression on Z, with regression 
              coefficients modeled as in the Ridge Regression plus and a random effect with co-variance structure G.
              }
    \item{intercept}{logical, if TRUE (default) an intercept is included.}
   	\item{nIter,burnIn, thin}{(integer) the number of iterations, burn-in and thinning.}
   	
   	\item{resCov}{A list used to define the co-variance matrix for model residuals (R). Four covariance strucures are 
   		    supported: i) Unstructured (\code{"UN"}), ii)Diagonal (\code{"DIAG"}), iii)Factor Analytic (\code{"FA"}) 
   		    and  iv)Recursive (\code{"REC"}), for example: 
   		    
   		    \preformatted{resCov=list(type="UN", df0=4, S0=V)} 
   		    
   		    specifies an UN-structured covariance matrix, 
   		    with an Inverse Whishart prior with degree of freedom df0 (scalar) and scale matrix (t x t) V.}
  	 \item{saveAt}{(string) this may include a path and a pre-fix that will be added to the name of the files that are saved as the program runs.}
  	 \item{R2}{(numeric, \code{0<R2<1}) The proportion of co-variance that one expects, a priori, to be explained by the regression. Only used if 
             the hyper-parameters are not specified; if that is the case, internaly, hyper-paramters are set so that the prior modes are consistent with the 
             co-variance partition specified by R2 and the prior distribution is relatively flat at the mode.}
     \item{verbose}{(logical) if TRUE the iteration history is printed, default TRUE.}
}

\details{
	
	\emph{Conditional distribution of the data}
	
	Model residuals are assumed to follow a multivariate normal distribution, with null mean and 
	covariance matrix   \ifelse{html}{\out{Cov((e'<sub>1</sub>,...,e'<sub>n</sub>)')=R<sub>0</sub> &otimes; I}}{\eqn{Cov((e'_1,...,e'_n)')=R_0 \otimes I}}
	where \ifelse{html}{\out{R<sub>0</sub>}}{\eqn{R_0}} is a t x t (within-subject) covariance matrix of model residuals
	and n-dimensional identity matrix. Therefore:
	
	\ifelse{html}{\out{p(y<sub>i1</sub>,...,y<sub>it</sub> | &theta; )=MN(&eta;<sub>i</sub>, R<sub>0</sub>)}}{\deqn{p(y_{i1},...,y_{it}|\theta)=MN(\eta_i,R_0)}}
	
	where MN(.,.), denotes the multivariate normal distribution with mean 
	\ifelse{html}{\out{&eta;<sub>i</sub>}}{\eqn{\eta_i}} and covariance matrix 
	\ifelse{html}{\out{R<sub>0</sub>}}{\eqn{R_0}}; here 
	\ifelse{html}{\out{&eta;<sub>i</sub>}}{\eqn{\eta_i}} is 
	a t-dimensional vector whose entries
	are the expected values of the response variable for the i-th individual. 
	
	\emph{Prior distribution}
	
	The prior distribution is structured hierarchically. The first level of the prior specifies the distribution 
	of the fixed and random effects given the codispersion parameters 
	(the covariance matrices of the random effects, see below). 
	The priors for the codispersion parameters are specified in a deeper level of the hierarchy.

	The intercepts and vectors of fixed effects are assigned flat prior 
	(each unknown is assigned a Gaussian prior with null mean and very large variance).
	
	The vectors of random effects \ifelse{html}{\out{u<sub>r</sub>}}{\eqn{u_r}} 
	are assigned independent multivariate normal priors with null mean and covariance matrices 
	\ifelse{html}{\out{Cov(u<sub>r</sub>)=G<sub>r</sub>&otimes;K<sub>r</sub>}}{\eqn{Cov(u_r)=G_r \otimes K_r}},
	\ifelse{html}{\out{u<sub>r</sub>}}{\eqn{u_r}} represent the vector of effects for the 
	r-th random effects (sorted by subject first and trait within subject),
	\ifelse{html}{\out{G<sub>r</sub>}}{\eqn{G_r}} is an t x t (within-subject) covariance matrix
	of the r-th random effect and \ifelse{html}{\out{K<sub>r</sub>}}{\eqn{K_r}} is 
	a user defined (between subjects) covariance matrix for the r-th
	random effect, for instance, may be a pedigree or marker-based 
	relationship matrix. The covariance matrix of random effects are assigned Inverse Wishart
	priors (for the case of unstructured options) or priors that are structured according to
	some model (diagonal, factor analytic or recursive).
	
	The vector or regression coefficients \ifelse{html}{\out{&beta;<sub>s</sub>}}{\eqn{\beta_s}}
	are assigned Gaussian  and Spike Slab priors whose 
	covariance matrixes depend on a \ifelse{html}{\out{&Omega;<sub>s</sub>}}{\eqn{\Omega_s}}
	covariance matrix of dimmensions t x t (within subject).
	The covariance matrix of random effects are assigned Inverse Wishart
	priors (for the case of unstructured options) or priors that are structured according to
	some model (diagonal, factor analytic or recursive).
	
	\emph{Algorithm}
	
	Internally, samples from all the model unknowns are drawn using a Gibbs sampler 
	(i.e., based on fully conditional distributions).
	
	
}

\value{
  List containing estimated posterior means and estimated posterior
  standard deviations.
}

\references{
  
  Burgueno, J., G. de los Campos, K. Weigel and J. Crossa. 2012. Genomic Prediction of Breeding Values when Modelling
  Genotype x Environment Interaction using Pedigree and Dense Molecular Markers. \emph{Crop Sci.}, \bold{52(2)}:707-719.
  
  de los Campos, G. and  D. Gianola. 2007. Factor analysis models for structuring covariance matrices of additive genetic 
  effects: a Bayesian implementation. \emph{Genet. Sel. Evol.}, \bold{39}:481-494.
  
  Cheng, H.,  K. Kadir, J., Zeng, D. Garrick and R. Fernando. 2018. Genomic Prediction from Multiple-Trait 
  Bayesian Regression Methods Using Mixture Priors. \emph{Genetics}, \bold{209(1)}: 89-103.
  
  Sorensen, D. and D. Gianola. 2002. Likelihood, Bayesian, and MCMC methods in quantitative genetics. 
  \emph{Springer-Verlag, New York}.
}

\author{
Gustavo de los Campos, Paulino Perez-Rodriguez
}

\examples{

\dontrun{
	
library(BGLR)
data(wheat)
X<-wheat.X
K<-wheat.A
y<-wheat.Y

#Example 1, Spike Slab regression
ETA1<-list(list(X=X,model="SpikeSlab"))

fm1<-Multitrait(y=y,ETA=ETA1,nIter=1000,burnIn=500)

#Example 2, Ridge Regression
ETA2<-list(list(X=X,model="BRR"))
fm2<-Multitrait(y=y,ETA=ETA2,nIter=1000,burnIn=500)

#Example 3, Random effects with user defined covariance structure
#for individuals derived from pedigree
ETA3<-list(list(K=K,model="RKHS"))
fm3<-Multitrait(y=y,ETA=ETA3,nIter=1000,burnIn=500)

#Example 4, Markers and pedigree
ETA4<-list(list(X=X,model="BRR"), list(K=K,model="RKHS"))

fm4<-Multitrait(y=y,ETA=ETA4,nIter=1000,burnIn=500)

#Example 5, recursive structures for within subject covariance matrix
M1 <- matrix(nrow = 4, ncol = 4, FALSE)
M1[3, 2] <- M1[4, 2] <- TRUE # Adding recursion from trait 2 onto traits 3 and 4
M1[4, 3] <- TRUE # Adding recursion from trait 3 on trait 4

ETA5<-list(list(K=K,model="RKHS",Cov=list(type="REC",M=M1)))
fm5<-Multitrait(y=y,ETA=ETA5,nIter=1000,burnIn=500)

#Example 6, diagonal residual covariance matrix with the predictor 
#used in example 5
residual1<-list(type="DIAG")
fm6<-Multitrait(y=y,ETA=ETA5,resCov=residual1,nIter=1000,burnIn=500)

}

}
\keyword{models}
