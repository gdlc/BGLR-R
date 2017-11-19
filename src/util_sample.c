#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rconfig.h>
#include <R_ext/Lapack.h>

/*
 * This is a generic function to sample betas in various models, including 
 * Bayesian LASSO, BayesA, Bayesian Ridge Regression, etc.
 
 * For example, in the Bayesian LASSO, we wish to draw samples from the full 
 * conditional distribution of each of the elements in the vector bL. The full conditional 
 * distribution is normal with mean and variance equal to the solution (inverse of the coefficient of the left hand side)
 * of the following equation (See suplementary materials in de los Campos et al., 2009 for details),
   
    (1/varE x_j' x_j + 1/(varE tau_j^2)) bL_j = 1/varE x_j' e
 
    or equivalently, 
    
    mean= (1/varE x_j' e)/ (1/varE x_j' x_j + 1/(varE tau_j^2))
    variance= 1/ (1/varE x_j' x_j + 1/(varE tau_j^2))
    
    xj= the jth column of the incidence matrix
    
 *The notation in the routine is as follows:
 
 n: Number of rows in X
 pL: Number of columns in X
 XL: the matrix X stacked by columns
 XL2: vector with x_j' x_j, j=1,...,p
 bL: vector of regression coefficients
 e: vector with residuals, e=y-yHat, yHat= predicted values
 varBj: vector with variances, 
	For Bayesian LASSO, varBj=tau_j^2 * varE, j=1,...,p
	For Ridge regression, varBj=varB, j=1,...,p, varB is the variance common to all betas.
	For BayesA, varBj=varB_j, j=1,...,p
	For BayesCpi, varBj=varB, j=1,...,p, varB is the variance common to all betas
	
 varE: residual variance
 minAbsBeta: in some cases values of betas near to zero can lead to numerical problems in BL, 
             so, instead of using this tiny number we assingn them minAbsBeta
 
 */


SEXP sample_beta(SEXP n, SEXP pL, SEXP XL, SEXP xL2, SEXP bL, SEXP e, SEXP varBj, SEXP varE, SEXP minAbsBeta)
{
    double *xj;
    double *pXL, *pxL2, *pbL, *pe, *pvarBj;
    double b;
    int inc=1;
    double rhs,c,sigma2e, smallBeta;
    int j, rows, cols;

    SEXP list;	

    GetRNGstate();
	
    rows=INTEGER_VALUE(n);
    cols=INTEGER_VALUE(pL);
    sigma2e=NUMERIC_VALUE(varE);
    smallBeta=NUMERIC_VALUE(minAbsBeta);
	
    PROTECT(XL=AS_NUMERIC(XL));
    pXL=NUMERIC_POINTER(XL);

    PROTECT(xL2=AS_NUMERIC(xL2));
    pxL2=NUMERIC_POINTER(xL2);

    PROTECT(bL=AS_NUMERIC(bL));
    pbL=NUMERIC_POINTER(bL);

    PROTECT(e=AS_NUMERIC(e));
    pe=NUMERIC_POINTER(e);

    PROTECT(varBj=AS_NUMERIC(varBj));
    pvarBj=NUMERIC_POINTER(varBj);

    for(j=0; j<cols;j++)
    {
          xj=pXL+j*rows;
          b=pbL[j];
          //F77_NAME(daxpy)(&rows, &b,xj,&inc, pe, &inc);
          rhs=F77_NAME(ddot)(&rows,xj,&inc,pe,&inc)/sigma2e;
          rhs+=pxL2[j]*b/sigma2e;
  	  c=pxL2[j]/sigma2e + 1.0/pvarBj[j];
	  pbL[j]=rhs/c + sqrt(1.0/c)*norm_rand();

          b-=pbL[j];
          //b=-pbL[j];        

          F77_NAME(daxpy)(&rows, &b,xj,&inc, pe,&inc);
          
          if(fabs(pbL[j])<smallBeta)
          {
             pbL[j]=smallBeta;
          }
      }
        
      // Creating a list with 2 vector elements:
      PROTECT(list = allocVector(VECSXP, 2));
      // attaching bL vector to list:
      SET_VECTOR_ELT(list, 0, bL);
      // attaching e vector to list:
      SET_VECTOR_ELT(list, 1, e);

      PutRNGstate();

      UNPROTECT(6);

      return(list);
}

SEXP sample_beta_lower_tri(SEXP n, SEXP pL, SEXP XL, SEXP xL2, SEXP bL, SEXP e, SEXP varB, SEXP varE, SEXP minAbsBeta)
{
    double *xj;
    double *pXL, *pxL2, *pbL, *pe;
    double b;
    int inc=1;
    double rhs,c, sigma2b, sigma2e, smallBeta;
    int j, rows, cols;

    SEXP list;	

    GetRNGstate();
	
    rows=INTEGER_VALUE(n);
    cols=INTEGER_VALUE(pL);
    sigma2b=NUMERIC_VALUE(varB);
    sigma2e=NUMERIC_VALUE(varE);
    smallBeta=NUMERIC_VALUE(minAbsBeta);
	
    PROTECT(XL=AS_NUMERIC(XL));
    pXL=NUMERIC_POINTER(XL);

    PROTECT(xL2=AS_NUMERIC(xL2));
    pxL2=NUMERIC_POINTER(xL2);

    PROTECT(bL=AS_NUMERIC(bL));
    pbL=NUMERIC_POINTER(bL);

    PROTECT(e=AS_NUMERIC(e));
    pe=NUMERIC_POINTER(e);

    xj=pXL;

    int r=rows;
    double *pe1;
    pe1=pe;

    for(j=0; j<cols;j++)
    {
          
          b=pbL[j];

          rhs=F77_NAME(ddot)(&r,xj,&inc,pe1,&inc)/sigma2e;
          rhs+=pxL2[j]*b/sigma2e;
  	  c=pxL2[j]/sigma2e + 1.0/sigma2b;
	  pbL[j]=rhs/c + sqrt(1.0/c)*norm_rand();

          b-=pbL[j];
          //b=-pbL[j];        

          F77_NAME(daxpy)(&r, &b,xj,&inc, pe1,&inc);
          
          if(fabs(pbL[j])<smallBeta)
          {
             pbL[j]=smallBeta;
          }
          
          //Update the pointer to covariates
          xj+=rows-j;
         
          //Update the pointer to the error
	  pe1+=1;

          r-=1;

      }
        
      // Creating a list with 2 vector elements:
      PROTECT(list = allocVector(VECSXP, 2));
      // attaching bL vector to list:
      SET_VECTOR_ELT(list, 0, bL);
      // attaching e vector to list:
      SET_VECTOR_ELT(list, 1, e);

      PutRNGstate();

      UNPROTECT(5);

      return(list);
}



/*
  Routine for sampling coefficients for BayesCpi and BayesB
*/

SEXP sample_beta_BB_BCp(SEXP n, SEXP p, SEXP X, SEXP x2, SEXP b, SEXP d, SEXP error, SEXP varBj, SEXP varE, SEXP minAbsBeta, SEXP probInside)
{
  int j,rows,cols;
  double sigma2e, probIn, logOdds,tmp,betaj;
  double logOddsPrior;
  double rhs,c;
  double dRSS, Xe;
  double *pX, *perror, *pb, *px2,*pvarBj;
  double *xj;
  int inc=1;
  double c1;
  int *pd;
  int change;
  SEXP list;

  cols=INTEGER_VALUE(p);
  rows=INTEGER_VALUE(n);
  sigma2e=NUMERIC_VALUE(varE);
  probIn=NUMERIC_VALUE(probInside);
  logOddsPrior=log(probIn/(1-probIn));

  PROTECT(X=AS_NUMERIC(X));
  pX=NUMERIC_POINTER(X);

  PROTECT(x2=AS_NUMERIC(x2));
  px2=NUMERIC_POINTER(x2);


  PROTECT(d=AS_INTEGER(d));
  pd=INTEGER_POINTER(d);

  PROTECT(b=AS_NUMERIC(b));
  pb=NUMERIC_POINTER(b);

  PROTECT(error=AS_NUMERIC(error));
  perror=NUMERIC_POINTER(error);

  PROTECT(varBj=AS_NUMERIC(varBj));
  pvarBj=NUMERIC_POINTER(varBj);

  GetRNGstate();

  // for clarity, maybe we need to use -1/2/sigma2e
  c1=-0.5/sigma2e;  

  for(j=0; j<cols; j++)
  {
     xj=pX+j*rows;
     Xe=F77_NAME(ddot)(&rows,perror,&inc,xj,&inc);

     if(pd[j])
     {
	//Indicator variable equal to one   [4]
       dRSS=-1*pb[j]*pb[j]*px2[j]-2*pb[j]*Xe;
     }else{
       //Indicator variable equal to zero [3]
       dRSS=pb[j]*pb[j]*px2[j]-2*pb[j]*Xe;
     }

     logOdds=logOddsPrior+c1*(dRSS);
 
     tmp=1.0/(1.0+exp(-logOdds));// [1]

     change=pd[j];

     pd[j]=unif_rand()<tmp ? 1 : 0;

     //Update residuals
     if(change!=pd[j])
     {
        if(pd[j]>change)// d=0 => d=1
        {
                betaj=-pb[j];
                F77_NAME(daxpy)(&rows, &betaj,xj,&inc, perror, &inc);
                Xe=F77_NAME(ddot)(&rows,perror,&inc,xj,&inc);
        }else{ // d=1 => d=0
                betaj=pb[j];
                F77_NAME(daxpy)(&rows, &betaj,xj,&inc, perror, &inc);
        }
     }

     //Sample the coefficients
     if(pd[j]==0)
     {
        //Sample from the prior
        pb[j]=sqrt(pvarBj[j])*norm_rand();
     }else{
	   //Sampling from the conditional
           rhs=(px2[j]*pb[j] + Xe)/sigma2e;
           c=px2[j]/sigma2e + 1.0/pvarBj[j];
           tmp=rhs/c + sqrt(1.0/c)*norm_rand();

           betaj=pb[j]-tmp;
           F77_NAME(daxpy)(&rows, &betaj,xj,&inc, perror, &inc);
           pb[j]=tmp;
     }

  }

  // Creating a list with 3 vector elements
  PROTECT(list = allocVector(VECSXP,3));

  // attaching b vector to list
  SET_VECTOR_ELT(list, 0,d);

  // attaching error vector to list
  SET_VECTOR_ELT(list, 1, error);

  // attaching b to the list
  SET_VECTOR_ELT(list,2,b);

  PutRNGstate();

  UNPROTECT(7);

  return(list);

}


/*
  Routine for sampling coefficients for BayesCpi and BayesB with different variances by groups
  FIXME: Work in progress, it is not yet finished neither well tested
*/

SEXP sample_beta_BB_BCp_groups(SEXP n, SEXP p, SEXP X, SEXP x2, SEXP b, SEXP d, SEXP error, SEXP varBj, SEXP varE, SEXP minAbsBeta, SEXP probInside, SEXP groups, SEXP nGroups)
{
  int i, j, k, rows, cols;
  int ngroups;  //Number of groups
  int *g;       //pointer for holding groups
  double *sigma2e; //now sigma2e is a vector
  double probIn, logOdds,tmp,betaj;
  double logOddsPrior;
  double sum_rhs,c;
  double sum_dRSS;
  double *dRSS;
  double *Xe;
  double *pX, *perror, *pb, *pX2,*pvarBj;
  double *xj;
  double *xj2;
  int inc=1;
  double *c1;
  int *pd;
  int change;
  SEXP list;

  cols=INTEGER_VALUE(p);
  rows=INTEGER_VALUE(n);
  sigma2e=NUMERIC_POINTER(varE);
  probIn=NUMERIC_VALUE(probInside);
  logOddsPrior=log(probIn/(1-probIn));

  PROTECT(X=AS_NUMERIC(X));
  pX=NUMERIC_POINTER(X);

  PROTECT(x2=AS_NUMERIC(x2));
  pX2=NUMERIC_POINTER(x2);

  PROTECT(d=AS_INTEGER(d));
  pd=INTEGER_POINTER(d);

  PROTECT(b=AS_NUMERIC(b));
  pb=NUMERIC_POINTER(b);

  PROTECT(error=AS_NUMERIC(error));
  perror=NUMERIC_POINTER(error);

  PROTECT(varBj=AS_NUMERIC(varBj));
  pvarBj=NUMERIC_POINTER(varBj);

  ngroups=INTEGER_VALUE(nGroups);
  g=INTEGER_POINTER(groups);

  //c1 =(double *) malloc(ngroups); 
  c1=(double *) R_alloc(ngroups, sizeof(double));

  //Xe =(double *) malloc(ngroups);
  Xe=(double *) R_alloc(ngroups, sizeof(double));
 
  //dRSS=(double *) malloc(ngroups);
  dRSS=(double *) R_alloc(ngroups,sizeof(double));

  GetRNGstate();

  for(k=0;k<ngroups;k++)
  {
     c1[k]=-0.5/sigma2e[k];
  }

  for(j=0; j<cols; j++)
  {
     xj=pX+j*rows;   //Pointer to the j-th column
     xj2=pX2+j*ngroups; //Pointer to the sum of squares of the j-th column
     

     //Clean initialization
     sum_dRSS=0;
     for(k=0; k<ngroups;k++)
     {
        dRSS[k]=0;
        Xe[k]=0;
     }
     
     //Compute Xe for groups
     for(i=0; i<rows;i++)
     {
        Xe[g[i]]+=xj[i]*perror[i];
     }

     if(pd[j])
     {
	//Indicator variable equal to one
        //dRSS=-1*pb[j]*pb[j]*pX2[j]-2*pb[j]*Xe;
        for(k=0; k<ngroups;k++)
        {
           sum_dRSS+=(-1*pb[j]*pb[j]*xj2[k]-2*pb[j]*Xe[k])*c1[k];
        }
        
     }else{
       //Indicator variable equal to zero
       //dRSS=pb[j]*pb[j]*pX2[j]-2*pb[j]*Xe;
       for(k=0; k<ngroups;k++)
       {
          sum_dRSS+=(pb[j]*pb[j]*xj2[k]-2*pb[j]*Xe[k])*c1[k];
       }
     }

     //logOdds=logOddsPrior+c1*(dRSS);
     logOdds=logOddsPrior+sum_dRSS;

     tmp=1.0/(1.0+exp(-logOdds));

     change=pd[j];

     pd[j]=unif_rand()<tmp ? 1 : 0;

     //Update residuals
     if(change!=pd[j])
     {
        if(pd[j]>change)
        {
                betaj=-pb[j];
                F77_NAME(daxpy)(&rows, &betaj,xj,&inc, perror, &inc);

                //Xe=F77_NAME(ddot)(&rows,perror,&inc,xj,&inc);
                for(k=0; k<ngroups;k++) 
                {
                  Xe[k]=0;
                }
                
                for(i=0; i<rows;i++)
                {
                   Xe[g[i]]+=xj[i]*perror[i];  
                }              
        }else{
                betaj=pb[j];
                F77_NAME(daxpy)(&rows, &betaj,xj,&inc, perror, &inc);
        }
     }

     //Sample the coefficients
     if(pd[j]==0)
     {
        //Sample from the prior
        pb[j]=sqrt(pvarBj[j])*norm_rand();
     }else{
	   //Sampling from the conditional
           //rhs=(px2[j]*pb[j] + Xe)/sigma2e;
           //c=px2[j]/sigma2e + 1.0/pvarBj[j];
           //tmp=rhs/c + sqrt(1.0/c)*norm_rand();

           sum_rhs=0;
           c=0;

           for(k=0;k<ngroups;k++)
           {
             sum_rhs+=(xj2[k]*pb[j]+Xe[k])/sigma2e[k];
             c+=xj2[k]/sigma2e[k];
           }
           c+=1.0/pvarBj[j];
           tmp=sum_rhs/c + sqrt(1.0/c)*norm_rand();

           betaj=pb[j]-tmp;
           F77_NAME(daxpy)(&rows, &betaj,xj,&inc, perror, &inc);
           pb[j]=tmp;
     }

  }

  // Creating a list with 3 vector elements
  PROTECT(list = allocVector(VECSXP,3));

  // attaching b vector to list
  SET_VECTOR_ELT(list, 0,d);

  // attaching error vector to list
  SET_VECTOR_ELT(list, 1, error);

  // attaching b to the list
  SET_VECTOR_ELT(list,2,b);

  PutRNGstate();

  UNPROTECT(7);

  return(list);
}

/*
ddot for groups use for unrolling
*/


void weighted_ddot(int n, double *dx, double *dy, int *groups, double *rhs)
{

   int m, i;

  /* Clean-up loop so remaining vector length is a multiple of 5.  */

      m = n % 5;
      if ( m != 0 ) {
         for ( i = 0 ; i < m; i++ )
            rhs[groups[i]]+= dx[i] * dy[i];
         if ( n < 5 )
            return;
      }
      for (i = m; i < n; i = i + 5)
      {
         rhs[groups[i]]+= dx[i] * dy[i]; 
         rhs[groups[i+1]]+= dx[i+1] * dy[i+1];
         rhs[groups[i+2]]+= dx[i+2] * dy[i+2]; 
         rhs[groups[i+3]]+= dx[i+3] * dy[i+3];
         rhs[groups[i+4]]+= dx[i+4] * dy[i+4];
      }
}



/*
 * This is a generic function to sample betas in various models, including 
 * Bayesian LASSO, BayesA, Bayesian Ridge Regression, etc.
 
 * For example, in the Bayesian LASSO, we wish to draw samples from the full 
 * conditional distribution of each of the elements in the vector bL. The full conditional 
 * distribution is normal with mean and variance equal to the solution (inverse of the coefficient of the left hand side)
 * of the following equation (See suplementary materials in de los Campos et al., 2009 for details),
   
    (1/varE x_j' x_j + 1/(varE tau_j^2)) bL_j = 1/varE x_j' e
 
    or equivalently, 
    
    mean= (1/varE x_j' e)/ (1/varE x_j' x_j + 1/(varE tau_j^2))
    variance= 1/ (1/varE x_j' x_j + 1/(varE tau_j^2))
    
    xj= the jth column of the incidence matrix
    
 *The notation in the routine is as follows:
 
 n: Number of rows in X
 pL: Number of columns in X
 XL: the matrix X stacked by columns
 XL2: vector with x_j' x_j, j=1,...,p
 bL: vector of regression coefficients
 e: vector with residuals, e=y-yHat, yHat= predicted values
 varBj: vector with variances, 
	For Bayesian LASSO, varBj=tau_j^2 * varE, j=1,...,p
	For Ridge regression, varBj=varB, j=1,...,p, varB is the variance common to all betas.
	For BayesA, varBj=varB_j, j=1,...,p
	For BayesCpi, varBj=varB, j=1,...,p, varB is the variance common to all betas
	
 varE: residual variance
 minAbsBeta: in some cases values of betas near to zero can lead to numerical problems in BL, 
             so, instead of using this tiny number we assingn them minAbsBeta
 
 */


SEXP sample_beta_groups(SEXP n, SEXP pL, SEXP XL, SEXP xL2, SEXP bL, SEXP e, SEXP varBj, SEXP varE, SEXP minAbsBeta, SEXP groups, SEXP nGroups)
{
    double *xj;
    double *xj2;
    double *pXL, *pXL2, *pbL, *pe, *pvarBj;
    double b;
    int inc=1;
    double smallBeta;
    int i, j, k, rows, cols;
    int ngroups;  //number of groups
    int *g;       //pointer for holding groups
    double *rhs, c, *sigma2e; //now the rhs, c and sigma2e are vectors, we have to take the sum over the corresponding elements.
    double sum_rhs;

    SEXP list;	

    GetRNGstate();
	
    rows=INTEGER_VALUE(n);
    cols=INTEGER_VALUE(pL);
    smallBeta=NUMERIC_VALUE(minAbsBeta);
	
    PROTECT(XL=AS_NUMERIC(XL));
    pXL=NUMERIC_POINTER(XL);

    PROTECT(xL2=AS_NUMERIC(xL2));
    pXL2=NUMERIC_POINTER(xL2);

    PROTECT(bL=AS_NUMERIC(bL));
    pbL=NUMERIC_POINTER(bL);

    PROTECT(e=AS_NUMERIC(e));
    pe=NUMERIC_POINTER(e);

    PROTECT(varBj=AS_NUMERIC(varBj));
    pvarBj=NUMERIC_POINTER(varBj);

    sigma2e=NUMERIC_POINTER(varE);
    ngroups=INTEGER_VALUE(nGroups);
    g=INTEGER_POINTER(groups);

    //rhs =(double *) malloc(ngroups);
    rhs=(double *) R_alloc(ngroups, sizeof(double));
    
    for(j=0; j<cols;j++)
    {
	  sum_rhs=0;

          for(k=0;k<ngroups;k++) rhs[k]=0;

          c=0;

          xj=pXL+j*rows;
          b=pbL[j];

          //rhs=F77_NAME(ddot)(&rows,xj,&inc,pe,&inc)/sigma2e;
          //rhs+=pxL2[j]*b/sigma2e;
  	  //c=pxL2[j]/sigma2e + 1.0/pvarBj[j];
	  //pbL[j]=rhs/c + sqrt(1.0/c)*norm_rand();

          /*
          This code works well but it is not efficient
	  for(i=0;i<rows;i++)
	  {
		rhs+=xj[i]*pe[i]/sigma2e[g[i]-1]+b*xj[i]*xj[i]/sigma2e[g[i]-1];
                c+=xj[i]*xj[i]/sigma2e[g[i]-1]; 
          }

          c+=1.0/pvarBj[j];
          pbL[j]=rhs/c + sqrt(1.0/c)*norm_rand();
          */
          
          //This code works well but perhaps we can do it more efficiently with 
          //a call to ddot 
          //
          
          for(i=0;i<rows;i++)
	  {
              /*rhs+=xj[i]*pe[i]/sigma2e[g[i]-1];*/
              rhs[g[i]]+=xj[i]*pe[i];
          }

          //weighted_ddot(rows, xj, pe, g, rhs);

          xj2=pXL2+j*ngroups;

          for(k=0;k<ngroups;k++)
	  {
             //rhs+=b*xj2[k]/sigma2e[k];
             sum_rhs+=(rhs[k]+b*xj2[k])/sigma2e[k];
             c+=xj2[k]/sigma2e[k];
          }
          c+=1.0/pvarBj[j];
          pbL[j]=sum_rhs/c + sqrt(1.0/c)*norm_rand();
          
          //End of code for groups

          b-=pbL[j];  

          F77_NAME(daxpy)(&rows, &b,xj,&inc, pe,&inc);
          
          if(fabs(pbL[j])<smallBeta)
          {
             pbL[j]=smallBeta;
          }
      }
        
      // Creating a list with 2 vector elements:
      PROTECT(list = allocVector(VECSXP, 2));
      // attaching bL vector to list:
      SET_VECTOR_ELT(list, 0, bL);
      // attaching e vector to list:
      SET_VECTOR_ELT(list, 1, e);

      PutRNGstate();

      UNPROTECT(6);

      //free(rhs);

      return(list);
}
