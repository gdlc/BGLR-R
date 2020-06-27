#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rconfig.h>
#include <R_ext/Lapack.h>
 
/*

Dirac Spikes and Slab

*/

SEXP sampler_DiracSS_mt(SEXP trait,
                        SEXP lpo,
                        SEXP n,
                        SEXP nCol,
                        SEXP nTraits,
                        SEXP Rinv,
                        SEXP X,
                        SEXP e,
                        SEXP beta,
                        SEXP b,
                        SEXP d,
                        SEXP x2,
                        SEXP tmp12,
                        SEXP Sigma2,
                        SEXP rOmegainv,
                        SEXP Oikk)
{
	double logPriorOdds, logOdds, s1, s2, s3, s4, u, theta, betaOld, shift, sigma2;
    double *pRinv;          //Pointer to Rinv
    double *cRinv;          //Pointer to a column of Rinv
    double *pX;             //Pointer to X
    double *xj;             //Pointer to a columns of X
    double *perror;         //Pointer to error
    double *cerror;         //Pointer to column of error
    double *pbeta;          //Pointer to beta matrix
    double *cbeta;          //Pointer to a column of beta
    double *pb;             //Pointer to b matrix
    double *cb;             //Pointer to a column of b matrix
    double *pd;             //Pointer to d
    double *cd;             //Pointer to a column of d
    double *px2;            //Pointer to x2
    double *ptmp12;         //Pointer to tmp12
    double *prOmegainv;     //Pointer to row of Omegainv
    double mu;
    double v;
    double Omegainversekk;
    
    int j, p, k, traits, t, rows, indicator, m;
    
    int inc=1;
	
	logPriorOdds=NUMERIC_VALUE(lpo);
    p=INTEGER_VALUE(nCol);
    traits=INTEGER_VALUE(nTraits);
    k=INTEGER_VALUE(trait)-1;       //In C we begin to count in 0
    rows=INTEGER_VALUE(n);
    sigma2=NUMERIC_VALUE(Sigma2);
    Omegainversekk=NUMERIC_VALUE(Oikk);
    
    PROTECT(Rinv=AS_NUMERIC(Rinv));
    pRinv=NUMERIC_POINTER(Rinv);
    
    PROTECT(X=AS_NUMERIC(X));
    pX=NUMERIC_POINTER(X);
    
    PROTECT(e=AS_NUMERIC(e));
    perror=NUMERIC_POINTER(e);
    
    PROTECT(beta=AS_NUMERIC(beta));
    pbeta=NUMERIC_POINTER(beta);
    
    PROTECT(b=AS_NUMERIC(b));
    pb=NUMERIC_POINTER(b);
    
    PROTECT(d=AS_NUMERIC(d));
    pd=NUMERIC_POINTER(d);
    
    PROTECT(x2=AS_NUMERIC(x2));
    px2=NUMERIC_POINTER(x2);
    
    PROTECT(tmp12=AS_NUMERIC(tmp12));
    ptmp12=NUMERIC_POINTER(tmp12);
    
    PROTECT(rOmegainv=AS_NUMERIC(rOmegainv));
    prOmegainv=NUMERIC_POINTER(rOmegainv);
    
    //k-th column of beta
    cbeta=pbeta+k*p;
    
    //k-th column of b
    cb=pb+k*p;
    
    //k-th column of d
    cd=pd+k*p;
    
    GetRNGstate();
    
    for(j=0; j<p;j++)
    {
        
        //j-th column of X
        xj=pX+j*rows;
        
        s1=0;
        for(t=0; t<traits;t++)
        {
            //t-th column of Rinv and error
            cRinv=pRinv+t*traits;
            cerror=perror+t*rows;
            
            s1+=cRinv[k]*F77_NAME(ddot)(&rows,xj,&inc,cerror,&inc);
        }
        
        //k-th column of Rinv and error
        cRinv=pRinv+k*traits;
        cerror=perror+k*rows;
        
        //e = y-X*beta = e* - x_j*beta_j, so e* = error + x_j * beta_j
        //we can not use  e* = e + x_j * b_j, because error was computed based on beta
        //in previous iteration, and if d_j was 0 beta_j = b_j * d_j = 0
        
        s2=cb[j]*s1 + cb[j]*cbeta[j]*cRinv[k]*px2[j] - 0.5*cRinv[k]*cb[j]*cb[j]*px2[j];
        
        logOdds=s2+logPriorOdds;
             
        u=unif_rand();
        theta=1.0/(1.0+exp(-logOdds));
             
        indicator = u<theta ? 1 : 0;
        
        cd[j]=indicator;
        
        if(indicator)
        {
            //Sample from the conditional
            s4=0.0;
            m=0;
            for(t=0;t<traits;t++)
            {
                if(t!=k)
                {
                    cb=pb+t*p;
                    s4+=cb[j]*prOmegainv[m];
                    m+=1;
                }
            }
            
            //Reset pointer to k-th column of cb
            cb=pb + k*p;
            s3=s1 + cbeta[j]*cRinv[k]*px2[j] - s4;
            v=cRinv[k]*px2[j] + Omegainversekk;
            mu=s3/v;
                
            cb[j]=mu+sqrt(1.0/v)*norm_rand();
            
        }else{
            
            //Sample from the prior
            
            mu=0.0;
            m=0;
            for(t=0;t<traits;t++)
            {
                if(t!=k)
                {
                    cb=pb+t*p;
                    mu+=cb[j]*ptmp12[m];
                    m+=1;
                }
            }
                       
            //Reset pointer to k-th column of b
            cb=pb+k*p;
            cb[j]=mu+sqrt(sigma2)*norm_rand();
            
        }
        
        betaOld=cbeta[j];
        
        cbeta[j]=cd[j]*cb[j];
        
        //Update the error
        shift=betaOld-cbeta[j];
        F77_NAME(daxpy)(&rows, &shift,xj,&inc, cerror, &inc);
        
    }
    
    PutRNGstate();
    
    UNPROTECT(9);
    
    return(R_NilValue);
}


/*

Bayesian Ridge Regression

*/

SEXP sampler_BRR_mt(SEXP trait,
                 SEXP n,
                 SEXP nCol,
                 SEXP nTraits,
                 SEXP Rinv,
                 SEXP X,
                 SEXP e,
                 SEXP beta,
                 SEXP x2,
                 SEXP rOmegainv,
                 SEXP Oikk)
{
	double s1, s3, s4, betaOld, shift;
    double *pRinv;          //Pointer to Rinv
    double *cRinv;          //Pointer to a column of Rinv
    double *pX;             //Pointer to X
    double *xj;             //Pointer to a columns of X
    double *perror;         //Pointer to error
    double *cerror;         //Pointer to column of error
    double *pbeta;          //Pointer to beta matrix
    double *cbeta;          //Pointer to a column of beta
    double *px2;            //Pointer to x2
    double *prOmegainv;     //Pointer to row of Omegainv
    double mu;
    double v;
    double Omegainversekk;
    
    int j, p, k, traits, t, rows, m;
    
    int inc=1;
    
    p=INTEGER_VALUE(nCol);
    traits=INTEGER_VALUE(nTraits);
    k=INTEGER_VALUE(trait)-1;       //In C we begin to count in 0
    rows=INTEGER_VALUE(n);
    Omegainversekk=NUMERIC_VALUE(Oikk);
    
    PROTECT(Rinv=AS_NUMERIC(Rinv));
    pRinv=NUMERIC_POINTER(Rinv);
    
    PROTECT(X=AS_NUMERIC(X));
    pX=NUMERIC_POINTER(X);
    
    PROTECT(e=AS_NUMERIC(e));
    perror=NUMERIC_POINTER(e);
    
    PROTECT(beta=AS_NUMERIC(beta));
    pbeta=NUMERIC_POINTER(beta);
        
    PROTECT(x2=AS_NUMERIC(x2));
    px2=NUMERIC_POINTER(x2);
    
    PROTECT(rOmegainv=AS_NUMERIC(rOmegainv));
    prOmegainv=NUMERIC_POINTER(rOmegainv);
    
    //k-th column of beta
    cbeta=pbeta+k*p;
    
    GetRNGstate();
    
    for(j=0; j<p;j++)
    {
        
        //j-th column of X
        xj=pX+j*rows;
        
        s1=0;
        for(t=0; t<traits;t++)
        {
            //t-th column of Rinv and error
            cRinv=pRinv+t*traits;
            cerror=perror+t*rows;
            s1+=cRinv[k]*F77_NAME(ddot)(&rows,xj,&inc,cerror,&inc);
        }
        
        //k-th column of Rinv and error
        cRinv=pRinv+k*traits;
        cerror=perror+k*rows;
            
        //Sample from the conditional
        betaOld=cbeta[j];
        s4=0.0;
        m=0;
        
        for(t=0;t<traits;t++)
        {
                if(t!=k)
                {
                    cbeta=pbeta+t*p;
                    s4+=cbeta[j]*prOmegainv[m];
                    m+=1;
                }
        }
            
        //Reset pointer to k-th column of cbeta
        cbeta=pbeta+k*p;
        s3=s1+cRinv[k]*cbeta[j]*px2[j]-s4;
        v=cRinv[k]*px2[j]+Omegainversekk;
        mu=s3/v;
            
        cbeta[j]=mu+sqrt(1.0/v)*norm_rand();
        
        //Update the error
        
        shift=betaOld-cbeta[j];
        F77_NAME(daxpy)(&rows, &shift,xj,&inc, cerror, &inc);
        
    }
    
    PutRNGstate();
    
    UNPROTECT(6);
    
    return(R_NilValue);
}
