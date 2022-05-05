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
        xj=pX+(long long)j*rows;
        
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
        xj=pX+(long long)j*rows;
        
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

/*

Bayesian Ridge Regression version 2

*/

SEXP sampler_BRR_mt_v2(SEXP n,
                       SEXP nCol,
                       SEXP nTraits,
                       SEXP Rinv,
                       SEXP X,
                       SEXP e,
                       SEXP beta,
                       SEXP x2,
                       SEXP Omegainv)
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
    double *pOmegainv;      //Pointer to Omegainv
    double *cOmegainv;      //Pointer to column of Omegainv
    double *xe;             //Pointer to [x_j'*e_1,...,x_j'*e_traits]
    double mu;
    double v;
    double Omegainversekk;
    
    int j, p, k, traits, t, rows;
    
    int inc=1;
    
    p=INTEGER_VALUE(nCol);
    traits=INTEGER_VALUE(nTraits);
    rows=INTEGER_VALUE(n);
    
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
    
    PROTECT(Omegainv=AS_NUMERIC(Omegainv));
    pOmegainv=NUMERIC_POINTER(Omegainv);
    
    //This memory is taken from the heap, and released at the end of the .Call
    xe=(double *) R_alloc(traits, sizeof(double));
    
    GetRNGstate();
    
    //Main loop, level 1
    for(j=0; j<p;j++)
    {
        
        //j-th column of X
        xj=pX+(long long)j*rows;
        
        //Compute [x_j'*e_1,...,x_j'*e_traits]
        for(t=0; t<traits;t++)
        {
        	//t-th column of error
        	cerror=perror+t*rows;
        	xe[t]=F77_NAME(ddot)(&rows,xj,&inc,cerror,&inc);
        }
        
        //loop, level 2
        for(k=0;k<traits;k++)
        {
        	 s1=0;
        	 for(t=0; t<traits;t++)
        	 {
        	 	 //t-th column of Rinv
            	 	cRinv=pRinv+t*traits;
            	 	s1+=cRinv[k]*xe[t];
        	 }
        	 
        	 //k-th column of Rinv, error, Omegainv and beta
        	 cRinv=pRinv+k*traits;
        	 cerror=perror+k*rows;
    		 cOmegainv=pOmegainv+k*traits;
    		 cbeta=pbeta+k*p;
    		 
    		 //Sample from the conditional
        	 betaOld=cbeta[j];
        	 s4=0.0;
        	  
       		 for(t=0;t<traits;t++)
       		 {
			if(t!=k)
			{
				cbeta=pbeta+t*p;
				s4+=cbeta[j]*cOmegainv[t];
			}
       		 }
       		  
       		 //Reset pointer to k-th column of cbeta
        	 cbeta=pbeta+k*p;
        	
        	 s3=s1+cRinv[k]*cbeta[j]*px2[j]-s4;
        	 Omegainversekk=cOmegainv[k];
        	 v=cRinv[k]*px2[j]+Omegainversekk;
        	 mu=s3/v;
            
        	 cbeta[j]=mu+sqrt(1.0/v)*norm_rand();
        	
		 //Update the error, e_k and xe_k
                 //e_k = shift*xj + e_k
                 //once we update e_k, we can compute xe_k = xj' * e_k, but it is more efficient
                 //if we update first xe_k = xj' * (shift * xj + e_k)=shift*xj'*xj + xj'*e_k
                 //so that we avoid a call to ddot

                 shift=betaOld-cbeta[j];
                 xe[k]=shift*px2[j]+xe[k];
                 F77_NAME(daxpy)(&rows, &shift,xj,&inc, cerror, &inc);
 
        	 //Update the error
        	 //shift=betaOld-cbeta[j];
        	 //F77_NAME(daxpy)(&rows, &shift,xj,&inc, cerror, &inc);
        	 
        	 //update xe
        	 //xe[k]=F77_NAME(ddot)(&rows,xj,&inc,cerror,&inc);
        	 
        }
        
    }
    
    PutRNGstate();
    
    UNPROTECT(6);
    
    return(R_NilValue);
}

/*

Auxiliary function to compute:

S11 <- ETA[[j]]$Cov$Omega[k,k,drop=FALSE]
S22 <- ETA[[j]]$Cov$Omega[-k,-k,drop=FALSE]
S12 <- ETA[[j]]$Cov$Omega[k,-k,drop=FALSE]
tmp12 <- as.vector(S12%*%solve(S22))
tmp11 <- tmp12%*%t(S12)
sigma2 <- as.numeric(S11-tmp11)

once that it finishes, returns sigma2 and tmp12 has the information that is needed

*/

double tmp12_sigma2(double *pOmega, int traits, int k, double *tmp12)
{
	int neq;
	int i,j,m;
	int *ipiv;
	int info;
	double *cOmega; 	//Pointer to column of Omega
	
	neq=traits-1;
    
    	double S11;
    	double *S22;
    	double *S12;
    	double tmp11;
    	double *Identity;
    	double *cIdentity;   //Pointer to column of Identity   
    
    	double sigma2;
    	double s;
    
    	//This memory is taken from the heap, and released at the end of the .Call
    	S22=(double *) R_alloc((traits-1)*(traits-1), sizeof(double));
    	S12=(double *) R_alloc((traits-1),sizeof(double));
    	Identity=(double *) R_alloc((traits-1)*(traits-1), sizeof(double));  
    	ipiv = (int *) R_alloc(neq, sizeof(int)); 
    
    
    	//k-th column of Omega
    	cOmega=pOmega+k*traits;
    	S11=cOmega[k];
    	//Rprintf("************************************************\n");
    	//Rprintf("k=%d\n",k);
    	//Rprintf("S11=%f\n",S11);
    	
    	//fill S12
    	//Rprintf("S12=\n");
    	m=0;
    	for(j=0;j<traits;j++)
    	{
    		if(j!=k)
    		{
    			S12[m]=cOmega[j];
    			//Rprintf("%f\n",S12[m]);
    			m+=1;
    		}
    	}
    	
    	//fill S22
    	//Rprintf("S22=\n");
    	m=0;
    	for(i=0;i<traits;i++)
    	{	
    		//i-th column of Omega
    		cOmega=pOmega+i*traits;
    		for(j=0;j<traits;j++)
    		{
    			if(i!=k && j!=k)
    			{
    				S22[m]=cOmega[j];
    				m+=1;	
    			}
    		}
    	}
    	
    	//for(m=0; m<(traits-1)*(traits-1);m++)
    	//{
    	//	Rprintf("%f\n",S22[m]);
    	//}
    	
    	// Obtain the inverse of S22
    	// DGESV - compute the solution to a real system of linear 
	// equations  A * X = B
			
	//F77_NAME(dgesv)(const int* n, const int* nrhs, double* a, const int* lda,
	//int* ipiv, double* b, const int* ldb, int* info);
		
	//Fill the identity matrix
	m=0;
    	for(i=0;i<traits-1;i++)
    	{
    		for(j=0;j<traits-1;j++)
    		{
    			if(i==j)
    			{
    				Identity[m]=1.0;
    			}else{
    				Identity[m]=0.0;
    			}
    			m+=1;
    		}
    	}
    	
    	//Inverse is stored in the Identity matrix after successful solution...
    	F77_NAME(dgesv)(&neq,&neq,S22,&neq,ipiv,Identity,&neq,&info);
    	
    	if(info < 0){
    		error("argument %d of Lapack routine %s had invalid value",-info);
    	}
    	
	if(info > 0){
		error("Lapack routine %s: system is exactly singular: U[%d,%d] = 0", info, info);	
	}
	
    	//for(m=0; m<(traits-1)*(traits-1);m++)
    	//{
    	//	Rprintf("%f\n",Identity[m]);
    	//}
    	
    	//S12*S22inv, S22inv is stored in the Identity matrix
    	//tmp11=tmp12*S12'
    	tmp11=0;
    	for(i=0;i<traits-1;i++)
    	{
    	    //i-th column of Identity/S22inv
    	    cIdentity=Identity+i*(traits-1);
    	    s=0;
    	 	   
    	    for(j=0; j<traits-1; j++)
    	    {
    	    	s+=S12[j]*cIdentity[j];
    	    }
    	    
    	    tmp12[i] = s;
    	    
    	    tmp11+=tmp12[i]*S12[i];
    	}
    	
    	//Rprintf("tmp12=\n");
    	//for(i=0;i<traits-1;i++)
    	//{
    	//	Rprintf("%f\n",tmp12[i]);
    	//}
    	
    	//Rprintf("tmp11=%f\n",tmp11);
    	
    	//sigma2=S11-tmp11
    	sigma2=S11-tmp11;
    	//Rprintf("sigma2=%f\n",sigma2);
    	
    	return(sigma2);
}

/*

Dirac Spikes and Slab version 2

*/

SEXP sampler_DiracSS_mt_v2(SEXP lpo,
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
                           SEXP Omega,
                           SEXP Omegainv)
{
    double logOdds, s1, s2, s3, s4, u, theta, betaOld, shift, sigma2;
    double *plogPriorOdds;  //Pointer to logPriorOdds;
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
    double *pOmega;	    //Pointer to Omega
    double *pOmegainv;      //Pointer to Omegainv
    double *cOmegainv;	    //Pointer to column of Omegainv
    double *xe;             //Pointer to [x_j'*e_1,...,x_j'*e_traits]
    double *tmp12;	    //Pointer to tmp12=Omega[k-k]*(Omega[-k,-k])^{-1}
    double *Mtmp12;         //Pointer to vectors that store values in tmp12
    double *vsigma2;        //Pointer to vectors that store sigma2 values
    
    
    double Omegainversekk;
    double mu;
    double v;
    
    
    int j, p, k, traits, t, rows, indicator, m;
    
    int inc=1;
	
	
    PROTECT(lpo=AS_NUMERIC(lpo));
    plogPriorOdds=NUMERIC_POINTER(lpo);
    
    p=INTEGER_VALUE(nCol);
    traits=INTEGER_VALUE(nTraits);
    rows=INTEGER_VALUE(n);
    
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
    
    PROTECT(Omega=AS_NUMERIC(Omega));
    pOmega=NUMERIC_POINTER(Omega);
    
    PROTECT(Omegainv=AS_NUMERIC(Omegainv));
    pOmegainv=NUMERIC_POINTER(Omegainv);
    
    //This memory is taken from the heap, and released at the end of the .Call
    xe=(double *) R_alloc(traits, sizeof(double));
    
    tmp12 = (double *) R_alloc((traits-1),sizeof(double));
    Mtmp12= (double *) R_alloc((traits-1)*traits,sizeof(double));
    vsigma2= (double *) R_alloc(traits,sizeof(double));
    
    
    //Precompute sigma2 and tmp12 for all the traits
    m=0;    
    for(k=0;k<traits;k++)
    {
    	sigma2=tmp12_sigma2(pOmega,traits, k, tmp12);
    	//copy values
    	vsigma2[k]=sigma2;
    	for(j=0;j<traits-1;j++)
    	{
    		Mtmp12[m]=tmp12[j];
    		m+=1;
    	}
    }
        
    GetRNGstate();
    
    
    //Main loop, level 1
    for(j=0; j<p;j++)
    {
    
    	//j-th column of X
        xj=pX+(long long)j*rows;
        
        //Compute [x_j'*e_1,...,x_j'*e_traits]
        for(t=0; t<traits;t++)
        {
        	//t-th column of error
        	cerror=perror+t*rows;
        	xe[t]=F77_NAME(ddot)(&rows,xj,&inc,cerror,&inc);
        }
        
        //loop, level 2
        for(k=0;k<traits;k++)
        {
        	//k-th row of Mtmp12
        	tmp12=Mtmp12+k*(traits-1);
        	sigma2=vsigma2[k];
        	
		s1=0;
        	for(t=0; t<traits;t++)
        	{
            		//t-th column of Rinv
            		cRinv=pRinv+t*traits;
            		s1+=cRinv[k]*xe[t];
        	}
        	
        	//k-th column of Rinv, error, beta, b, d
        	cRinv=pRinv+k*traits;
        	cerror=perror+k*rows;
        	cbeta=pbeta+k*p;
        	cb=pb+k*p;
        	cd=pd+k*p;
        	
        	//k-th column of Omegainv
        	cOmegainv=pOmegainv+k*traits;
        	
        	//e = y-X*beta = e* - x_j*beta_j, so e* = error + x_j * beta_j
        	//we can not use  e* = e + x_j * b_j, because error was computed based on beta
        	//in previous iteration, and if d_j was 0 beta_j = b_j * d_j = 0
        
        	s2=cb[j]*s1 + cb[j]*cbeta[j]*cRinv[k]*px2[j] - 0.5*cRinv[k]*cb[j]*cb[j]*px2[j];
        
        	logOdds=s2+plogPriorOdds[k];
             
        	u=unif_rand();
        	theta=1.0/(1.0+exp(-logOdds));
             
        	indicator = u<theta ? 1 : 0;
        
        	cd[j]=indicator;
        
        	if(indicator)
        	{
            		//Sample from the conditional
            		s4=0.0;
            		for(t=0;t<traits;t++)
            		{
                		if(t!=k)
                		{
                    			cb=pb+t*p;
                    			s4+=cb[j]*cOmegainv[t];
                		}
            		}
            
            		//Reset pointer to k-th column of cb
            		cb=pb + k*p;
            		s3=s1 + cbeta[j]*cRinv[k]*px2[j] - s4;
            		Omegainversekk=cOmegainv[k];
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
                    			mu+=cb[j]*tmp12[m];
                    			m+=1;
                		}
            		}
                       
            		//Reset pointer to k-th column of b
            		cb=pb+k*p;
            		cb[j]=mu+sqrt(sigma2)*norm_rand();
            
        	}
        
        	betaOld=cbeta[j];
        
        	cbeta[j]=cd[j]*cb[j];

		//Update the error, e_k and xe_k
                //e_k = shift*xj + e_k
                //once we update e_k, we can compute xe_k = xj' * e_k, but it is more efficient
                //if we update first xe_k = xj' * (shift * xj + e_k)=shift*xj'*xj + xj'*e_k
                //so that we avoid a call to ddot
                
                shift=betaOld-cbeta[j];
                xe[k]=shift*px2[j]+xe[k];
                F77_NAME(daxpy)(&rows, &shift,xj,&inc, cerror, &inc);

        	//Update the error
        	//shift=betaOld-cbeta[j];
        	//F77_NAME(daxpy)(&rows, &shift,xj,&inc, cerror, &inc);
        	
        	//Update xe
        	//xe[k]=F77_NAME(ddot)(&rows,xj,&inc,cerror,&inc);
        	
        } //loop for k
        
    } //loop for j
    
    PutRNGstate();
    
    UNPROTECT(10);
    
    return(R_NilValue);
}

/*
See tmp12_sigma2 routine, this routine was developed for testing purposes.
*/

/*
SEXP test_computations(SEXP Omega, SEXP nTraits)
{
	
    int traits;
    int neq;
    int i,j,k,m;
    int *ipiv;
    int info;
    double *pOmega;	//Pointer to matrix
    double *cOmega; 	//Pointer to column of Omega
	
    traits=INTEGER_VALUE(nTraits);
    neq=traits-1;
    PROTECT(Omega=AS_NUMERIC(Omega));
    pOmega=NUMERIC_POINTER(Omega);
    
    double S11;
    double *S22;
    double *S12;
    double *tmp12;
    double tmp11;
    double *Identity;
    double *cIdentity;   //Pointer to column of Identity   
    
    double sigma2;
    double s;
    
    //This memory is taken from the heap, and released at the end of the .Call
    S22=(double *) R_alloc((traits-1)*(traits-1), sizeof(double));
    S12=(double *) R_alloc((traits-1),sizeof(double));
    Identity=(double *) R_alloc((traits-1)*(traits-1), sizeof(double));  
    ipiv = (int *) R_alloc(neq, sizeof(int)); 
    tmp12 = (double *) R_alloc((traits-1),sizeof(double));
    
    for(k=0; k<traits;k++)
    {
    	//k-th column of Omega
    	cOmega=pOmega+k*traits;
    	S11=cOmega[k];
    	Rprintf("************************************************\n");
    	Rprintf("k=%d\n",k);
    	Rprintf("S11=%f\n",S11);
    	
    	//fill S12
    	//Rprintf("S12=\n");
    	m=0;
    	for(j=0;j<traits;j++)
    	{
    		if(j!=k)
    		{
    			S12[m]=cOmega[j];
    			//Rprintf("%f\n",S12[m]);
    			m+=1;
    		}
    	}
    	
    	//fill S22
    	//Rprintf("S22=\n");
    	m=0;
    	for(i=0;i<traits;i++)
    	{	
    		//i-th column of Omega
    		cOmega=pOmega+i*traits;
    		for(j=0;j<traits;j++)
    		{
    			if(i!=k && j!=k)
    			{
    				S22[m]=cOmega[j];
    				m+=1;	
    			}
    		}
    	}
    	
    	//for(m=0; m<(traits-1)*(traits-1);m++)
    	//{
    	//	Rprintf("%f\n",S22[m]);
    	//}
    	
    	// Obtain the inverse of S22
    	// DGESV - compute the solution to a real system of linear 
	// equations  A * X = B
			
	//F77_NAME(dgesv)(const int* n, const int* nrhs, double* a, const int* lda,
	//int* ipiv, double* b, const int* ldb, int* info);
		
	//Fill the identity matrix
	m=0;
    	for(i=0;i<traits-1;i++)
    	{
    		for(j=0;j<traits-1;j++)
    		{
    			if(i==j)
    			{
    				Identity[m]=1.0;
    			}else{
    				Identity[m]=0.0;
    			}
    			m+=1;
    		}
    	}
    	
    	//Inverse is stored in the Identity matrix after successful solution...
    	F77_NAME(dgesv)(&neq,&neq,S22,&neq,ipiv,Identity,&neq,&info);
    	
    	if(info < 0){
    		error("argument %d of Lapack routine %s had invalid value",-info);
    	}
    	
	if(info > 0){
		error("Lapack routine %s: system is exactly singular: U[%d,%d] = 0", info, info);	
	}
	
    	//for(m=0; m<(traits-1)*(traits-1);m++)
    	//{
    	//	Rprintf("%f\n",Identity[m]);
    	//}
    	
    	//S12*S22inv, S22inv is stored in the Identity matrix
    	//tmp11=tmp12*S12'
    	tmp11=0;
    	for(i=0;i<traits-1;i++)
    	{
    	    //i-th column of Identity/S22inv
    	    cIdentity=Identity+i*(traits-1);
    	    s=0;
    	 	   
    	    for(j=0; j<traits-1; j++)
    	    {
    	    	s+=S12[j]*cIdentity[j];
    	    }
    	    
    	    tmp12[i] = s;
    	    
    	    tmp11+=tmp12[i]*S12[i];
    	}
    	
    	//Rprintf("tmp12=\n");
    	//for(i=0;i<traits-1;i++)
    	//{
    	//	Rprintf("%f\n",tmp12[i]);
    	//}
    	
    	//Rprintf("tmp11=%f\n",tmp11);
    	
    	//sigma2=S11-tmp11
    	sigma2=S11-tmp11;
    	Rprintf("sigma2=%f\n",sigma2);
    	
    }
    
    UNPROTECT(1);
    
    return(R_NilValue);
}

*/

/*
SEXP test_computations2(SEXP trait,
                 		SEXP nTraits,
                 		SEXP rOmegainv)
{
	
    double *prOmegainv;     //Pointer to row of Omegainv
    
    int k,traits, t, m;
    
    k=INTEGER_VALUE(trait)-1;       //In C we begin to count in 0
    traits=INTEGER_VALUE(nTraits);
    
    PROTECT(rOmegainv=AS_NUMERIC(rOmegainv));
    prOmegainv=NUMERIC_POINTER(rOmegainv);
    
    m=0;
        
    for(t=0;t<traits;t++)
    {
        if(t!=k)
        {
            
        	Rprintf("%f\n",prOmegainv[m]);
            m+=1;
        }
    }
    
    UNPROTECT(1);
    
    return(R_NilValue);
}

*/

/*

SEXP test_computations3(SEXP nTraits,
                        SEXP Omegainv)
{
	double *pOmegainv;      //Pointer to Omegainv
    double *cOmegainv;      //Pointer to column of Omegainv
    
    int k, traits, t;
    
    traits=INTEGER_VALUE(nTraits);
    
    PROTECT(Omegainv=AS_NUMERIC(Omegainv));
    pOmegainv=NUMERIC_POINTER(Omegainv);
    
    for(k=0;k<traits;k++)
    {
    	cOmegainv=pOmegainv+k*traits;
    	
    	for(t=0;t<traits;t++)
       	{
			if(t!=k)
			{
						
				Rprintf("%f\n",cOmegainv[t]);
			}
       	}
    }
    
    
    UNPROTECT(1);
    
    return(R_NilValue);
}

*/
