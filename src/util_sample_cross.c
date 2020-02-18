#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rconfig.h>
#include <R_ext/Lapack.h>
 
/*

Dirac Spikes and Slab

nCol the number of columns in X'X
XX matrix X'X
XY matrix X'y
idColumns A subset of columns of X'X used for updating b, beta and d
length the length of vector idColumns

*/

SEXP sampler_DiracSS(SEXP nCol, 
         SEXP XX, 
         SEXP XY,
         SEXP idColumns,
         SEXP length, 
         SEXP b, 
         SEXP beta, 
         SEXP d, 
         SEXP varB, 
         SEXP varE, 
         SEXP probIn, 
         SEXP RSS)
{

    double Cjj, lhs, rhs,  sol, offset, vare, z, u, logOdds, logPriorOdds,logP, probin, old_beta;
    double *pXX, *pXY, *pb, *pd, *pvarB, *pbeta, *pRSS;
    int *pidColumns;
 
    int inc=1;
    int p;
    int q;
    int j;
    ptrdiff_t j_global;
 
    SEXP list;
 
    GetRNGstate();
 
    p=INTEGER_VALUE(nCol);
    //Rprintf("p=%d\n",p);
    q=INTEGER_VALUE(length);
    //Rprintf("q=%d\n",q);
 
    vare=NUMERIC_VALUE(varE);
    probin=NUMERIC_VALUE(probIn);
    
    
    PROTECT(XX=AS_NUMERIC(XX));
    pXX=NUMERIC_POINTER(XX);
 
    PROTECT(XY=AS_NUMERIC(XY));
    pXY=NUMERIC_POINTER(XY);
 
    PROTECT(b=AS_NUMERIC(b));
    pb=NUMERIC_POINTER(b);
 
    PROTECT(d=AS_NUMERIC(d));
    pd=NUMERIC_POINTER(d);

    PROTECT(beta=AS_NUMERIC(beta)); // beta=b*D
    pbeta=NUMERIC_POINTER(beta);

    PROTECT(varB=AS_NUMERIC(varB));
    pvarB=NUMERIC_POINTER(varB);

    PROTECT(RSS=AS_NUMERIC(RSS));
    pRSS=NUMERIC_POINTER(RSS);
    
    PROTECT(idColumns=AS_INTEGER(idColumns));
    pidColumns=INTEGER_POINTER(idColumns);
    
    logPriorOdds=log(probin/(1-probin));    

	/*
	
	This loop is for the idsColumns
	
	*/   
	 
     for(j=0; j<q;j++)
     {
          
          j_global=pidColumns[j]-1;  //Substract one because in C we count in zero
          Cjj=pXX[j_global*(p+1)];
          rhs=pXY[j_global];
          
          old_beta=pbeta[j_global];
          
          offset=F77_NAME(ddot)(&p, pXX+j_global*p, &inc, pbeta, &inc);

          offset-=Cjj*pbeta[j_global];
          
          lhs=Cjj+vare/pvarB[j];
          
          // random variables
          z=norm_rand();
          u=unif_rand();
          logOdds=log(u/(1-u));
              	
              	
         // include the effect in the model?
	 	logP= -(0.5/vare)*( Cjj*pow(pb[j],2)-2*pb[j]*(rhs-offset) ) +  logPriorOdds ;
   
         if( logP>logOdds? 1:0){
            pd[j]=1;
            sol=(rhs-offset)/lhs;
            pb[j]=sol+sqrt(vare/lhs)*z;
            pbeta[j_global]=pb[j]  ;       	
         }else{    
         	pd[j]=0;
         	pb[j]=z*sqrt(pvarB[j]);
         	pbeta[j_global]=0;
         }
         
         pRSS[0]+=(pow(pbeta[j_global],2) - pow(old_beta,2))*Cjj  -2*(pbeta[j_global]-old_beta)*(rhs-offset);

    }
 
      // Creating a list with 1 vector elements:
      PROTECT(list = allocVector(VECSXP, 4));
      SET_VECTOR_ELT(list, 0, b);
      SET_VECTOR_ELT(list, 1, d);      
      SET_VECTOR_ELT(list, 2, beta);
      SET_VECTOR_ELT(list, 3, RSS);
       
      PutRNGstate();
 
      UNPROTECT(9);
 
      return(list);
}


SEXP sampler_others(SEXP nCol, 
                    SEXP XX, 
                    SEXP XY,
                    SEXP idColumns,
                    SEXP length, 
                    SEXP beta, 
                    SEXP varB, 
                    SEXP varE, 
                    SEXP RSS)
{

    double Cjj, lhs, rhs,  sol, offset, vare, z, old_beta;
    double *pXX, *pXY, *pvarB, *pbeta, *pRSS;
    int *pidColumns;
 
    int inc=1;
    int p;
    int q;
    int j;
    ptrdiff_t j_global;
 
    SEXP list;
 
    GetRNGstate();
 
    p=INTEGER_VALUE(nCol);
    q=INTEGER_VALUE(length);
 
    vare=NUMERIC_VALUE(varE);
       
    PROTECT(XX=AS_NUMERIC(XX));
    pXX=NUMERIC_POINTER(XX);
 
    PROTECT(XY=AS_NUMERIC(XY));
    pXY=NUMERIC_POINTER(XY);
 
    PROTECT(beta=AS_NUMERIC(beta)); 
    pbeta=NUMERIC_POINTER(beta);

    PROTECT(varB=AS_NUMERIC(varB));
    pvarB=NUMERIC_POINTER(varB);

    PROTECT(RSS=AS_NUMERIC(RSS));
    pRSS=NUMERIC_POINTER(RSS);
    
    PROTECT(idColumns=AS_INTEGER(idColumns));
    pidColumns=INTEGER_POINTER(idColumns);
    
 
	/*
	
	This loop is for the idsColumns
	
	*/   
	 
    for(j=0; j<q;j++)
    {
    	  /*
    	  Substract one because in C we begin to count in 0
    	  */
    	  
    	  j_global=pidColumns[j]-1;
    	  
          Cjj=pXX[j_global*(p+1)];					
          rhs=pXY[j_global];
          
          old_beta=pbeta[j_global];
          
          offset=F77_NAME(ddot)(&p, pXX+j_global*p, &inc, pbeta, &inc);

          offset-=Cjj*pbeta[j_global];
          
          lhs=Cjj+vare/pvarB[j];
          
          // random variables
          z=norm_rand();
          sol=(rhs-offset)/lhs;
          pbeta[j_global]=sol+sqrt(vare/lhs)*z;
         
          pRSS[0]+=(pow(pbeta[j_global],2) - pow(old_beta,2))*Cjj  -2*(pbeta[j_global]-old_beta)*(rhs-offset);

    }
 
      // Creating a list with 1 vector elements:
      PROTECT(list = allocVector(VECSXP, 2));      
      SET_VECTOR_ELT(list, 0, beta);
      SET_VECTOR_ELT(list, 1, RSS);
       
      PutRNGstate();
 
      UNPROTECT(7);
 
      return(list);
}

/*

Absolutely continuous Spikes and Slab

nCol the number of columns in X'X
XX matrix X'X
XY matrix X'y
idColumns A subset of columns of X'X used for updating a, beta and d
d here is the indicator variable, a=1 if d==1 and a=c of d==0 
length the length of vector idColumns

*/

SEXP sampler_ACSS(SEXP nCol, 
         SEXP XX, 
         SEXP XY,
         SEXP idColumns,
         SEXP length, 
         SEXP a, 
         SEXP beta, 
         SEXP d, 
         SEXP varB, 
         SEXP varE, 
         SEXP probIn, 
         SEXP RSS,
         SEXP c)
{
	double Cjj, lhs, rhs,  sol, offset, vare, z, u, logOdds, logPriorOdds,logP, probin, old_beta, cs;
    double *pXY, *pXX, *pa, *pd, *pvarB, *pbeta, *pRSS;
    
    int *pidColumns;
 
    int inc=1;
    int p;
    int q;
    int j;
    int j_global;
 
    SEXP list;
 
    GetRNGstate();
 
    p=INTEGER_VALUE(nCol);
    q=INTEGER_VALUE(length);
 
    vare=NUMERIC_VALUE(varE);
    probin=NUMERIC_VALUE(probIn);
    cs=NUMERIC_VALUE(c);
    
	PROTECT(XX=AS_NUMERIC(XX));
    pXX=NUMERIC_POINTER(XX);
     
    PROTECT(XY=AS_NUMERIC(XY));
    pXY=NUMERIC_POINTER(XY);
 
    PROTECT(a=AS_NUMERIC(a));
    pa=NUMERIC_POINTER(a);
 
    PROTECT(d=AS_NUMERIC(d));
    pd=NUMERIC_POINTER(d);

    PROTECT(beta=AS_NUMERIC(beta)); 
    pbeta=NUMERIC_POINTER(beta);

    PROTECT(varB=AS_NUMERIC(varB));
    pvarB=NUMERIC_POINTER(varB);

    PROTECT(RSS=AS_NUMERIC(RSS));
    pRSS=NUMERIC_POINTER(RSS);
    
    PROTECT(idColumns=AS_INTEGER(idColumns));
    pidColumns=INTEGER_POINTER(idColumns);
    
    logPriorOdds=log(probin/(1-probin));    

	/*
	This loop is for the idsColumns
	*/   
	 
     for(j=0; j<q;j++)
     {
     
     	  /*
    	  Substract one because in C we begin to count in 0
    	  */
    	  
    	  j_global=pidColumns[j]-1;
    	  
          Cjj=pXX[j_global*(p+1)];					
          rhs=pXY[j_global];
          
          old_beta=pbeta[j_global];
          
          offset=F77_NAME(ddot)(&p, pXX+j_global*p, &inc, pbeta, &inc);
          
          offset-=Cjj*pbeta[j_global];
          
          // random variables
          z=norm_rand();
          u=unif_rand();
          logOdds=log(u/(1-u));
              	
          // include the effect in the model?
          logP= 1.0/(2*pvarB[j])*pow(old_beta,2.0)*(1.0/pow(cs,2.0)-1.0)+log(cs) + logPriorOdds;
          
          if(logP>logOdds){
            pd[j]=1;
            pa[j]=1;       	
          }else{    
         	pd[j]=0;
         	pa[j]=cs;
          }
         
          lhs=Cjj+vare/(pvarB[j]*pow(pa[j],2.0));
         
          sol=(rhs-offset)/lhs;
          pbeta[j_global]=sol+sqrt(vare/lhs)*z;
         
          pRSS[0]+=(pow(pbeta[j_global],2) - pow(old_beta,2))*Cjj  -2*(pbeta[j_global]-old_beta)*(rhs-offset);
         
    }
 
      // Creating a list with 1 vector elements:
      PROTECT(list = Rf_allocVector(VECSXP, 4));
      SET_VECTOR_ELT(list, 0, a);
      SET_VECTOR_ELT(list, 1, d);      
      SET_VECTOR_ELT(list, 2, beta);
      SET_VECTOR_ELT(list, 3, RSS);
       
      PutRNGstate();
 
      UNPROTECT(9);
 
      return(list);

}
