#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <Rinternals.h> // for SXPs


void kenscore(double *kenscore,
	      double *x, double *y, int *n)
{
    int i,j; double prod=0,kst=0;
    
    for ( i=0; i<*n; i++ ){
	for ( j=i+1; j<*n; j++ ){
	    prod = (x[i] - x[j]) * (y[i] - y[j]);
	    kst  += (prod > 0) - (prod < 0);
	}
    }
    
    *kenscore = kst;
}


static R_NativePrimitiveArgType kenscore_t[] = {
  REALSXP, REALSXP, REALSXP, INTSXP
};

static R_CMethodDef cMethods[] = {
  {"kenscore", (DL_FUNC) &kenscore, 4, kenscore_t},
   {NULL, NULL, 0}
};


#include <Rversion.h>
void R_init_meta(DllInfo *dll)
{
    R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
#if defined(R_VERSION) && R_VERSION >= R_Version(2, 16, 0)
    R_forceSymbols(dll, TRUE);
#endif
}
