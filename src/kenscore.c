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
