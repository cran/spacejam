#include "stdlib.h"
#include "math.h"
#include "R.h"

void grpsel( double *b0, double *bigX, double *ssx, double *lambda, double *res, int *Active, int *n, int *p, int *k, int *q){
  double cond, rl1[*n], rl2[*n], ang[2*(*k)], lambdaeff;
  int i, j;

  for(int a = 0; a < *q; a++){
    i = Active[a] - 1;
    j = Active[*q + a] - 1;
    lambdaeff = 0;

	// Remove group effects from residuals
	for(int r=0; r < *n; r++){

	  rl1[r] = res[r*(*p) + i];
	  rl2[r] = res[r*(*p) + j];

	  for(int c=0; c < *k; c++){
	    rl1[r] += b0[i*(*k)*(*p) + j +c*(*p)]*bigX[r*(*k)*(*p) + j + (*p)*c];
	    rl2[r] += b0[j*(*k)*(*p) + i +c*(*p)]*bigX[r*(*k)*(*p) + i + (*p)*c];
	  }
	}
      
	// Calculate angle
	for(int d=0; d < 2*(*k); d++){
	  ang[d] = 0;
	}
	for(int r=0; r < *n; r++){
	  for(int v =0; v < *k; v++){
	    ang[v] += rl1[r]*(bigX[r*(*k)*(*p) + j + (*p)*v])/(*n-1);
	    ang[v+ (*k)] += rl2[r]*(bigX[r*(*k)*(*p) + i + (*p)*v])/(*n-1);
	  }
	}
	// Threshold 
	cond= 0;
	for(int c=0; c < 2*(*k); c++){
	  cond += ang[c]*ang[c];
	}
	
	if(sqrt(cond) > lambda[0]){
	  for(int l=0; l < *k; l++){
	    b0[i*(*k)*(*p) + j +l*(*p)] = (1-lambda[0]/sqrt(cond))*ang[l];
	    b0[j*(*k)*(*p) + i +l*(*p)] = (1-lambda[0]/sqrt(cond))*ang[l + *k];
	    //Rprintf("%f \n",  ssx[1] );
	  }
      
	} else {
	  for(int l=0; l < *k; l++){
	    b0[i*(*k)*(*p) + j +l*(*p)] = 0;
	    b0[j*(*k)*(*p) + i +l*(*p)] = 0;
	  }
	}

	// Update Residuals
      	for(int r=0; r < *n; r++){

	  res[r*(*p) + i] =  rl1[r];
	  res[r*(*p) + j] =  rl2[r];

	  for(int c=0; c < *k; c++){
	    res[r*(*p) + i] -= b0[i*(*k)*(*p) + j +c*(*p)]*bigX[r*(*k)*(*p) + j + (*p)*c];
	    res[r*(*p) + j] -= b0[j*(*k)*(*p) + i +c*(*p)]*bigX[r*(*k)*(*p) + i + (*p)*c];
	  }
	}
  }
 return;
}

void grpselDAG( double *b0, double *bigX, double *ssx, double *lambda, double *res, int *Active, int *n, int *p, int *k, int *q){
  double cond, rl1[*n], ang[*k];
  int i, j;

  for(int a = 0; a < *q; a++){
    i = Active[a] - 1; //outcome (top row)
    j = Active[*q + a] - 1; //predictor (bottom row)

	// Remove group effects from residuals
	for(int r=0; r < *n; r++){
	  rl1[r] = res[r*(*p) + i];

	  for(int c=0; c < *k; c++){
	    rl1[r] += b0[i*(*k)*(*p) + j +c*(*p)]*bigX[r*(*k)*(*p) + j + (*p)*c];
	  }
	}
	// Calculate angle
	for(int d=0; d < *k; d++){
	  ang[d] = 0;
	}
	for(int r=0; r < *n; r++){
	  for(int v =0; v < *k; v++){
	    ang[v] += rl1[r]*(bigX[r*(*k)*(*p) + j + (*p)*v])/(*n-1);
	  }
	}
	// Threshold 
	cond= 0;
	for(int c=0; c < *k; c++){
	  cond += ang[c]*ang[c];
	}
	if(sqrt(cond) > *lambda){
	  for(int l=0; l < *k; l++){
	    b0[i*(*k)*(*p) + j +l*(*p)] = (1-*lambda/sqrt(cond))*ang[l];
	  }
	} else {
	  for(int l=0; l < *k; l++){
	    b0[i*(*k)*(*p) + j +l*(*p)] = 0;
	  }
	}

	// Update Residuals
      	for(int r=0; r < *n; r++){
	  res[r*(*p) + i] =  rl1[r];
	  for(int c=0; c < *k; c++){
	    res[r*(*p) + i] -= b0[i*(*k)*(*p) + j +c*(*p)]*bigX[r*(*k)*(*p) + j + (*p)*c];
	  }
	}
  }
 return;
}
