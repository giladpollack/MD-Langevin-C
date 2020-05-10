#include "declarations.h"

#include <stdio.h>
#include <math.h>
/****************************************************
-----------------------------------------------
        Cholesky decomposition.

        input    NParticles*2  size of matrix
        input    D  Symmetric positive def. matrix
        output   C  lower deomposed matrix
----------------------------------------------- */

int cholesky(int NParticles, double** D, double** C){
  
  int ii, jj, kk;
  double s;
      
  //copy matrix
   
        for (ii = 0; ii < 2*NParticles; ii++) {
	    for (jj = 0; jj < 2*NParticles; jj++){ 
	      C[ii][jj] = D[ii][jj];
	      	    }
	  }

	//do decomposition
	for (ii = 0; ii<2*NParticles; ii++){
	  for (jj = 0; jj<ii+1; jj++){
	    s=0;
	    for (kk = 0; kk<jj; kk++){
	      s+=C[ii][kk]*C[jj][kk];
	    }
	    if (ii==jj){
	      if (C[ii][ii]-s<=0){
		printf("matrix is not positive definite\n");
		throw;
		return 0;
	      }
	      C[ii][jj]=sqrt(C[ii][ii]-s);
	    }
	    else {
	      C[ii][jj]=(1.0/C[jj][jj])*(C[ii][jj]-s);
	    }
	  }
	}
	 
	//make lower triangle	
	for (ii = 0; ii < 2*NParticles; ii++) {
            for (jj = ii + 1; jj < 2*NParticles; jj++) {
	      C[ii][jj]  = 0;
	     }
	 }

      return 1;
}




