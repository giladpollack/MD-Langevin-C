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

int cholesky(int NParticles, double *D, double *C, double *P){
  
  int ii, jj, kk;
  double s;
      
  //copy matrix
   
        for (ii = 0; ii < 2*NParticles; ii++) {
	    for (jj = 0; jj < 2*NParticles; jj++){ 
	      C[Index(ii,jj,NParticles)] = D[Index(ii,jj,NParticles)];
	      	    }
	  }

	//do decomposition
	for (ii = 0; ii<2*NParticles; ii++){
	  for (jj = 0; jj<ii+1; jj++){
	    s=0;
	    for (kk = 0; kk<jj; kk++){
	      s+=C[Index(ii,kk,NParticles)]*C[Index(jj,kk,NParticles)];
	    }
	    if (ii==jj){
	      if (C[Index(ii,ii,NParticles)]-s<=0){
		printf("matrix is not positive definite\n");
		return 0;
	      }
	      C[Index(ii,jj,NParticles)]=sqrt(C[Index(ii,ii,NParticles)]-s);
	    }
	    else {
	      C[Index(ii,jj,NParticles)]=(1.0/C[Index(jj,jj,NParticles)])*(C[Index(ii,jj,NParticles)]-s);
	    }
	  }
	}
	 
	//make lower triangle	
	for (ii = 0; ii < 2*NParticles; ii++) {
            for (jj = ii + 1; jj < 2*NParticles; jj++) {
	      C[Index(ii,jj,NParticles)]  = 0;
	     }
	 }

      return 1;
}




