#include "declarations.h"

/* this function calculate the rotne prager  tensor for n particles.
x and y are vector of x, y positions of n particles.
x and y must have the same length.
we create a matix of 2*length(x) X 2*length(x)
eta is the viscosity of the fluid.
a is the particle's radius.*/

void MobRot(Cartesian_Point *PosC, double *D, double *C, double *P, double aa, int NParticles)

{
  int ii, jj, result,n,m,mm;
  double  dx, dy, r, g, h;

  g  = (3.0*aa)/4.0;
  h  = pow(aa,3.0)/2.0;
  
    
  for(ii=0 ; ii < 2*NParticles ; ii++){
    for(jj=0 ; jj < 2*NParticles ; jj++){ 
      n=floor(ii/2);
      m=floor(jj/2);
      dx=PosC[n].x-PosC[m].x;
      dy=PosC[n].y-PosC[m].y;
      r=sqrt(dx*dx+dy*dy); 
      // if (r<2){  
      //printf("%f %f %f\n", dx,dy,r);
      //}
      if(ii == jj)   
	D[Index(ii,ii,NParticles)] = 1;
      
        else if (jj == (ii+1) && (jj % 2) == 1)
	D[Index(ii,jj,NParticles)] = 0;
      
      else if (ii == (jj+1) && (ii % 2) == 1)
	D[Index(ii,jj,NParticles)] = 0;
      
      else if ((ii % 2) != (jj % 2)){
	
      	D[Index(ii,jj,NParticles)] = (g/r)*((dx*dy)/(r*r))    - (h/pow(r,3))*(3*dx*dy)/(r*r);
      }
      else if ((ii % 2) == (jj % 2) && (ii % 2) == 1){
	
        D[Index(ii,jj,NParticles)] = (g/r)*(1.+(dy*dy)/(r*r)) + (h/(pow(r,3)))*(1-(3.*dy*dy)/(r*r));
     }
      else if ((ii % 2) == (jj % 2) && (ii % 2) == 0){

	D[Index(ii,jj,NParticles)] = (g/r)*(1.+(dx*dx)/(r*r)) + (h/(pow(r,3)))*(1-(3.*dx*dx)/(r*r));
	}
      }
    }
   result = cholesky(NParticles, D, C, P); 
    /*  for (mm=0 ;mm<4 ;mm++){
    printf("%f %f %f %f\n", D[Index(mm,0,2)],D[Index(mm,1,2)],D[Index(mm,2,2)],D[Index(mm,3,2)]);
  }
for (mm=0;mm<4;mm++){
    printf("%f %f %f %f\n", C[Index(mm,0,2)],C[Index(mm,1,2)],C[Index(mm,2,2)],C[Index(mm,3,2)]);
    }*/

  if(result==0){
    printf("cholesky failed\n");
    exit(0);
  }	    
}

