#include "declarations.h"

// This function calculates the force acting on a particle in a vortex

void VorForce(double *F, Cartesian_Point *PosC, Polar_Point *PosP, double aa, double K, double Kt, double Eps, double Rad, int NParticles)
{
  int ii,jj;
  double dx, dy,r2,rmax2;

  for(ii=0 ; ii < NParticles ; ii++){
    PosP[ii].Theta = atan2(PosC[ii].y,PosC[ii].x);
    PosP[ii].R = sqrt(sqr(PosC[ii].x) + sqr(PosC[ii].y));
  }
  
  // F rotation and trapping
  
  for(ii=0 ; ii < NParticles ; ii++){
    F[2*ii]   = PosP[ii].R * Kt * (-sin(PosP[ii].Theta)) - K * ( PosP[ii].R - Rad) * cos(PosP[ii].Theta);
    F[2*ii+1] = PosP[ii].R * Kt * ( cos(PosP[ii].Theta)) - K * ( PosP[ii].R - Rad) * sin(PosP[ii].Theta);    
    if(NParticles > 1){
      for(jj=0 ; jj < NParticles ; jj++){ 
	if (jj!=ii){
	dx=PosC[ii].x-PosC[jj].x;
	dy=PosC[ii].y-PosC[jj].y;
	r2=dx*dx+dy*dy;
	rmax2=pow(pow(2,1./6.)*(2*aa),2);
	if(r2 < rmax2){
	  //printf("Fx=%f, Fy=%f\n",Eps,F[2*ii+1]);
	  F[2*ii]     += 24.0*(Eps/(r2)) * (2*pow((2*aa/r2),6) - pow((2*aa/r2),3))*dx;
	  F[2*ii+1]   += 24.0*(Eps/(r2)) * (2*pow((2*aa/r2),6) - pow((2*aa/r2),3))*dy;
	  //printf("Fx=%f, Fy=%f\n",F[2*ii],F[2*ii+1]);
	  //exit(0);
	}
	}
	}
	}
  }
}
   
