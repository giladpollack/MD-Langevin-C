#include "declarations.h"

/*****************************************************/

void HydroRot(Cartesian_Point *PosC, Polar_Point *PosP, double *D, double *C, double *R, double *F, double aa, double dt, double K, double Kt, double Eps, double Rad, double *P, int NParticles, double T, long &idum)
{
  int ii, jj;
  double S, Fs;

  MobRot(PosC,D,C,P,aa,NParticles);

  VorForce(F,PosC,PosP,aa,K,Kt,Eps,Rad,NParticles);

  //#pragma omp parallel
  //{
  //#pragma omp for
    for(jj=0 ; jj < 2*NParticles ; jj++) R[jj] = Randn(idum);
    // }

    //#pragma omp parallel
      // {
    //#pragma omp for private(jj,S,Fs)
    for(ii=0 ; ii < 2*NParticles ; ii++){
      S=Fs=0.0;

      if(!(ii%2)){

	for(jj=0 ; jj < 2*NParticles ; jj++){
	  S  += C[Index(ii,jj,NParticles)]*R[jj];

	  Fs += D[Index(ii,jj,NParticles)]*F[jj];
	  //printf("R %f, F %f\n",R[jj],F[jj]);
	}


	PosC[ii/2].x += dt*Fs + (sqrt(dt*T)*S);
	//  printf("Fs %f, S %f\n",Fs,S);
      }
      else{

	for(jj=0 ; jj < 2*NParticles ; jj++){
	  S  += C[Index(ii,jj,NParticles)]*R[jj];
	  Fs += D[Index(ii,jj,NParticles)]*F[jj];
	}

	PosC[ii/2].y += dt*Fs + (sqrt(dt*T)*S);
	//printf("Fs %f, S %f\n",Fs,S);
      }
    }
    //}
}
