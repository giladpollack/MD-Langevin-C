#ifndef declerations_h
#define declerations_h

// ********************** Libraries **********************

//#include <iostream.h>
//#include <iomanip.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <malloc.h>
#include <assert.h>
#include "complex.h"
#include <sys/time.h>
#include <omp.h>
#include <string.h>
#include <time.h>
#include <omp.h>

// ********************** constants  **********************

#define PIE        (3.141592653589793)
#define hbar       (1.054571726e-34)    // Joule*Sec (MKS).
#define e          (-1.6021765653e-19)  // Coulomb (MKS).
#define EPS        (1.0e-12)
#define TINY       (1.0e-16)
#define kB         (3.167e-6)           // Boltzmann factor in a.u. (in kB*T T is give in Kelvin units).
#define eV2au      (0.03675)            // 0.03675 a.u.      = 1 ev.
#define au2sec     (2.418884326502e-17) // 1 a.u. of time    = 2.4188843265e-17 sec.
#define au2ampere  (0.00662361763)      // 1 a.u. of current = 0.00662361763 ampere
#define au2Joule   (4.3597482e-18)

#define sqr(x)     ((x)*(x))

// *********************** structures  **********************

typedef struct
{
  double x, y;
  
} Cartesian_Point;


typedef struct
{
  double Theta, R;
  
} Polar_Point;

// ********************** Functions  **********************

void InitSim(int NParticles, double Rad, double Err, long &idum);
void HydroRot(Cartesian_Point *PosC, Polar_Point *PosP, double *D, double *C, double *R, double *F, double aa, double dt, double K, double Kt, double Eps, double Rad, double *P, int NParticles,double T, long &idum);
void MobRot(Cartesian_Point *PosC, double *D, double *C, double *P, double aa, int NParticles);
void VorForce(double *F, Cartesian_Point *PosC, Polar_Point *PosP, double aa, double K, double Kt, double Eps, double Rad, int NParticles);
int cholesky(int NParticles, double *D, double *C, double *P);
int Index(int ii,int jj, int NParticles);
void Randomize();
double Randn(long &idum);
double ran_nrc(long *idum);

#endif