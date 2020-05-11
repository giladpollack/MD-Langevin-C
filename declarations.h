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
#include <time.h>
#include <omp.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include "RandGen.h"


// ********************** constants  **********************

#define PI         (3.141592653589793)
#define hbar       (1.054571726e-34)    // Joule*Sec (MKS).
#define e          (-1.6021765653e-19)  // Coulomb (MKS).
#define EPS        (1.0e-12)
#define TINY       (1.0e-16)
#define kB_AU      (3.167e-6)           // Boltzmann factor in a.u. (in kB*T T is give in Kelvin units).
#define kB         (1.3806e-23)         // Boltzmann factor in a.u. (in kB*T T is give in Kelvin units).
#define eV2au      (0.03675)            // 0.03675 a.u.      = 1 ev.
#define au2sec     (2.418884326502e-17) // 1 a.u. of time    = 2.4188843265e-17 sec.
#define au2ampere  (0.00662361763)      // 1 a.u. of current = 0.00662361763 ampere
#define au2Joule   (4.3597482e-18)

#define sqr(x)     ((x)*(x))

// *********************** structures  **********************

typedef struct
{
  double x, y;
  
} Point;


typedef struct
{
  double Theta, R;
  
} Polar_Point;
// ********************** Functions  **********************

void InitSim(int NParticles, double Rad, double Err, long &idum);
void HydroRot(Point *PosC, Polar_Point *PosP, double *D, double *C, double *R, double *F, double aa, double dt, double K, double Kt, double Eps, double *P, int NParticles,double T, long &idum);
void MobRot(Point *PosC, double *D, double *C, double *P, double aa, int NParticles);
void VorForce(double *F, Point *PosC, Polar_Point *PosP, double aa, double K, double Kt, double Eps, double Rad, int NParticles);
int cholesky(int NParticles, double** D, double** C);
int Index(int ii,int jj, int NParticles);
void Randomize(long* idum);
void RandomizePositions(int NumOfParticles, double* WallPositionsX, double* WallPositionsY, double R, Point* Positions, RandGen rng);
void CopyPositions(Point* TargetArray, Point* SourceArray, int NumOfParticles);
bool doesDirExist(char* Path);
void GetPositionsString(Point* Positions, int NumOfParticles, char* OutString);
void GetSingleAxisSavedSteps(Point* ParticlePositions, int NumOfParticles, char axis, char* OutString);
void SetToZero(double* Array, int ArrayLen);
void InfoChamber(int N, double Dt, double SampleRate,
                 double R,double T, double Eta,
                 double Lx, double Ly,double WallShrink,
                 int NumOfParticles, char* SaveFoldername, bool DisplayLive, RandGen rng);
void RotnePrager(   Point* ParticlePositions,
                    int NumOfParticles,
                    double R,
                    double D,
                    double** DMat,
                    double** AMat,
                    double* Dx,
                    double* Dy,
                    double* Ax,
                    double* Ay);                 
#endif
