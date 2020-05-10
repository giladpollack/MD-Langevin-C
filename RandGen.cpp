#include "declarations.h"
#include "RandGen.h"
#include <random>

#define IM1 2147483563
#define IM2 2147483399
#define AM  (1.0 / (double)IM1)
#define IMM1 (IM1 - 1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1 + IMM1 / NTAB)
#define REPS 1.2e-7
#define RNMX (1.0 - REPS)

/****************************************************************************/

// Members
long RandGen::idum;

// Functions

void RandGen::Init(long* idum)
{
  RandGen::ran_nrc(idum);
  RandGen::idum = *idum;
}

/**
 * Returns a uniformly distributed random number between 0 and 1
 * Based on Numerical Recipes' ran2() algorithm page 282
 *
 * @idum The seed of the random generator
 * @return The random number
 */
float RandGen::ran_nrc(long *idum)
{
  int    j;
  long   k;
  static long idum2 = 123456789;
  static long iy = 0;
  static long iv[NTAB];
  double temp;

  if (*idum <= 0){  /* Initialize */
    *idum += time(0); // Randomizing the seed a bit more
    if (-(*idum) < 1) *idum = 1;    /* Be sure to prevent idum = 0 */
    else *idum = -(*idum);
    idum2 = (*idum);

    for (j = NTAB+7; j >= 0; j--){/* Load the shuffle table after 8 warm-ups */
      k = (*idum) / IQ1;
      *idum = IA1 * (*idum - k * IQ1) - k * IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy = iv[0];
  }

  k = (*idum) / IQ1;      /* Start Here when not initializing */
  *idum = IA1 * (*idum - k * IQ1) - k * IR1; /* Compute idum = IA1*idum % IM1*/
  if (*idum < 0) *idum += IM1;               /* without overflow by Scharge  */

  k = idum2 / IQ2;                           /* method.                      */
  idum2 = IA2 * (idum2 - k * IQ2) - k * IR2; /* Compute idum2=IA2*idum % IM2 */
  if (idum2 < 0) idum2 += IM2;               /* likewise.                    */

  j = iy / NDIV;          /* Will be in the range 0..NTAB-1 */
  iy = iv[j] - idum2;     /* Here idum is shuffled, idum and idum2 are */
  iv[j] = *idum;          /* combined to generate output               */

  if (iy < 1) iy += IMM1;
  if ((temp = AM * iy) > RNMX) return RNMX; /* Because users don't expect */
  else return temp;                         /* endpoint value             */
}

/****************************************************************************/

/**
 * Moves the random seed some steps ahead according to the current time
 */
void RandGen::Randomize()
{
  long* idum = &(RandGen::idum);
  double tt = (double)time(0);
  long   rn = (long)tt % 1000;
  int    i;

  for (i = 0; i < rn; i++) ran_nrc(idum);
  return;
}

/****************************************************************************/
/**
 * Returns a normally distributed random number with mean 0 and sigma 1
 * Based on numerical recipes page 289
 *
 * @return The normal random number
 */
double RandGen::Randn()
{
  long* idum = &(RandGen::idum);
  double x, y, s=1.0;

  while(s >= 1){
    // srand((unsigned)time(NULL));
    x = 2.0*(ran_nrc(idum)) - 1;
    y = 2.0*(ran_nrc(idum)) - 1;
    s = ((x*x) + (y*y));
  }
  double fac = sqrt(-2.0*log(s)/s);

  return(x*fac);
}

double RandGen::Randu(double LBound, double UBound)
{
  double Len = abs(UBound - LBound);
  float x = ran_nrc(&(RandGen::idum));
  double res = LBound + (Len*ran_nrc(&(RandGen::idum)));
  return res;
}