#include "declarations.h"
#include "Utility.h"
#include <iostream>
#include <iomanip>
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

void PrintMat(double** Mat, int MatSize, double DivSize)
{
    // std::cout << std::fixed;
    // std::cout << std::setprecision(3);
    for (int CurrRow = 0; CurrRow < MatSize; CurrRow++)
    {
        for (int CurrCol = 0; CurrCol < MatSize; CurrCol++)
        {
            std::cout << Mat[CurrRow][CurrCol] / DivSize << '\t';
        }

        std::cout << std::endl;
    }
    std::cout << std::endl << std::endl << std::endl;
}

/**
 * Computes the hydrodynamic interactions using the Rotne Prager method (Following Nagar-Roichman 2018)
 *
 * Forthe particles, rather than using a tensor, we use a 2D matrix of the form:
 * 
 * ******************************************
 * * x1   y1  x2  y2  ...                   *
 * * y1                                     *
 * * x2                                     *
 * * y2                                     *
 * * .                                      *
 * * .                                      *
 * * .                                      *
 * ******************************************
 * 
 * @ParticlePositions The positions of the particles
 * @NumOfParticles The number of particles in the simulation
 * @R The radius of the particles (denoted as "a" in the article)
 * @D The diffusion coefficient of particles in the system without hydrodynamic interactions
 * @DMat, AMat The matrices to use, allocated once to save on memory allocations TODO: Move to a static class member
 * @Dx, Dy Output parameters for the corrected diffusion coefficients for the particles
 * @Ax, Ay Output parameters fot the square root of the corrected diffusion coefficients (using Cholesky decomposition)
 */
void RotnePrager(   Point* ParticlePositions,
                    int NumOfParticles,
                    double R,
                    double D,
                    double** DMat,
                    double** AMat,
                    double* Dx,
                    double* Dy,
                    double* Ax,
                    double* Ay,
                    double* Psi)
{
    // The code works for only 2 dimensions, thus:
    int d = 2;
    // Constants for the matrix values
    double C1 = 0.75*R;
    double C2 = 0.5*pow(R,3);

    // Setting the D matrix to be a unit matrix (multiplied by D) and initializing the A matrix
    for (int CurrRow = 0; CurrRow < d*NumOfParticles; CurrRow++)
    {
        for (int CurrCol = 0; CurrCol < d*NumOfParticles; CurrCol++)
        {
            DMat[CurrRow][CurrCol] = 0;
            AMat[CurrRow][CurrCol] = 0;

        }
        
        DMat[CurrRow][CurrRow] = D;
    }  

    // Computing the lower triangle of the symmetrical D matrix
    for (int CheckedParticle = 0; CheckedParticle < NumOfParticles; CheckedParticle++)
    {
        for (int PairedParticle = 0; PairedParticle < CheckedParticle; PairedParticle++)
        {
            // Computing the indexes in the matrix corresponding to x and y coordinates
            int CheckedXInd = CheckedParticle * d;
            int CheckedYInd = CheckedParticle * d + 1;
            int PairedXInd = PairedParticle * d;
            int PairedYInd = PairedParticle * d + 1;

            double XDiff = ParticlePositions[CheckedParticle].x - ParticlePositions[PairedParticle].x;
            double YDiff = ParticlePositions[CheckedParticle].y - ParticlePositions[PairedParticle].y;
            double r = sqrt(sqr(XDiff) + sqr(YDiff));

            // For the x-x matrix position, setting (C1/r)(1 + x^2/r^2) + (C2/r^3)(1 - 3*(x^2/r^2))
            DMat[CheckedXInd][PairedXInd] = D*((C1 / r)*(1 + sqr(XDiff) / sqr(r)) + (C2 / pow(r,3))*(1 - 3*(sqr(XDiff) / sqr(r))));

            // For the y-y matrix position, setting (C1/r)(1 + y^2/r^2) + (C2/r^3)(1 - 3*(x^2/r^2))
            DMat[CheckedYInd][PairedYInd] = D*((C1 / r)*(1 + sqr(YDiff) / sqr(r)) + (C2 / pow(r,3))*(1 - 3*(sqr(YDiff) / sqr(r))));

            // For the x-y and y-x positons, setting (C1/r)(xy/r^2) - 3(C2/r^3)(xy/r^2)
            DMat[CheckedXInd][PairedYInd] = D*((C1 / r)*(XDiff*YDiff/sqr(r)) - 3*(C2/pow(r,3))*(XDiff*YDiff/sqr(r)));
            DMat[CheckedYInd][PairedXInd] = DMat[CheckedXInd][PairedYInd];

            
        }
    }

    // PrintMat(DMat, d*NumOfParticles, D);
    // Using the matrix's symmetry to fill the upper triangle
    for (int i = 0; i < NumOfParticles*d; i++)
    {
        for (int j = i + 1; j < NumOfParticles*d; j++)
        {
            DMat[i][j] = DMat[j][i];
        }
    }

    // PrintMat(DMat, d*NumOfParticles, D);
    // Getting the cholesky lower triangle decomposition of the matrix
    cholesky(NumOfParticles, DMat, AMat);

    // PrintMat(AMat, d*NumOfParticles, sqrt(D));

    // Getting the DVector as defined by the D matrix times the unit vector of the proper length
    SetToZero(Dx, NumOfParticles);
    SetToZero(Dy, NumOfParticles);
    SetToZero(Ax, NumOfParticles);
    SetToZero(Ay, NumOfParticles);
    for (int CurrRow = 0; CurrRow < d*NumOfParticles; CurrRow++)
    {
        for (int CurrCol = 0; CurrCol < d*NumOfParticles; CurrCol++)
        {
            // Checking if the current value beongs in the x vector or the y vector
            int VectorInd = CurrRow / 2; // Automatically takes the floor of the result by c definitions
            if (CurrRow % 2 == 0)
            {
                Dx[VectorInd] += DMat[CurrRow][CurrCol];
                Ax[VectorInd] += sqrt(2) * AMat[CurrRow][CurrCol] * Psi[CurrCol];
            }
            else
            {
                Dy[VectorInd] += DMat[CurrRow][CurrCol];
                Ay[VectorInd] += sqrt(2) * AMat[CurrRow][CurrCol] * Psi[CurrCol];
            }
        }
        
    }
}