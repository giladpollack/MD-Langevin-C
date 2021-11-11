#include "declarations.h"
void GetWCAParticleForces(Point* Positions, int NumOfParticles, double R, double Eps, double* Fx, double* Fy);
void GetWCAWallForces(  Point* ParticlePositions,
                        int NumOfParticles,
                        double R,
                        double Eps,
                        double* WallPositionsX,
                        double* WallPositionsY,
                        double* Fx,
                        double* Fy);
void GetHarmonicWallForces( Point* ParticlePositions,
                            int NumOfParticles,
                            double R,
                            double K,
                            double* WallPositionsX,
                            double* WallPositionsY,
                            double* Fx,
                            double* Fy);
                            
void GetGaussianWallForces( Point* ParticlePositions,
                            int NumOfParticles,
                            double A,
                            double sigmaSq,
                            double* WallPositionsX,
                            double* WallPositionsY,
                            double* Fx,
                            double* Fy,
                            AdditionalData& AddedData, int SampleInd);