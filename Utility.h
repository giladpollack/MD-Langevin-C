#ifndef Utility_h
#define Utility_h
#include "declarations.h"

// Functions
void RandomizePositions(int NumOfParticles, double* WallPositionsX, double* WallPositionsY, double R, Point* Positions, RandGen rng);
void CopyPositions(Point* TargetArray, Point* SourceArray, int NumOfParticles);
void GetPositionsString(Point* Positions, int NumOfParticles, char* OutString);
void GetSingleAxisSavedSteps(Point* ParticlePositions, int NumOfParticles, char axis, char* OutString);
bool doesDirExist(char* Path);
void SetToZero(double* Array, int ArrayLen);

#endif