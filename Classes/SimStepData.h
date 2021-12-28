#ifndef SimStepData_h
#define SimStepData_h
#include "../declarations.h"
#include "Vector.h"
class SimStepData
{
    public:
    // General
    int StepNum;
    Point* ParticlePositions;
    Vector<Point> ParticlePosVec;
    double* ParticleOrientations;
    Vector<double> ParticleOrientVec;
    Vector<double> FxVec;
    Vector<double> FyVec;
    double* Fx;
    double* Fy;

    // Trap related
    Point* TrapPositions;
    double TrapsA;
    double TrapsS;

    // Wall Related
    double WallPositionsX[2];
    double WallPositionsY[2];

    // Ctor
    SimStepData() = delete;
    SimStepData(int NumOfParticles, bool UseWalls = false, bool UseTraps = false, int NumOfTraps = 0);
};

#endif