#include "declarations.h"
class SimStepData
{
    public:
    // General
    int StepNum;
    Point* ParticlePositions;
    double* Fx;
    double* Fy;

    // Trap related
    Point* TrapPositions;
    double TrapsA;
    double TrapsS;

    // Wall Related
    double* WallPositionsX;
    double* WallPositionsY;
};