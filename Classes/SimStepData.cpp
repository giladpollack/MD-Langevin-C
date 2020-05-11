#include "SimStepData.h"

SimStepData::SimStepData(int NumOfParticles, bool UseWalls, bool UseTraps, int NumOfTraps)
:ParticlePosVec(NumOfParticles),FxVec(NumOfParticles),FyVec(NumOfParticles)
{
    this->ParticlePositions = ParticlePosVec.ptr;
    this->Fx = FxVec.ptr;
    this->Fy = FyVec.ptr;
}