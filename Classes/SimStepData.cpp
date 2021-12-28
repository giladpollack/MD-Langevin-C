#include "SimStepData.h"

SimStepData::SimStepData(int NumOfParticles, bool UseWalls, bool UseTraps, int NumOfTraps)
:ParticlePosVec(NumOfParticles),FxVec(NumOfParticles),FyVec(NumOfParticles),ParticleOrientVec(NumOfParticles)
{
    this->ParticlePositions = ParticlePosVec.ptr;
    this->ParticleOrientations = ParticleOrientVec.ptr;
    this->Fx = FxVec.ptr;
    this->Fy = FyVec.ptr;
}

