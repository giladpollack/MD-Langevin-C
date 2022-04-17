#include <iostream>
#include "SimConfig.h"
#include "../declarations.h"
#include "../Utility.h"

void SimConfig::ToString(char* CfgString)
{
    CfgString[0] = '\0';
    sprintf(CfgString, 
            "N = %d\n"
            "Dt = %f\n"
            "R = %f\n"
            "Eta = %f\n"
            "T=%f\n"
            "SampleRate = %f\n"
            "SavePeriod = %d\n"
            "NumOfParticles = %d\n"
            "UseParticleRepulsion = %d\n"
            "UseHydro = %d\n"
            "UseWalls = %d\n"
            "WallGaussianA = %d\n"
            "WallGaussianS = %d\n"
            "IsActive = %d\n"
            "ActiveV = %f\n"
            "ActiveChirality = %f\n"
            "FluctForce = %f\n"
            "FluctSwitchFreq = %f\n",
            this->N,
            this->Dt,
            this->R,
            this->Eta,
            this->T,
            this->SampleRate,
            this->SavePeriod,
            this->NumOfParticles,
            this->UseParticleRepulsion,
            this->UseHydro,
            this->UseWalls,
            this->WallGaussianA,
            this->WallGaussianS,
            this->IsActive,
            this->ActiveV,
            this->ActiveChirality,
            this->FluctForce,
            this->FluctSwitchFreq);
    strcat(CfgString, "SaveFoldername = ");        
    strcat(CfgString, this->SaveFoldername);
    strcat(CfgString, "\n");
    strcat(CfgString, "Positions:\n");
    char PositionsString[500];
    GetPositionsString(this->InitPositions,this->NumOfParticles,PositionsString);
    strcat(CfgString, PositionsString);
}