#include "SimConfig.h"
#include "declarations.h"
#include <iostream>

char* SimConfig::ToString(char* CfgString)
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
            "UseWalls = %d\n",
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
            this->UseWalls);
    strcat(CfgString, "SaveFoldername = ");        
    strcat(CfgString, this->SaveFoldername);
    strcat(CfgString, "\n");
    strcat(CfgString, "Positions:\n");
    for (int i = 0; i < this->NumOfParticles; i++)
    {
        char CurrPos[50];
        sprintf(CurrPos, "%e,%e\n", this->InitPositions[i].x, this->InitPositions[i].y);
        strcat(CfgString, CurrPos);
    }
}