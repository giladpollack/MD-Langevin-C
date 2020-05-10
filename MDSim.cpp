#include "declarations.h"
#include "MDSim.h"
#include "ForceFunctions.h"
#include <iostream>
#include <direct.h>

const char COORD_FILE_X[] = "pos_x.csv";
const char COORD_FILE_Y[] = "pos_y.csv";
const char CFG_FILE[] = "Cfg.txt";

MDSim::MDSim(SimConfig Cfg)
{
    this->Cfg = Cfg;
};

bool MDSim::CheckFeedbackFunc(SimConfig Cfg, SimStepData CurrStepData, AdditionalData AddedData)
{
    return false;
}

void MDSim::FeedbackFunc(SimConfig Cfg, SimStepData CurrStepData, AdditionalData AddedData)
{
    return;
}

void MDSim::PrintFunc(SimConfig Cfg, SimStepData CurrStepData, AdditionalData AddedData)
{
    std::cout << Cfg.SaveFoldername << (CurrStepData.StepNum / Cfg.N) * 100 << "%" << std::endl;
}

void MDSim::ForcesFunc(SimConfig Cfg, SimStepData CurrStepData, AdditionalData AddedData)
{
    if (Cfg.UseParticleRepulsion)
    {
        GetWCAParticleForces(CurrStepData.ParticlePositions,
                             Cfg.NumOfParticles,
                             Cfg.R, Cfg.WCAEpsilon,
                             CurrStepData.Fx,
                             CurrStepData.Fy);
    }
}

void MDSim::RunSim()
{
    this->RunSim(this->Cfg, this->addedData);
}
void MDSim::RunSim(SimConfig Cfg, AdditionalData AddedData)
{
    // Basic definitions
    double Gamma = 6*PI*Cfg.R*Cfg.Eta;
    double D = kB*Cfg.T/Gamma;
    int d = 2;
    int SamplePeriod = 1 / (Cfg.Dt * Cfg.SampleRate);
    char PosFileX[50];
    char PosFileY[50];
    char CfgFile[50];
    strcpy(PosFileX, Cfg.SaveFoldername);
    strcpy(PosFileY, Cfg.SaveFoldername);
    strcpy(CfgFile, Cfg.SaveFoldername);

    strcat(strcat(PosFileX, "/"), COORD_FILE_X);
    strcat(strcat(PosFileY, "/"), COORD_FILE_Y);
    strcat(strcat(CfgFile, "/"), CFG_FILE);

    // Checking if the save directory already exists
    if (!doesDirExist(Cfg.SaveFoldername))
    {
        mkdir(Cfg.SaveFoldername);
    }
    else
    {
        char response;
        std::cout << "The directory exists. Overwrite? (insert y to overwrite, anything else to abort)" <<  std::endl;
        std::cin >> response;
        if (response != 'y')
        {
            std::cout << "Aborting Sim";
            throw;
        }
        remove(PosFileX);
        remove(PosFileY);
        remove(CfgFile);
    }
    // Saving the configuration to a file

    char CfgString[10000];
    Cfg.ToString(CfgString);
    std::cout << CfgString << std::endl;
    FILE* CfgFilestream = fopen(CfgFile,"w");
    fprintf(CfgFilestream, CfgString);
    fclose(CfgFilestream);
    

    // Setting up the run variables
    SimStepData CurrStepData;
    CurrStepData.ParticlePositions = (Point*) malloc (Cfg.NumOfParticles*sizeof(Point));
    CopyPositions(CurrStepData.ParticlePositions, Cfg.InitPositions, Cfg.NumOfParticles);

    free(CurrStepData.ParticlePositions);
}