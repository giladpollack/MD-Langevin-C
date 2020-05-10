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
    if (Cfg.DisplayLive)
    {
        char PositionsString[500];
        GetPositionsString(CurrStepData.ParticlePositions,Cfg.NumOfParticles,PositionsString);
        std::cout << PositionsString;
    }
    
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

    // Allocating the diffusion coefficient vectors
    double* Dx = (double*) malloc(Cfg.NumOfParticles*sizeof(double));
    double* Dy = (double*) malloc(Cfg.NumOfParticles*sizeof(double));
    double* Ax = (double*) malloc(Cfg.NumOfParticles*sizeof(double));
    double* Ay = (double*) malloc(Cfg.NumOfParticles*sizeof(double));

    // Allocating the memory buffer of particle positions to use between saves
    int NumOfSavedSteps = Cfg.SavePeriod / SamplePeriod;
    Point** ParticlePositions = (Point**) malloc(NumOfSavedSteps*sizeof(Point*));
    for (int i = 0; i < NumOfSavedSteps; i++)
    {
        ParticlePositions[i] = (Point*) malloc(Cfg.NumOfParticles*sizeof(Point));
    }
    
    // Creating the path strings for the save files
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
    CopyPositions(ParticlePositions[0], CurrStepData.ParticlePositions, Cfg.NumOfParticles);

    // Checking whether to initialize wall-related variables
    if (Cfg.UseWalls)
    {
        // Copying the wall positions from SimConfig to SimStepData
        double CurrWallPositionsX[2];
        double CurrWallPositionsY[2];
        CurrWallPositionsX[0] = Cfg.WallPositionsX[0];
        CurrWallPositionsX[1] = Cfg.WallPositionsX[1];
        CurrWallPositionsY[0] = Cfg.WallPositionsY[0];
        CurrWallPositionsY[1] = Cfg.WallPositionsY[1];
        CurrStepData.WallPositionsX = CurrWallPositionsX;
        CurrStepData.WallPositionsY = CurrWallPositionsY;
    }

    // Printing the initial placements
    PrintFunc(Cfg,CurrStepData, AddedData);

    // Checking whether to use hydrodynamic interactions
    if (!Cfg.UseHydro)
    {
        // Setting the diffusion vectors to standard diffusion without hydrodynamic corrections
        for (int i = 0; i < Cfg.NumOfParticles; i++)
        {
            Dx[i] = D;
            Dy[i] = D;
            Ax[i] = sqrt(2*D);
            Ay[i] = sqrt(2*D);
        }
    }


    
    // Freeing the allocated variables
    for (int i = 0; i < NumOfSavedSteps; i++)
    {
        free(ParticlePositions[i]);
    }
    free(ParticlePositions);
    free(CurrStepData.ParticlePositions);
    free(ParticlePositions);
    free(Dx);
    free(Dy);
    free(Ax);
    free(Ay);
}