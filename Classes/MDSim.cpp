#include <iostream>
#include <direct.h>
#include "MDSim.h"
#include "RandGen.h"
#include "Matrix.h"
#include "Vector.h"
#include "../declarations.h"
#include "../ForceFunctions.h"
#include "../Utility.h"

const char COORD_FILE_X[] = "pos_x.csv";
const char COORD_FILE_Y[] = "pos_y.csv";
const char CFG_FILE[] = "Cfg.txt";
const char ADDED_DATA_FILE[] = "AddedData.csv";

// Declarations
void SaveDataToFiles(char* PosFileX, char* PosFileY, char* AddedDataFile, Point** ParticlePositions, int SampleInd, SimConfig& Cfg, AdditionalData& addedData);
void SampleAddedData(SimConfig& Cfg, SimStepData& CurrStepData, AdditionalData& AddedData, int SampleInd);

MDSim::MDSim(SimConfig Cfg, RandGen* rng)
{
    this->Cfg = Cfg;
    this->rng = rng;
};

bool MDSim::CheckFeedbackFunc(SimConfig& Cfg, SimStepData& CurrStepData, AdditionalData& AddedData)
{
    return false;
}

void MDSim::FeedbackFunc(SimConfig& Cfg, SimStepData& CurrStepData, AdditionalData& AddedData)
{
    return;
}

void MDSim::PrintFunc(SimConfig& Cfg, SimStepData& CurrStepData, AdditionalData& AddedData)
{
    //std::cout << Cfg.SaveFoldername << " - " << ((CurrStepData.StepNum + 1) / (Cfg.N / 100)) << "%" << std::endl;
    if (Cfg.DisplayLive)
    {
        char PositionsString[500];
        GetPositionsString(CurrStepData.ParticlePositions,Cfg.NumOfParticles,PositionsString);
        std::cout << PositionsString;
    }
    
}

void MDSim::ForcesFunc(SimConfig& Cfg, SimStepData& CurrStepData, AdditionalData& AddedData)
{
    // Initializing the forces at 0
    SetToZero(CurrStepData.Fx, Cfg.NumOfParticles);
    SetToZero(CurrStepData.Fy, Cfg.NumOfParticles);

    // Computing particle repulsion forces if necessary
    if (Cfg.UseParticleRepulsion)
    {
        GetWCAParticleForces(CurrStepData.ParticlePositions,
                             Cfg.NumOfParticles,
                             Cfg.R, Cfg.WCAEpsilon,
                             CurrStepData.Fx,
                             CurrStepData.Fy);
    }

    // Computing the wall repulsion forces if necessary
    if (Cfg.UseWalls)
    {
        if (strcmp(Cfg.WallRepulsionType, "WCA") == 0)
        {
            GetWCAWallForces(   CurrStepData.ParticlePositions,
                                Cfg.NumOfParticles,
                                Cfg.R,
                                Cfg.WCAEpsilon,
                                CurrStepData.WallPositionsX,
                                CurrStepData.WallPositionsY,
                                CurrStepData.Fx,
                                CurrStepData.Fy);
        }
        else if (strcmp(Cfg.WallRepulsionType, "Harmonic") == 0)
        {
            GetHarmonicWallForces(  CurrStepData.ParticlePositions,
                                    Cfg.NumOfParticles,
                                    Cfg.R,
                                    Cfg.WallHarmonicK,
                                    CurrStepData.WallPositionsX,
                                    CurrStepData.WallPositionsY,
                                    CurrStepData.Fx,
                                    CurrStepData.Fy);
        }
        else
        {
            throw;
        }
    }
    
}

void MDSim::RunSim()
{
    this->RunSim(this->Cfg, this->AddedData);
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
    char AddedDataFile[50];
    char CfgFile[50];

    // Allocating the diffusion coefficient vectors
    Vector<double> DxVec(Cfg.NumOfParticles);
    Vector<double> DyVec(Cfg.NumOfParticles);
    Vector<double> AxVec(Cfg.NumOfParticles);
    Vector<double> AyVec(Cfg.NumOfParticles);
    Vector<double> PsiVec(d*Cfg.NumOfParticles); // A vector of Gaussian distributed random variables
    double* Psi = PsiVec.ptr;
    double* Dx = DxVec.ptr;
    double* Dy = DyVec.ptr;
    double* Ax = AxVec.ptr;
    double* Ay = AyVec.ptr;

    // Allocating the memory buffer of particle positions to use between saves
    int NumOfSavedSteps = Cfg.SavePeriod / SamplePeriod + 2;
    Matrix<Point> ParticlePosMat(NumOfSavedSteps, Cfg.NumOfParticles);
    Point** ParticlePositions = ParticlePosMat.Mat;

    // Allocating memory for the distance from wall added variable, if necessary
    Vector<double> ClosestPosVec(NumOfSavedSteps);
    AddedData.ClosestParticlePositions = ClosestPosVec.ptr;

    // Allocating matrices to compute the hydrodynamic interactions if needed
    Matrix<double> DMat(Cfg.NumOfParticles*d, Cfg.NumOfParticles*d);
    Matrix<double> AMat(Cfg.NumOfParticles*d, Cfg.NumOfParticles*d);
    
    // Creating the path strings for the save files
    strcpy(PosFileX, Cfg.SaveFoldername);
    strcpy(PosFileY, Cfg.SaveFoldername);
    strcpy(CfgFile, Cfg.SaveFoldername);
    strcpy(AddedDataFile, Cfg.SaveFoldername);


    strcat(strcat(PosFileX, "/"), COORD_FILE_X);
    strcat(strcat(PosFileY, "/"), COORD_FILE_Y);
    strcat(strcat(CfgFile, "/"), CFG_FILE);
    strcat(strcat(AddedDataFile, "/"), ADDED_DATA_FILE);


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
        remove(AddedDataFile);
    }

    // Saving the configuration to a file
    char CfgString[10000];
    Cfg.ToString(CfgString);
    //std::cout << CfgString << std::endl;
    FILE* CfgFilestream = fopen(CfgFile,"w");
    fprintf(CfgFilestream, CfgString);
    fclose(CfgFilestream);

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

    // Setting up the run variables
    SimStepData CurrStepData(Cfg.NumOfParticles);
    CopyPositions(CurrStepData.ParticlePositions, Cfg.InitPositions, Cfg.NumOfParticles);
    CopyPositions(ParticlePositions[0], CurrStepData.ParticlePositions, Cfg.NumOfParticles);
    SampleAddedData(Cfg,CurrStepData,AddedData,0);

    // Checking whether to initialize wall-related variables
    if (Cfg.UseWalls)
    {
        // Copying the wall positions from SimConfig to SimStepData
        CurrStepData.WallPositionsX[0] = Cfg.WallPositionsX[0];
        CurrStepData.WallPositionsX[1] = Cfg.WallPositionsX[1];
        CurrStepData.WallPositionsY[0] = Cfg.WallPositionsY[0];
        CurrStepData.WallPositionsY[1] = Cfg.WallPositionsY[1];
    }

    // Printing the initial placements
    PrintFunc(Cfg,CurrStepData, AddedData);

    // Running the simulation
    int SampleInd = 1;
    for (int i = 1; i < Cfg.N; i++)
    {
        CurrStepData.StepNum = i;

        // Checking whether to sample and print
        if (i % SamplePeriod == 0)
        {
            CopyPositions(ParticlePositions[SampleInd], CurrStepData.ParticlePositions, Cfg.NumOfParticles);
            SampleAddedData(Cfg, CurrStepData, AddedData, SampleInd);
            SampleInd++;
        }

        // Checking whether to save to file
        if ((i % Cfg.SavePeriod == 0) || (i == Cfg.N - 1))
        {
            PrintFunc(Cfg, CurrStepData, AddedData);
            SaveDataToFiles(PosFileX, PosFileY, AddedDataFile, ParticlePositions, SampleInd, Cfg, AddedData);
            SampleInd = 0;
        }

        // Checking whether to reseed the random number generator
        if (i % Cfg.ReseedPeriod == 0)
        {
            (*this->rng).Reseed((*this->rng).Randu(-1e7, 1e7));
        }
        

        // Forces Computation
        ForcesFunc(Cfg, CurrStepData, AddedData);

        // Computing the vector of random gaussian variables for the brownian motion
        for (int CurrVar = 0; CurrVar < d*Cfg.NumOfParticles; CurrVar++)
        {
            Psi[CurrVar] = (*this->rng).Randn();
        }

        // Hydrodynamic interactions computation
        if (Cfg.UseHydro)
        {
            RotnePrager(    CurrStepData.ParticlePositions,
                            Cfg.NumOfParticles,
                            Cfg.R,
                            D,
                            DMat.Mat,
                            AMat.Mat,
                            Dx,
                            Dy,
                            Ax,
                            Ay,
                            Psi);
        }
        else
        {
            for (int CurrA = 0; CurrA < d*Cfg.NumOfParticles; CurrA++)
            {
                if (CurrA % 2 == 0)
                {
                    Ax[CurrA / 2] = Psi[CurrA] * sqrt(2*D);
                }
                else
                {
                    Ay[CurrA / 2] = Psi[CurrA] * sqrt(2*D);
                }
            }
            
        }
        
        
        // Running the step
        for (int CurrParticle = 0; CurrParticle < Cfg.NumOfParticles; CurrParticle++)
        {
            double First = Ax[CurrParticle]*sqrt(Cfg.Dt)*(*this->rng).Randn();
            double Second = (Dx[CurrParticle]/(kB*Cfg.T))*CurrStepData.Fx[CurrParticle]*Cfg.Dt;
            CurrStepData.ParticlePositions[CurrParticle].x += Ax[CurrParticle]*sqrt(Cfg.Dt) +
                                                              (Dx[CurrParticle]/(kB*Cfg.T))*CurrStepData.Fx[CurrParticle]*Cfg.Dt;
            CurrStepData.ParticlePositions[CurrParticle].y += Ay[CurrParticle]*sqrt(Cfg.Dt) +
                                                              (Dy[CurrParticle]/(kB*Cfg.T))*CurrStepData.Fy[CurrParticle]*Cfg.Dt;
        }
        
        // Checking whether to perform feedback
        if (CheckFeedbackFunc(Cfg,CurrStepData,AddedData))
        {
            FeedbackFunc(Cfg,CurrStepData,AddedData);
        }
    }
}

void SaveDataToFiles(char* PosFileX, char* PosFileY, char* AddedDataFile, Point** ParticlePositions, int SampleInd, SimConfig& Cfg, AdditionalData& AddedData)
{
            char StepString[500];

            // Saving the X positions
            FILE* XFilestream = fopen(PosFileX,"a");
            for (int SavedStepInd = 0; SavedStepInd < SampleInd - 1; SavedStepInd++)
            {
                GetSingleAxisSavedSteps(ParticlePositions[SavedStepInd], Cfg.NumOfParticles, 'x', StepString);
                fprintf(XFilestream, StepString); 
            }

            fclose(XFilestream);

            // Saving the Y Positions
            FILE* YFilestream = fopen(PosFileY,"a");
            for (int SavedStepInd = 0; SavedStepInd < SampleInd - 1; SavedStepInd++)
            {
                GetSingleAxisSavedSteps(ParticlePositions[SavedStepInd], Cfg.NumOfParticles, 'y', StepString);
                fprintf(YFilestream, StepString); 
            }
            
            fclose(YFilestream);

            // Saving the additional tracked data
            FILE* DataFilestream = fopen(AddedDataFile,"a");
            for (int SavedStepInd = 0; SavedStepInd < SampleInd - 1; SavedStepInd++)
            {
                fprintf(DataFilestream, "%e\n", AddedData.ClosestParticlePositions[SavedStepInd]); 
            }
            fclose(DataFilestream);
}

void SampleAddedData(SimConfig& Cfg, SimStepData& CurrStepData, AdditionalData& AddedData, int SampleInd)
{
    AddedData.ClosestParticlePositions[SampleInd] = CurrStepData.WallPositionsX[1] - CurrStepData.ParticlePositions[0].x - Cfg.R;
}