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
const char ORIENT_FILE[] = "orientations.csv";
const char FLUCT_DIR_FILE[] = "FluctDirections.csv";
const char CFG_FILE[] = "Cfg.txt";
const char ADDED_DATA_FILE[] = "AddedData.csv";

// Declarations
void SaveDataToFiles(char* PosFileX, char* PosFileY, char* AddedDataFile, Point** ParticlePositions, int SampleInd, SimConfig& Cfg, AdditionalData& addedData);
void SaveDataToFiles(char* PosFileX, char* PosFileY, char* OrientationsFile, char* AddedDataFile,  Point** ParticlePositions, double** ParticleOrientations, int SampleInd, SimConfig& Cfg, AdditionalData& AddedData);
void SampleAddedData(SimConfig& Cfg, SimStepData& CurrStepData, AdditionalData& AddedData, int SampleInd);

MDSim::MDSim(SimConfig Cfg, RandGen* rng)
{
    this->Cfg = Cfg;
    this->rng = rng;
};

bool MDSim::CheckFeedbackFunc(SimConfig& Cfg, SimStepData& CurrStepData, AdditionalData& AddedData)
{
    if (!Cfg.IsAthermal)
    {
        return false;
    }
    else
    {
        // TODO: Optimize by computing this once outside of the function
        int SwitchPeriod = 1 / (Cfg.Dt * Cfg.FluctSwitchFreq);
        if (CurrStepData.StepNum % SwitchPeriod == 0)
        {
            return true;
        }
        
    }
    
    return false;
}

void MDSim::FeedbackFunc(SimConfig& Cfg, SimStepData& CurrStepData, AdditionalData& AddedData)
{
    if (Cfg.IsAthermal)
    {
        for (int CurrParticle = 0; CurrParticle < Cfg.NumOfParticles; CurrParticle++)
        {
            CurrStepData.FluctDirections[CurrParticle] = (*this->rng).Randu(0, 2*PI);
        }
    }
    
}

void MDSim::PrintFunc(SimConfig& Cfg, SimStepData& CurrStepData, AdditionalData& AddedData)
{
    // std::cout << Cfg.SaveFoldername << " - " << ((CurrStepData.StepNum + 1) / (Cfg.N / 100)) << "%" << std::endl;
    if (Cfg.DisplayLive)
    {
        char PositionsString[500];
        GetPositionsString(CurrStepData.ParticlePositions,Cfg.NumOfParticles,PositionsString);
        std::cout << PositionsString;
    }
    
}

void MDSim::ForcesFunc(SimConfig& Cfg, SimStepData& CurrStepData, AdditionalData& AddedData, int SampleInd)
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
        else if (strcmp(Cfg.WallRepulsionType, "Gaussian") == 0)
        {
            GetGaussianWallForces(  CurrStepData.ParticlePositions,
                                    Cfg.NumOfParticles,
                                    Cfg.WallGaussianA,
                                    pow(Cfg.WallGaussianS,2),
                                    CurrStepData.WallPositionsX,
                                    CurrStepData.WallPositionsY,
                                    CurrStepData.Fx,
                                    CurrStepData.Fy,
                                    AddedData, SampleInd);
        }        
        else
        {
            throw;
        }
    }

    if (Cfg.IsAthermal)
    {
        GetAthermalFluctForces(CurrStepData.ParticlePositions,
                                    Cfg.NumOfParticles,
                                    Cfg.FluctForce,
                                    CurrStepData.FluctDirections,
                                    CurrStepData.Fx,
                                    CurrStepData.Fy,
                                    AddedData, SampleInd);
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
    double DRotate = kB*Cfg.T/(8*PI*Cfg.Eta*pow(Cfg.R,3));
    int d = 2;
    int SamplePeriod = 1 / (Cfg.Dt * Cfg.SampleRate);
    char PosFileX[100];
    char PosFileY[100];
    char OrientFile[100];
    char FluctDirFile[100];
    char AddedDataFile[100];
    char CfgFile[100];

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

    // Allocating memory for saving the particle orientations to use between saves
    Matrix<double> ParticleOrientMat(NumOfSavedSteps, Cfg.NumOfParticles);
    double** ParticleOrientations = ParticleOrientMat.Mat;

    // Allocating memory for saving the particle Fluctuation directions to use between saves
    Matrix<double> FluctDirectionsMat(NumOfSavedSteps, Cfg.NumOfParticles);
    double** FluctDirections = ParticleOrientMat.Mat;    

    // Allocating memory for the distance from wall added variable, if necessary
    Vector<double> ClosestPosVec(NumOfSavedSteps);
    AddedData.ClosestParticlePositions = ClosestPosVec.ptr;

    // Allocating memory for the force measurement variables, if necessary
    Vector<double> FRightVec(NumOfSavedSteps);
    Vector<double> FLeftVec(NumOfSavedSteps);
    Vector<double> FUpVec(NumOfSavedSteps);
    Vector<double> FDownVec(NumOfSavedSteps);
    AddedData.ForcesRight = FRightVec.ptr;
    AddedData.ForcesLeft = FLeftVec.ptr;
    AddedData.ForcesUp = FUpVec.ptr;
    AddedData.ForcesDown = FDownVec.ptr;

    // Allocating matrices to compute the hydrodynamic interactions if needed
    Matrix<double> DMat(Cfg.NumOfParticles*d, Cfg.NumOfParticles*d);
    Matrix<double> AMat(Cfg.NumOfParticles*d, Cfg.NumOfParticles*d);
    
    // Creating the path strings for the save files
    strcpy(PosFileX, Cfg.SaveFoldername);
    strcpy(PosFileY, Cfg.SaveFoldername);
    strcpy(OrientFile, Cfg.SaveFoldername);
    strcpy(FluctDirFile, Cfg.SaveFoldername);
    strcpy(CfgFile, Cfg.SaveFoldername);
    strcpy(AddedDataFile, Cfg.SaveFoldername);


    strcat(strcat(PosFileX, "/"), COORD_FILE_X);
    strcat(strcat(PosFileY, "/"), COORD_FILE_Y);
    strcat(strcat(OrientFile, "/"), ORIENT_FILE);
    strcat(strcat(FluctDirFile, "/"), FLUCT_DIR_FILE);
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
    if (Cfg.IsActive)
    {
        CopyDoubleArr(CurrStepData.ParticleOrientations, Cfg.InitOrientations, Cfg.NumOfParticles);
        CopyDoubleArr(ParticleOrientations[0], CurrStepData.ParticleOrientations, Cfg.NumOfParticles);
    }

    SampleAddedData(Cfg,CurrStepData,AddedData,0);

    if (Cfg.IsAthermal)
    {
        // Randomizing the random fluctuation directions
        FeedbackFunc(Cfg, CurrStepData, AddedData);
        CopyDoubleArr(FluctDirections[0], CurrStepData.FluctDirections, Cfg.NumOfParticles);
    }

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
            CopyDoubleArr(ParticleOrientations[SampleInd], CurrStepData.ParticleOrientations, Cfg.NumOfParticles);
            CopyDoubleArr(FluctDirections[SampleInd], CurrStepData.FluctDirections, Cfg.NumOfParticles);
            SampleAddedData(Cfg, CurrStepData, AddedData, SampleInd);
            AddedData.ForcesLeft[SampleInd] /= SamplePeriod;
            AddedData.ForcesRight[SampleInd] /= SamplePeriod;
            AddedData.ForcesUp[SampleInd] /= SamplePeriod;
            AddedData.ForcesDown[SampleInd] /= SamplePeriod;
            SampleInd++;
        }

        // Checking whether to save to file
        if ((i % Cfg.SavePeriod == 0) || (i == Cfg.N - 1))
        {
            PrintFunc(Cfg, CurrStepData, AddedData);
            if (Cfg.IsActive && Cfg.IsAthermal)
            {
                // FUTURE: Add save for active particles with athermal fluctuations
            }
            else if (!Cfg.IsActive && Cfg.IsAthermal)
            {
                SaveDataToFiles(PosFileX, PosFileY, FluctDirFile, AddedDataFile, ParticlePositions, FluctDirections, SampleInd, Cfg, AddedData);
            }
            else if (Cfg.IsActive && !Cfg.IsAthermal)
            {
                SaveDataToFiles(PosFileX, PosFileY, OrientFile, AddedDataFile, ParticlePositions, ParticleOrientations, SampleInd, Cfg, AddedData);
            }
            else
            {
                SaveDataToFiles(PosFileX, PosFileY, AddedDataFile, ParticlePositions, SampleInd, Cfg, AddedData);
            }
            
            
            SampleInd = 0;
        }

        // Checking whether to reseed the random number generator
        if (i % Cfg.ReseedPeriod == 0)
        {
            (*this->rng).Reseed((*this->rng).Randu(-1e7, 1e7));
        }
        

        // Forces Computation
        ForcesFunc(Cfg, CurrStepData, AddedData, SampleInd);

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
            CurrStepData.ParticlePositions[CurrParticle].x += Ax[CurrParticle]*sqrt(Cfg.Dt) +
                                                            (Dx[CurrParticle]/(kB*Cfg.T))*CurrStepData.Fx[CurrParticle]*Cfg.Dt;
            CurrStepData.ParticlePositions[CurrParticle].y += Ay[CurrParticle]*sqrt(Cfg.Dt) +
                                                            (Dy[CurrParticle]/(kB*Cfg.T))*CurrStepData.Fy[CurrParticle]*Cfg.Dt;            
            if (Cfg.IsActive)
            {
                // Adding the active motion component
                CurrStepData.ParticlePositions[CurrParticle].x += Cfg.ActiveV * cos(CurrStepData.ParticleOrientations[CurrParticle]) * Cfg.Dt;
                CurrStepData.ParticlePositions[CurrParticle].y += Cfg.ActiveV * sin(CurrStepData.ParticleOrientations[CurrParticle]) * Cfg.Dt;

                CurrStepData.ParticleOrientations[CurrParticle] += Cfg.ActiveChirality * Cfg.Dt + sqrt(2*DRotate*Cfg.Dt)*(*this->rng).Randn();
            }

        }
        
        // Checking whether to perform feedback
        if (CheckFeedbackFunc(Cfg,CurrStepData,AddedData))
        {
            FeedbackFunc(Cfg,CurrStepData,AddedData);
        }
    }
}

void SaveDataToFiles(char* PosFileX, char* PosFileY, char* OrientationsFile, char* AddedDataFile,  Point** ParticlePositions, double** ParticleOrientations, int SampleInd, SimConfig& Cfg, AdditionalData& AddedData)
{
            char StepString[500];

            // Saving the X positions
            FILE* OrientationsFileStream = fopen(OrientationsFile,"a");
            for (int SavedStepInd = 0; SavedStepInd < SampleInd - 1; SavedStepInd++)
            {
                GetRotationSavedSteps(ParticleOrientations[SavedStepInd], Cfg.NumOfParticles, StepString);
                fprintf(OrientationsFileStream, StepString); 
            }

            fclose(OrientationsFileStream);

            SaveDataToFiles(PosFileX, PosFileY, AddedDataFile, ParticlePositions, SampleInd, Cfg, AddedData);
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
                fprintf(DataFilestream, "%e,%e,%e,%e,%e\n",
                        AddedData.ForcesRight[SavedStepInd], AddedData.ForcesLeft[SavedStepInd],
                        AddedData.ForcesUp[SavedStepInd], AddedData.ForcesDown[SavedStepInd],
                        AddedData.ClosestParticlePositions[SavedStepInd]);
            }
            fclose(DataFilestream);
}

void SampleAddedData(SimConfig& Cfg, SimStepData& CurrStepData, AdditionalData& AddedData, int SampleInd)
{
    double MaxPosition = -1;
    for (int CurrParticleInd = 0; CurrParticleInd < Cfg.NumOfParticles; CurrParticleInd++)
    {
        if (CurrStepData.ParticlePositions[CurrParticleInd].x > MaxPosition)
        {
            MaxPosition = CurrStepData.ParticlePositions[CurrParticleInd].x;
        }
    }
    
    AddedData.ClosestParticlePositions[SampleInd] = CurrStepData.WallPositionsX[1] - MaxPosition;
    // for (int CurrParticle = 0; CurrParticle < Cfg.NumOfParticles; CurrParticle++)
    // {
    //     double CurrFx = CurrStepData.Fx[CurrParticle];
    //     double CurrFy = CurrStepData.Fy[CurrParticle];
    //     if (CurrFx > 0)
    //     {
    //         AddedData.ForcesRight[SampleInd] += abs(CurrFx);
    //     }
    //     else
    //     {
    //         AddedData.ForcesLeft[SampleInd] += abs(CurrFx);
    //     }
        
    //     if (CurrFy > 0)
    //     {
    //         AddedData.ForcesUp[SampleInd] += abs(CurrFy);
    //     }
    //     else
    //     {
    //         AddedData.ForcesDown[SampleInd] += abs(CurrFy);
    //     }
    // }
    
}