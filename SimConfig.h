#ifndef SimConfig_h
#define SimConfig_h
#include "declarations.h"
#include "RandGen.h"


class SimConfig
{
    // Access Modifier
    public:

    // Data members:
    // Basic settings
    int N;
    int NumOfParticles;
    double R;
    double Eta;
    double T;
    double Dt;
    double SampleRate;
    Point* InitPositions;

    // Display Settings
    bool DisplayLive;
    double* Xlimits;
    double* Ylimits;

    // Save Settings
    char* SaveFoldername;
    int SavePeriod;

    // Settings for forces
    bool UseHydro;
    bool UseWalls;
    char* WallRepulsionType;
    double WallHarmonicK;
    double WallGaussianA;
    double WallGaussianS;
    double* WallPositionsX;
    double* WallPositionsY;

    bool UseParticleRepulsion;
    double WCAEpsilon;

    bool UseTraps;
    Point* InitTrapPositions;
    double TrapsA;
    double TrapsS;

    // Random Generation
    RandGen* rng;

    // Methods
    void ToString(char* CfgStirng);

    // // Default Constructor
    SimConfig()
    {
        this->NumOfParticles = 0;
        this->R = 0;
        this->InitPositions = nullptr;
        this->Eta = 0;
        this->T = 0;
        this->N = 0;
        this->Dt = 0;
        this->SampleRate = 0;
        this->SaveFoldername = nullptr;
        this->DisplayLive = false;
        this->Xlimits = nullptr;
        this->Ylimits = nullptr;
        this->SavePeriod = 0;
        this->UseHydro = false;
        this->UseWalls = false;
        this->UseTraps = false;
        this->WallRepulsionType = nullptr;
        this->WallGaussianA = 0;
        this->WallGaussianS = 0;
        this->WallHarmonicK = 0;
        this->WallPositionsX = nullptr;
        this->WallPositionsY = nullptr;
        this->UseParticleRepulsion = false;
        this->WCAEpsilon = 0;
        this->InitTrapPositions = nullptr;
        this->TrapsA = 0;
        this->TrapsS = 0;
    }
    // SimConfig() = delete;
    // SimConfig(int NumOfParticles, double** R, double** InitPositions, double Eta,
    //           double T, double N, double Dt, double sampleRate, double SavePeriod, char** SaveFoldername)
    // {
    //     this->NumOfParticles = NumOfParticles;
    //     this->R = *R;
    //     this->InitPositions = *InitPositions;
    //     this->Eta = Eta;
    //     this->T = T;
    //     this->N = N;
    //     this->Dt = Dt;
    //     this->SampleRate = sampleRate;
    //     this->SaveFoldername = *SaveFoldername;
    //     this->DisplayLive = false;
    //     this->Xlimits = nullptr;
    //     this->Ylimits = nullptr;
    //     this->SavePeriod = SavePeriod;
    //     this->UseHydro = false;
    //     this->UseWalls = false;
    //     this->UseTraps = false;
    //     this->WallRepulsionType = nullptr;
    //     this->WallGaussianA = 0;
    //     this->WallGaussianS = 0;
    //     this->WallHarmonicK = 0;
    //     this->WallPositionsX = nullptr;
    //     this->WallPositionsY = nullptr;
    //     this->UseParticleRepulsion = false;
    //     this->WCAEpsilon = 0;
    //     this->InitTrapPositions = nullptr;
    //     this->TrapsA = 0;
    //     this->TrapsS = 0;
    // }
};
#endif