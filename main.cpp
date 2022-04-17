#include <iostream>
#include <time.h>
#include <omp.h>
#include <assert.h>
#include <random>
#include <direct.h>
#include "declarations.h"
#include "Classes/RandGen.h"



const int MAX_PARALLEL_SIMS = 10;

void RunActiveSims( RandGen rng,
                    const int NumOfSims,
                    long StepsInSim,
                    double Dt,
                    int SampleRate,
                    int T,
                    double ParticleDiameter,
                    int NumOfParticles,
                    bool UseParticleInteractions,
                    double ChamberLen,
                    double WallShrinkPerc,
                    double WallShrink,
                    double SolutionViscosity);

void RunAthermalSims( RandGen rng,
                      const int NumOfSims,
                      long StepsInSim,
                      double Dt,
                      int SampleRate,
                      int T,
                      double ParticleDiameter,
                      int NumOfParticles,
                      bool UseParticleInteractions,
                      double ChamberLen,
                      double WallShrinkPerc,
                      double WallShrink,
                      double SolutionViscosity);                    

/*this function is simulating hydrodynamic interaction using the rotne
prager tensor and not the simple mobility function [X,Y]=full_hydro(x,
y, eta, a, dt, simtime, Kr, Kt, eps, rad, path) x,y are the starting
position nof the particles in the x, y axes.  eta is the viscosity, a
is the radius and dt is the time between collisions, simtime is the
number of steps for the simulation.  the output is x and y matrixs for
the positions */
int main(int argc, char *argv[])
{  
  RandGen rng(time(NULL));
  bool IsActive = false;

  for (int ShrinkInd = 0; ShrinkInd < 5; ShrinkInd++)
  {
    for (int CurrParticleNum = 8; CurrParticleNum < 9; CurrParticleNum++)
    {

      const int NumOfSims = 300;
      long StepsInSim = 1e7;

      double Dt = 1e-4; // sec
      int SampleRate = 20; // 1/sec
      int T = 300;
      double ParticleDiameter = 2e-6; // meter
      int NumOfParticles = CurrParticleNum;
      bool UseParticleInteractions;
      double ChamberLen = 16e-6; // meter
      double WallShrinkPerc = 10 * ShrinkInd;
      double WallShrink = (WallShrinkPerc / 100) * ChamberLen;

      double WaterViscosity = 0.001; // Pascal-second
      double SolutionViscosity = 2*WaterViscosity;



      if (IsActive)
      {
        RunActiveSims(rng, NumOfSims, StepsInSim, Dt, SampleRate, T, ParticleDiameter,
                      NumOfParticles, UseParticleInteractions, ChamberLen, WallShrinkPerc,
                      WallShrink, SolutionViscosity);
      }
      else
      {
        RunAthermalSims(rng, NumOfSims, StepsInSim, Dt, SampleRate, T, ParticleDiameter,
                        NumOfParticles, UseParticleInteractions, ChamberLen, WallShrinkPerc,
                        WallShrink, SolutionViscosity);
      }
      
        
      
    }
  }
  // Non-Parallel run
  // char SaveFoldername[20] = "Single system rand";
  // InfoChamber(1e8, 1e-5, 20, 1e-6, 300, 0.001, 11e-6, 11e-6, 0, 1, SaveFoldername, false, rng);
}

void RunActiveSims( RandGen rng,
                    const int NumOfSims,
                    long StepsInSim,
                    double Dt,
                    int SampleRate,
                    int T,
                    double ParticleDiameter,
                    int NumOfParticles,
                    bool UseParticleInteractions,
                    double ChamberLen,
                    double WallShrinkPerc,
                    double WallShrink,
                    double SolutionViscosity)
{
  char Foldernames[NumOfSims][50];
  char BaseFoldername[200];
  for (int currVInd = 1; currVInd < 2; currVInd++)
  {
    
    double ActiveV = 0.2e-6 + currVInd*0.2e-6; // meter/sec

    for (int InteractionsInd = 1; InteractionsInd < 2; InteractionsInd++)
    { 
      UseParticleInteractions = (InteractionsInd == 1);
      
      


      if (!UseParticleInteractions)
      {
        sprintf(BaseFoldername, "%d particles Gaussian shrunk %d active v=%0.3f no interactions", NumOfParticles, (int)WallShrinkPerc, ActiveV*1e6);
      }
      else
      {
        sprintf(BaseFoldername, "%d particles Gaussian shrunk %d active v=%0.3f", NumOfParticles, (int)WallShrinkPerc, ActiveV*1e6);
      }
      
      
      // char BaseFoldername[] = "5 Particles With Interactions Shrunk 6 micron";
      int NumOfSimsRun = 0;
      
      // Creating the base folder to save the simulations in
      mkdir(BaseFoldername);

      // Running the parallel simulations
      omp_set_num_threads(8);
      #pragma omp parallel for
      for (int i = 0; i < NumOfSims; i++)
      {
        sprintf(Foldernames[i],"%s/%d",BaseFoldername, i + 1);
        RandGen currGen(rng.Randu(-1e7, 1e7));
        InfoChamber(StepsInSim, Dt, SampleRate, ParticleDiameter / 2, T, SolutionViscosity,
                    ChamberLen, ChamberLen, WallShrink, NumOfParticles, Foldernames[i],
                    false, UseParticleInteractions, ActiveV,
                    currGen);
        NumOfSimsRun++;
        std::cout << "Run " << NumOfSimsRun << " out of " << NumOfSims << std::endl;
      }
    }
  }
}

void RunAthermalSims( RandGen rng,
                      const int NumOfSims,
                      long StepsInSim,
                      double Dt,
                      int SampleRate,
                      int T,
                      double ParticleDiameter,
                      int NumOfParticles,
                      bool UseParticleInteractions,
                      double ChamberLen,
                      double WallShrinkPerc,
                      double WallShrink,
                      double SolutionViscosity)
{
  char Foldernames[NumOfSims][50];
  char BaseFoldername[200];
  for (int currFluctForceInd = 0; currFluctForceInd < 5; currFluctForceInd++)
  {
    double FluctForceNum = 10 + currFluctForceInd*10;
    double FluctForce = kB*T*1e6*FluctForceNum; // Newton
    double FluctSwitchFreq = 20;

    for (int InteractionsInd = 1; InteractionsInd < 2; InteractionsInd++)
    { 
      UseParticleInteractions = (InteractionsInd == 1);
      
      


      if (!UseParticleInteractions)
      {
        sprintf(BaseFoldername, "%d particles Gaussian shrunk %d athermal force=%0.1f no interactions", NumOfParticles, (int)WallShrinkPerc,  FluctForceNum);
      }
      else
      {
        sprintf(BaseFoldername, "%d particles Gaussian shrunk %d athermal force=%0.1f", NumOfParticles, (int)WallShrinkPerc, FluctForceNum);
      }
      
      int NumOfSimsRun = 0;
      
      // Creating the base folder to save the simulations in
      mkdir(BaseFoldername);

      // Running the parallel simulations
      omp_set_num_threads(8);
      #pragma omp parallel for
      for (int i = 0; i < NumOfSims; i++)
      {
        sprintf(Foldernames[i],"%s/%d",BaseFoldername, i + 1);
        RandGen currGen(rng.Randu(-1e7, 1e7));
        InfoChamberAthermal(StepsInSim, Dt, SampleRate, ParticleDiameter / 2, T, SolutionViscosity,
                    ChamberLen, ChamberLen, WallShrink, NumOfParticles, Foldernames[i],
                    false, UseParticleInteractions, FluctForce, FluctSwitchFreq,
                    currGen);
        NumOfSimsRun++;
        std::cout << "Shrunk " << (int)WallShrinkPerc << " f=" << (int)FluctForceNum  << " Run " << NumOfSimsRun << " out of " << NumOfSims << std::endl;
      }
    }
  }
}