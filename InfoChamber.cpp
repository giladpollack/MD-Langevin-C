#include "declarations.h"
#include "Classes/SimConfig.h"
#include "Classes/MDSim.h"
#include "Classes/RandGen.h"
#include "Utility.h"
#include "Classes/Vector.h"

void InfoChamber(int N, double Dt, double SampleRate,
                 double R,double T, double Eta,
                 double Lx, double Ly,double WallShrink,
                 int NumOfParticles, char* SaveFoldername,
                 bool DisplayLive, bool UseParticleInteractions,
                 double ActiveV, RandGen rng)
{
  double WallPositionsX[2];
  double WallPositionsY[2];
  SimConfig Cfg;
  Cfg.N = N;
  Cfg.Dt = Dt;
  Cfg.SampleRate = SampleRate;
  Cfg.NumOfParticles = NumOfParticles;
  Cfg.Eta = Eta;
  Cfg.R = R;
  Cfg.T = T;
  Cfg.SaveFoldername = SaveFoldername;
  Cfg.SavePeriod = Cfg.N / 10;
  Cfg.UseWalls = true;
  char WallRepulsionType[20] = "Gaussian";
  Cfg.WallRepulsionType = WallRepulsionType;
  Cfg.WallHarmonicK =  2e6*kB*Cfg.T /1e-6;
  Cfg.WallGaussianA = 50*kB*Cfg.T;
  Cfg.WallGaussianS = 2.5e-7;
  Cfg.UseHydro = false;
  Cfg.UseParticleRepulsion = UseParticleInteractions;
  Cfg.UseTraps = false;
  Cfg.WCAEpsilon = 0.2*kB*Cfg.T;
  Cfg.ReseedPeriod = 1e5;
  Cfg.IsActive = true;
  Cfg.ActiveV = ActiveV; // meter/sec
  Cfg.ActiveChirality = 0; // Rad/sec

  Cfg.IsAthermal = false;
  
  WallPositionsX[0] = -Lx / 2;
  WallPositionsX[1] = Lx /2 - WallShrink;
  WallPositionsY[0] = -Ly / 2;
  WallPositionsY[1] = Ly / 2;
  Cfg.WallPositionsX = WallPositionsX;
  Cfg.WallPositionsY = WallPositionsY;
  Vector<Point> InitPositions(Cfg.NumOfParticles);
  RandomizePositions(Cfg.NumOfParticles, WallPositionsX, WallPositionsY, R, InitPositions.ptr, rng);
  Cfg.InitPositions = InitPositions.ptr;

  // Randomizing the initial orientations
  Vector<double> InitOrientations(Cfg.NumOfParticles);
  for (int CurrParticleNum = 0; CurrParticleNum < Cfg.NumOfParticles; CurrParticleNum++)
  {
    InitOrientations.ptr[CurrParticleNum] = rng.Randu(0, 2*PI);
  }
  Cfg.InitOrientations = InitOrientations.ptr;

  // Initializing the simulation with the configuration
  MDSim currSim(Cfg, &rng);
  currSim.RunSim();
}

void InfoChamberAthermal(int N, double Dt, double SampleRate,
                 double R,double T, double Eta,
                 double Lx, double Ly,double WallShrink,
                 int NumOfParticles, char* SaveFoldername,
                 bool DisplayLive, bool UseParticleInteractions,
                 double FluctForce, double FluctFreq, RandGen rng)
{
  double WallPositionsX[2];
  double WallPositionsY[2];
  SimConfig Cfg;
  Cfg.N = N;
  Cfg.Dt = Dt;
  Cfg.SampleRate = SampleRate;
  Cfg.NumOfParticles = NumOfParticles;
  Cfg.Eta = Eta;
  Cfg.R = R;
  Cfg.T = T;
  Cfg.SaveFoldername = SaveFoldername;
  Cfg.SavePeriod = Cfg.N / 10;
  Cfg.UseWalls = true;
  char WallRepulsionType[20] = "Gaussian";
  Cfg.WallRepulsionType = WallRepulsionType;
  Cfg.WallHarmonicK =  2e6*kB*Cfg.T /1e-6;
  Cfg.WallGaussianA = 50*kB*Cfg.T;
  Cfg.WallGaussianS = 2.5e-7;
  Cfg.UseHydro = false;
  Cfg.UseParticleRepulsion = UseParticleInteractions;
  Cfg.UseTraps = false;
  Cfg.WCAEpsilon = 0.2*kB*Cfg.T;
  Cfg.ReseedPeriod = 1e5;
  Cfg.IsActive = false;
  Cfg.ActiveV = 0; // meter/sec
  Cfg.ActiveChirality = 0; // Rad/sec

  Cfg.IsAthermal = true;
  Cfg.FluctForce = FluctForce;
  Cfg.FluctSwitchFreq = FluctFreq;
  
  WallPositionsX[0] = -Lx / 2;
  WallPositionsX[1] = Lx /2 - WallShrink;
  WallPositionsY[0] = -Ly / 2;
  WallPositionsY[1] = Ly / 2;
  Cfg.WallPositionsX = WallPositionsX;
  Cfg.WallPositionsY = WallPositionsY;
  Vector<Point> InitPositions(Cfg.NumOfParticles);
  RandomizePositions(Cfg.NumOfParticles, WallPositionsX, WallPositionsY, R, InitPositions.ptr, rng);
  Cfg.InitPositions = InitPositions.ptr;

  // Initializing the simulation with the configuration
  MDSim currSim(Cfg, &rng);
  currSim.RunSim();
}