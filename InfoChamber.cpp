#include "declarations.h"
#include "Classes/SimConfig.h"
#include "Classes/MDSim.h"
#include "Classes/RandGen.h"
#include "Utility.h"
#include "Classes/Vector.h"

void InfoChamber(int N, double Dt, double SampleRate,
                 double R,double T, double Eta,
                 double Lx, double Ly,double WallShrink,
                 int NumOfParticles, char* SaveFoldername, bool DisplayLive, RandGen rng)
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
  char WallRepulsionType[20] = "WCA";
  Cfg.WallRepulsionType = WallRepulsionType;
  Cfg.UseHydro = true;
  Cfg.UseParticleRepulsion = true;
  Cfg.UseTraps = false;
  Cfg.WCAEpsilon = 0.2*kB*Cfg.T;
  
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