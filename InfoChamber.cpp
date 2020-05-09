#include "declarations.h"
#include "SimConfig.h"

void InfoChamber(int N, double Dt, double SampleRate,
                 double R,double T, double Eta,
                 double Lx, double Ly,double WallShrink,
                 int NumOfParticles, char* SaveFoldername, bool DisplayLive)
{
  double WallPositionsX[2];
  double WallPositionsY[2];
  SimConfig cfg;
  cfg.N = N;
  cfg.Dt = Dt;
  cfg.SampleRate = SampleRate;
  cfg.NumOfParticles = NumOfParticles;
  cfg.Eta = Eta;
  
//   WallPositionsX[0] = -Lx / 2;
//   WallPositionsX[1] = Lx /2 - WallShrink;
//   WallPositionsY[0] = -Ly / 2;
//   WallPositionsY[1] = Ly / 2;
//   cfg.WallPositionsX = WallPositionsX;
//   cfg.WallPositionsY = WallPositionsY;
//   Point* InitPositions = (Point*) malloc (cfg.NumOfParticles*sizeof(Point));
//   if(InitPositions==NULL) exit(10);
//   RandomizePositions(cfg.NumOfParticles, WallPositionsX, WallPositionsY, R, InitPositions);
//   cfg.InitPositions = InitPositions;

//   free(InitPositions);
}