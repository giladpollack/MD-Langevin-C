#include "declarations.h"
#include "RandGen.h"

// Declarations
bool checkCollision(double x1, double y1, double r1,  
                    double x2, double y2, double r2);

// Functions

/**
 * Randomizes the positions of particles within a given rectangle defined by the wall positions while preventing overlap
 *
 * @NumOfParticles The number of particles to randomize
 * @WallPositionsX The lower and upper bounds of the rectangle in X
 * @WallPositionsY The lower and upper bounds of the rectangle in Y
 * @R The radius of the particles
 * @Positions The output parameter in which to save the positions
 */
void RandomizePositions(int NumOfParticles, double* WallPositionsX, double* WallPositionsY, double R, Point* Positions)
{
  double LBoundX = WallPositionsX[0] + R;
  double UBoundX = WallPositionsX[1] - R;
  double LBoundY = WallPositionsY[0] + R;
  double UBoundY = WallPositionsY[1] - R;
  bool Collision;
  for (int CurrParticle = 0; CurrParticle < NumOfParticles; CurrParticle++)
  {
      Collision = true;
      // Running until there is no collision
      while (Collision)
      {
          Collision = false;
          Positions[CurrParticle].x = RandGen::Randu(LBoundX, UBoundX);
          Positions[CurrParticle].y = RandGen::Randu(LBoundY, UBoundY);
          for (int CollidedParticle = 0; CollidedParticle < CurrParticle; CollidedParticle++)
          {
              bool CurrCollision = checkCollision(Positions[CurrParticle].x, Positions[CurrParticle].y, R,
                                                  Positions[CollidedParticle].x, Positions[CollidedParticle].y, R);
              if (CurrCollision)
              {
                  Collision = true;
                  break;
              }
          }
          
      }
      
  }
}

/**
 * Checks if two circles overlap
 *
 * @x1 The x coordinate of the center of the first circle
 * @y1 The y coordinate of the center of the first circle
 * @r1 The radius of the first circle
 * @x2 The x coordinate of the center of the second circle
 * @y2 The y coordinate of the center of the second circle
 * @r2 The radius ot the second circle
 * @return true if the circles overlap, false if they don't
 */
bool checkCollision(double x1, double y1, double r1,  
                    double x2, double y2, double r2)
{ 
    double distSq = (x1 - x2) * (x1 - x2) + 
                 (y1 - y2) * (y1 - y2); 
    double radSumSq = (r1 + r2) * (r1 + r2); 
    if (distSq < radSumSq)
        return true; 
    else
        return false; 
} 