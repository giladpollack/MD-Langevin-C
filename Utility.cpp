#include "declarations.h"
#include "Classes/RandGen.h"
#include <sys/types.h>
#include <sys/stat.h>

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
void RandomizePositions(int NumOfParticles, double* WallPositionsX, double* WallPositionsY, double R, Point* Positions, RandGen rng)
{
  double LBoundX = WallPositionsX[0] + R;
  double UBoundX = WallPositionsX[1] - R;
  double LBoundY = WallPositionsY[0] + R;
  double UBoundY = WallPositionsY[1] - R;
  bool Collision;
  int MaxIterations = 1000;
  int CurrIteration = 0;
  for (int CurrParticle = 0; CurrParticle < NumOfParticles; CurrParticle++)
  {

      
    CurrIteration = 0;
    Collision = true;
    // Running until there is no collision
    while (Collision && CurrIteration < MaxIterations)
    {
        Collision = false;
        Positions[CurrParticle].x = rng.Randu(LBoundX, UBoundX);
        Positions[CurrParticle].y = rng.Randu(LBoundY, UBoundY);
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
        CurrIteration++;
    }
    if (CurrIteration >= MaxIterations - 1)
    {
        CurrParticle = -1;
    }
  }
}

/**
 * Copies values of points from the source array to the target array.
 * The target array MUST BE ALLOCATED BEFORE CALLING THIS FUNCTiON!
 *
 * @TargetArray The array of points into which values will be copied. MUST BE ALLOCATED!
 * @SourceArray The array from which the values are taken.
 * @NumOfParticles The number of particles to copy from the source to the target
 */
void CopyPositions(Point* TargetArray, Point* SourceArray, int NumOfParticles)
{
    for (int i = 0; i < NumOfParticles; i++)
    {
        TargetArray[i].x = SourceArray[i].x;
        TargetArray[i].y = SourceArray[i].y;
    }
}

/**
 * Creates a string of the particle positions. OutString MUST BE ALLOCATED BEFORE CALLING THE FUNCTION!
 *
 * @Positions The positions of the particles
 * @NumOfParticles The number of particles to get positions for
 * @Outstring The output parameter returning the string of positions
 */
void GetPositionsString(Point* Positions, int NumOfParticles, char* OutString)
{
    for (int i = 0; i < NumOfParticles; i++)
    {
        char CurrPos[50];
        sprintf(CurrPos, "%e,%e\n", Positions[i].x, Positions[i].y);
        strcat(OutString, CurrPos);
    }
}

/**
 * Creates a string of a single coordinate (x or y) of the particle position.
 * OutString MUST BE ALLOCATED BEFORE CALLING THE FUNCTION!
 *
 * @ParticlePositions The positions of the particles
 * @NumOfParticles The number of particles to get positions for
 * @axis Which axis to get
 * @MultFactor A factor to allow the user to save in desired units (micrometers for instance)
 * @Outstring The output parameter returning the string of positions
 */
void GetSingleAxisSavedSteps(Point* ParticlePositions, int NumOfParticles, char axis, char* OutString)
{
    OutString[0] = '\0';

    for (int i = 0; i < NumOfParticles; i++)
    {
        char TempString[20];
        if (axis == 'x')
        {
            sprintf(TempString, "%e,",ParticlePositions[i].x);
        }
        else if (axis == 'y')
        {
            sprintf(TempString, "%e,",ParticlePositions[i].y);
        }
        else
        {
            throw;
        }
        
        strcat(OutString,TempString);
    }

    OutString[strlen(OutString) - 1] = '\n';
    
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

/**
 * Checks if a directory exists
 *
 * @Path The path to the directory
 * @return true if the directory exists and is a directory, false otherwise
 */
bool doesDirExist(char* Path)
{
    struct stat info;
    bool res = false;
    if( stat( Path, &info ) != 0 )
    {
        //printf( "cannot access %s\n", Path );
    }
    else if( info.st_mode & S_IFDIR )
    {
        res = true;
    }
    
    return res;
}

/**
 * Sets the values of an array to be 0
 *
 * @Array The array to set
 * @ArrayLen The number of cells in the array
 */
void SetToZero(double* Array, int ArrayLen)
{
    for (int i = 0; i < ArrayLen; i++)
    {
        Array[i] = 0;
    }
}