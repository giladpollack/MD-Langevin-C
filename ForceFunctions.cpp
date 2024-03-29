#include "declarations.h"
#include "Classes/AdditionalData.h"

double GetLJForce(double R, double Eps, double Sigma)
{
    double Force = 48*Eps*((1/R)*(pow(Sigma/R,12)-0.5*pow(Sigma/R, 6)));
    return Force;
}

void GetWCAParticleForces(Point* Positions, int NumOfParticles, double R, double Eps, double* Fx, double* Fy)
{
    double Rc = 2.5*R;

    for (int CurrParticle = 0; CurrParticle < NumOfParticles; CurrParticle++)
    {
        Point CheckedParticlePos = Positions[CurrParticle];

        // For each particle, checking its interaction with all others
        for (int PairedParticle = 0; PairedParticle < CurrParticle; PairedParticle++)
        {
            Point PairedParticlePos = Positions[PairedParticle];
            double Ydiff = CheckedParticlePos.y - PairedParticlePos.y;
            double Xdiff = CheckedParticlePos.x - PairedParticlePos.x;

            double Distance = sqrt(pow(Xdiff,2) + pow(Ydiff,2)) - R;
            if (Distance <= Rc)
            {
                // Computing the forces between the particles
                double TotalForce = GetLJForce(Distance, Eps, R);
                double Angle = atan2(Ydiff, Xdiff);

                // Adding the force to both particles with opposite signs
                Fx[CurrParticle]    += TotalForce*cos(Angle);
                Fx[PairedParticle]  -= TotalForce*cos(Angle);
                Fy[CurrParticle]    += TotalForce*sin(Angle);
                Fy[PairedParticle]  -= TotalForce*sin(Angle);
            }
        }
        
    }
}

void GetAthermalFluctForces( Point* ParticlePositions,
                            int NumOfParticles,
                            double FluctForce,
                            double* FluctDirections,
                            double* Fx,
                            double* Fy,
                            AdditionalData& AddedData, int SampleInd)
{
    for (int CurrParticle = 0; CurrParticle < NumOfParticles; CurrParticle++)
    {
        Fx[CurrParticle] += FluctForce * cos(FluctDirections[CurrParticle]);
        Fy[CurrParticle] += FluctForce * sin(FluctDirections[CurrParticle]);
    }
}

void GetWCAWallForces(  Point* ParticlePositions,
                        int NumOfParticles,
                        double R,
                        double Eps,
                        double* WallPositionsX,
                        double* WallPositionsY,
                        double* Fx,
                        double* Fy)
{
    double Rc = 2.5*R;

    for (int CurrParticle = 0; CurrParticle < NumOfParticles; CurrParticle++)
    {
        Point CheckedParticlePosition = ParticlePositions[CurrParticle];
        double DistanceFromRightWall = WallPositionsX[1] - CheckedParticlePosition.x;
        double DistanceFromLeftWall = CheckedParticlePosition.x - WallPositionsX[0];
        double DistanceFromTopWall = WallPositionsY[1] - CheckedParticlePosition.y;
        double DistanceFromBottomWall = CheckedParticlePosition.y - WallPositionsY[0];

        if (DistanceFromRightWall <= Rc)
        {
            Fx[CurrParticle] -= GetLJForce(DistanceFromRightWall,Eps,R);
        }

        if (DistanceFromLeftWall <= Rc)
        {
            Fx[CurrParticle] += GetLJForce(DistanceFromLeftWall,Eps,R);
        }

        if (DistanceFromTopWall <= Rc)
        {
            Fy[CurrParticle] -= GetLJForce(DistanceFromTopWall, Eps, R);
        }
        
        if (DistanceFromBottomWall <= Rc)
        {
            Fy[CurrParticle] += GetLJForce(DistanceFromBottomWall, Eps, R);
        }
        
    }
    
}

void GetHarmonicWallForces( Point* ParticlePositions,
                            int NumOfParticles,
                            double R,
                            double K,
                            double* WallPositionsX,
                            double* WallPositionsY,
                            double* Fx,
                            double* Fy)
{
    for (int CurrParticle = 0; CurrParticle < NumOfParticles; CurrParticle++)
    {
        Point CheckedParticlePosition = ParticlePositions[CurrParticle];
        double XBeyondLeftWall = WallPositionsX[0] - (CheckedParticlePosition.x - R);
        double XBeyondRightWall = CheckedParticlePosition.x + R - WallPositionsX[1];
        double YBeyondTopWall = CheckedParticlePosition.y + R - WallPositionsY[1];
        double YBeyondBottomWall = WallPositionsY[0] - (CheckedParticlePosition.y - R);

        if (XBeyondRightWall > 0)
        {
            Fx[CurrParticle] -= K * XBeyondRightWall;
        }
        if (XBeyondLeftWall > 0)
        {
            Fx[CurrParticle] += K * XBeyondLeftWall;
        }
        if (YBeyondTopWall > 0)
        {
            Fy[CurrParticle] -= K * YBeyondTopWall;
        }
        if (YBeyondBottomWall > 0)
        {
            Fy[CurrParticle] += K * YBeyondBottomWall;
        }
    }
}

void GetGaussianWallForces( Point* ParticlePositions,
                            int NumOfParticles,
                            double A,
                            double sigmaSq,
                            double* WallPositionsX,
                            double* WallPositionsY,
                            double* Fx,
                            double* Fy,
                            AdditionalData& AddedData, int SampleInd)
{
    for (int CurrParticle = 0; CurrParticle < NumOfParticles; CurrParticle++)
    {
        Point CheckedParticlePosition = ParticlePositions[CurrParticle];
        double LeftDist =  CheckedParticlePosition.x - WallPositionsX[0];
        double RightDist = WallPositionsX[1] - CheckedParticlePosition.x;
        double TopDist = WallPositionsY[1] - CheckedParticlePosition.y;
        double BottomDist = CheckedParticlePosition.y - WallPositionsY[0];

        double ForceRight = -(A/sigmaSq)*(LeftDist*exp(-pow(LeftDist,2)/(2*sigmaSq)));
        double ForceLeft = -(A/sigmaSq)*(RightDist*exp(-pow(RightDist,2)/(2*sigmaSq)));
        double ForceUp = -(A/sigmaSq)*(BottomDist*exp(-pow(BottomDist,2)/(2*sigmaSq)));
        double ForceDown = -(A/sigmaSq)*(TopDist*exp(-pow(TopDist,2)/(2*sigmaSq)));
        ForceRight = abs(ForceRight);
        ForceLeft = -abs(ForceLeft);
        ForceDown = -abs(ForceDown);
        ForceUp = abs(ForceUp);

        Fx[CurrParticle] += ForceLeft + ForceRight;
        Fy[CurrParticle] += ForceUp + ForceDown;

        AddedData.ForcesRight[SampleInd] += abs(ForceRight);
        AddedData.ForcesLeft[SampleInd] += abs(ForceLeft);
        AddedData.ForcesDown[SampleInd] += abs(ForceDown);
        AddedData.ForcesUp[SampleInd] += abs(ForceUp);
    }
}