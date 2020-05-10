#include "declarations.h"

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

            double Distance = sqrt(pow(Xdiff,2) + pow(Ydiff,2));
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