#ifndef RandGen_h
#define RandGen_h

#include <random>
class RandGen
{
    public:
    
    void Reseed(long idum);
    void Randomize();
    double Randn();
    double Randu(double LBound, double UBound);

    long idum;
    float ran_nrc(long* idum);
    std::default_random_engine generator;

    RandGen() = delete;
    RandGen(long idum);
};
#endif