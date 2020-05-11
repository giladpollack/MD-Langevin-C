#ifndef RandGen_h
#define RandGen_h
class RandGen
{
    public:
    
    void Init(long* idum);
    void Randomize();
    double Randn();
    double Randu(double LBound, double UBound);

    long idum;
    float ran_nrc(long* idum);

    RandGen() = delete;
    RandGen(long* idum);
};
#endif