class RandGen
{
    public:
    static void Init(long* idum);
    static void Randomize();
    static double Randn();
    static double Randu(double LBound, double UBound);

    // private:
    static long idum;
    static float ran_nrc(long* idum);
};