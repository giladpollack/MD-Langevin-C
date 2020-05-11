class Matrix
{
    public:
    // Data Members
    double** Mat;
    int MatSize;

    // Ctors/Dtors
    Matrix() = delete;
    Matrix(int Size, bool InitToZero);

    ~Matrix();
};