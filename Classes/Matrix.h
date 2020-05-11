#ifndef Matrix_h
#define Matrix_h
#include "../declarations.h"

template <typename T> 
class Matrix
{
    public:
    // Data Members
    T** Mat;
    int NumRows;
    int NumCols;

    // Ctors/Dtors
    Matrix() = delete;
    Matrix(int NumRows, int NumCols);

    ~Matrix();
};

template <typename T> 
Matrix<T>::Matrix(int NumRows, int NumCols)
{
    this->NumRows = NumRows;
    this->NumCols = NumCols;
    this->Mat = new T*[NumRows];
    for (int CurrRow = 0; CurrRow < NumRows; CurrRow++)
    {
        this->Mat[CurrRow] = new T[NumCols];
        
    }
}

template <typename T> 
Matrix<T>::~Matrix()
{
    for (int CurrRow = 0; CurrRow < this->NumRows; CurrRow++)
    {
        delete[] this->Mat[CurrRow];
    }

    delete[] this->Mat;
}
#endif