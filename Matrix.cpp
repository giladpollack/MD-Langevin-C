#include "Matrix.h"
#include "declarations.h"

Matrix::Matrix(int Size, bool InitToZero)
{
    this->MatSize = Size;
    this->Mat = (double**) malloc(Size*sizeof(double*));
    for (int CurrRow = 0; CurrRow < Size; CurrRow++)
    {
        this->Mat[CurrRow] = (double*) malloc(Size*sizeof(double));
        if (InitToZero)
        {
            for (int CurrCol = 0; CurrCol < Size; CurrCol++)
            {
                this->Mat[CurrRow][CurrCol] = 0;
            }
            
        }
        
    }
}

Matrix::~Matrix()
{
    for (int CurrRow = 0; CurrRow < this->MatSize; CurrRow++)
    {
        free(this->Mat[CurrRow]);
    }
    free(this->Mat);
}