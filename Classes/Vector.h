#ifndef Vector_h
#define Vector_h
#include <iostream> 
using namespace std; 

/**
 * A class responsible for wrapping memory allocations for dynamical vectors of any type
 *
 */
template <typename T> 
class Vector { 
private: 
    int size; 
  
public:
    T* ptr;
    Vector() = delete;
    Vector(int Len); 
    void print();
    ~Vector<T>();
};


template <typename T> 
Vector<T>::Vector(int Len) 
{
    ptr = new T[Len];
    if (ptr == NULL)
    {
        throw;
    }
    
}
template <typename T>
Vector<T>::~Vector()
{
    delete[] ptr;
}
#endif