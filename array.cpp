#include <new>
#include "array.hpp"
#include <algorithm> //std::fill
#include <iostream>
#include "init.hpp"
#include <omp.h>
#include <cstring>

void Array::Init()
{
    n0 = 0; n1 = 0; n2 = 0;
    ptr1 = nullptr; ptr2 = nullptr; ptr3 = nullptr;
}

void Array::Destroy()
{
    delete [] ptr1;
    delete [] ptr2;
    delete [] ptr3;
}

void Array::Allocate(int n0)
{
    try{
        ptr1 = new real [n0];
    }
    catch (std::bad_alloc &e)
    {
        std::cerr << "ERROR: Memory Allocation failed: " << e.what() <<std::endl;
        std::cerr << "File [" << __FILE__ << "], Line [" << __LINE__ << "]." << std::endl;
        Destroy();
        exit(1);
    }
    
    std::fill(ptr1,ptr1+n0,0.0);
}

void Array::Allocate(int n0, int n1)
{
    try{
        ptr1 = new real [n1*n0];
        ptr2 = new real* [n1];
    }
    catch (std::bad_alloc &e)
    {
        std::cerr << "ERROR: Memory Allocation failed: " << e.what() <<std::endl;
        std::cerr << "File [" << __FILE__ << "], Line [" << __LINE__ << "]." << std::endl;
        Destroy();
        exit(1);
    }
    
    for(int i1 = 0; i1 < n1; i1++)
        ptr2[i1] = &ptr1[i1*n0];
    
    std::fill(ptr1,ptr1+n0*n1,0.0);
}

void Array::Allocate(int n0, int n1, int n2)
{
    try{
        ptr1 = new real [n0*n1*n2];
        ptr2 = new real* [n1*n2];
        ptr3 = new real** [n2];
    }
    catch (std::bad_alloc &e)
    {
        std::cerr << "ERROR: Memory Allocation failed: " << e.what() <<std::endl;
        std::cerr << "File [" << __FILE__ << "], Line [" << __LINE__ << "]." << std::endl;
        Destroy();
        exit(1);
    }
    
    //TODO usigned int
    for(int i1 = 0; i1 < n2*n1; i1++)
        ptr2[i1] = &ptr1[i1*n0];
    for(int i2 = 0; i2 < n2; i2++)
        ptr3[i2] = &ptr2[i2*n1];
    
    std::fill(ptr1,ptr1+n0*n1*n2,0.0);
}


Array::Array()
{
    Init();
}

Array::Array(int n0)
{
    Init();
    Allocate(n0);
    this->n0 = n0;
}


Array::Array(int n0, int n1)
{
    Init();
    Allocate(n0, n1);
    this->n1 = n1;
    this->n0 = n0;
}

Array::Array(int n0, int n1, int n2)
{
    Init();
    Allocate(n0, n1, n2);
    this->n2 = n2;
    this->n1 = n1;
    this->n0 = n0;
}

Array::~Array()
{
    Destroy();
}

void Array::size(int *n0) const
{
    *n0 = this->n0;
}

void Array::size(int *n0, int *n1) const
{
    *n0 = this->n0;
    *n1 = this->n1;
}

void Array::size(int *n0, int *n1, int *n2) const
{
    *n0 = this->n0;
    *n1 = this->n1;
    *n2 = this->n2;
}


Array& Array::operator= (const Array &rhs)
{
    int rhs_n0, rhs_n1, rhs_n2;
    rhs.size(&rhs_n0, &rhs_n1, &rhs_n2);
    
    bool flag = false;
    if(n0!=rhs_n0 || n1!=rhs_n1 || n2!=rhs_n2)
    {
        flag = true;
        Destroy();
        n0 = rhs_n0;
        n1 = rhs_n1;
        n2 = rhs_n2;
    }
    
    if(n2 > 0)
    {
        if(flag) Allocate(n0,n1,n2);
        memcpy(ptr1, rhs.GetPtr1(), n0*n1*n2);
    }
    else if(n1 > 0)
    {
        if(flag) Allocate(n0,n1);
        memcpy(ptr1, rhs.GetPtr1(), n0*n1);
    }
    else
    {
        if(flag) Allocate(n0);
        memcpy(ptr1, rhs.GetPtr1(), n0);
    }
    
    return (*this);
}


real Array::AvgAbs()
{
    real AvgAbs=0.0;
    
    if(n2!=0)
    {
#ifdef _OPENMP
#pragma omp parallel for default(shared) reduction(+:AvgAbs)
#endif
        for(int i2 = 0; i2 < n2; i2++)
            for(int i1 = 0; i1 < n1; i1++)
                for(int i0 = 0; i0 < n0; i0++)
                    AvgAbs += fabs(GetVal(i0,i1,i2));
        
        AvgAbs /= n2*n1*n0;
    }
    else if(n1!=0)
    {
#ifdef _OPENMP
#pragma omp parallel for default(shared) reduction(+:AvgAbs)
#endif
        for(int i1 = 0; i1 < n1; i1++)
            for(int i0 = 0; i0 < n0; i0++)
                AvgAbs += fabs(GetVal(i0,i1));
        
        AvgAbs /= n1*n0;
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for default(shared) reduction(+:AvgAbs)
#endif
        for(int i0 = 0; i0 < n0; i0++)
            AvgAbs += fabs(GetVal(i0));
        
        AvgAbs /= n0;
    }
    
    return AvgAbs;
}