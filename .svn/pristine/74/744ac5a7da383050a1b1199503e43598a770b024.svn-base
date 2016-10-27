#include "init.hpp"

#ifndef _ARRAY_HPP
#define _ARRAY_HPP

class Array {
private:
    int n0,n1,n2;   //fastest/slower/slowest dim
    
    real ***ptr3;
    real **ptr2;
    real *ptr1;
    
    void Init();
    void Destroy();
    void Allocate(int n0);
    void Allocate(int n0, int n1);
    void Allocate(int n0, int n1, int n2);
    
    Array(Array &A); //disable copy constructor
    
public:
    Array();
    Array(int n0);
    Array(int n0, int n1);
    Array(int n0, int n1, int n2);
    ~Array();
    
    real GetVal(int i0) const {return ptr1[i0];};
    real GetVal(int i0, int i1) const {return ptr2[i1][i0];};
    real GetVal(int i0, int i1, int i2) const {return ptr3[i2][i1][i0];};
    real *GetPtr1() const {return ptr1;}
    
    real& SetVal(int i0) {return ptr1[i0];};
    real& SetVal(int i0, int i1) {return ptr2[i1][i0];};
    real& SetVal(int i0, int i1, int i2) {return ptr3[i2][i1][i0];};
    
    //return size of this array
    void size(int *n0) const;
    void size(int *n0, int *n1) const;
    void size(int *n0, int *n1, int *n2) const;
    
    //overload operator=
    Array& operator= (const Array &rhs);
    
    //return absolute averge value
    real AvgAbs();
};

#endif
