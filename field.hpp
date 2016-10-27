#include <vector>
#include <cmath>
#include <iostream>
//#include "mat.h"
#include "init.hpp"
#include "grid.hpp"
#include <rsf.hh>

#ifndef _FIELD_HPP
#define _FIELD_HPP

class Field {
private:
    std::string field_name;
    std::string grid_name;
    //staggered(1) or not(0) on i0/i1/i2 dim
    int pos[3];
    //start/end index of computational domain (include boundary) on i0/i1/i2 dim.
    int start[3], end[3], len[3];
    real *ptr1;
    real **ptr2;
    real ***ptr3;
    real ***ptr;
    bool allocated; //identifies if this Field has been correctly allocated and ready to work with operations.
    
    //disable init function from calling by the user
    void Init();
    
    /*
     *  Declare() declares attributes of a field, e.g., field_name, grid_name,
     *  pos, start, end, len; it does not perform allocation!
     *
     *  A CheckDeclareError() is performed at the end of the function. If something
     *  is wrong, Declare() returns FALSE, otherwise returns TRUE, which indicates
     *  this field is ready to allocate.
     */
    bool Declare(std::string string_name_source,
                 std::string grid_name_source,
                 int grid_pos[3], int grid_size[3],
                 int radius);
    /*
     *  Allocate() allocates memory for correctly declared field,
     *  so that each value in the field can be obtained using GetVal(i0,i1,i2);
     *  It also sets "allocated" to TRUE, initiali2es all field values to 0.0.
     *
     *  It returns TRUE if the field is correctly allocated, otherwise it returns
     *  FALSE, and the field is untouched.
     */
    bool Allocate();
    /*
     *  Deallocate() deallocates memories pointed by ptr1, ptr2, ptr3, and
     *  resets them and ptr to nullptr, "allocated" to false.
     */
    void Deallocate();
    
    //member function to check if this Field has been correctly declared
    bool CheckDeclarationError() const;

public:
    Field() {Init();}
    Field(std::string string_name_source,
          std::string grid_name_source,
          int grid_pos[3], int grid_size[3],
          int radius);
    Field(const Field & field_source);
    ~Field();
    
    
    //member functions for retrieving and setting private data members
    std::string GetFieldName() const {return field_name;}
    void SetFieldName(std::string field_name_source) {field_name=field_name_source;}
    std::string GetGridName() const {return grid_name;}
    int GetPos(int i) const {return pos[i];}
    int GetStart(int i) const {return start[i];}
    int GetEnd(int i) const {return end[i];}
    int GetLen(int i) const {return len[i];}
    real*** GetPtr() const {return ptr;}
    real* GetPtr1() const {return ptr1;}
    
    
    //member functions for accessing and operations on field values
    //(member functions defined inside the class are inlined by default,
    // so using a member function to access value will not hurt the performance!)
    real GetVal(int i0, int i1, int i2) const {return ptr[i2][i1][i0];}
    real& SetVal(int i0, int i1, int i2) {return ptr[i2][i1][i0];}
    void SetOneVal(real val);
    real MaxVal() const;
    real MinVal() const;
    //returns the average of the absolute field values
    real AvgAbs() const;
    
    //member functions to overload operators
    Field & operator= (const Field & rhs);
    const Field operator+ (const Field & rhs) const;
    const Field operator- (const Field & rhs) const;
    const Field operator* (const Field & rhs) const;
    const Field operator* (real rhs) const;
    
    //member function to check if this Field has been correctly allocated
    bool IsAllocated() const {return allocated;}
};

Field InvField(const Field& f);
//(This function assumes the Field is on primary grid for all dims
//, and it reshapes the Field to desired Grid positions.)
Field ReshapeField(const Field &f, int grid_pos[3], std::string new_field_name="");
//(This function is created for computing the energy.)
real TriMulFields(const Field &fa, const Field &fb, const Field &fc);

void PrintFieldAttr(const Field& f);
#endif