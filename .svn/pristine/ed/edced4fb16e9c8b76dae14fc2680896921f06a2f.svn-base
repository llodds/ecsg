#include <vector>
#include <cmath>
#include <iostream>
#include "init.hpp"
#include "grid.hpp"
#include "field.hpp"
#include <cassert>
#include <new>
#include <algorithm> //std::fill
#include <valarray>
#include <rsf.hh>
#include <omp.h>


void Field::Init()
{
    for(int i = 0; i < 3; i++)
    { pos[i] = -1; start[i] = 0; end[i] = 0; len[i] = -1; }
    ptr1 = nullptr; ptr2 = nullptr; ptr3 = nullptr; ptr = nullptr;
    allocated = false;
}

Field::Field(std::string string_name_source,
             std::string grid_name_source,
             int grid_pos[3], int grid_size[3],
             int radius)
{
    Init();
    
    if(!Declare(string_name_source,
                grid_name_source,
                grid_pos, grid_size,
                radius))
    {
        std::cerr << "ERROR: Failed to declare Field! "<<std::endl;
        std::cerr << "File [" << __FILE__ << "], Line [" << __LINE__ << "]." << std::endl;
        exit(-1);
    }
    if(!Allocate())
    {
        std::cerr << "ERROR: Failed to allocate Field! "<<std::endl;
        std::cerr << "File [" << __FILE__ << "], Line [" << __LINE__ << "]." << std::endl;
        exit(-1);
    }
}

Field::Field(const Field & field_source)
{
    Init();
    
    field_name = field_source.GetFieldName();
    grid_name = field_source.GetGridName();
    
    for(int i = 0; i < 3; i++)
    {
        pos[i] = field_source.GetPos(i);
        start[i] = field_source.GetStart(i);
        end[i] = field_source.GetEnd(i);
        len[i] = field_source.GetLen(i);
    }
    
    if(!Allocate())
    {
        std::cerr << "ERROR: Failed to allocate Field! "<<std::endl;
        std::cerr << "File [" << __FILE__ << "], Line [" << __LINE__ << "]." << std::endl;
        exit(-1);
    }
    
    memcpy(ptr1, field_source.GetPtr1(), sizeof(real)*len[2]*len[1]*len[0]);
}

Field::~Field()
{
    Deallocate();
}

bool Field::Declare(std::string string_name_source,
                    std::string grid_name_source,
                    int grid_pos[3], int grid_size[3],
                    int radius)
{
    field_name = string_name_source;
    grid_name = grid_name_source;
    
    for(int i = 0; i < 3; i++)
    {
        pos[i] = grid_pos[i];
        if(pos[i]==1)   //staggered dim
        {
            start[i] = 0;
            end[i] = grid_size[i]-1;
            len[i] = grid_size[i]+2*radius;
        }else{  //non-staggered dim
            start[i] = 0;
            end[i] = grid_size[i];
            len[i] = grid_size[i]+1+2*radius;
        }
    }
    
    if(CheckDeclarationError())
        return false;
    
    return true;
}

bool Field::Allocate()
{
    //make sure the Field has been correctly declared.
    if(CheckDeclarationError())
    {
        std::cerr << "ERROR: This field has not been correctly declared!" << std::endl;
        return false;
    }
    
    int radius[3];//radius MIGHT be different on each dim.
    for(int i = 0; i < 3; i++)
        radius[i] = (len[i] - (end[i]-start[i]+1))/2;
    
    //allocate memories
    try{
        ptr1 = new real [len[0]*len[1]*len[2]];
        ptr2 = new real* [len[1]*len[2]];
        ptr3 = new real** [len[2]];
    }catch(std::bad_alloc& e)
    {
        std::cerr << "ERROR: Memory Allocation failed: " << e.what() <<std::endl;
        std::cerr << "File [" << __FILE__ << "], Line [" << __LINE__ << "]." << std::endl;
        
        delete [] ptr1; ptr1 = nullptr;
        delete [] ptr2; ptr2 = nullptr;
        delete [] ptr3; ptr3 = nullptr;
        
        return false;
    }
    
    //set pointers
    for(int i = 0; i < len[1]*len[2]; i++)
        ptr2[i] = &ptr1[radius[0] + i*len[0]];
    for(int i = 0; i < len[2]; i++)
        ptr3[i] = &ptr2[radius[1] + i*len[1]];
    ptr = &ptr3[radius[2]];
    
    //initiali0e all values to 0
    std::fill(ptr1,ptr1+len[0]*len[1]*len[2],0.0);
    
    allocated = true;
    
    return true;
}

void Field::Deallocate()
{
    //c++ allows to delete a nullptr and it already performs the check
    delete [] ptr1; ptr1 = nullptr;
    delete [] ptr2; ptr2 = nullptr;
    delete [] ptr3; ptr3 = nullptr;
    ptr = nullptr;
    allocated = false;
}

void Field::SetOneVal(real val)
{
    assert(allocated);
    
#ifdef _OPENMP
#pragma omp parallel for default(shared) firstprivate(val)
#endif
    for(int i2 = start[2]; i2 <= end[2]; i2++)
        for(int i1 = start[1]; i1 <= end[1]; i1++)
            for(int i0 = start[0]; i0 <= end[0]; i0++)
                SetVal(i0,i1,i2) = val;
}

real Field::MaxVal() const
{
    assert(allocated);
    
    real max_val = GetVal(start[0], start[1], start[2]);
    
#ifdef _OPENMP
#pragma omp parallel for default(shared) reduction(max:max_val)
#endif
    for(int i2 = start[2]; i2 <= end[2]; i2++)
        for(int i1 = start[1]; i1 <= end[1]; i1++)
            for(int i0 = start[0]; i0 <= end[0]; i0++)
                if(max_val < GetVal(i0,i1,i2))
                    max_val = GetVal(i0,i1,i2);
    
    return max_val;
}

real Field::MinVal() const
{
    assert(allocated);
    
    real min_val = GetVal(start[0], start[1], start[2]);
    
#ifdef _OPENMP
#pragma omp parallel for default(shared) reduction(min:min_val)
#endif
    for(int i2 = start[2]; i2 <= end[2]; i2++)
        for(int i1 = start[1]; i1 <= end[1]; i1++)
            for(int i0 = start[0]; i0 <= end[0]; i0++)
                if(min_val > GetVal(i0,i1,i2))
                    min_val = GetVal(i0,i1,i2);
    
    return min_val;
}

real Field::AvgAbs() const
{
    assert(allocated);
    
    real AvgAbs = 0.0;
    
#ifdef _OPENMP
#pragma omp parallel for default(shared) reduction(+:AvgAbs)
#endif
    for(int i2 = start[2]; i2 <= end[2]; i2++)
        for(int i1 = start[1]; i1 <= end[1]; i1++)
            for(int i0 = start[0]; i0 <= end[0]; i0++)
                AvgAbs += fabs(GetVal(i0,i1,i2));
    
    AvgAbs /= (end[2]-start[2]+1);
    AvgAbs /= (end[1]-start[1]+1);
    AvgAbs /= (end[0]-start[0]+1);
    
    return AvgAbs;
}

Field & Field::operator= (const Field & rhs)
{
    //check if rhs exists
    assert(rhs.IsAllocated());
    
    bool flag = true;
    if(allocated && GetGridName() == rhs.GetGridName())
    {
        flag = false;
        for(int i = 0; i < 3; i++)
            if(GetPos(i) != rhs.GetPos(i))
            {flag = true; break;}
    }
    
    if(flag)
    {
        Deallocate();
        
        field_name = rhs.GetFieldName();
        grid_name = rhs.GetGridName();
        
        for(int i = 0; i < 3; i++)
        {
            pos[i] = rhs.GetPos(i);
            start[i] = rhs.GetStart(i);
            end[i] = rhs.GetEnd(i);
            len[i] = rhs.GetLen(i);
        }
        
        if(!Allocate())
        {
            std::cerr << "ERROR: Failed to allocate Field! "<<std::endl;
            std::cerr << "File [" << __FILE__ << "], Line [" << __LINE__ << "]." << std::endl;
            exit(-1);
        }
    }
    
    memcpy(ptr1, rhs.GetPtr1(), sizeof(real)*len[2]*len[1]*len[0]);
    
    return (*this);
}

const Field Field::operator+ (const Field & rhs) const
{
    //check if they are on the same grid
    assert(GetGridName() == rhs.GetGridName());
    
    //check if they are on the same grid position
    for(int i = 0; i < 3; i++)
        assert(GetPos(i)==rhs.GetPos(i));
    
    //check if they are correctly allocated
    assert(allocated);
    assert(rhs.IsAllocated());
    
    Field res((*this));
    
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int i2 = start[2]; i2 <= end[2]; i2++)
        for(int i1 = start[1]; i1 <= end[1]; i1++)
            for(int i0 = start[0]; i0 <= end[0]; i0++)
                res.SetVal(i0,i1,i2) += rhs.GetVal(i0,i1,i2);
    
    return res;
}

const Field Field::operator- (const Field & rhs) const
{
    //check if they are on the same grid
    assert(GetGridName() == rhs.GetGridName());
    
    //check if they are on the same grid position
    for(int i = 0; i < 3; i++)
        assert(GetPos(i)==rhs.GetPos(i));
    
    //check if they are correctly allocated
    assert(allocated);
    assert(rhs.IsAllocated());
    
    Field res((*this));
    
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int i2 = start[2]; i2 <= end[2]; i2++)
        for(int i1 = start[1]; i1 <= end[1]; i1++)
            for(int i0 = start[0]; i0 <= end[0]; i0++)
                res.SetVal(i0,i1,i2) -= rhs.GetVal(i0,i1,i2);
    
    return res;
}

const Field Field::operator* (const Field & rhs) const
{
    //check if they are on the same grid
    assert(GetGridName() == rhs.GetGridName());
    
    //check if they are on the same grid position
    for(int i = 0; i < 3; i++)
        assert(GetPos(i)==rhs.GetPos(i));
    
    //check if they are correctly allocated
    assert(allocated);
    assert(rhs.IsAllocated());
    
    Field res((*this));
    
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int i2 = start[2]; i2 <= end[2]; i2++)
        for(int i1 = start[1]; i1 <= end[1]; i1++)
            for(int i0 = start[0]; i0 <= end[0]; i0++)
                res.SetVal(i0,i1,i2) *= rhs.GetVal(i0,i1,i2);
    
    return res;
}

const Field Field::operator* (real rhs) const
{
    assert(allocated);
    
    Field res((*this));
    
#ifdef _OPENMP
#pragma omp parallel for default(shared) firstprivate(rhs)
#endif
    for(int i2 = start[2]; i2 <= end[2]; i2++)
        for(int i1 = start[1]; i1 <= end[1]; i1++)
            for(int i0 = start[0]; i0 <= end[0]; i0++)
                res.SetVal(i0,i1,i2) *= rhs;
    
    return res;
}

bool Field::CheckDeclarationError() const
{
    if(field_name == "")
    {
        std::cerr << "ERROR: Field::field_name==NULL."<< std::endl;
        std::cerr << "File [" << __FILE__ << "], Line [" << __LINE__ << "]." << std::endl;
        return true;
    }
    
    if(grid_name == "")
    {
        std::cerr << "ERROR: Field::grid_name==NULL."<< std::endl;
        std::cerr << "File [" << __FILE__ << "], Line [" << __LINE__ << "]." << std::endl;
        return true;
    }
    
    for(int i = 0; i < 3; i++)
    {
        assert(pos[i]==0||pos[i]==1);
        assert(start[i]<=end[i]);
        assert(len[i]>0);
    }
    
    return false;
}


Field InvField(const Field &f)
{
    //make sure there is something to invert
    assert(f.IsAllocated());
    
    Field res(f);
    
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int i2 = res.GetStart(2); i2 <= res.GetEnd(2); i2++)
        for(int i1 = res.GetStart(1); i1 <= res.GetEnd(1); i1++)
            for(int i0 = res.GetStart(0); i0 <= res.GetEnd(0); i0++)
            {
                //make sure each value is big enough to invert
                if(std::abs(res.GetVal(i0,i1,i2))<EPS)
                {
#ifdef _OPENMP
#pragma omp critical
#endif
                    {
                        std::cerr << "ERROR: can't invert ["<<res.GetFieldName()<<"], GetVal["
                        <<i2<<"]["<<i1<<"]["<<i0<<"]="<<res.GetVal(i0,i1,i2)<<"."<<std::endl;
                        std::cerr << "File [" << __FILE__ << "], Line [" << __LINE__ << "]." << std::endl;
                        exit(-1);
                    }
                }
                res.SetVal(i0,i1,i2) = 1/res.GetVal(i0,i1,i2);
            }
    
    return res;
}

Field ReshapeField(const Field &f, int grid_pos[3], std::string new_field_name)
{
    //check if reshaping is needed
    if(grid_pos[0]==f.GetPos(0)&&grid_pos[1]==f.GetPos(1)&&grid_pos[2]==f.GetPos(2))
        return f;
    
    //make sure grid_pos[] is a valid position on the grid
    for(int i = 0; i < 3; i++)
        assert(grid_pos[i]==0||grid_pos[i]==1);
    
    //make sure the field is on the primary grid on all dims
    for(int i = 0; i < 3; i++)
        assert(f.GetPos(i)==0);
    
    //make sure this field f has been correctly allocated
    assert(f.IsAllocated());
    
    //reshape this field!
    Field res;
    
    int grid_size[3];
    for(int i = 0; i < 3; i++) grid_size[i] = f.GetEnd(i)-f.GetStart(i);
    
    std::string field_name = (new_field_name=="")? f.GetFieldName():new_field_name;
    res = Field(field_name,
                f.GetGridName(),
                grid_pos, grid_size,
                0);
    
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int i2 = res.GetStart(2); i2 <= res.GetEnd(2); i2++)
        for(int i1 = res.GetStart(1); i1 <= res.GetEnd(1); i1++)
            for(int i0 = res.GetStart(0); i0 <= res.GetEnd(0); i0++)
            {
                
                res.SetVal(i0,i1,i2) = (f.GetVal(i0,i1,i2)
                                        +f.GetVal(i0,i1,i2+grid_pos[2])
                                        +f.GetVal(i0,i1+grid_pos[1],i2)
                                        +f.GetVal(i0,i1+grid_pos[1],i2+grid_pos[2])
                                        +f.GetVal(i0+grid_pos[0],i1,i2)
                                        +f.GetVal(i0+grid_pos[0],i1,i2+grid_pos[2])
                                        +f.GetVal(i0+grid_pos[0],i1+grid_pos[1],i2)
                                        +f.GetVal(i0+grid_pos[0],i1+grid_pos[1],i2+grid_pos[2]))/8;
            }
    
    return res;
}


/* Pointwise multiplication among fa, fb, fc on their computational domain.
 * This multiplication takes care of special boundary treatments.
 */
real TriMulFields(const Field &fa, const Field &fb, const Field &fc)
{
    //make sure they are on the same grid
    assert(fa.GetGridName()==fb.GetGridName()&&fa.GetGridName()==fc.GetGridName());
    
    //make sure they are on the same grid position
    for(int i=0; i<3; i++)
        assert(fa.GetPos(i)==fb.GetPos(i)&&fa.GetPos(i)==fc.GetPos(i));
    
    //make sure they have been correctly allocated
    assert(fa.IsAllocated());
    assert(fb.IsAllocated());
    assert(fc.IsAllocated());
    
    //Do the TriMul!
    real res = 0;
    
#ifdef _OPENMP
#pragma omp parallel for default(shared) reduction(+:res)
#endif
    for(int i2 = fa.GetStart(2); i2 <= fa.GetEnd(2); i2++)
    {
        real mask2 = 1-(fa.GetPos(2) ? 0:0.5*((i2==fa.GetStart(2))||(i2==fa.GetEnd(2))));
        for(int i1 = fa.GetStart(1); i1 <= fa.GetEnd(1); i1++)
        {
            real mask1 = 1-(fa.GetPos(1) ? 0:0.5*((i1==fa.GetStart(1))||(i1==fa.GetEnd(1))));
            for(int i0 = fa.GetStart(0); i0 <= fa.GetEnd(0); i0++)
            {
                //TODO might not be vectori0ed
                real mask0 = 1-(fa.GetPos(0) ? 0:0.5*((i0==fa.GetStart(0))||(i0==fa.GetEnd(0))));
                res += mask2*mask1*mask0*fa.GetVal(i0,i1,i2)*fb.GetVal(i0,i1,i2)*fc.GetVal(i0,i1,i2);
            }
        }
    }
    
    return res;
}

void PrintFieldAttr(const Field& f)
{
    std::cout << "====================" <<std::endl;
    
    std::cout << "Field Name: ["<<f.GetFieldName()<<"]"<<std::endl;
    std::cout << "Grid Name: ["<<f.GetGridName()<<"]"<<std::endl;
    
    std::cout << "Pos: ["<<f.GetPos(0)<<","
    <<f.GetPos(1)<<","<<f.GetPos(2)<<"]"<<std::endl;
    
    std::cout << "Start: ["<<f.GetStart(0)<<","
    <<f.GetStart(1)<<","<<f.GetStart(2)<<"]"<<std::endl;
    
    std::cout << "End: ["<<f.GetEnd(0)<<","
    <<f.GetEnd(1)<<","<<f.GetEnd(2)<<"]"<<std::endl;
    
    std::cout << "Len: ["<<f.GetLen(0)<<","
    <<f.GetLen(1)<<","<<f.GetLen(2)<<"]"<<std::endl;
    
    std::cout << "Min: ["<<f.MinVal()<<"]"<<std::endl;
    
    std::cout << "Max: ["<<f.MaxVal()<<"]"<<std::endl;
    
    std::cout << "====================" <<std::endl;
}









