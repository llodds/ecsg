#include "init.hpp"

#ifndef _GRID_HPP
#define _GRID_HPP

struct Grid {
    std::string grid_name;
    int size[3];
    real step[3];
    real origin[3];
};


void PrintGrid(Grid grid); //print.cpp
#endif