#ifndef gmx2pqr_my_structs_hpp
#define gmx2pqr_my_structs_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cmath>
#include <cstring>



#ifndef MPI
#define MPI  4.*atan(1.)
#endif
//Vacuum permittivity C^2 s^2 /(kg m^3)
#define eps0  8.8541878e-12
// Coulombs/e-
#define ec  1.6021773e-19
// Avogadro's number, molecules/mol
#define NA  6.0221367e23
// J/K
#define boltzmann  1.3806581e-23
#define zscale  (ec / (4 * MPI * eps0 * 1.e-10))
#define zmagic  (ec * NA * 1.e-3)
#define kbT  (boltzmann * 300. * NA * 1.e-3)
#define cfac  (zscale * zmagic / kbT)



typedef struct gmx2amb gmx2amb;

struct gmx2amb {
    char resname[256];
    char gmxname[256];
    char ambername[256];
    float charge, radius;
};




#endif
