#ifndef gmx2pqr_electrostatics_hpp
#define gmx2pqr_electrostatics_hpp
#include "my_structs.hpp"
#include "linalg.hpp"

#ifdef __cplusplus
extern "C"
#endif
#include <gromacs/copyrite.h>
#include <gromacs/gmx_fatal.h>
#include <gromacs/vec.h>

typedef struct
{
    double field_proj_total;
    double field_exclude[3];
    std::vector<double> field_proj_each_exclude;
    std::vector<std::vector<double> > field_each_exclude;
    double field_proj_exclude;
    double field_mid[3];
    double field_a1[3];
    double field_a2[3];
    std::vector<double> protein_site;
    std::vector<double> water_site;
} t_electrostatics;

double calculate_potential(const double &r, const int &j, const t_topology &top, const t_trxframe &fr);
double calculate_potential(const float &r, const int &j, const t_topology &top, const t_trxframe &fr);

void calculate_field(const rvec &point, const int &i, const t_topology &top, const t_trxframe &fr, double *field, double alpha );
void calculate_field(const double *point, const int &i, const t_topology &top, const t_trxframe &fr, double *field, double alpha );
void calculate_field(const float *point, const int &i, const t_topology &top, const t_trxframe &fr, double *field, double alpha );

double project_field(const rvec &point, const double *field);
double project_field(const double *point, const double *field);
double project_field(const float *point, const double *field);

#endif