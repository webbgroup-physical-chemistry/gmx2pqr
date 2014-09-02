#ifndef gmx2pqr_linalg_hpp
#define gmx2pqr_linalg_hpp
#include "my_structs.hpp"

#ifdef __cplusplus
extern "C"
#endif
#include <gromacs/copyrite.h>
#include <gromacs/gmx_fatal.h>
#include <gromacs/vec.h>

#endif

void cross( const rvec &a, const rvec &b, double *ab);
void cross( const double *a, const double *b, double *ab);
void cross( const float *a, const float *b, float *ab);

double dot( const rvec &a, const rvec &b);
double dot( const double *a, const double *b);
double dot( const float *a, const float *b);

double calculate_r(const int &i, const int &j, const t_topology &top, const t_trxframe &fr);
double calculate_r(const double *point, const int &j, const t_topology &top, const t_trxframe &fr);
double calculate_r(const float *point, const int &j, const t_topology &top, const t_trxframe &fr);

void ring_points( const rvec point, const double *xvec, const double *yvec, const double *zvec, const float &ring_dist, std::vector<std::vector<double> > &each_point, int &j );

void write_dummy_atom( FILE *pqr, const double *coord, const std::string &atomname );