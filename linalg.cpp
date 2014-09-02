#include "linalg.hpp"

void cross( const rvec &a, const rvec &b, double *ab) {
    ab[0] = a[1] * b[2] - a[2] * b[1];
    ab[1] = a[2] * b[0] - a[0] * b[2];
    ab[2] = a[0] * b[1] - a[1] * b[0];
    return;
}
void cross( const double *a, const double *b, double *ab) {
    ab[0] = a[1] * b[2] - a[2] * b[1];
    ab[1] = a[2] * b[0] - a[0] * b[2];
    ab[2] = a[0] * b[1] - a[1] * b[0];
    return;
}
void cross( const float *a, const float *b, float *ab) {
    ab[0] = a[1] * b[2] - a[2] * b[1];
    ab[1] = a[2] * b[0] - a[0] * b[2];
    ab[2] = a[0] * b[1] - a[1] * b[0];
    return;
}

double dot( const rvec &a, const rvec &b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
double dot( const double *a, const double *b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
double dot( const float *a, const float *b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

double calculate_r(const int &i, const int &j, const t_topology &top, const t_trxframe &fr) {
    double xyz[3];
    double r2, r;
    for (int k=0; k<3; k++){
        xyz[k] = 10 * (fr.x[i][k] - fr.x[j][k]);
    }
    r2 = dot(xyz,xyz);
    r = pow(r2, 0.5);
    return r;
}
double calculate_r(const double *point, const int &j, const t_topology &top, const t_trxframe &fr) {
    double xyz[3];
    double r2, r;
    for (int k=0; k<3; k++){
        xyz[k] = 10 * (point[k] - fr.x[j][k]);
    }
    r2 = dot(xyz,xyz);
    r = pow(r2, 0.5);
    return r;
}
double calculate_r(const float *point, const int &j, const t_topology &top, const t_trxframe &fr){
    double xyz[3];
    double r2, r;
    for (int k=0; k<3; k++){
        xyz[k] = 10 * (point[k] - fr.x[j][k]);
    }
    r2 = dot(xyz,xyz);
    r = pow(r2, 0.5);
    return r;
}

void ring_points( const rvec point, const double *xvec, const double *yvec, const double *zvec, const float &ring_dist, std::vector<std::vector<double> > &each_point, int &j )
{
    // plus x
    for (int i=0; i<3; i++){
        each_point[j][i] = point[i] + ring_dist * xvec[i];
    }
    j++;
    // minus x
    for (int i=0; i<3; i++){
        each_point[j][i] = point[i] - ring_dist * xvec[i];
    }
    j++;
    // plus y
    for (int i=0; i<3; i++){
        each_point[j][i] = point[i] + ring_dist * yvec[i];
    }
    j++;
    // minus y
    for (int i=0; i<3; i++){
        each_point[j][i] = point[i] - ring_dist * yvec[i];
    }
    j++;
    double pos[3];
    // normalized[ (plus x) + (plus y) ] -> j-4, j-2
    for (int i=0; i<3; i++)
    {
        each_point[j][i] = point[i] + ring_dist * (xvec[i] + yvec[i]) / sqrt(2);
    }
    j++;
    // normalized[ (plus x) + (minus y) ] -> j-4, j-1
    for (int i=0; i<3; i++)
    {
        each_point[j][i] = point[i] + ring_dist * (xvec[i] - yvec[i]) / sqrt(2);
    }
    j++;
    // normalized[ (minus x) + (plus y) ] -> j-3, j-2
    for (int i=0; i<3; i++)
    {
        each_point[j][i] = point[i] - ring_dist * (xvec[i] + yvec[i]) / sqrt(2);
    }
    j++;
    // normalized[ (minus x) + (minus y) ] -> j-3, j-1
    for (int i=0; i<3; i++)
    {
        each_point[j][i] = point[i] - ring_dist * (xvec[i] - yvec[i]) / sqrt(2);
    }
    j++;
    return;
}

void write_dummy_atom( FILE *pqr, const double* coord, const std::string &atomname ) {
    fprintf(pqr, "ATOM %6i %4s %4s %5i %11.3f %7.3f %7.3f %7.4f %7.3f\n", 0, atomname.c_str(), "DUM", 0, coord[0]*10, coord[1]*10, coord[2]*10, 0.0, 0.0);
    return;
}
