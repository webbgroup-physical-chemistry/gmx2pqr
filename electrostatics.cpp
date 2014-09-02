#include "electrostatics.hpp"

double calculate_potential(const double &r, const int &j, const t_topology &top, const t_trxframe &fr) {
    return top.atoms.atom[j].q * cfac / r;
}
double calculate_potential(const float &r, const int &j, const t_topology &top, const t_trxframe &fr) {
    return top.atoms.atom[j].q * cfac / r;
}

void calculate_field(const rvec &point, const int &i, const t_topology &top, const t_trxframe &fr, double *field, double alpha = 1.0) {
    double xyz[3];
    double r, r2, rr3;
    // don't forget to convert from nm to Angstroms
    for (int j=0;j<3;j++)
    {
        xyz[j] = 10 * (point[j] - fr.x[i][j]);
    }
    r2 = dot(xyz,xyz);
    r = pow(r2, 0.5);
    rr3 = 1.0 / (r * r * r);
    for (int j=0; j<3; j++)
    {
        field[j] += alpha * rr3 * top.atoms.atom[i].q * xyz[j] * cfac;
    }
    return;
}
void calculate_field(const double *point, const int &i, const t_topology &top, const t_trxframe &fr, double *field, double alpha = 1.0) {
    double xyz[3];
    double r, r2, rr3;
    // don't forget to convert from nm to Angstroms
    for (int j=0;j<3;j++)
    {
        xyz[j] = 10 * (point[j] - fr.x[i][j]);
    }
    r2 = dot(xyz,xyz);
    r = pow(r2, 0.5);
    rr3 = 1.0 / (r * r * r);
    for (int j=0; j<3; j++)
    {
        field[j] += alpha * rr3 * top.atoms.atom[i].q * xyz[j] * cfac;
    }
    return;
}
void calculate_field(const float *point, const int &i, const t_topology &top, const t_trxframe &fr, double *field, double alpha = 1.0) {
    double xyz[3];
    double r, r2, rr3;
    // don't forget to convert from nm to Angstroms
    for (int j=0;j<3;j++)
    {
        xyz[j] = 10 * (point[j] - fr.x[i][j]);
    }
    r2 = dot(xyz,xyz);
    r = pow(r2, 0.5);
    rr3 = 1.0 / (r * r * r);
    for (int j=0; j<3; j++)
    {
        field[j] += alpha * rr3 * top.atoms.atom[i].q * xyz[j] * cfac;
    }
    return;
}

double project_field(const rvec &vector, const double *field) {
    double projection = 0;
    double vlen = pow(iprod(vector,vector),.5);
    for (int i=0; i<3; i++)
    {
        projection += field[i] * vector[i] / vlen;
    }
    if (! std::isnan(projection)) {
        return projection;
    }
    else
    {
        std::cerr << "\nERROR: NAN is a result for a projection, there's a bug somewhere!\n";
        //std::cout << vector[0] << " " << vector[1] << " " << vector[2] << std::endl;
        //std::exit(1);
        return 0;
    }
}
double project_field(const double *vector, const double *field) {
    double projection = 0;
    double vlen = pow(dot(vector,vector),.5);
    for (int i=0; i<3; i++)
    {
        projection += field[i] * vector[i] / vlen;
    }
    if (! std::isnan(projection)) {
        return projection;
    }
    else
    {
        std::cerr << "\nERROR: NAN is a result for a projection, there's a bug somewhere!\n";
        std::cerr << vector[0] << " " << vector[1] << " " << vector[2] << std::endl;
        std::cerr << field[0] << " " << field[1] << " " << field[2] << std::endl;

        std::exit(1);
        return 0;
    }
}
double project_field(const float *vector, const double *field) {
    double projection = 0;
    double vlen = pow(dot(vector,vector),.5);
    for (int i=0; i<3; i++)
    {
        projection += field[i] * vector[i] / vlen;
    }
    if (! std::isnan(projection)) {
        return projection;
    }
    else
    {
        std::cerr << "\nERROR: NAN is a result for a projection, there's a bug somewhere!\n";
        //std::exit(1);
        return 0;
    }
}