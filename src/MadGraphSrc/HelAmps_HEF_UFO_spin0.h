//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph 5 v. 1.5.3, 2012-11-01
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#ifndef HelAmps_HEF_UFO_spin0_H
#define HelAmps_HEF_UFO_spin0_H

#include <cmath> 
#include <complex> 

using namespace std; 

namespace MG5_HEF_UFO_spin0
{
double Sgn(double e, double f); 

void ixxxxx(double p[4], double fmass, int nhel, int nsf, std::complex<double>
    fi[6]);

void oxxxxx(double p[4], double fmass, int nhel, int nsf, std::complex<double>
    fo[6]);

void sxxxxx(double p[4], int nss, std::complex<double> sc[3]); 

void vxxxxx(double p[4], double vmass, int nhel, int nsv, std::complex<double>
    v[6]);

void VVS11_0(complex<double> V1[], complex<double> V2[], complex<double> S3[],
    complex<double> COUP, complex<double> & vertex);
void VVS11_12_13_0(complex<double> V1[], complex<double> V2[], complex<double>
    S3[], complex<double> COUP1, complex<double> COUP2, complex<double> COUP3,
    complex<double> & vertex);

void VVS12_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> S3[]);

void VVS12_0(complex<double> V1[], complex<double> V2[], complex<double> S3[],
    complex<double> COUP, complex<double> & vertex);

void FFV44_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[]);

void VVS9_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> S3[]);

void FFV42_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[]);
void FFV42_44_3(complex<double> F1[], complex<double> F2[], complex<double>
    COUP1, complex<double> COUP2, double M3, double W3, complex<double> V3[]);

void VVS13_0(complex<double> V1[], complex<double> V2[], complex<double> S3[],
    complex<double> COUP, complex<double> & vertex);

void VVS13_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> S3[]);

void VVS11_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> S3[]);
void VVS11_12_13_3(complex<double> V1[], complex<double> V2[], complex<double>
    COUP1, complex<double> COUP2, complex<double> COUP3, double M3, double W3,
    complex<double> S3[]);

void VVS10_0(complex<double> V1[], complex<double> V2[], complex<double> S3[],
    complex<double> COUP, complex<double> & vertex);

}  // end namespace MG5_HEF_UF

#endif  // HelAmps_HEF_UFO_spin0_H
