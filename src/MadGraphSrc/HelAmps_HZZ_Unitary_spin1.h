//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph 5 v. 1.5.5, 2012-11-18
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#ifndef HelAmps_HZZ_Unitary_H
#define HelAmps_HZZ_Unitary_H

#include <cmath> 
#include <complex> 

using namespace std; 

namespace MG5_HZZ_Unitary 
{
void vxxxxx(double p[4], double vmass, int nhel, int nsv, std::complex<double>
    v[6]);

void sxxxxx(double p[4], int nss, std::complex<double> sc[3]); 

void oxxxxx(double p[4], double fmass, int nhel, int nsf, std::complex<double>
    fo[6]);

void ixxxxx(double p[4], double fmass, int nhel, int nsf, std::complex<double>
    fi[6]);

double Sgn(double e, double f); 

void VVV3_0(complex<double> V1[], complex<double> V2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex);

void VVV3_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[]);

void FFV4_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[]);

void VVV2_0(complex<double> V1[], complex<double> V2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex);
void VVV2_3_0(complex<double> V1[], complex<double> V2[], complex<double> V3[],
    complex<double> COUP1, complex<double> COUP2, complex<double> & vertex);

void FFV2_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[]);
void FFV2_4_3(complex<double> F1[], complex<double> F2[], complex<double>
    COUP1, complex<double> COUP2, double M3, double W3, complex<double> V3[]);

void VVV1_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[]);
void VVV1_3_3(complex<double> V1[], complex<double> V2[], complex<double>
    COUP1, complex<double> COUP2, double M3, double W3, complex<double> V3[]);

}  // end namespace MG5_HZZ_Unitar

#endif  // HelAmps_HZZ_Unitary_H
