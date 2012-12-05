//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph 5 v. 1.5.3, 2012-11-01
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#ifndef HelAmps_HEF_UFO_spin2_H
#define HelAmps_HEF_UFO_spin2_H

#include <cmath> 
#include <complex> 

using namespace std; 

namespace MG5_HEF_UFO_spin2
{
double Sgn(double e, double f); 

void ixxxxx(double p[4], double fmass, int nhel, int nsf, std::complex<double>
    fi[6]);

void oxxxxx(double p[4], double fmass, int nhel, int nsf, std::complex<double>
    fo[6]);

void sxxxxx(double p[4], int nss, std::complex<double> sc[3]); 

void vxxxxx(double p[4], double vmass, int nhel, int nsv, std::complex<double>
    v[6]);

void VVT16_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void FFV42_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[]);
void FFV42_44_3(complex<double> F1[], complex<double> F2[], complex<double>
    COUP1, complex<double> COUP2, double M3, double W3, complex<double> V3[]);

void VVT19_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex);

void VVT18_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex);

void VVT12_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex);

void VVT11_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex);
void VVT11_12_15_16_17_18_19_20_21_22_0(complex<double> V1[], complex<double>
    V2[], complex<double> T3[], complex<double> COUP1, complex<double> COUP2,
    complex<double> COUP3, complex<double> COUP4, complex<double> COUP5,
    complex<double> COUP6, complex<double> COUP7, complex<double> COUP8,
    complex<double> COUP9, complex<double> COUP10, complex<double> & vertex);

void VVT22_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex);

void VVT20_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void VVT18_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void VVT21_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void VVT10_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);
void VVT10_12_14_16_17_18_19_20_21_22_3(complex<double> V1[], complex<double>
    V2[], complex<double> COUP1, complex<double> COUP2, complex<double> COUP3,
    complex<double> COUP4, complex<double> COUP5, complex<double> COUP6,
    complex<double> COUP7, complex<double> COUP8, complex<double> COUP9,
    complex<double> COUP10, double M3, double W3, complex<double> T3[]);

void VVT21_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex);

void VVT12_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void VVT19_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void VVT17_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex);

void VVT22_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void FFV44_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[]);

void VVT16_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex);

void VVT17_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void VVT13_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void VVT20_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex);

void VVT15_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex);

void VVT14_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

}  // end namespace MG5_HEF_UF

#endif  // HelAmps_HEF_UFO_spin2_H
