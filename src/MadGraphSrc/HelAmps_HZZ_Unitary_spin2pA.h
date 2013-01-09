//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph 5 v. 1.5.5, 2012-11-18
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#ifndef HelAmps_HZZ_Unitary_spin2pA_H
#define HelAmps_HZZ_Unitary_spin2pA_H

#include <cmath> 
#include <complex> 

using namespace std; 

namespace MG5_HZZ_Unitary_spin2pA
{
void vxxxxx(double p[4], double vmass, int nhel, int nsv, std::complex<double>
    v[6]);

void sxxxxx(double p[4], int nss, std::complex<double> sc[3]); 

void oxxxxx(double p[4], double fmass, int nhel, int nsf, std::complex<double>
    fo[6]);

void ixxxxx(double p[4], double fmass, int nhel, int nsf, std::complex<double>
    fi[6]);

double Sgn(double e, double f); 

void VVT13_1(complex<double> V2[], complex<double> T3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[]);

void FFV2_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex);
void FFV2_4_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP1, complex<double> COUP2, complex<double> & vertex);

void VVT13_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void VVT3_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void FFV1_2(complex<double> F1[], complex<double> V3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[]);

void VVT5_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void FFV2_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[]);
void FFV2_4_3(complex<double> F1[], complex<double> F2[], complex<double>
    COUP1, complex<double> COUP2, double M3, double W3, complex<double> V3[]);

void VVT3_1(complex<double> V2[], complex<double> T3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[]);

void VVT10_1(complex<double> V2[], complex<double> T3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[]);
void VVT10_11_12_13_2_3_6_7_8_9_1(complex<double> V2[], complex<double> T3[],
    complex<double> COUP1, complex<double> COUP2, complex<double> COUP3,
    complex<double> COUP4, complex<double> COUP5, complex<double> COUP6,
    complex<double> COUP7, complex<double> COUP8, complex<double> COUP9,
    complex<double> COUP10, double M1, double W1, complex<double> V1[]);

void FFV4_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex);

void VVT6_1(complex<double> V2[], complex<double> T3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[]);

void VVT10_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void VVT9_1(complex<double> V2[], complex<double> T3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[]);

void VVT12_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void VVT9_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void VVT4_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void VVT12_1(complex<double> V2[], complex<double> T3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[]);

void VVT1_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);
void VVT1_10_11_12_13_3_5_7_8_9_3(complex<double> V1[], complex<double> V2[],
    complex<double> COUP1, complex<double> COUP2, complex<double> COUP3,
    complex<double> COUP4, complex<double> COUP5, complex<double> COUP6,
    complex<double> COUP7, complex<double> COUP8, complex<double> COUP9,
    complex<double> COUP10, double M3, double W3, complex<double> T3[]);

void VVT11_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void VVT7_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void VVT11_1(complex<double> V2[], complex<double> T3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[]);

void VVT2_1(complex<double> V2[], complex<double> T3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[]);

void FFV1_1(complex<double> F2[], complex<double> V3[], complex<double> COUP,
    double M1, double W1, complex<double> F1[]);

void VVT8_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void FFV4_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[]);

void VVT8_1(complex<double> V2[], complex<double> T3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[]);

void VVT7_1(complex<double> V2[], complex<double> T3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[]);

}  // end namespace MG5_HZZ_Unitar

#endif  // HelAmps_HZZ_Unitary_H
