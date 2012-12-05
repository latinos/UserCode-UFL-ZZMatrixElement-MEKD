//==========================================================================
// This file has been automatically generated for C++ by
// MadGraph 5 v. 1.5.2, 2012-10-15
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#ifndef Parameters_HEF_UFO_CC
#define Parameters_HEF_UFO_CC

#include <iostream> 
#include <iomanip> 
#include "Parameters_HEF_UFO.h"

// Initialize static instance
Parameters_HEF_UFO * Parameters_HEF_UFO::instance = 0; 

// Function to get static instance - only one instance per program
Parameters_HEF_UFO * Parameters_HEF_UFO::getInstance()
{
  if (instance == 0)
    instance = new Parameters_HEF_UFO(); 

  return instance; 
}

void Parameters_HEF_UFO::setIndependentParameters(SLHAReader& slha)
{
  // Define "zero"
  zero = 0; 
  ZERO = 0; 
  // Prepare a vector for indices
  vector<int> indices(2, 0); 
  WXG = slha.get_block_entry("decay", 9000007, 5.753088e-03); 
  WH = slha.get_block_entry("decay", 9000006, 5.753088e-03); 
  WW = slha.get_block_entry("decay", 24, 2.085000e+00); 
  WZ = slha.get_block_entry("decay", 23, 2.495200e+00); 
  WT = slha.get_block_entry("decay", 6, 1.508336e+00); 
  ymtau = slha.get_block_entry("yukawa", 15, 1.777000e+00); 
  ymm = slha.get_block_entry("yukawa", 13, 1.056600e-01); 
  yme = slha.get_block_entry("yukawa", 11, 5.110000e-04); 
  ymt = slha.get_block_entry("yukawa", 6, 1.720000e+02); 
  ymb = slha.get_block_entry("yukawa", 5, 4.700000e+00); 
  ymc = slha.get_block_entry("yukawa", 4, 1.270000e+00); 
  yms = slha.get_block_entry("yukawa", 3, 1.010000e-01); 
  ymup = slha.get_block_entry("yukawa", 2, 2.550000e-03); 
  ymdo = slha.get_block_entry("yukawa", 1, 5.040000e-03); 
  aS = slha.get_block_entry("sminputs", 3, 1.184000e-01); 
  Gf = slha.get_block_entry("sminputs", 2, 1.166370e-05); 
  aEWM1 = slha.get_block_entry("sminputs", 1, 1.279000e+02); 
  MXG = slha.get_block_entry("mass", 9000007, 1.250000e+02); 
  MH = slha.get_block_entry("mass", 9000006, 1.250000e+02); 
  MZ = slha.get_block_entry("mass", 23, 9.118760e+01); 
  MTA = slha.get_block_entry("mass", 15, 1.777000e+00); 
  MM = slha.get_block_entry("mass", 13, 1.056600e-01); 
  Me = slha.get_block_entry("mass", 11, 5.110000e-04); 
  MT = slha.get_block_entry("mass", 6, 1.720000e+02); 
  MB = slha.get_block_entry("mass", 5, 4.700000e+00); 
  MC = slha.get_block_entry("mass", 4, 1.270000e+00); 
  MS = slha.get_block_entry("mass", 3, 1.010000e-01); 
  MU = slha.get_block_entry("mass", 2, 2.550000e-03); 
  MD = slha.get_block_entry("mass", 1, 5.040000e-03); 
  g4z = slha.get_block_entry("heff", 8, 1.000000e-01); 
  g3z = slha.get_block_entry("heff", 7, 1.000000e-01); 
  g2z = slha.get_block_entry("heff", 6, 1.000000e-01); 
  g1z = slha.get_block_entry("heff", 5, 1.000000e-01); 
  g4g = slha.get_block_entry("heff", 4, 1.000000e-01); 
  g3g = slha.get_block_entry("heff", 3, 1.000000e-01); 
  g2g = slha.get_block_entry("heff", 2, 1.000000e-01); 
  g1g = slha.get_block_entry("heff", 1, 1.000000e-01); 
  k10g = slha.get_block_entry("gravity", 10, 1.000000e-01); 
  k9g = slha.get_block_entry("gravity", 9, 1.000000e-01); 
  k8g = slha.get_block_entry("gravity", 8, 1.000000e-01); 
  k7g = slha.get_block_entry("gravity", 7, 1.000000e-01); 
  k6g = slha.get_block_entry("gravity", 6, 1.000000e-01); 
  k5g = slha.get_block_entry("gravity", 5, 1.000000e-01); 
  k4g = slha.get_block_entry("gravity", 4, 1.000000e-01); 
  k3g = slha.get_block_entry("gravity", 3, 1.000000e-01); 
  k2g = slha.get_block_entry("gravity", 2, 1.000000e-01); 
  k1g = slha.get_block_entry("gravity", 1, 1.000000e-01); 
  k10z = slha.get_block_entry("gravity", 20, 1.000000e-01); 
  k9z = slha.get_block_entry("gravity", 19, 1.000000e-01); 
  k8z = slha.get_block_entry("gravity", 18, 1.000000e-01); 
  k7z = slha.get_block_entry("gravity", 17, 1.000000e-01); 
  k6z = slha.get_block_entry("gravity", 16, 1.000000e-01); 
  k5z = slha.get_block_entry("gravity", 15, 1.000000e-01); 
  k4z = slha.get_block_entry("gravity", 14, 1.000000e-01); 
  k3z = slha.get_block_entry("gravity", 13, 1.000000e-01); 
  k2z = slha.get_block_entry("gravity", 12, 1.000000e-01); 
  k1z = slha.get_block_entry("gravity", 11, 1.000000e-01); 
  cabi = slha.get_block_entry("ckmblock", 1, 2.277360e-01); 
  gw = 1.; 
  g1 = 1.; 
  cos__cabi = cos(cabi); 
  CKM11 = cos__cabi; 
  sin__cabi = sin(cabi); 
  CKM12 = sin__cabi; 
  CKM13 = 0.; 
  CKM21 = -sin__cabi; 
  CKM22 = cos__cabi; 
  CKM23 = 0.; 
  CKM31 = 0.; 
  CKM32 = 0.; 
  CKM33 = 1.; 
  MZ__exp__2 = pow(MZ, 2.); 
  MZ__exp__4 = pow(MZ, 4.); 
  sqrt__2 = sqrt(2.); 
  MH__exp__2 = pow(MH, 2.); 
  conjg__CKM11 = conj(CKM11); 
  conjg__CKM21 = conj(CKM21); 
  conjg__CKM31 = conj(CKM31); 
  conjg__CKM12 = conj(CKM12); 
  conjg__CKM22 = conj(CKM22); 
  conjg__CKM32 = conj(CKM32); 
  conjg__CKM13 = conj(CKM13); 
  conjg__CKM23 = conj(CKM23); 
  conjg__CKM33 = conj(CKM33); 
  complexi = std::complex<double> (0., 1.); 
  aEW = 1./aEWM1; 
  MW = sqrt(MZ__exp__2/2. + sqrt(MZ__exp__4/4. - (aEW * M_PI * MZ__exp__2)/(Gf
      * sqrt__2)));
  sqrt__aEW = sqrt(aEW); 
  ee = 2. * sqrt__aEW * sqrt(M_PI); 
  MW__exp__2 = pow(MW, 2.); 
  sw2 = 1. - MW__exp__2/MZ__exp__2; 
  cw = sqrt(1. - sw2); 
  sqrt__sw2 = sqrt(sw2); 
  sw = sqrt__sw2; 
  vev = (2. * MW * sw)/ee; 
  vev__exp__2 = pow(vev, 2.); 
  lam = MH__exp__2/(2. * vev__exp__2); 
  yb = (ymb * sqrt__2)/vev; 
  yc = (ymc * sqrt__2)/vev; 
  ydo = (ymdo * sqrt__2)/vev; 
  ye = (yme * sqrt__2)/vev; 
  ym = (ymm * sqrt__2)/vev; 
  ys = (yms * sqrt__2)/vev; 
  yt = (ymt * sqrt__2)/vev; 
  ytau = (ymtau * sqrt__2)/vev; 
  yup = (ymup * sqrt__2)/vev; 
  muH = sqrt(lam * vev__exp__2); 
  I1x11 = ydo * conjg__CKM11; 
  I1x12 = ydo * conjg__CKM21; 
  I1x13 = ydo * conjg__CKM31; 
  I1x21 = ys * conjg__CKM12; 
  I1x22 = ys * conjg__CKM22; 
  I1x23 = ys * conjg__CKM32; 
  I1x31 = yb * conjg__CKM13; 
  I1x32 = yb * conjg__CKM23; 
  I1x33 = yb * conjg__CKM33; 
  I2x11 = yup * conjg__CKM11; 
  I2x12 = yc * conjg__CKM21; 
  I2x13 = yt * conjg__CKM31; 
  I2x21 = yup * conjg__CKM12; 
  I2x22 = yc * conjg__CKM22; 
  I2x23 = yt * conjg__CKM32; 
  I2x31 = yup * conjg__CKM13; 
  I2x32 = yc * conjg__CKM23; 
  I2x33 = yt * conjg__CKM33; 
  I3x11 = CKM11 * yup; 
  I3x12 = CKM21 * yc; 
  I3x13 = CKM31 * yt; 
  I3x21 = CKM12 * yup; 
  I3x22 = CKM22 * yc; 
  I3x23 = CKM32 * yt; 
  I3x31 = CKM13 * yup; 
  I3x32 = CKM23 * yc; 
  I3x33 = CKM33 * yt; 
  I4x11 = CKM11 * ydo; 
  I4x12 = CKM21 * ydo; 
  I4x13 = CKM31 * ydo; 
  I4x21 = CKM12 * ys; 
  I4x22 = CKM22 * ys; 
  I4x23 = CKM32 * ys; 
  I4x31 = CKM13 * yb; 
  I4x32 = CKM23 * yb; 
  I4x33 = CKM33 * yb; 
  ee__exp__2 = pow(ee, 2.); 
  sw__exp__2 = pow(sw, 2.); 
  cw__exp__2 = pow(cw, 2.); 
}

void Parameters_HEF_UFO::setIndependentCouplings()
{
    GC_1 = -(ee * complexi)/3.;
    GC_2 = (2. * ee * complexi)/3.;
    GC_3 = -(ee * complexi);
    GC_11 = -(complexi * g1g);
    GC_12 = -(complexi * g1z);
    GC_13 = -(complexi * g2g);
    GC_16 = -(complexi * g2z);
    GC_17 = -(complexi * g3g);
    GC_20 = -(complexi * g3z);
    GC_21 = (complexi * g4g)/8.;
    GC_23 = (complexi * g4z)/2.;
    GC_110 = (cw * ee * complexi)/(2. * sw);
    GC_115 = -(ee * complexi * sw)/(6. * cw); 
    
  GC_109 = -(cw * ee * complexi)/(2. * sw); 
  GC_116 = (ee * complexi * sw)/(2. * cw); 
  GC_60 = -(complexi * k10g)/2.; 
  GC_61 = -(complexi * k10z)/2.; 
  GC_62 = -(complexi * k1g); 
  GC_65 = -(complexi * k1z); 
  GC_66 = complexi * k2g; 
  GC_69 = complexi * k2z; 
  GC_70 = complexi * k3g; 
  GC_73 = complexi * k3z; 
  GC_74 = -2. * complexi * k4g; 
  GC_77 = -2. * complexi * k4z; 
  GC_78 = complexi * k5g; 
  GC_79 = complexi * k5z; 
  GC_80 = -(complexi * k6g)/2.; 
  GC_81 = -(complexi * k6z)/2.; 
  GC_82 = -(complexi * k7g); 
  GC_83 = -(complexi * k7z); 
  GC_84 = (complexi * k8g)/4.; 
  GC_88 = complexi * k8z; 
  GC_89 = -(complexi * k9g)/2.; 
  GC_90 = -(complexi * k9z)/2.; 
}

void Parameters_HEF_UFO::setDependentParameters()
{
  sqrt__aS = sqrt(aS); 
  G = 2. * sqrt__aS * sqrt(M_PI); 
  G__exp__2 = pow(G, 2.); 
}

void Parameters_HEF_UFO::setDependentCouplings()
{

}



#endif
