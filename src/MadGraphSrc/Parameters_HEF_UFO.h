//==========================================================================
// This file has been automatically generated for C++
// MadGraph 5 v. 1.5.2, 2012-10-15
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#ifndef Parameters_HEF_UFO_H
#define Parameters_HEF_UFO_H

#include <complex> 

#include "read_slha.h"
using namespace std;

 
class Parameters_HEF_UFO
{
  public:

    static Parameters_HEF_UFO * getInstance(); 

    // Define "zero"
    double zero, ZERO; 
    // Model parameters independent of aS
    double WXG, WH, WW, WZ, WT, ymtau, ymm, yme, ymt, ymb, ymc, yms, ymup,
        ymdo, aS, Gf, aEWM1, MXG, MH, MZ, MTA, MM, Me, MT, MB, MC, MS, MU, MD,
        g4z, g3z, g2z, g1z, g4g, g3g, g2g, g1g, k10g, k9g, k8g, k7g, k6g, k5g,
        k4g, k3g, k2g, k1g, k10z, k9z, k8z, k7z, k6z, k5z, k4z, k3z, k2z, k1z,
        cabi, gw, g1, cos__cabi, sin__cabi, MZ__exp__2, MZ__exp__4, sqrt__2,
        MH__exp__2, aEW, MW, sqrt__aEW, ee, MW__exp__2, sw2, cw, sqrt__sw2, sw,
        vev, vev__exp__2, lam, yb, yc, ydo, ye, ym, ys, yt, ytau, yup, muH,
        ee__exp__2, sw__exp__2, cw__exp__2;
    std::complex<double> CKM11, CKM12, CKM13, CKM21, CKM22, CKM23, CKM31,
        CKM32, CKM33, conjg__CKM11, conjg__CKM21, conjg__CKM31, conjg__CKM12,
        conjg__CKM22, conjg__CKM32, conjg__CKM13, conjg__CKM23, conjg__CKM33,
        complexi, I1x11, I1x12, I1x13, I1x21, I1x22, I1x23, I1x31, I1x32,
        I1x33, I2x11, I2x12, I2x13, I2x21, I2x22, I2x23, I2x31, I2x32, I2x33,
        I3x11, I3x12, I3x13, I3x21, I3x22, I3x23, I3x31, I3x32, I3x33, I4x11,
        I4x12, I4x13, I4x21, I4x22, I4x23, I4x31, I4x32, I4x33;
    // Model parameters dependent on aS
    double sqrt__aS, G, G__exp__2; 
    // Model couplings independent of aS
    std::complex<double> GC_109, GC_116, GC_60, GC_61, GC_62, GC_65, GC_66,
        GC_69, GC_70, GC_73, GC_74, GC_77, GC_78, GC_79, GC_80, GC_81, GC_82,
        GC_83, GC_84, GC_88, GC_89, GC_90, GC_1, GC_2, GC_3, GC_11, GC_12, GC_13, GC_16, GC_17, GC_20, GC_21, GC_23, GC_110, GC_115;
    // Model couplings dependent on aS


    // Set parameters that are unchanged during the run
    void setIndependentParameters(SLHAReader& slha); 
    // Set couplings that are unchanged during the run
    void setIndependentCouplings(); 
    // Set parameters that are changed event by event
    void setDependentParameters(); 
    // Set couplings that are changed event by event
    void setDependentCouplings(); 


  private:
    static Parameters_HEF_UFO * instance; 
}; 

#endif  // Parameters_HEF_UFO_H

