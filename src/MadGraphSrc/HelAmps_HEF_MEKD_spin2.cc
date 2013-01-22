//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph 5 v. 1.5.7, 2013-01-15
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#include "HelAmps_HEF_MEKD_spin2.h"

namespace MG5_HEF_MEKD_spin2 
{


double Sgn(double a, double b)
{
  return (b < 0)? - abs(a):abs(a); 
}

void oxxxxx(double p[4], double fmass, int nhel, int nsf, complex<double> fo[6])
{
  complex<double> chi[2]; 
  double sf[2], sfomeg[2], omega[2], pp, pp3, sqp0p3, sqm[2]; 
  int nh, ip, im; 
  fo[0] = complex<double> (p[0] * nsf, p[3] * nsf); 
  fo[1] = complex<double> (p[1] * nsf, p[2] * nsf); 
  nh = nhel * nsf; 
  if (fmass != 0.000)
  {
    pp = min(p[0], pow(pow(p[1], 2) + pow(p[2], 2) + pow(p[3], 2), 0.5)); 
    if (pp == 0.000)
    {
      sqm[0] = pow(abs(fmass), 0.5); 
      sqm[1] = Sgn(sqm[0], fmass); 
      ip = -((1 - nh)/2) * nhel; 
      im = (1 + nh)/2 * nhel; 
      fo[2] = im * sqm[im]; 
      fo[3] = ip * nsf * sqm[im]; 
      fo[4] = im * nsf * sqm[ - ip]; 
      fo[5] = ip * sqm[ - ip]; 
    }
    else
    {
      pp = min(p[0], pow(pow(p[1], 2) + pow(p[2], 2) + pow(p[3], 2), 0.5)); 
      sf[0] = double(1 + nsf + (1 - nsf) * nh) * 0.5; 
      sf[1] = double(1 + nsf - (1 - nsf) * nh) * 0.5; 
      omega[0] = pow(p[0] + pp, 0.5); 
      omega[1] = fmass/omega[0]; 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      sfomeg[0] = sf[0] * omega[ip]; 
      sfomeg[1] = sf[1] * omega[im]; 
      pp3 = max(pp + p[3], 0.00); 
      chi[0] = complex<double> (pow(pp3 * 0.5/pp, 0.5), 0.00); 
      if (pp3 == 0.00)
      {
        chi[1] = complex<double> (-nh, 0.00); 
      }
      else
      {
        chi[1] = complex<double> (nh * p[1], -p[2])/pow(2.0 * pp * pp3, 0.5); 
      }
      fo[2] = sfomeg[1] * chi[im]; 
      fo[3] = sfomeg[1] * chi[ip]; 
      fo[4] = sfomeg[0] * chi[im]; 
      fo[5] = sfomeg[0] * chi[ip]; 
    }
  }
  else
  {
    if((p[1] == 0.00) and (p[2] == 0.00) and (p[3] < 0.00))
    {
      sqp0p3 = 0.00; 
    }
    else
    {
      sqp0p3 = pow(max(p[0] + p[3], 0.00), 0.5) * nsf; 
    }
    chi[0] = complex<double> (sqp0p3, 0.00); 
    if(sqp0p3 == 0.000)
    {
      chi[1] = complex<double> (-nhel, 0.00) * pow(2.0 * p[0], 0.5); 
    }
    else
    {
      chi[1] = complex<double> (nh * p[1], -p[2])/sqp0p3; 
    }
    if(nh == 1)
    {
      fo[2] = chi[0]; 
      fo[3] = chi[1]; 
      fo[4] = complex<double> (0.00, 0.00); 
      fo[5] = complex<double> (0.00, 0.00); 
    }
    else
    {
      fo[2] = complex<double> (0.00, 0.00); 
      fo[3] = complex<double> (0.00, 0.00); 
      fo[4] = chi[1]; 
      fo[5] = chi[0]; 
    }
  }
  return; 
}

void ixxxxx(double p[4], double fmass, int nhel, int nsf, complex<double> fi[6])
{
  complex<double> chi[2]; 
  double sf[2], sfomega[2], omega[2], pp, pp3, sqp0p3, sqm[2]; 
  int ip, im, nh; 
  fi[0] = complex<double> (-p[0] * nsf, -p[3] * nsf); 
  fi[1] = complex<double> (-p[1] * nsf, -p[2] * nsf); 
  nh = nhel * nsf; 
  if (fmass != 0.0)
  {
    pp = min(p[0], pow((pow(p[1], 2) + pow(p[2], 2) + pow(p[3], 2)), 0.5)); 
    if (pp == 0.0)
    {
      sqm[0] = pow(abs(fmass), 0.5); 
      sqm[1] = Sgn(sqm[0], fmass); 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      fi[2] = ip * sqm[ip]; 
      fi[3] = im * nsf * sqm[ip]; 
      fi[4] = ip * nsf * sqm[im]; 
      fi[5] = im * sqm[im]; 
    }
    else
    {
      sf[0] = (1 + nsf + (1 - nsf) * nh) * 0.5; 
      sf[1] = (1 + nsf - (1 - nsf) * nh) * 0.5; 
      omega[0] = pow(p[0] + pp, 0.5); 
      omega[1] = fmass/omega[0]; 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      sfomega[0] = sf[0] * omega[ip]; 
      sfomega[1] = sf[1] * omega[im]; 
      pp3 = max(pp + p[3], 0.0); 
      chi[0] = complex<double> (pow(pp3 * 0.5/pp, 0.5), 0); 
      if (pp3 == 0.0)
      {
        chi[1] = complex<double> (-nh, 0); 
      }
      else
      {
        chi[1] = complex<double> (nh * p[1], p[2])/pow(2.0 * pp * pp3, 0.5); 
      }
      fi[2] = sfomega[0] * chi[im]; 
      fi[3] = sfomega[0] * chi[ip]; 
      fi[4] = sfomega[1] * chi[im]; 
      fi[5] = sfomega[1] * chi[ip]; 
    }
  }
  else
  {
    if (p[1] == 0.0 and p[2] == 0.0 and p[3] < 0.0)
    {
      sqp0p3 = 0.0; 
    }
    else
    {
      sqp0p3 = pow(max(p[0] + p[3], 0.0), 0.5) * nsf; 
    }
    chi[0] = complex<double> (sqp0p3, 0.0); 
    if (sqp0p3 == 0.0)
    {
      chi[1] = complex<double> (-nhel * pow(2.0 * p[0], 0.5), 0.0); 
    }
    else
    {
      chi[1] = complex<double> (nh * p[1], p[2])/sqp0p3; 
    }
    if (nh == 1)
    {
      fi[2] = complex<double> (0.0, 0.0); 
      fi[3] = complex<double> (0.0, 0.0); 
      fi[4] = chi[0]; 
      fi[5] = chi[1]; 
    }
    else
    {
      fi[2] = chi[1]; 
      fi[3] = chi[0]; 
      fi[4] = complex<double> (0.0, 0.0); 
      fi[5] = complex<double> (0.0, 0.0); 
    }
  }
  return; 
}

void sxxxxx(double p[4], int nss, complex<double> sc[3])
{
  sc[2] = complex<double> (1.00, 0.00); 
  sc[0] = complex<double> (p[0] * nss, p[3] * nss); 
  sc[1] = complex<double> (p[1] * nss, p[2] * nss); 
  return; 
}

void vxxxxx(double p[4], double vmass, int nhel, int nsv, complex<double> vc[6])
{
  double hel, hel0, pt, pt2, pp, pzpt, emp, sqh; 
  int nsvahl; 
  sqh = pow(0.5, 0.5); 
  hel = double(nhel); 
  nsvahl = nsv * abs(hel); 
  pt2 = pow(p[1], 2) + pow(p[2], 2); 
  pp = min(p[0], pow(pt2 + pow(p[3], 2), 0.5)); 
  pt = min(pp, pow(pt2, 0.5)); 
  vc[0] = complex<double> (p[0] * nsv, p[3] * nsv); 
  vc[1] = complex<double> (p[1] * nsv, p[2] * nsv); 
  if (vmass != 0.0)
  {
    hel0 = 1.0 - abs(hel); 
    if(pp == 0.0)
    {
      vc[2] = complex<double> (0.0, 0.0); 
      vc[3] = complex<double> (-hel * sqh, 0.0); 
      vc[4] = complex<double> (0.0, nsvahl * sqh); 
      vc[5] = complex<double> (hel0, 0.0); 
    }
    else
    {
      emp = p[0]/(vmass * pp); 
      vc[2] = complex<double> (hel0 * pp/vmass, 0.0); 
      vc[5] = complex<double> (hel0 * p[3] * emp + hel * pt/pp * sqh, 0.0); 
      if (pt != 0.0)
      {
        pzpt = p[3]/(pp * pt) * sqh * hel; 
        vc[3] = complex<double> (hel0 * p[1] * emp - p[1] * pzpt, -nsvahl *
            p[2]/pt * sqh);
        vc[4] = complex<double> (hel0 * p[2] * emp - p[2] * pzpt, nsvahl *
            p[1]/pt * sqh);
      }
      else
      {
        vc[3] = complex<double> (-hel * sqh, 0.0); 
        vc[4] = complex<double> (0.0, nsvahl * Sgn(sqh, p[3])); 
      }
    }
  }
  else
  {
    pp = p[0]; 
    pt = pow(pow(p[1], 2) + pow(p[2], 2), 0.5); 
    vc[2] = complex<double> (0.0, 0.0); 
    vc[5] = complex<double> (hel * pt/pp * sqh, 0.0); 
    if (pt != 0.0)
    {
      pzpt = p[3]/(pp * pt) * sqh * hel; 
      vc[3] = complex<double> (-p[1] * pzpt, -nsv * p[2]/pt * sqh); 
      vc[4] = complex<double> (-p[2] * pzpt, nsv * p[1]/pt * sqh); 
    }
    else
    {
      vc[3] = complex<double> (-hel * sqh, 0.0); 
      vc[4] = complex<double> (0.0, nsv * Sgn(sqh, p[3])); 
    }
  }
  return; 
}

void VVT12_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP2; 
  complex<double> TMP1; 
  double P1[4]; 
  complex<double> TMP10; 
  double P2[4]; 
  complex<double> TMP23; 
  complex<double> TMP7; 
  double P3[4]; 
  double OM3; 
  complex<double> denom; 
  complex<double> TMP26; 
  complex<double> TMP4; 
  complex<double> TMP9; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +V1[0] + V2[0]; 
  T3[1] = +V1[1] + V2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP26 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP23 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP9 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP4 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP7 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP1 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP10 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP2 = (P3[0] * P2[0] - P3[1] * P2[1] - P3[2] * P2[2] - P3[3] * P2[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * 2. * (OM3 * (P3[0] * (P3[0] * (OM3 * (TMP1 * 0.666666667 *
      (-cI * (TMP9 * TMP23) + cI * (TMP2 * TMP4)) + 0.666666667 * (TMP26 * (-cI
      * (TMP2 * TMP7) + cI * (TMP10 * TMP23)))) + (-0.666666667 * cI * (TMP7 *
      TMP9) + 0.666666667 * cI * (TMP4 * TMP10))) + (P1[0] * (-cI * (TMP2 *
      TMP4) + cI * (TMP9 * TMP23)) + (P2[0] * (-cI * (TMP1 * TMP4) + cI * (TMP7
      * TMP26)) + (TMP10 * - 1. * (+cI * (V2[2] * TMP26 + V1[2] * TMP23)) +
      (+cI * (TMP2 * V1[2] * TMP7 + TMP1 * V2[2] * TMP9)))))) + (TMP1 *
      0.333333333 * (-cI * (TMP9 * TMP23) + cI * (TMP2 * TMP4)) + 0.333333333 *
      (TMP26 * (-cI * (TMP2 * TMP7) + cI * (TMP10 * TMP23))))) + (TMP10 *
      (-0.666666667 * cI * (TMP4) + cI * (V2[2] * V1[2])) + (TMP7 * (-cI *
      (P2[0] * V1[2]) + 0.666666667 * cI * (TMP9)) + P1[0] * (-cI * (V2[2] *
      TMP9) + cI * (P2[0] * TMP4)))));
  T3[3] = denom * (OM3 * (P3[0] * (P3[1] * (OM3 * (TMP1 * 1.333333333 * (-cI *
      (TMP9 * TMP23) + cI * (TMP2 * TMP4)) + 1.333333333 * (TMP26 * (-cI *
      (TMP2 * TMP7) + cI * (TMP10 * TMP23)))) + (-1.333333333 * cI * (TMP7 *
      TMP9) + 1.333333333 * cI * (TMP4 * TMP10))) + (P1[1] * (-cI * (TMP2 *
      TMP4) + cI * (TMP9 * TMP23)) + (P2[1] * (-cI * (TMP1 * TMP4) + cI * (TMP7
      * TMP26)) + (TMP10 * - 1. * (+cI * (V1[3] * TMP23 + V2[3] * TMP26)) +
      (+cI * (TMP2 * V1[3] * TMP7 + TMP1 * V2[3] * TMP9)))))) + P3[1] * (P1[0]
      * (-cI * (TMP2 * TMP4) + cI * (TMP9 * TMP23)) + (P2[0] * (-cI * (TMP1 *
      TMP4) + cI * (TMP7 * TMP26)) + (TMP10 * - 1. * (+cI * (V2[2] * TMP26 +
      V1[2] * TMP23)) + (+cI * (TMP2 * V1[2] * TMP7 + TMP1 * V2[2] * TMP9))))))
      + (P1[0] * (-cI * (V2[3] * TMP9) + cI * (P2[1] * TMP4)) + (P1[1] * (-cI *
      (V2[2] * TMP9) + cI * (P2[0] * TMP4)) + (TMP10 * (+cI * (V2[2] * V1[3] +
      V2[3] * V1[2])) - TMP7 * (+cI * (P2[1] * V1[2] + P2[0] * V1[3]))))));
  T3[4] = denom * (OM3 * (P3[0] * (P3[2] * (OM3 * (TMP1 * 1.333333333 * (-cI *
      (TMP9 * TMP23) + cI * (TMP2 * TMP4)) + 1.333333333 * (TMP26 * (-cI *
      (TMP2 * TMP7) + cI * (TMP10 * TMP23)))) + (-1.333333333 * cI * (TMP7 *
      TMP9) + 1.333333333 * cI * (TMP4 * TMP10))) + (P1[2] * (-cI * (TMP2 *
      TMP4) + cI * (TMP9 * TMP23)) + (P2[2] * (-cI * (TMP1 * TMP4) + cI * (TMP7
      * TMP26)) + (TMP10 * - 1. * (+cI * (V1[4] * TMP23 + V2[4] * TMP26)) +
      (+cI * (TMP2 * V1[4] * TMP7 + TMP1 * V2[4] * TMP9)))))) + P3[2] * (P1[0]
      * (-cI * (TMP2 * TMP4) + cI * (TMP9 * TMP23)) + (P2[0] * (-cI * (TMP1 *
      TMP4) + cI * (TMP7 * TMP26)) + (TMP10 * - 1. * (+cI * (V2[2] * TMP26 +
      V1[2] * TMP23)) + (+cI * (TMP2 * V1[2] * TMP7 + TMP1 * V2[2] * TMP9))))))
      + (P1[0] * (-cI * (V2[4] * TMP9) + cI * (P2[2] * TMP4)) + (P1[2] * (-cI *
      (V2[2] * TMP9) + cI * (P2[0] * TMP4)) + (TMP10 * (+cI * (V2[2] * V1[4] +
      V2[4] * V1[2])) - TMP7 * (+cI * (P2[2] * V1[2] + P2[0] * V1[4]))))));
  T3[5] = denom * (OM3 * (P3[0] * (P3[3] * (OM3 * (TMP1 * 1.333333333 * (-cI *
      (TMP9 * TMP23) + cI * (TMP2 * TMP4)) + 1.333333333 * (TMP26 * (-cI *
      (TMP2 * TMP7) + cI * (TMP10 * TMP23)))) + (-1.333333333 * cI * (TMP7 *
      TMP9) + 1.333333333 * cI * (TMP4 * TMP10))) + (P1[3] * (-cI * (TMP2 *
      TMP4) + cI * (TMP9 * TMP23)) + (P2[3] * (-cI * (TMP1 * TMP4) + cI * (TMP7
      * TMP26)) + (TMP10 * - 1. * (+cI * (V1[5] * TMP23 + V2[5] * TMP26)) +
      (+cI * (TMP2 * V1[5] * TMP7 + TMP1 * V2[5] * TMP9)))))) + P3[3] * (P1[0]
      * (-cI * (TMP2 * TMP4) + cI * (TMP9 * TMP23)) + (P2[0] * (-cI * (TMP1 *
      TMP4) + cI * (TMP7 * TMP26)) + (TMP10 * - 1. * (+cI * (V2[2] * TMP26 +
      V1[2] * TMP23)) + (+cI * (TMP2 * V1[2] * TMP7 + TMP1 * V2[2] * TMP9))))))
      + (P1[0] * (-cI * (V2[5] * TMP9) + cI * (P2[3] * TMP4)) + (P1[3] * (-cI *
      (V2[2] * TMP9) + cI * (P2[0] * TMP4)) + (TMP10 * (+cI * (V2[2] * V1[5] +
      V2[5] * V1[2])) - TMP7 * (+cI * (P2[3] * V1[2] + P2[0] * V1[5]))))));
  T3[6] = denom * (OM3 * (P3[0] * (P3[1] * (OM3 * (TMP1 * 1.333333333 * (-cI *
      (TMP9 * TMP23) + cI * (TMP2 * TMP4)) + 1.333333333 * (TMP26 * (-cI *
      (TMP2 * TMP7) + cI * (TMP10 * TMP23)))) + (-1.333333333 * cI * (TMP7 *
      TMP9) + 1.333333333 * cI * (TMP4 * TMP10))) + (P1[1] * (-cI * (TMP2 *
      TMP4) + cI * (TMP9 * TMP23)) + (P2[1] * (-cI * (TMP1 * TMP4) + cI * (TMP7
      * TMP26)) + (TMP10 * - 1. * (+cI * (V2[3] * TMP26 + V1[3] * TMP23)) +
      (+cI * (TMP2 * V1[3] * TMP7 + TMP1 * V2[3] * TMP9)))))) + P3[1] * (P1[0]
      * (-cI * (TMP2 * TMP4) + cI * (TMP9 * TMP23)) + (P2[0] * (-cI * (TMP1 *
      TMP4) + cI * (TMP7 * TMP26)) + (TMP10 * - 1. * (+cI * (V1[2] * TMP23 +
      V2[2] * TMP26)) + (+cI * (TMP2 * V1[2] * TMP7 + TMP1 * V2[2] * TMP9))))))
      + (P1[0] * (-cI * (V2[3] * TMP9) + cI * (P2[1] * TMP4)) + (P1[1] * (-cI *
      (V2[2] * TMP9) + cI * (P2[0] * TMP4)) + (TMP10 * (+cI * (V2[3] * V1[2] +
      V2[2] * V1[3])) - TMP7 * (+cI * (P2[0] * V1[3] + P2[1] * V1[2]))))));
  T3[7] = denom * 2. * (OM3 * (P3[1] * (P3[1] * (OM3 * (TMP1 * 0.666666667 *
      (-cI * (TMP9 * TMP23) + cI * (TMP2 * TMP4)) + 0.666666667 * (TMP26 * (-cI
      * (TMP2 * TMP7) + cI * (TMP10 * TMP23)))) + (-0.666666667 * cI * (TMP7 *
      TMP9) + 0.666666667 * cI * (TMP4 * TMP10))) + (P1[1] * (-cI * (TMP2 *
      TMP4) + cI * (TMP9 * TMP23)) + (P2[1] * (-cI * (TMP1 * TMP4) + cI * (TMP7
      * TMP26)) + (TMP10 * - 1. * (+cI * (V2[3] * TMP26 + V1[3] * TMP23)) +
      (+cI * (TMP2 * V1[3] * TMP7 + TMP1 * V2[3] * TMP9)))))) + (TMP1 *
      0.333333333 * (-cI * (TMP2 * TMP4) + cI * (TMP9 * TMP23)) + 0.333333333 *
      (TMP26 * (-cI * (TMP10 * TMP23) + cI * (TMP2 * TMP7))))) + (TMP10 * (+cI
      * (V2[3] * V1[3]) + 0.666666667 * cI * (TMP4)) + (TMP7 * - 1. * (+cI *
      (P2[1] * V1[3]) + 0.666666667 * cI * (TMP9)) + P1[1] * (-cI * (V2[3] *
      TMP9) + cI * (P2[1] * TMP4)))));
  T3[8] = denom * (OM3 * (P3[1] * (P3[2] * (OM3 * (TMP1 * 1.333333333 * (-cI *
      (TMP9 * TMP23) + cI * (TMP2 * TMP4)) + 1.333333333 * (TMP26 * (-cI *
      (TMP2 * TMP7) + cI * (TMP10 * TMP23)))) + (-1.333333333 * cI * (TMP7 *
      TMP9) + 1.333333333 * cI * (TMP4 * TMP10))) + (P1[2] * (-cI * (TMP2 *
      TMP4) + cI * (TMP9 * TMP23)) + (P2[2] * (-cI * (TMP1 * TMP4) + cI * (TMP7
      * TMP26)) + (TMP10 * - 1. * (+cI * (V1[4] * TMP23 + V2[4] * TMP26)) +
      (+cI * (TMP2 * V1[4] * TMP7 + TMP1 * V2[4] * TMP9)))))) + P3[2] * (P1[1]
      * (-cI * (TMP2 * TMP4) + cI * (TMP9 * TMP23)) + (P2[1] * (-cI * (TMP1 *
      TMP4) + cI * (TMP7 * TMP26)) + (TMP10 * - 1. * (+cI * (V2[3] * TMP26 +
      V1[3] * TMP23)) + (+cI * (TMP2 * V1[3] * TMP7 + TMP1 * V2[3] * TMP9))))))
      + (P1[1] * (-cI * (V2[4] * TMP9) + cI * (P2[2] * TMP4)) + (P1[2] * (-cI *
      (V2[3] * TMP9) + cI * (P2[1] * TMP4)) + (TMP10 * (+cI * (V2[3] * V1[4] +
      V2[4] * V1[3])) - TMP7 * (+cI * (P2[2] * V1[3] + P2[1] * V1[4]))))));
  T3[9] = denom * (OM3 * (P3[1] * (P3[3] * (OM3 * (TMP1 * 1.333333333 * (-cI *
      (TMP9 * TMP23) + cI * (TMP2 * TMP4)) + 1.333333333 * (TMP26 * (-cI *
      (TMP2 * TMP7) + cI * (TMP10 * TMP23)))) + (-1.333333333 * cI * (TMP7 *
      TMP9) + 1.333333333 * cI * (TMP4 * TMP10))) + (P1[3] * (-cI * (TMP2 *
      TMP4) + cI * (TMP9 * TMP23)) + (P2[3] * (-cI * (TMP1 * TMP4) + cI * (TMP7
      * TMP26)) + (TMP10 * - 1. * (+cI * (V1[5] * TMP23 + V2[5] * TMP26)) +
      (+cI * (TMP2 * V1[5] * TMP7 + TMP1 * V2[5] * TMP9)))))) + P3[3] * (P1[1]
      * (-cI * (TMP2 * TMP4) + cI * (TMP9 * TMP23)) + (P2[1] * (-cI * (TMP1 *
      TMP4) + cI * (TMP7 * TMP26)) + (TMP10 * - 1. * (+cI * (V2[3] * TMP26 +
      V1[3] * TMP23)) + (+cI * (TMP2 * V1[3] * TMP7 + TMP1 * V2[3] * TMP9))))))
      + (P1[1] * (-cI * (V2[5] * TMP9) + cI * (P2[3] * TMP4)) + (P1[3] * (-cI *
      (V2[3] * TMP9) + cI * (P2[1] * TMP4)) + (TMP10 * (+cI * (V2[3] * V1[5] +
      V2[5] * V1[3])) - TMP7 * (+cI * (P2[3] * V1[3] + P2[1] * V1[5]))))));
  T3[10] = denom * (OM3 * (P3[0] * (P3[2] * (OM3 * (TMP1 * 1.333333333 * (-cI *
      (TMP9 * TMP23) + cI * (TMP2 * TMP4)) + 1.333333333 * (TMP26 * (-cI *
      (TMP2 * TMP7) + cI * (TMP10 * TMP23)))) + (-1.333333333 * cI * (TMP7 *
      TMP9) + 1.333333333 * cI * (TMP4 * TMP10))) + (P1[2] * (-cI * (TMP2 *
      TMP4) + cI * (TMP9 * TMP23)) + (P2[2] * (-cI * (TMP1 * TMP4) + cI * (TMP7
      * TMP26)) + (TMP10 * - 1. * (+cI * (V2[4] * TMP26 + V1[4] * TMP23)) +
      (+cI * (TMP2 * V1[4] * TMP7 + TMP1 * V2[4] * TMP9)))))) + P3[2] * (P1[0]
      * (-cI * (TMP2 * TMP4) + cI * (TMP9 * TMP23)) + (P2[0] * (-cI * (TMP1 *
      TMP4) + cI * (TMP7 * TMP26)) + (TMP10 * - 1. * (+cI * (V1[2] * TMP23 +
      V2[2] * TMP26)) + (+cI * (TMP2 * V1[2] * TMP7 + TMP1 * V2[2] * TMP9))))))
      + (P1[0] * (-cI * (V2[4] * TMP9) + cI * (P2[2] * TMP4)) + (P1[2] * (-cI *
      (V2[2] * TMP9) + cI * (P2[0] * TMP4)) + (TMP10 * (+cI * (V2[4] * V1[2] +
      V2[2] * V1[4])) - TMP7 * (+cI * (P2[0] * V1[4] + P2[2] * V1[2]))))));
  T3[11] = denom * (OM3 * (P3[1] * (P3[2] * (OM3 * (TMP1 * 1.333333333 * (-cI *
      (TMP9 * TMP23) + cI * (TMP2 * TMP4)) + 1.333333333 * (TMP26 * (-cI *
      (TMP2 * TMP7) + cI * (TMP10 * TMP23)))) + (-1.333333333 * cI * (TMP7 *
      TMP9) + 1.333333333 * cI * (TMP4 * TMP10))) + (P1[2] * (-cI * (TMP2 *
      TMP4) + cI * (TMP9 * TMP23)) + (P2[2] * (-cI * (TMP1 * TMP4) + cI * (TMP7
      * TMP26)) + (TMP10 * - 1. * (+cI * (V2[4] * TMP26 + V1[4] * TMP23)) +
      (+cI * (TMP2 * V1[4] * TMP7 + TMP1 * V2[4] * TMP9)))))) + P3[2] * (P1[1]
      * (-cI * (TMP2 * TMP4) + cI * (TMP9 * TMP23)) + (P2[1] * (-cI * (TMP1 *
      TMP4) + cI * (TMP7 * TMP26)) + (TMP10 * - 1. * (+cI * (V1[3] * TMP23 +
      V2[3] * TMP26)) + (+cI * (TMP2 * V1[3] * TMP7 + TMP1 * V2[3] * TMP9))))))
      + (P1[1] * (-cI * (V2[4] * TMP9) + cI * (P2[2] * TMP4)) + (P1[2] * (-cI *
      (V2[3] * TMP9) + cI * (P2[1] * TMP4)) + (TMP10 * (+cI * (V2[4] * V1[3] +
      V2[3] * V1[4])) - TMP7 * (+cI * (P2[1] * V1[4] + P2[2] * V1[3]))))));
  T3[12] = denom * 2. * (OM3 * (P3[2] * (P3[2] * (OM3 * (TMP1 * 0.666666667 *
      (-cI * (TMP9 * TMP23) + cI * (TMP2 * TMP4)) + 0.666666667 * (TMP26 * (-cI
      * (TMP2 * TMP7) + cI * (TMP10 * TMP23)))) + (-0.666666667 * cI * (TMP7 *
      TMP9) + 0.666666667 * cI * (TMP4 * TMP10))) + (P1[2] * (-cI * (TMP2 *
      TMP4) + cI * (TMP9 * TMP23)) + (P2[2] * (-cI * (TMP1 * TMP4) + cI * (TMP7
      * TMP26)) + (TMP10 * - 1. * (+cI * (V2[4] * TMP26 + V1[4] * TMP23)) +
      (+cI * (TMP2 * V1[4] * TMP7 + TMP1 * V2[4] * TMP9)))))) + (TMP1 *
      0.333333333 * (-cI * (TMP2 * TMP4) + cI * (TMP9 * TMP23)) + 0.333333333 *
      (TMP26 * (-cI * (TMP10 * TMP23) + cI * (TMP2 * TMP7))))) + (TMP10 * (+cI
      * (V2[4] * V1[4]) + 0.666666667 * cI * (TMP4)) + (TMP7 * - 1. * (+cI *
      (P2[2] * V1[4]) + 0.666666667 * cI * (TMP9)) + P1[2] * (-cI * (V2[4] *
      TMP9) + cI * (P2[2] * TMP4)))));
  T3[13] = denom * (OM3 * (P3[2] * (P3[3] * (OM3 * (TMP1 * 1.333333333 * (-cI *
      (TMP9 * TMP23) + cI * (TMP2 * TMP4)) + 1.333333333 * (TMP26 * (-cI *
      (TMP2 * TMP7) + cI * (TMP10 * TMP23)))) + (-1.333333333 * cI * (TMP7 *
      TMP9) + 1.333333333 * cI * (TMP4 * TMP10))) + (P1[3] * (-cI * (TMP2 *
      TMP4) + cI * (TMP9 * TMP23)) + (P2[3] * (-cI * (TMP1 * TMP4) + cI * (TMP7
      * TMP26)) + (TMP10 * - 1. * (+cI * (V1[5] * TMP23 + V2[5] * TMP26)) +
      (+cI * (TMP2 * V1[5] * TMP7 + TMP1 * V2[5] * TMP9)))))) + P3[3] * (P1[2]
      * (-cI * (TMP2 * TMP4) + cI * (TMP9 * TMP23)) + (P2[2] * (-cI * (TMP1 *
      TMP4) + cI * (TMP7 * TMP26)) + (TMP10 * - 1. * (+cI * (V2[4] * TMP26 +
      V1[4] * TMP23)) + (+cI * (TMP2 * V1[4] * TMP7 + TMP1 * V2[4] * TMP9))))))
      + (P1[2] * (-cI * (V2[5] * TMP9) + cI * (P2[3] * TMP4)) + (P1[3] * (-cI *
      (V2[4] * TMP9) + cI * (P2[2] * TMP4)) + (TMP10 * (+cI * (V2[4] * V1[5] +
      V2[5] * V1[4])) - TMP7 * (+cI * (P2[3] * V1[4] + P2[2] * V1[5]))))));
  T3[14] = denom * (OM3 * (P3[0] * (P3[3] * (OM3 * (TMP1 * 1.333333333 * (-cI *
      (TMP9 * TMP23) + cI * (TMP2 * TMP4)) + 1.333333333 * (TMP26 * (-cI *
      (TMP2 * TMP7) + cI * (TMP10 * TMP23)))) + (-1.333333333 * cI * (TMP7 *
      TMP9) + 1.333333333 * cI * (TMP4 * TMP10))) + (P1[3] * (-cI * (TMP2 *
      TMP4) + cI * (TMP9 * TMP23)) + (P2[3] * (-cI * (TMP1 * TMP4) + cI * (TMP7
      * TMP26)) + (TMP10 * - 1. * (+cI * (V2[5] * TMP26 + V1[5] * TMP23)) +
      (+cI * (TMP2 * V1[5] * TMP7 + TMP1 * V2[5] * TMP9)))))) + P3[3] * (P1[0]
      * (-cI * (TMP2 * TMP4) + cI * (TMP9 * TMP23)) + (P2[0] * (-cI * (TMP1 *
      TMP4) + cI * (TMP7 * TMP26)) + (TMP10 * - 1. * (+cI * (V1[2] * TMP23 +
      V2[2] * TMP26)) + (+cI * (TMP2 * V1[2] * TMP7 + TMP1 * V2[2] * TMP9))))))
      + (P1[0] * (-cI * (V2[5] * TMP9) + cI * (P2[3] * TMP4)) + (P1[3] * (-cI *
      (V2[2] * TMP9) + cI * (P2[0] * TMP4)) + (TMP10 * (+cI * (V2[5] * V1[2] +
      V2[2] * V1[5])) - TMP7 * (+cI * (P2[0] * V1[5] + P2[3] * V1[2]))))));
  T3[15] = denom * (OM3 * (P3[1] * (P3[3] * (OM3 * (TMP1 * 1.333333333 * (-cI *
      (TMP9 * TMP23) + cI * (TMP2 * TMP4)) + 1.333333333 * (TMP26 * (-cI *
      (TMP2 * TMP7) + cI * (TMP10 * TMP23)))) + (-1.333333333 * cI * (TMP7 *
      TMP9) + 1.333333333 * cI * (TMP4 * TMP10))) + (P1[3] * (-cI * (TMP2 *
      TMP4) + cI * (TMP9 * TMP23)) + (P2[3] * (-cI * (TMP1 * TMP4) + cI * (TMP7
      * TMP26)) + (TMP10 * - 1. * (+cI * (V2[5] * TMP26 + V1[5] * TMP23)) +
      (+cI * (TMP2 * V1[5] * TMP7 + TMP1 * V2[5] * TMP9)))))) + P3[3] * (P1[1]
      * (-cI * (TMP2 * TMP4) + cI * (TMP9 * TMP23)) + (P2[1] * (-cI * (TMP1 *
      TMP4) + cI * (TMP7 * TMP26)) + (TMP10 * - 1. * (+cI * (V1[3] * TMP23 +
      V2[3] * TMP26)) + (+cI * (TMP2 * V1[3] * TMP7 + TMP1 * V2[3] * TMP9))))))
      + (P1[1] * (-cI * (V2[5] * TMP9) + cI * (P2[3] * TMP4)) + (P1[3] * (-cI *
      (V2[3] * TMP9) + cI * (P2[1] * TMP4)) + (TMP10 * (+cI * (V2[5] * V1[3] +
      V2[3] * V1[5])) - TMP7 * (+cI * (P2[1] * V1[5] + P2[3] * V1[3]))))));
  T3[16] = denom * (OM3 * (P3[2] * (P3[3] * (OM3 * (TMP1 * 1.333333333 * (-cI *
      (TMP9 * TMP23) + cI * (TMP2 * TMP4)) + 1.333333333 * (TMP26 * (-cI *
      (TMP2 * TMP7) + cI * (TMP10 * TMP23)))) + (-1.333333333 * cI * (TMP7 *
      TMP9) + 1.333333333 * cI * (TMP4 * TMP10))) + (P1[3] * (-cI * (TMP2 *
      TMP4) + cI * (TMP9 * TMP23)) + (P2[3] * (-cI * (TMP1 * TMP4) + cI * (TMP7
      * TMP26)) + (TMP10 * - 1. * (+cI * (V2[5] * TMP26 + V1[5] * TMP23)) +
      (+cI * (TMP2 * V1[5] * TMP7 + TMP1 * V2[5] * TMP9)))))) + P3[3] * (P1[2]
      * (-cI * (TMP2 * TMP4) + cI * (TMP9 * TMP23)) + (P2[2] * (-cI * (TMP1 *
      TMP4) + cI * (TMP7 * TMP26)) + (TMP10 * - 1. * (+cI * (V1[4] * TMP23 +
      V2[4] * TMP26)) + (+cI * (TMP2 * V1[4] * TMP7 + TMP1 * V2[4] * TMP9))))))
      + (P1[2] * (-cI * (V2[5] * TMP9) + cI * (P2[3] * TMP4)) + (P1[3] * (-cI *
      (V2[4] * TMP9) + cI * (P2[2] * TMP4)) + (TMP10 * (+cI * (V2[5] * V1[4] +
      V2[4] * V1[5])) - TMP7 * (+cI * (P2[2] * V1[5] + P2[3] * V1[4]))))));
  T3[17] = denom * 2. * (OM3 * (P3[3] * (P3[3] * (OM3 * (TMP1 * 0.666666667 *
      (-cI * (TMP9 * TMP23) + cI * (TMP2 * TMP4)) + 0.666666667 * (TMP26 * (-cI
      * (TMP2 * TMP7) + cI * (TMP10 * TMP23)))) + (-0.666666667 * cI * (TMP7 *
      TMP9) + 0.666666667 * cI * (TMP4 * TMP10))) + (P1[3] * (-cI * (TMP2 *
      TMP4) + cI * (TMP9 * TMP23)) + (P2[3] * (-cI * (TMP1 * TMP4) + cI * (TMP7
      * TMP26)) + (TMP10 * - 1. * (+cI * (V2[5] * TMP26 + V1[5] * TMP23)) +
      (+cI * (TMP2 * V1[5] * TMP7 + TMP1 * V2[5] * TMP9)))))) + (TMP1 *
      0.333333333 * (-cI * (TMP2 * TMP4) + cI * (TMP9 * TMP23)) + 0.333333333 *
      (TMP26 * (-cI * (TMP10 * TMP23) + cI * (TMP2 * TMP7))))) + (TMP10 * (+cI
      * (V2[5] * V1[5]) + 0.666666667 * cI * (TMP4)) + (TMP7 * - 1. * (+cI *
      (P2[3] * V1[5]) + 0.666666667 * cI * (TMP9)) + P1[3] * (-cI * (V2[5] *
      TMP9) + cI * (P2[3] * TMP4)))));
}


void VVT9_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP23; 
  double P3[4]; 
  complex<double> TMP26; 
  complex<double> TMP28; 
  complex<double> TMP27; 
  complex<double> TMP24; 
  complex<double> TMP25; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = T3[0].real(); 
  P3[1] = T3[1].real(); 
  P3[2] = T3[1].imag(); 
  P3[3] = T3[0].imag(); 
  TMP24 = (P2[0] * - 1. * (V1[3] * T3[6] + V1[4] * T3[10] + V1[5] * T3[14] -
      V1[2] * T3[2]) + (P2[1] * (V1[3] * T3[7] + V1[4] * T3[11] + V1[5] *
      T3[15] - V1[2] * T3[3]) + (P2[2] * (V1[3] * T3[8] + V1[4] * T3[12] +
      V1[5] * T3[16] - V1[2] * T3[4]) + P2[3] * (V1[3] * T3[9] + V1[4] * T3[13]
      + V1[5] * T3[17] - V1[2] * T3[5]))));
  TMP25 = (P2[0] * - 1. * (V1[3] * T3[3] + V1[4] * T3[4] + V1[5] * T3[5] -
      V1[2] * T3[2]) + (P2[1] * (V1[3] * T3[7] + V1[4] * T3[8] + V1[5] * T3[9]
      - V1[2] * T3[6]) + (P2[2] * (V1[3] * T3[11] + V1[4] * T3[12] + V1[5] *
      T3[13] - V1[2] * T3[10]) + P2[3] * (V1[3] * T3[15] + V1[4] * T3[16] +
      V1[5] * T3[17] - V1[2] * T3[14]))));
  TMP26 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP27 = (P1[0] * - 1. * (V2[3] * T3[6] + V2[4] * T3[10] + V2[5] * T3[14] -
      V2[2] * T3[2]) + (P1[1] * (V2[3] * T3[7] + V2[4] * T3[11] + V2[5] *
      T3[15] - V2[2] * T3[3]) + (P1[2] * (V2[3] * T3[8] + V2[4] * T3[12] +
      V2[5] * T3[16] - V2[2] * T3[4]) + P1[3] * (V2[3] * T3[9] + V2[4] * T3[13]
      + V2[5] * T3[17] - V2[2] * T3[5]))));
  TMP23 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP28 = (P1[0] * - 1. * (V2[3] * T3[3] + V2[4] * T3[4] + V2[5] * T3[5] -
      V2[2] * T3[2]) + (P1[1] * (V2[3] * T3[7] + V2[4] * T3[8] + V2[5] * T3[9]
      - V2[2] * T3[6]) + (P1[2] * (V2[3] * T3[11] + V2[4] * T3[12] + V2[5] *
      T3[13] - V2[2] * T3[10]) + P1[3] * (V2[3] * T3[15] + V2[4] * T3[16] +
      V2[5] * T3[17] - V2[2] * T3[14]))));
  vertex = COUP * - 1. * (TMP23 * (+cI * (TMP24 + TMP25)) + TMP26 * (+cI *
      (TMP27 + TMP28)));
}


void VVT12_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> TMP22; 
  complex<double> TMP10; 
  double P2[4]; 
  complex<double> TMP7; 
  complex<double> TMP36; 
  complex<double> TMP21; 
  complex<double> TMP4; 
  complex<double> TMP28; 
  complex<double> TMP27; 
  complex<double> TMP24; 
  complex<double> TMP25; 
  complex<double> TMP9; 
  complex<double> TMP35; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  TMP24 = (P2[0] * - 1. * (V1[3] * T3[6] + V1[4] * T3[10] + V1[5] * T3[14] -
      V1[2] * T3[2]) + (P2[1] * (V1[3] * T3[7] + V1[4] * T3[11] + V1[5] *
      T3[15] - V1[2] * T3[3]) + (P2[2] * (V1[3] * T3[8] + V1[4] * T3[12] +
      V1[5] * T3[16] - V1[2] * T3[4]) + P2[3] * (V1[3] * T3[9] + V1[4] * T3[13]
      + V1[5] * T3[17] - V1[2] * T3[5]))));
  TMP25 = (P2[0] * - 1. * (V1[3] * T3[3] + V1[4] * T3[4] + V1[5] * T3[5] -
      V1[2] * T3[2]) + (P2[1] * (V1[3] * T3[7] + V1[4] * T3[8] + V1[5] * T3[9]
      - V1[2] * T3[6]) + (P2[2] * (V1[3] * T3[11] + V1[4] * T3[12] + V1[5] *
      T3[13] - V1[2] * T3[10]) + P2[3] * (V1[3] * T3[15] + V1[4] * T3[16] +
      V1[5] * T3[17] - V1[2] * T3[14]))));
  TMP27 = (P1[0] * - 1. * (V2[3] * T3[6] + V2[4] * T3[10] + V2[5] * T3[14] -
      V2[2] * T3[2]) + (P1[1] * (V2[3] * T3[7] + V2[4] * T3[11] + V2[5] *
      T3[15] - V2[2] * T3[3]) + (P1[2] * (V2[3] * T3[8] + V2[4] * T3[12] +
      V2[5] * T3[16] - V2[2] * T3[4]) + P1[3] * (V2[3] * T3[9] + V2[4] * T3[13]
      + V2[5] * T3[17] - V2[2] * T3[5]))));
  TMP21 = (P1[0] * - 1. * (P2[1] * T3[6] + P2[2] * T3[10] + P2[3] * T3[14] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[11] + P2[3] *
      T3[15] - P2[0] * T3[3]) + (P1[2] * (P2[1] * T3[8] + P2[2] * T3[12] +
      P2[3] * T3[16] - P2[0] * T3[4]) + P1[3] * (P2[1] * T3[9] + P2[2] * T3[13]
      + P2[3] * T3[17] - P2[0] * T3[5]))));
  TMP22 = (P1[0] * - 1. * (P2[1] * T3[3] + P2[2] * T3[4] + P2[3] * T3[5] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[8] + P2[3] * T3[9]
      - P2[0] * T3[6]) + (P1[2] * (P2[1] * T3[11] + P2[2] * T3[12] + P2[3] *
      T3[13] - P2[0] * T3[10]) + P1[3] * (P2[1] * T3[15] + P2[2] * T3[16] +
      P2[3] * T3[17] - P2[0] * T3[14]))));
  TMP28 = (P1[0] * - 1. * (V2[3] * T3[3] + V2[4] * T3[4] + V2[5] * T3[5] -
      V2[2] * T3[2]) + (P1[1] * (V2[3] * T3[7] + V2[4] * T3[8] + V2[5] * T3[9]
      - V2[2] * T3[6]) + (P1[2] * (V2[3] * T3[11] + V2[4] * T3[12] + V2[5] *
      T3[13] - V2[2] * T3[10]) + P1[3] * (V2[3] * T3[15] + V2[4] * T3[16] +
      V2[5] * T3[17] - V2[2] * T3[14]))));
  TMP10 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP9 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP4 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP7 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP36 = (V1[2] * - 1. * (V2[3] * T3[3] + V2[4] * T3[4] + V2[5] * T3[5] -
      V2[2] * T3[2]) + (V1[3] * (V2[3] * T3[7] + V2[4] * T3[8] + V2[5] * T3[9]
      - V2[2] * T3[6]) + (V1[4] * (V2[3] * T3[11] + V2[4] * T3[12] + V2[5] *
      T3[13] - V2[2] * T3[10]) + V1[5] * (V2[3] * T3[15] + V2[4] * T3[16] +
      V2[5] * T3[17] - V2[2] * T3[14]))));
  TMP35 = (V1[2] * - 1. * (V2[3] * T3[6] + V2[4] * T3[10] + V2[5] * T3[14] -
      V2[2] * T3[2]) + (V1[3] * (V2[3] * T3[7] + V2[4] * T3[11] + V2[5] *
      T3[15] - V2[2] * T3[3]) + (V1[4] * (V2[3] * T3[8] + V2[4] * T3[12] +
      V2[5] * T3[16] - V2[2] * T3[4]) + V1[5] * (V2[3] * T3[9] + V2[4] * T3[13]
      + V2[5] * T3[17] - V2[2] * T3[5]))));
  vertex = COUP * (TMP10 * - 1. * (+cI * (TMP35 + TMP36)) + (TMP4 * - 1. * (+cI
      * (TMP21 + TMP22)) + (TMP7 * (+cI * (TMP24 + TMP25)) + TMP9 * (+cI *
      (TMP27 + TMP28)))));
}


void VVT3_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP2; 
  complex<double> TMP1; 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP56; 
  double P3[4]; 
  complex<double> TMP55; 
  complex<double> denom; 
  double OM3; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +V1[0] + V2[0]; 
  T3[1] = +V1[1] + V2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP55 = -1. * (P1[0] * (P3[1] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] *
      (V2[5] * V1[3] - V2[3] * V1[5]) + P3[3] * (V2[3] * V1[4] - V2[4] *
      V1[3]))) + (P1[1] * (P3[0] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] * (V2[4] * V1[2] - V2[2] *
      V1[4]))) + (P1[2] * (P3[0] * (V2[3] * V1[5] - V2[5] * V1[3]) + (P3[1] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] * (V2[2] * V1[3] - V2[3] *
      V1[2]))) + P1[3] * (P3[0] * (V2[4] * V1[3] - V2[3] * V1[4]) + (P3[1] *
      (V2[2] * V1[4] - V2[4] * V1[2]) + P3[2] * (V2[3] * V1[2] - V2[2] *
      V1[3]))))));
  TMP56 = -1. * (P2[0] * (P3[1] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] *
      (V2[5] * V1[3] - V2[3] * V1[5]) + P3[3] * (V2[3] * V1[4] - V2[4] *
      V1[3]))) + (P2[1] * (P3[0] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] * (V2[4] * V1[2] - V2[2] *
      V1[4]))) + (P2[2] * (P3[0] * (V2[3] * V1[5] - V2[5] * V1[3]) + (P3[1] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] * (V2[2] * V1[3] - V2[3] *
      V1[2]))) + P2[3] * (P3[0] * (V2[4] * V1[3] - V2[3] * V1[4]) + (P3[1] *
      (V2[2] * V1[4] - V2[4] * V1[2]) + P3[2] * (V2[3] * V1[2] - V2[2] *
      V1[3]))))));
  TMP1 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP2 = (P3[0] * P2[0] - P3[1] * P2[1] - P3[2] * P2[2] - P3[3] * P2[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * - 2. * cI * (OM3 * P3[0] * (TMP1 * (P3[1] * (V2[5] * V1[4] -
      V2[4] * V1[5]) + (P3[2] * (V2[3] * V1[5] - V2[5] * V1[3]) + P3[3] *
      (V2[4] * V1[3] - V2[3] * V1[4]))) + (TMP2 * (P3[1] * (V2[4] * V1[5] -
      V2[5] * V1[4]) + (P3[2] * (V2[5] * V1[3] - V2[3] * V1[5]) + P3[3] *
      (V2[3] * V1[4] - V2[4] * V1[3]))) + 0.333333333 * (P3[0] * (TMP56 -
      TMP55)))) + (P1[0] * (P3[1] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] *
      (V2[5] * V1[3] - V2[3] * V1[5]) + P3[3] * (V2[3] * V1[4] - V2[4] *
      V1[3]))) + (P2[0] * (P3[1] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] *
      (V2[3] * V1[5] - V2[5] * V1[3]) + P3[3] * (V2[4] * V1[3] - V2[3] *
      V1[4]))) + (-0.333333333 * (TMP56) + 0.333333333 * (TMP55)))));
  T3[3] = denom * cI * (OM3 * (P3[0] * (TMP1 * (P3[0] * (V2[4] * V1[5] - V2[5]
      * V1[4]) + (P3[2] * (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] * (V2[2] *
      V1[4] - V2[4] * V1[2]))) + (TMP2 * (P3[0] * (V2[5] * V1[4] - V2[4] *
      V1[5]) + (P3[2] * (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] * (V2[4] *
      V1[2] - V2[2] * V1[4]))) + 0.666666667 * (P3[1] * (TMP55 - TMP56)))) +
      P3[1] * (TMP1 * (P3[1] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] *
      (V2[5] * V1[3] - V2[3] * V1[5]) + P3[3] * (V2[3] * V1[4] - V2[4] *
      V1[3]))) + TMP2 * (P3[1] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] *
      (V2[3] * V1[5] - V2[5] * V1[3]) + P3[3] * (V2[4] * V1[3] - V2[3] *
      V1[4]))))) + (P3[2] * (V1[5] * (V2[2] * (P1[0] - P2[0]) + V2[3] * (P1[1]
      - P2[1])) + V2[5] * (V1[2] * (P2[0] - P1[0]) + V1[3] * (P2[1] - P1[1])))
      + (P3[3] * (V1[4] * (V2[2] * (P2[0] - P1[0]) + V2[3] * (P2[1] - P1[1])) +
      V2[4] * (V1[2] * (P1[0] - P2[0]) + V1[3] * (P1[1] - P2[1]))) + (P3[0] *
      (V1[4] * V2[5] * (P1[0] - P2[0]) + V1[5] * V2[4] * (P2[0] - P1[0])) +
      P3[1] * (V1[4] * V2[5] * (P1[1] - P2[1]) + V1[5] * V2[4] * (P2[1] -
      P1[1]))))));
  T3[4] = denom * cI * (OM3 * (P3[0] * (TMP1 * (P3[0] * (V2[5] * V1[3] - V2[3]
      * V1[5]) + (P3[1] * (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] * (V2[3] *
      V1[2] - V2[2] * V1[3]))) + (TMP2 * (P3[0] * (V2[3] * V1[5] - V2[5] *
      V1[3]) + (P3[1] * (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] * (V2[2] *
      V1[3] - V2[3] * V1[2]))) + 0.666666667 * (P3[2] * (TMP55 - TMP56)))) +
      P3[2] * (TMP1 * (P3[1] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] *
      (V2[5] * V1[3] - V2[3] * V1[5]) + P3[3] * (V2[3] * V1[4] - V2[4] *
      V1[3]))) + TMP2 * (P3[1] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] *
      (V2[3] * V1[5] - V2[5] * V1[3]) + P3[3] * (V2[4] * V1[3] - V2[3] *
      V1[4]))))) + (P3[1] * (V1[5] * (V2[2] * (P2[0] - P1[0]) + V2[4] * (P2[2]
      - P1[2])) + V2[5] * (V1[2] * (P1[0] - P2[0]) + V1[4] * (P1[2] - P2[2])))
      + (P3[3] * (V1[3] * (V2[2] * (P1[0] - P2[0]) + V2[4] * (P1[2] - P2[2])) +
      V2[3] * (V1[2] * (P2[0] - P1[0]) + V1[4] * (P2[2] - P1[2]))) + (P3[0] *
      (V1[3] * V2[5] * (P2[0] - P1[0]) + V1[5] * V2[3] * (P1[0] - P2[0])) +
      P3[2] * (V1[3] * V2[5] * (P2[2] - P1[2]) + V1[5] * V2[3] * (P1[2] -
      P2[2]))))));
  T3[5] = denom * cI * (OM3 * (P3[0] * (TMP1 * (P3[0] * (V2[3] * V1[4] - V2[4]
      * V1[3]) + (P3[1] * (V2[4] * V1[2] - V2[2] * V1[4]) + P3[2] * (V2[2] *
      V1[3] - V2[3] * V1[2]))) + (TMP2 * (P3[0] * (V2[4] * V1[3] - V2[3] *
      V1[4]) + (P3[1] * (V2[2] * V1[4] - V2[4] * V1[2]) + P3[2] * (V2[3] *
      V1[2] - V2[2] * V1[3]))) + 0.666666667 * (P3[3] * (TMP55 - TMP56)))) +
      P3[3] * (TMP1 * (P3[1] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] *
      (V2[5] * V1[3] - V2[3] * V1[5]) + P3[3] * (V2[3] * V1[4] - V2[4] *
      V1[3]))) + TMP2 * (P3[1] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] *
      (V2[3] * V1[5] - V2[5] * V1[3]) + P3[3] * (V2[4] * V1[3] - V2[3] *
      V1[4]))))) + (P3[1] * (V1[4] * (V2[2] * (P1[0] - P2[0]) + V2[5] * (P1[3]
      - P2[3])) + V2[4] * (V1[2] * (P2[0] - P1[0]) + V1[5] * (P2[3] - P1[3])))
      + (P3[2] * (V1[3] * (V2[2] * (P2[0] - P1[0]) + V2[5] * (P2[3] - P1[3])) +
      V2[3] * (V1[2] * (P1[0] - P2[0]) + V1[5] * (P1[3] - P2[3]))) + (P3[0] *
      (V1[3] * V2[4] * (P1[0] - P2[0]) + V1[4] * V2[3] * (P2[0] - P1[0])) +
      P3[3] * (V1[3] * V2[4] * (P1[3] - P2[3]) + V1[4] * V2[3] * (P2[3] -
      P1[3]))))));
  T3[6] = denom * - cI * (OM3 * (P3[0] * (TMP1 * (P3[0] * (V2[5] * V1[4] -
      V2[4] * V1[5]) + (P3[2] * (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] *
      (V2[4] * V1[2] - V2[2] * V1[4]))) + (TMP2 * (P3[0] * (V2[4] * V1[5] -
      V2[5] * V1[4]) + (P3[2] * (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] *
      (V2[2] * V1[4] - V2[4] * V1[2]))) + 0.666666667 * (P3[1] * (TMP56 -
      TMP55)))) + P3[1] * (TMP1 * (P3[1] * (V2[5] * V1[4] - V2[4] * V1[5]) +
      (P3[2] * (V2[3] * V1[5] - V2[5] * V1[3]) + P3[3] * (V2[4] * V1[3] - V2[3]
      * V1[4]))) + TMP2 * (P3[1] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] *
      (V2[5] * V1[3] - V2[3] * V1[5]) + P3[3] * (V2[3] * V1[4] - V2[4] *
      V1[3]))))) + (P3[2] * (V1[5] * (V2[2] * (P2[0] - P1[0]) + V2[3] * (P2[1]
      - P1[1])) + V2[5] * (V1[2] * (P1[0] - P2[0]) + V1[3] * (P1[1] - P2[1])))
      + (P3[3] * (V1[4] * (V2[2] * (P1[0] - P2[0]) + V2[3] * (P1[1] - P2[1])) +
      V2[4] * (V1[2] * (P2[0] - P1[0]) + V1[3] * (P2[1] - P1[1]))) + (P3[0] *
      (V1[4] * V2[5] * (P2[0] - P1[0]) + V1[5] * V2[4] * (P1[0] - P2[0])) +
      P3[1] * (V1[4] * V2[5] * (P2[1] - P1[1]) + V1[5] * V2[4] * (P1[1] -
      P2[1]))))));
  T3[7] = denom * 2. * cI * (OM3 * P3[1] * (TMP1 * (P3[0] * (V2[4] * V1[5] -
      V2[5] * V1[4]) + (P3[2] * (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] *
      (V2[2] * V1[4] - V2[4] * V1[2]))) + (TMP2 * (P3[0] * (V2[5] * V1[4] -
      V2[4] * V1[5]) + (P3[2] * (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] *
      (V2[4] * V1[2] - V2[2] * V1[4]))) + 0.333333333 * (P3[1] * (TMP55 -
      TMP56)))) + (P1[1] * (P3[0] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] * (V2[4] * V1[2] - V2[2] *
      V1[4]))) + (P2[1] * (P3[0] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] * (V2[2] * V1[4] - V2[4] *
      V1[2]))) + (-0.333333333 * (TMP56) + 0.333333333 * (TMP55)))));
  T3[8] = denom * cI * (OM3 * (P3[1] * (TMP1 * (P3[0] * (V2[5] * V1[3] - V2[3]
      * V1[5]) + (P3[1] * (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] * (V2[3] *
      V1[2] - V2[2] * V1[3]))) + (TMP2 * (P3[0] * (V2[3] * V1[5] - V2[5] *
      V1[3]) + (P3[1] * (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] * (V2[2] *
      V1[3] - V2[3] * V1[2]))) + 0.666666667 * (P3[2] * (TMP55 - TMP56)))) +
      P3[2] * (TMP1 * (P3[0] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] * (V2[2] * V1[4] - V2[4] *
      V1[2]))) + TMP2 * (P3[0] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] * (V2[4] * V1[2] - V2[2] *
      V1[4]))))) + (P3[0] * (V1[5] * (V2[3] * (P1[1] - P2[1]) + V2[4] * (P2[2]
      - P1[2])) + V2[5] * (V1[3] * (P2[1] - P1[1]) + V1[4] * (P1[2] - P2[2])))
      + (P3[3] * (V1[2] * (V2[3] * (P2[1] - P1[1]) + V2[4] * (P1[2] - P2[2])) +
      V2[2] * (V1[3] * (P1[1] - P2[1]) + V1[4] * (P2[2] - P1[2]))) + (P3[1] *
      (V1[2] * V2[5] * (P1[1] - P2[1]) + V1[5] * V2[2] * (P2[1] - P1[1])) +
      P3[2] * (V1[2] * V2[5] * (P2[2] - P1[2]) + V1[5] * V2[2] * (P1[2] -
      P2[2]))))));
  T3[9] = denom * cI * (OM3 * (P3[1] * (TMP1 * (P3[0] * (V2[3] * V1[4] - V2[4]
      * V1[3]) + (P3[1] * (V2[4] * V1[2] - V2[2] * V1[4]) + P3[2] * (V2[2] *
      V1[3] - V2[3] * V1[2]))) + (TMP2 * (P3[0] * (V2[4] * V1[3] - V2[3] *
      V1[4]) + (P3[1] * (V2[2] * V1[4] - V2[4] * V1[2]) + P3[2] * (V2[3] *
      V1[2] - V2[2] * V1[3]))) + 0.666666667 * (P3[3] * (TMP55 - TMP56)))) +
      P3[3] * (TMP1 * (P3[0] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] * (V2[2] * V1[4] - V2[4] *
      V1[2]))) + TMP2 * (P3[0] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] * (V2[4] * V1[2] - V2[2] *
      V1[4]))))) + (P3[0] * (V1[4] * (V2[3] * (P2[1] - P1[1]) + V2[5] * (P1[3]
      - P2[3])) + V2[4] * (V1[3] * (P1[1] - P2[1]) + V1[5] * (P2[3] - P1[3])))
      + (P3[2] * (V1[2] * (V2[3] * (P1[1] - P2[1]) + V2[5] * (P2[3] - P1[3])) +
      V2[2] * (V1[3] * (P2[1] - P1[1]) + V1[5] * (P1[3] - P2[3]))) + (P3[1] *
      (V1[2] * V2[4] * (P2[1] - P1[1]) + V1[4] * V2[2] * (P1[1] - P2[1])) +
      P3[3] * (V1[2] * V2[4] * (P1[3] - P2[3]) + V1[4] * V2[2] * (P2[3] -
      P1[3]))))));
  T3[10] = denom * - cI * (OM3 * (P3[0] * (TMP1 * (P3[0] * (V2[3] * V1[5] -
      V2[5] * V1[3]) + (P3[1] * (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] *
      (V2[2] * V1[3] - V2[3] * V1[2]))) + (TMP2 * (P3[0] * (V2[5] * V1[3] -
      V2[3] * V1[5]) + (P3[1] * (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] *
      (V2[3] * V1[2] - V2[2] * V1[3]))) + 0.666666667 * (P3[2] * (TMP56 -
      TMP55)))) + P3[2] * (TMP1 * (P3[1] * (V2[5] * V1[4] - V2[4] * V1[5]) +
      (P3[2] * (V2[3] * V1[5] - V2[5] * V1[3]) + P3[3] * (V2[4] * V1[3] - V2[3]
      * V1[4]))) + TMP2 * (P3[1] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] *
      (V2[5] * V1[3] - V2[3] * V1[5]) + P3[3] * (V2[3] * V1[4] - V2[4] *
      V1[3]))))) + (P3[1] * (V1[5] * (V2[2] * (P1[0] - P2[0]) + V2[4] * (P1[2]
      - P2[2])) + V2[5] * (V1[2] * (P2[0] - P1[0]) + V1[4] * (P2[2] - P1[2])))
      + (P3[3] * (V1[3] * (V2[2] * (P2[0] - P1[0]) + V2[4] * (P2[2] - P1[2])) +
      V2[3] * (V1[2] * (P1[0] - P2[0]) + V1[4] * (P1[2] - P2[2]))) + (P3[0] *
      (V1[3] * V2[5] * (P1[0] - P2[0]) + V1[5] * V2[3] * (P2[0] - P1[0])) +
      P3[2] * (V1[3] * V2[5] * (P1[2] - P2[2]) + V1[5] * V2[3] * (P2[2] -
      P1[2]))))));
  T3[11] = denom * cI * (OM3 * (P3[1] * (TMP1 * (P3[0] * (V2[5] * V1[3] - V2[3]
      * V1[5]) + (P3[1] * (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] * (V2[3] *
      V1[2] - V2[2] * V1[3]))) + (TMP2 * (P3[0] * (V2[3] * V1[5] - V2[5] *
      V1[3]) + (P3[1] * (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] * (V2[2] *
      V1[3] - V2[3] * V1[2]))) + 0.666666667 * (P3[2] * (TMP55 - TMP56)))) +
      P3[2] * (TMP1 * (P3[0] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] * (V2[2] * V1[4] - V2[4] *
      V1[2]))) + TMP2 * (P3[0] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] * (V2[4] * V1[2] - V2[2] *
      V1[4]))))) + (P3[0] * (V1[5] * (V2[3] * (P1[1] - P2[1]) + V2[4] * (P2[2]
      - P1[2])) + V2[5] * (V1[3] * (P2[1] - P1[1]) + V1[4] * (P1[2] - P2[2])))
      + (P3[3] * (V1[2] * (V2[3] * (P2[1] - P1[1]) + V2[4] * (P1[2] - P2[2])) +
      V2[2] * (V1[3] * (P1[1] - P2[1]) + V1[4] * (P2[2] - P1[2]))) + (P3[1] *
      (V1[2] * V2[5] * (P1[1] - P2[1]) + V1[5] * V2[2] * (P2[1] - P1[1])) +
      P3[2] * (V1[2] * V2[5] * (P2[2] - P1[2]) + V1[5] * V2[2] * (P1[2] -
      P2[2]))))));
  T3[12] = denom * 2. * cI * (OM3 * P3[2] * (TMP1 * (P3[0] * (V2[5] * V1[3] -
      V2[3] * V1[5]) + (P3[1] * (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] *
      (V2[3] * V1[2] - V2[2] * V1[3]))) + (TMP2 * (P3[0] * (V2[3] * V1[5] -
      V2[5] * V1[3]) + (P3[1] * (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] *
      (V2[2] * V1[3] - V2[3] * V1[2]))) + 0.333333333 * (P3[2] * (TMP55 -
      TMP56)))) + (P1[2] * (P3[0] * (V2[3] * V1[5] - V2[5] * V1[3]) + (P3[1] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] * (V2[2] * V1[3] - V2[3] *
      V1[2]))) + (P2[2] * (P3[0] * (V2[5] * V1[3] - V2[3] * V1[5]) + (P3[1] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] * (V2[3] * V1[2] - V2[2] *
      V1[3]))) + (-0.333333333 * (TMP56) + 0.333333333 * (TMP55)))));
  T3[13] = denom * cI * (OM3 * (P3[2] * (TMP1 * (P3[0] * (V2[3] * V1[4] - V2[4]
      * V1[3]) + (P3[1] * (V2[4] * V1[2] - V2[2] * V1[4]) + P3[2] * (V2[2] *
      V1[3] - V2[3] * V1[2]))) + (TMP2 * (P3[0] * (V2[4] * V1[3] - V2[3] *
      V1[4]) + (P3[1] * (V2[2] * V1[4] - V2[4] * V1[2]) + P3[2] * (V2[3] *
      V1[2] - V2[2] * V1[3]))) + 0.666666667 * (P3[3] * (TMP55 - TMP56)))) +
      P3[3] * (TMP1 * (P3[0] * (V2[5] * V1[3] - V2[3] * V1[5]) + (P3[1] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] * (V2[3] * V1[2] - V2[2] *
      V1[3]))) + TMP2 * (P3[0] * (V2[3] * V1[5] - V2[5] * V1[3]) + (P3[1] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] * (V2[2] * V1[3] - V2[3] *
      V1[2]))))) + (P3[0] * (V1[3] * (V2[4] * (P1[2] - P2[2]) + V2[5] * (P2[3]
      - P1[3])) + V2[3] * (V1[4] * (P2[2] - P1[2]) + V1[5] * (P1[3] - P2[3])))
      + (P3[1] * (V1[2] * (V2[4] * (P2[2] - P1[2]) + V2[5] * (P1[3] - P2[3])) +
      V2[2] * (V1[4] * (P1[2] - P2[2]) + V1[5] * (P2[3] - P1[3]))) + (P3[2] *
      (V1[2] * V2[3] * (P1[2] - P2[2]) + V1[3] * V2[2] * (P2[2] - P1[2])) +
      P3[3] * (V1[2] * V2[3] * (P2[3] - P1[3]) + V1[3] * V2[2] * (P1[3] -
      P2[3]))))));
  T3[14] = denom * - cI * (OM3 * (P3[0] * (TMP1 * (P3[0] * (V2[4] * V1[3] -
      V2[3] * V1[4]) + (P3[1] * (V2[2] * V1[4] - V2[4] * V1[2]) + P3[2] *
      (V2[3] * V1[2] - V2[2] * V1[3]))) + (TMP2 * (P3[0] * (V2[3] * V1[4] -
      V2[4] * V1[3]) + (P3[1] * (V2[4] * V1[2] - V2[2] * V1[4]) + P3[2] *
      (V2[2] * V1[3] - V2[3] * V1[2]))) + 0.666666667 * (P3[3] * (TMP56 -
      TMP55)))) + P3[3] * (TMP1 * (P3[1] * (V2[5] * V1[4] - V2[4] * V1[5]) +
      (P3[2] * (V2[3] * V1[5] - V2[5] * V1[3]) + P3[3] * (V2[4] * V1[3] - V2[3]
      * V1[4]))) + TMP2 * (P3[1] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] *
      (V2[5] * V1[3] - V2[3] * V1[5]) + P3[3] * (V2[3] * V1[4] - V2[4] *
      V1[3]))))) + (P3[1] * (V1[4] * (V2[2] * (P2[0] - P1[0]) + V2[5] * (P2[3]
      - P1[3])) + V2[4] * (V1[2] * (P1[0] - P2[0]) + V1[5] * (P1[3] - P2[3])))
      + (P3[2] * (V1[3] * (V2[2] * (P1[0] - P2[0]) + V2[5] * (P1[3] - P2[3])) +
      V2[3] * (V1[2] * (P2[0] - P1[0]) + V1[5] * (P2[3] - P1[3]))) + (P3[0] *
      (V1[3] * V2[4] * (P2[0] - P1[0]) + V1[4] * V2[3] * (P1[0] - P2[0])) +
      P3[3] * (V1[3] * V2[4] * (P2[3] - P1[3]) + V1[4] * V2[3] * (P1[3] -
      P2[3]))))));
  T3[15] = denom * cI * (OM3 * (P3[1] * (TMP1 * (P3[0] * (V2[3] * V1[4] - V2[4]
      * V1[3]) + (P3[1] * (V2[4] * V1[2] - V2[2] * V1[4]) + P3[2] * (V2[2] *
      V1[3] - V2[3] * V1[2]))) + (TMP2 * (P3[0] * (V2[4] * V1[3] - V2[3] *
      V1[4]) + (P3[1] * (V2[2] * V1[4] - V2[4] * V1[2]) + P3[2] * (V2[3] *
      V1[2] - V2[2] * V1[3]))) + 0.666666667 * (P3[3] * (TMP55 - TMP56)))) +
      P3[3] * (TMP1 * (P3[0] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] * (V2[2] * V1[4] - V2[4] *
      V1[2]))) + TMP2 * (P3[0] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] * (V2[4] * V1[2] - V2[2] *
      V1[4]))))) + (P3[0] * (V1[4] * (V2[3] * (P2[1] - P1[1]) + V2[5] * (P1[3]
      - P2[3])) + V2[4] * (V1[3] * (P1[1] - P2[1]) + V1[5] * (P2[3] - P1[3])))
      + (P3[2] * (V1[2] * (V2[3] * (P1[1] - P2[1]) + V2[5] * (P2[3] - P1[3])) +
      V2[2] * (V1[3] * (P2[1] - P1[1]) + V1[5] * (P1[3] - P2[3]))) + (P3[1] *
      (V1[2] * V2[4] * (P2[1] - P1[1]) + V1[4] * V2[2] * (P1[1] - P2[1])) +
      P3[3] * (V1[2] * V2[4] * (P1[3] - P2[3]) + V1[4] * V2[2] * (P2[3] -
      P1[3]))))));
  T3[16] = denom * cI * (OM3 * (P3[2] * (TMP1 * (P3[0] * (V2[3] * V1[4] - V2[4]
      * V1[3]) + (P3[1] * (V2[4] * V1[2] - V2[2] * V1[4]) + P3[2] * (V2[2] *
      V1[3] - V2[3] * V1[2]))) + (TMP2 * (P3[0] * (V2[4] * V1[3] - V2[3] *
      V1[4]) + (P3[1] * (V2[2] * V1[4] - V2[4] * V1[2]) + P3[2] * (V2[3] *
      V1[2] - V2[2] * V1[3]))) + 0.666666667 * (P3[3] * (TMP55 - TMP56)))) +
      P3[3] * (TMP1 * (P3[0] * (V2[5] * V1[3] - V2[3] * V1[5]) + (P3[1] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] * (V2[3] * V1[2] - V2[2] *
      V1[3]))) + TMP2 * (P3[0] * (V2[3] * V1[5] - V2[5] * V1[3]) + (P3[1] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] * (V2[2] * V1[3] - V2[3] *
      V1[2]))))) + (P3[0] * (V1[3] * (V2[4] * (P1[2] - P2[2]) + V2[5] * (P2[3]
      - P1[3])) + V2[3] * (V1[4] * (P2[2] - P1[2]) + V1[5] * (P1[3] - P2[3])))
      + (P3[1] * (V1[2] * (V2[4] * (P2[2] - P1[2]) + V2[5] * (P1[3] - P2[3])) +
      V2[2] * (V1[4] * (P1[2] - P2[2]) + V1[5] * (P2[3] - P1[3]))) + (P3[2] *
      (V1[2] * V2[3] * (P1[2] - P2[2]) + V1[3] * V2[2] * (P2[2] - P1[2])) +
      P3[3] * (V1[2] * V2[3] * (P2[3] - P1[3]) + V1[3] * V2[2] * (P1[3] -
      P2[3]))))));
  T3[17] = denom * 2. * cI * (OM3 * P3[3] * (TMP1 * (P3[0] * (V2[3] * V1[4] -
      V2[4] * V1[3]) + (P3[1] * (V2[4] * V1[2] - V2[2] * V1[4]) + P3[2] *
      (V2[2] * V1[3] - V2[3] * V1[2]))) + (TMP2 * (P3[0] * (V2[4] * V1[3] -
      V2[3] * V1[4]) + (P3[1] * (V2[2] * V1[4] - V2[4] * V1[2]) + P3[2] *
      (V2[3] * V1[2] - V2[2] * V1[3]))) + 0.333333333 * (P3[3] * (TMP55 -
      TMP56)))) + (P1[3] * (P3[0] * (V2[4] * V1[3] - V2[3] * V1[4]) + (P3[1] *
      (V2[2] * V1[4] - V2[4] * V1[2]) + P3[2] * (V2[3] * V1[2] - V2[2] *
      V1[3]))) + (P2[3] * (P3[0] * (V2[3] * V1[4] - V2[4] * V1[3]) + (P3[1] *
      (V2[4] * V1[2] - V2[2] * V1[4]) + P3[2] * (V2[2] * V1[3] - V2[3] *
      V1[2]))) + (-0.333333333 * (TMP56) + 0.333333333 * (TMP55)))));
}


void VVT5_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP34; 
  complex<double> TMP1; 
  double P1[4]; 
  complex<double> TMP2; 
  double P2[4]; 
  complex<double> TMP23; 
  double P3[4]; 
  double OM3; 
  complex<double> TMP33; 
  complex<double> TMP26; 
  complex<double> denom; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +V1[0] + V2[0]; 
  T3[1] = +V1[1] + V2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP26 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP23 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP33 = -1. * (P1[0] * (P2[1] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] *
      (P3[1] * V1[5] - P3[3] * V1[3]) + P2[3] * (P3[2] * V1[3] - P3[1] *
      V1[4]))) + (P1[1] * (P2[0] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] *
      (P3[3] * V1[2] - P3[0] * V1[5]) + P2[3] * (P3[0] * V1[4] - P3[2] *
      V1[2]))) + (P1[2] * (P2[0] * (P3[3] * V1[3] - P3[1] * V1[5]) + (P2[1] *
      (P3[0] * V1[5] - P3[3] * V1[2]) + P2[3] * (P3[1] * V1[2] - P3[0] *
      V1[3]))) + P1[3] * (P2[0] * (P3[1] * V1[4] - P3[2] * V1[3]) + (P2[1] *
      (P3[2] * V1[2] - P3[0] * V1[4]) + P2[2] * (P3[0] * V1[3] - P3[1] *
      V1[2]))))));
  TMP34 = -1. * (P1[0] * (P2[1] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P2[2] *
      (P3[3] * V2[3] - P3[1] * V2[5]) + P2[3] * (P3[1] * V2[4] - P3[2] *
      V2[3]))) + (P1[1] * (P2[0] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P2[2] *
      (P3[0] * V2[5] - P3[3] * V2[2]) + P2[3] * (P3[2] * V2[2] - P3[0] *
      V2[4]))) + (P1[2] * (P2[0] * (P3[1] * V2[5] - P3[3] * V2[3]) + (P2[1] *
      (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3] * (P3[0] * V2[3] - P3[1] *
      V2[2]))) + P1[3] * (P2[0] * (P3[2] * V2[3] - P3[1] * V2[4]) + (P2[1] *
      (P3[0] * V2[4] - P3[2] * V2[2]) + P2[2] * (P3[1] * V2[2] - P3[0] *
      V2[3]))))));
  TMP1 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP2 = (P3[0] * P2[0] - P3[1] * P2[1] - P3[2] * P2[2] - P3[3] * P2[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * - 2. * cI * (TMP23 * (OM3 * P3[0] * (TMP1 * (P2[1] * (P3[2] *
      V1[5] - P3[3] * V1[4]) + (P2[2] * (P3[3] * V1[3] - P3[1] * V1[5]) + P2[3]
      * (P3[1] * V1[4] - P3[2] * V1[3]))) - 0.333333333 * (P3[0] * TMP33)) +
      (P1[0] * (P2[1] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[1] *
      V1[5] - P3[3] * V1[3]) + P2[3] * (P3[2] * V1[3] - P3[1] * V1[4]))) +
      0.333333333 * (TMP33))) + TMP26 * (OM3 * P3[0] * (TMP2 * (P1[1] * (P3[2]
      * V2[5] - P3[3] * V2[4]) + (P1[2] * (P3[3] * V2[3] - P3[1] * V2[5]) +
      P1[3] * (P3[1] * V2[4] - P3[2] * V2[3]))) - 0.333333333 * (P3[0] *
      TMP34)) + (P2[0] * (P1[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P1[2] *
      (P3[1] * V2[5] - P3[3] * V2[3]) + P1[3] * (P3[2] * V2[3] - P3[1] *
      V2[4]))) + 0.333333333 * (TMP34))));
  T3[3] = denom * cI * (OM3 * (P3[0] * (TMP23 * (TMP1 * (P2[0] * (P3[3] * V1[4]
      - P3[2] * V1[5]) + (P2[2] * (P3[0] * V1[5] - P3[3] * V1[2]) + P2[3] *
      (P3[2] * V1[2] - P3[0] * V1[4]))) + 0.666666667 * (P3[1] * TMP33)) +
      TMP26 * (TMP2 * (P1[0] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P1[2] *
      (P3[0] * V2[5] - P3[3] * V2[2]) + P1[3] * (P3[2] * V2[2] - P3[0] *
      V2[4]))) + 0.666666667 * (P3[1] * TMP34))) + P3[1] * (TMP1 * TMP23 *
      (P2[1] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[1] * V1[5] -
      P3[3] * V1[3]) + P2[3] * (P3[2] * V1[3] - P3[1] * V1[4]))) + TMP2 * TMP26
      * (P1[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[1] * V2[5] -
      P3[3] * V2[3]) + P1[3] * (P3[2] * V2[3] - P3[1] * V2[4]))))) + (TMP23 *
      (P1[0] * (P2[0] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] * (P3[3] *
      V1[2] - P3[0] * V1[5]) + P2[3] * (P3[0] * V1[4] - P3[2] * V1[2]))) +
      P1[1] * (P2[1] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] * (P3[3] *
      V1[3] - P3[1] * V1[5]) + P2[3] * (P3[1] * V1[4] - P3[2] * V1[3])))) +
      TMP26 * (P2[0] * (P1[0] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P1[2] *
      (P3[3] * V2[2] - P3[0] * V2[5]) + P1[3] * (P3[0] * V2[4] - P3[2] *
      V2[2]))) + P2[1] * (P1[1] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P1[2] *
      (P3[3] * V2[3] - P3[1] * V2[5]) + P1[3] * (P3[1] * V2[4] - P3[2] *
      V2[3]))))));
  T3[4] = denom * cI * (OM3 * (P3[0] * (TMP23 * (TMP1 * (P2[0] * (P3[1] * V1[5]
      - P3[3] * V1[3]) + (P2[1] * (P3[3] * V1[2] - P3[0] * V1[5]) + P2[3] *
      (P3[0] * V1[3] - P3[1] * V1[2]))) + 0.666666667 * (P3[2] * TMP33)) +
      TMP26 * (TMP2 * (P1[0] * (P3[1] * V2[5] - P3[3] * V2[3]) + (P1[1] *
      (P3[3] * V2[2] - P3[0] * V2[5]) + P1[3] * (P3[0] * V2[3] - P3[1] *
      V2[2]))) + 0.666666667 * (P3[2] * TMP34))) + P3[2] * (TMP1 * TMP23 *
      (P2[1] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[1] * V1[5] -
      P3[3] * V1[3]) + P2[3] * (P3[2] * V1[3] - P3[1] * V1[4]))) + TMP2 * TMP26
      * (P1[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[1] * V2[5] -
      P3[3] * V2[3]) + P1[3] * (P3[2] * V2[3] - P3[1] * V2[4]))))) + (TMP23 *
      (P1[0] * (P2[0] * (P3[3] * V1[3] - P3[1] * V1[5]) + (P2[1] * (P3[0] *
      V1[5] - P3[3] * V1[2]) + P2[3] * (P3[1] * V1[2] - P3[0] * V1[3]))) +
      P1[2] * (P2[1] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] * (P3[3] *
      V1[3] - P3[1] * V1[5]) + P2[3] * (P3[1] * V1[4] - P3[2] * V1[3])))) +
      TMP26 * (P2[0] * (P1[0] * (P3[3] * V2[3] - P3[1] * V2[5]) + (P1[1] *
      (P3[0] * V2[5] - P3[3] * V2[2]) + P1[3] * (P3[1] * V2[2] - P3[0] *
      V2[3]))) + P2[2] * (P1[1] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P1[2] *
      (P3[3] * V2[3] - P3[1] * V2[5]) + P1[3] * (P3[1] * V2[4] - P3[2] *
      V2[3]))))));
  T3[5] = denom * cI * (OM3 * (P3[0] * (TMP23 * (TMP1 * (P2[0] * (P3[2] * V1[3]
      - P3[1] * V1[4]) + (P2[1] * (P3[0] * V1[4] - P3[2] * V1[2]) + P2[2] *
      (P3[1] * V1[2] - P3[0] * V1[3]))) + 0.666666667 * (P3[3] * TMP33)) +
      TMP26 * (TMP2 * (P1[0] * (P3[2] * V2[3] - P3[1] * V2[4]) + (P1[1] *
      (P3[0] * V2[4] - P3[2] * V2[2]) + P1[2] * (P3[1] * V2[2] - P3[0] *
      V2[3]))) + 0.666666667 * (P3[3] * TMP34))) + P3[3] * (TMP1 * TMP23 *
      (P2[1] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[1] * V1[5] -
      P3[3] * V1[3]) + P2[3] * (P3[2] * V1[3] - P3[1] * V1[4]))) + TMP2 * TMP26
      * (P1[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[1] * V2[5] -
      P3[3] * V2[3]) + P1[3] * (P3[2] * V2[3] - P3[1] * V2[4]))))) + (TMP23 *
      (P1[0] * (P2[0] * (P3[1] * V1[4] - P3[2] * V1[3]) + (P2[1] * (P3[2] *
      V1[2] - P3[0] * V1[4]) + P2[2] * (P3[0] * V1[3] - P3[1] * V1[2]))) +
      P1[3] * (P2[1] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] * (P3[3] *
      V1[3] - P3[1] * V1[5]) + P2[3] * (P3[1] * V1[4] - P3[2] * V1[3])))) +
      TMP26 * (P2[0] * (P1[0] * (P3[1] * V2[4] - P3[2] * V2[3]) + (P1[1] *
      (P3[2] * V2[2] - P3[0] * V2[4]) + P1[2] * (P3[0] * V2[3] - P3[1] *
      V2[2]))) + P2[3] * (P1[1] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P1[2] *
      (P3[3] * V2[3] - P3[1] * V2[5]) + P1[3] * (P3[1] * V2[4] - P3[2] *
      V2[3]))))));
  T3[6] = denom * cI * (OM3 * (P3[0] * (TMP23 * (TMP1 * (P2[0] * (P3[3] * V1[4]
      - P3[2] * V1[5]) + (P2[2] * (P3[0] * V1[5] - P3[3] * V1[2]) + P2[3] *
      (P3[2] * V1[2] - P3[0] * V1[4]))) + 0.666666667 * (P3[1] * TMP33)) +
      TMP26 * (TMP2 * (P1[0] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P1[2] *
      (P3[0] * V2[5] - P3[3] * V2[2]) + P1[3] * (P3[2] * V2[2] - P3[0] *
      V2[4]))) + 0.666666667 * (P3[1] * TMP34))) + P3[1] * (TMP1 * TMP23 *
      (P2[1] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[1] * V1[5] -
      P3[3] * V1[3]) + P2[3] * (P3[2] * V1[3] - P3[1] * V1[4]))) + TMP2 * TMP26
      * (P1[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[1] * V2[5] -
      P3[3] * V2[3]) + P1[3] * (P3[2] * V2[3] - P3[1] * V2[4]))))) + (TMP23 *
      (P1[0] * (P2[0] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] * (P3[3] *
      V1[2] - P3[0] * V1[5]) + P2[3] * (P3[0] * V1[4] - P3[2] * V1[2]))) +
      P1[1] * (P2[1] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] * (P3[3] *
      V1[3] - P3[1] * V1[5]) + P2[3] * (P3[1] * V1[4] - P3[2] * V1[3])))) +
      TMP26 * (P2[0] * (P1[0] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P1[2] *
      (P3[3] * V2[2] - P3[0] * V2[5]) + P1[3] * (P3[0] * V2[4] - P3[2] *
      V2[2]))) + P2[1] * (P1[1] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P1[2] *
      (P3[3] * V2[3] - P3[1] * V2[5]) + P1[3] * (P3[1] * V2[4] - P3[2] *
      V2[3]))))));
  T3[7] = denom * 2. * cI * (TMP23 * (OM3 * P3[1] * (TMP1 * (P2[0] * (P3[3] *
      V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[0] * V1[5] - P3[3] * V1[2]) + P2[3]
      * (P3[2] * V1[2] - P3[0] * V1[4]))) + 0.333333333 * (P3[1] * TMP33)) +
      (P1[1] * (P2[0] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] * (P3[3] *
      V1[2] - P3[0] * V1[5]) + P2[3] * (P3[0] * V1[4] - P3[2] * V1[2]))) +
      0.333333333 * (TMP33))) + TMP26 * (OM3 * P3[1] * (TMP2 * (P1[0] * (P3[3]
      * V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[0] * V2[5] - P3[3] * V2[2]) +
      P1[3] * (P3[2] * V2[2] - P3[0] * V2[4]))) + 0.333333333 * (P3[1] *
      TMP34)) + (P2[1] * (P1[0] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P1[2] *
      (P3[3] * V2[2] - P3[0] * V2[5]) + P1[3] * (P3[0] * V2[4] - P3[2] *
      V2[2]))) + 0.333333333 * (TMP34))));
  T3[8] = denom * cI * (OM3 * (P3[1] * (TMP23 * (TMP1 * (P2[0] * (P3[1] * V1[5]
      - P3[3] * V1[3]) + (P2[1] * (P3[3] * V1[2] - P3[0] * V1[5]) + P2[3] *
      (P3[0] * V1[3] - P3[1] * V1[2]))) + 0.666666667 * (P3[2] * TMP33)) +
      TMP26 * (TMP2 * (P1[0] * (P3[1] * V2[5] - P3[3] * V2[3]) + (P1[1] *
      (P3[3] * V2[2] - P3[0] * V2[5]) + P1[3] * (P3[0] * V2[3] - P3[1] *
      V2[2]))) + 0.666666667 * (P3[2] * TMP34))) + P3[2] * (TMP1 * TMP23 *
      (P2[0] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[0] * V1[5] -
      P3[3] * V1[2]) + P2[3] * (P3[2] * V1[2] - P3[0] * V1[4]))) + TMP2 * TMP26
      * (P1[0] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[0] * V2[5] -
      P3[3] * V2[2]) + P1[3] * (P3[2] * V2[2] - P3[0] * V2[4]))))) + (TMP23 *
      (P1[1] * (P2[0] * (P3[3] * V1[3] - P3[1] * V1[5]) + (P2[1] * (P3[0] *
      V1[5] - P3[3] * V1[2]) + P2[3] * (P3[1] * V1[2] - P3[0] * V1[3]))) +
      P1[2] * (P2[0] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] * (P3[3] *
      V1[2] - P3[0] * V1[5]) + P2[3] * (P3[0] * V1[4] - P3[2] * V1[2])))) +
      TMP26 * (P2[1] * (P1[0] * (P3[3] * V2[3] - P3[1] * V2[5]) + (P1[1] *
      (P3[0] * V2[5] - P3[3] * V2[2]) + P1[3] * (P3[1] * V2[2] - P3[0] *
      V2[3]))) + P2[2] * (P1[0] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P1[2] *
      (P3[3] * V2[2] - P3[0] * V2[5]) + P1[3] * (P3[0] * V2[4] - P3[2] *
      V2[2]))))));
  T3[9] = denom * cI * (OM3 * (P3[1] * (TMP23 * (TMP1 * (P2[0] * (P3[2] * V1[3]
      - P3[1] * V1[4]) + (P2[1] * (P3[0] * V1[4] - P3[2] * V1[2]) + P2[2] *
      (P3[1] * V1[2] - P3[0] * V1[3]))) + 0.666666667 * (P3[3] * TMP33)) +
      TMP26 * (TMP2 * (P1[0] * (P3[2] * V2[3] - P3[1] * V2[4]) + (P1[1] *
      (P3[0] * V2[4] - P3[2] * V2[2]) + P1[2] * (P3[1] * V2[2] - P3[0] *
      V2[3]))) + 0.666666667 * (P3[3] * TMP34))) + P3[3] * (TMP1 * TMP23 *
      (P2[0] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[0] * V1[5] -
      P3[3] * V1[2]) + P2[3] * (P3[2] * V1[2] - P3[0] * V1[4]))) + TMP2 * TMP26
      * (P1[0] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[0] * V2[5] -
      P3[3] * V2[2]) + P1[3] * (P3[2] * V2[2] - P3[0] * V2[4]))))) + (TMP23 *
      (P1[1] * (P2[0] * (P3[1] * V1[4] - P3[2] * V1[3]) + (P2[1] * (P3[2] *
      V1[2] - P3[0] * V1[4]) + P2[2] * (P3[0] * V1[3] - P3[1] * V1[2]))) +
      P1[3] * (P2[0] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] * (P3[3] *
      V1[2] - P3[0] * V1[5]) + P2[3] * (P3[0] * V1[4] - P3[2] * V1[2])))) +
      TMP26 * (P2[1] * (P1[0] * (P3[1] * V2[4] - P3[2] * V2[3]) + (P1[1] *
      (P3[2] * V2[2] - P3[0] * V2[4]) + P1[2] * (P3[0] * V2[3] - P3[1] *
      V2[2]))) + P2[3] * (P1[0] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P1[2] *
      (P3[3] * V2[2] - P3[0] * V2[5]) + P1[3] * (P3[0] * V2[4] - P3[2] *
      V2[2]))))));
  T3[10] = denom * cI * (OM3 * (P3[0] * (TMP23 * (TMP1 * (P2[0] * (P3[1] *
      V1[5] - P3[3] * V1[3]) + (P2[1] * (P3[3] * V1[2] - P3[0] * V1[5]) + P2[3]
      * (P3[0] * V1[3] - P3[1] * V1[2]))) + 0.666666667 * (P3[2] * TMP33)) +
      TMP26 * (TMP2 * (P1[0] * (P3[1] * V2[5] - P3[3] * V2[3]) + (P1[1] *
      (P3[3] * V2[2] - P3[0] * V2[5]) + P1[3] * (P3[0] * V2[3] - P3[1] *
      V2[2]))) + 0.666666667 * (P3[2] * TMP34))) + P3[2] * (TMP1 * TMP23 *
      (P2[1] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[1] * V1[5] -
      P3[3] * V1[3]) + P2[3] * (P3[2] * V1[3] - P3[1] * V1[4]))) + TMP2 * TMP26
      * (P1[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[1] * V2[5] -
      P3[3] * V2[3]) + P1[3] * (P3[2] * V2[3] - P3[1] * V2[4]))))) + (TMP23 *
      (P1[0] * (P2[0] * (P3[3] * V1[3] - P3[1] * V1[5]) + (P2[1] * (P3[0] *
      V1[5] - P3[3] * V1[2]) + P2[3] * (P3[1] * V1[2] - P3[0] * V1[3]))) +
      P1[2] * (P2[1] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] * (P3[3] *
      V1[3] - P3[1] * V1[5]) + P2[3] * (P3[1] * V1[4] - P3[2] * V1[3])))) +
      TMP26 * (P2[0] * (P1[0] * (P3[3] * V2[3] - P3[1] * V2[5]) + (P1[1] *
      (P3[0] * V2[5] - P3[3] * V2[2]) + P1[3] * (P3[1] * V2[2] - P3[0] *
      V2[3]))) + P2[2] * (P1[1] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P1[2] *
      (P3[3] * V2[3] - P3[1] * V2[5]) + P1[3] * (P3[1] * V2[4] - P3[2] *
      V2[3]))))));
  T3[11] = denom * cI * (OM3 * (P3[1] * (TMP23 * (TMP1 * (P2[0] * (P3[1] *
      V1[5] - P3[3] * V1[3]) + (P2[1] * (P3[3] * V1[2] - P3[0] * V1[5]) + P2[3]
      * (P3[0] * V1[3] - P3[1] * V1[2]))) + 0.666666667 * (P3[2] * TMP33)) +
      TMP26 * (TMP2 * (P1[0] * (P3[1] * V2[5] - P3[3] * V2[3]) + (P1[1] *
      (P3[3] * V2[2] - P3[0] * V2[5]) + P1[3] * (P3[0] * V2[3] - P3[1] *
      V2[2]))) + 0.666666667 * (P3[2] * TMP34))) + P3[2] * (TMP1 * TMP23 *
      (P2[0] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[0] * V1[5] -
      P3[3] * V1[2]) + P2[3] * (P3[2] * V1[2] - P3[0] * V1[4]))) + TMP2 * TMP26
      * (P1[0] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[0] * V2[5] -
      P3[3] * V2[2]) + P1[3] * (P3[2] * V2[2] - P3[0] * V2[4]))))) + (TMP23 *
      (P1[1] * (P2[0] * (P3[3] * V1[3] - P3[1] * V1[5]) + (P2[1] * (P3[0] *
      V1[5] - P3[3] * V1[2]) + P2[3] * (P3[1] * V1[2] - P3[0] * V1[3]))) +
      P1[2] * (P2[0] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] * (P3[3] *
      V1[2] - P3[0] * V1[5]) + P2[3] * (P3[0] * V1[4] - P3[2] * V1[2])))) +
      TMP26 * (P2[1] * (P1[0] * (P3[3] * V2[3] - P3[1] * V2[5]) + (P1[1] *
      (P3[0] * V2[5] - P3[3] * V2[2]) + P1[3] * (P3[1] * V2[2] - P3[0] *
      V2[3]))) + P2[2] * (P1[0] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P1[2] *
      (P3[3] * V2[2] - P3[0] * V2[5]) + P1[3] * (P3[0] * V2[4] - P3[2] *
      V2[2]))))));
  T3[12] = denom * 2. * cI * (TMP23 * (OM3 * P3[2] * (TMP1 * (P2[0] * (P3[1] *
      V1[5] - P3[3] * V1[3]) + (P2[1] * (P3[3] * V1[2] - P3[0] * V1[5]) + P2[3]
      * (P3[0] * V1[3] - P3[1] * V1[2]))) + 0.333333333 * (P3[2] * TMP33)) +
      (P1[2] * (P2[0] * (P3[3] * V1[3] - P3[1] * V1[5]) + (P2[1] * (P3[0] *
      V1[5] - P3[3] * V1[2]) + P2[3] * (P3[1] * V1[2] - P3[0] * V1[3]))) +
      0.333333333 * (TMP33))) + TMP26 * (OM3 * P3[2] * (TMP2 * (P1[0] * (P3[1]
      * V2[5] - P3[3] * V2[3]) + (P1[1] * (P3[3] * V2[2] - P3[0] * V2[5]) +
      P1[3] * (P3[0] * V2[3] - P3[1] * V2[2]))) + 0.333333333 * (P3[2] *
      TMP34)) + (P2[2] * (P1[0] * (P3[3] * V2[3] - P3[1] * V2[5]) + (P1[1] *
      (P3[0] * V2[5] - P3[3] * V2[2]) + P1[3] * (P3[1] * V2[2] - P3[0] *
      V2[3]))) + 0.333333333 * (TMP34))));
  T3[13] = denom * cI * (OM3 * (P3[2] * (TMP23 * (TMP1 * (P2[0] * (P3[2] *
      V1[3] - P3[1] * V1[4]) + (P2[1] * (P3[0] * V1[4] - P3[2] * V1[2]) + P2[2]
      * (P3[1] * V1[2] - P3[0] * V1[3]))) + 0.666666667 * (P3[3] * TMP33)) +
      TMP26 * (TMP2 * (P1[0] * (P3[2] * V2[3] - P3[1] * V2[4]) + (P1[1] *
      (P3[0] * V2[4] - P3[2] * V2[2]) + P1[2] * (P3[1] * V2[2] - P3[0] *
      V2[3]))) + 0.666666667 * (P3[3] * TMP34))) + P3[3] * (TMP1 * TMP23 *
      (P2[0] * (P3[1] * V1[5] - P3[3] * V1[3]) + (P2[1] * (P3[3] * V1[2] -
      P3[0] * V1[5]) + P2[3] * (P3[0] * V1[3] - P3[1] * V1[2]))) + TMP2 * TMP26
      * (P1[0] * (P3[1] * V2[5] - P3[3] * V2[3]) + (P1[1] * (P3[3] * V2[2] -
      P3[0] * V2[5]) + P1[3] * (P3[0] * V2[3] - P3[1] * V2[2]))))) + (TMP23 *
      (P1[2] * (P2[0] * (P3[1] * V1[4] - P3[2] * V1[3]) + (P2[1] * (P3[2] *
      V1[2] - P3[0] * V1[4]) + P2[2] * (P3[0] * V1[3] - P3[1] * V1[2]))) +
      P1[3] * (P2[0] * (P3[3] * V1[3] - P3[1] * V1[5]) + (P2[1] * (P3[0] *
      V1[5] - P3[3] * V1[2]) + P2[3] * (P3[1] * V1[2] - P3[0] * V1[3])))) +
      TMP26 * (P2[2] * (P1[0] * (P3[1] * V2[4] - P3[2] * V2[3]) + (P1[1] *
      (P3[2] * V2[2] - P3[0] * V2[4]) + P1[2] * (P3[0] * V2[3] - P3[1] *
      V2[2]))) + P2[3] * (P1[0] * (P3[3] * V2[3] - P3[1] * V2[5]) + (P1[1] *
      (P3[0] * V2[5] - P3[3] * V2[2]) + P1[3] * (P3[1] * V2[2] - P3[0] *
      V2[3]))))));
  T3[14] = denom * cI * (OM3 * (P3[0] * (TMP23 * (TMP1 * (P2[0] * (P3[2] *
      V1[3] - P3[1] * V1[4]) + (P2[1] * (P3[0] * V1[4] - P3[2] * V1[2]) + P2[2]
      * (P3[1] * V1[2] - P3[0] * V1[3]))) + 0.666666667 * (P3[3] * TMP33)) +
      TMP26 * (TMP2 * (P1[0] * (P3[2] * V2[3] - P3[1] * V2[4]) + (P1[1] *
      (P3[0] * V2[4] - P3[2] * V2[2]) + P1[2] * (P3[1] * V2[2] - P3[0] *
      V2[3]))) + 0.666666667 * (P3[3] * TMP34))) + P3[3] * (TMP1 * TMP23 *
      (P2[1] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[1] * V1[5] -
      P3[3] * V1[3]) + P2[3] * (P3[2] * V1[3] - P3[1] * V1[4]))) + TMP2 * TMP26
      * (P1[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[1] * V2[5] -
      P3[3] * V2[3]) + P1[3] * (P3[2] * V2[3] - P3[1] * V2[4]))))) + (TMP23 *
      (P1[0] * (P2[0] * (P3[1] * V1[4] - P3[2] * V1[3]) + (P2[1] * (P3[2] *
      V1[2] - P3[0] * V1[4]) + P2[2] * (P3[0] * V1[3] - P3[1] * V1[2]))) +
      P1[3] * (P2[1] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] * (P3[3] *
      V1[3] - P3[1] * V1[5]) + P2[3] * (P3[1] * V1[4] - P3[2] * V1[3])))) +
      TMP26 * (P2[0] * (P1[0] * (P3[1] * V2[4] - P3[2] * V2[3]) + (P1[1] *
      (P3[2] * V2[2] - P3[0] * V2[4]) + P1[2] * (P3[0] * V2[3] - P3[1] *
      V2[2]))) + P2[3] * (P1[1] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P1[2] *
      (P3[3] * V2[3] - P3[1] * V2[5]) + P1[3] * (P3[1] * V2[4] - P3[2] *
      V2[3]))))));
  T3[15] = denom * cI * (OM3 * (P3[1] * (TMP23 * (TMP1 * (P2[0] * (P3[2] *
      V1[3] - P3[1] * V1[4]) + (P2[1] * (P3[0] * V1[4] - P3[2] * V1[2]) + P2[2]
      * (P3[1] * V1[2] - P3[0] * V1[3]))) + 0.666666667 * (P3[3] * TMP33)) +
      TMP26 * (TMP2 * (P1[0] * (P3[2] * V2[3] - P3[1] * V2[4]) + (P1[1] *
      (P3[0] * V2[4] - P3[2] * V2[2]) + P1[2] * (P3[1] * V2[2] - P3[0] *
      V2[3]))) + 0.666666667 * (P3[3] * TMP34))) + P3[3] * (TMP1 * TMP23 *
      (P2[0] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[0] * V1[5] -
      P3[3] * V1[2]) + P2[3] * (P3[2] * V1[2] - P3[0] * V1[4]))) + TMP2 * TMP26
      * (P1[0] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[0] * V2[5] -
      P3[3] * V2[2]) + P1[3] * (P3[2] * V2[2] - P3[0] * V2[4]))))) + (TMP23 *
      (P1[1] * (P2[0] * (P3[1] * V1[4] - P3[2] * V1[3]) + (P2[1] * (P3[2] *
      V1[2] - P3[0] * V1[4]) + P2[2] * (P3[0] * V1[3] - P3[1] * V1[2]))) +
      P1[3] * (P2[0] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] * (P3[3] *
      V1[2] - P3[0] * V1[5]) + P2[3] * (P3[0] * V1[4] - P3[2] * V1[2])))) +
      TMP26 * (P2[1] * (P1[0] * (P3[1] * V2[4] - P3[2] * V2[3]) + (P1[1] *
      (P3[2] * V2[2] - P3[0] * V2[4]) + P1[2] * (P3[0] * V2[3] - P3[1] *
      V2[2]))) + P2[3] * (P1[0] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P1[2] *
      (P3[3] * V2[2] - P3[0] * V2[5]) + P1[3] * (P3[0] * V2[4] - P3[2] *
      V2[2]))))));
  T3[16] = denom * cI * (OM3 * (P3[2] * (TMP23 * (TMP1 * (P2[0] * (P3[2] *
      V1[3] - P3[1] * V1[4]) + (P2[1] * (P3[0] * V1[4] - P3[2] * V1[2]) + P2[2]
      * (P3[1] * V1[2] - P3[0] * V1[3]))) + 0.666666667 * (P3[3] * TMP33)) +
      TMP26 * (TMP2 * (P1[0] * (P3[2] * V2[3] - P3[1] * V2[4]) + (P1[1] *
      (P3[0] * V2[4] - P3[2] * V2[2]) + P1[2] * (P3[1] * V2[2] - P3[0] *
      V2[3]))) + 0.666666667 * (P3[3] * TMP34))) + P3[3] * (TMP1 * TMP23 *
      (P2[0] * (P3[1] * V1[5] - P3[3] * V1[3]) + (P2[1] * (P3[3] * V1[2] -
      P3[0] * V1[5]) + P2[3] * (P3[0] * V1[3] - P3[1] * V1[2]))) + TMP2 * TMP26
      * (P1[0] * (P3[1] * V2[5] - P3[3] * V2[3]) + (P1[1] * (P3[3] * V2[2] -
      P3[0] * V2[5]) + P1[3] * (P3[0] * V2[3] - P3[1] * V2[2]))))) + (TMP23 *
      (P1[2] * (P2[0] * (P3[1] * V1[4] - P3[2] * V1[3]) + (P2[1] * (P3[2] *
      V1[2] - P3[0] * V1[4]) + P2[2] * (P3[0] * V1[3] - P3[1] * V1[2]))) +
      P1[3] * (P2[0] * (P3[3] * V1[3] - P3[1] * V1[5]) + (P2[1] * (P3[0] *
      V1[5] - P3[3] * V1[2]) + P2[3] * (P3[1] * V1[2] - P3[0] * V1[3])))) +
      TMP26 * (P2[2] * (P1[0] * (P3[1] * V2[4] - P3[2] * V2[3]) + (P1[1] *
      (P3[2] * V2[2] - P3[0] * V2[4]) + P1[2] * (P3[0] * V2[3] - P3[1] *
      V2[2]))) + P2[3] * (P1[0] * (P3[3] * V2[3] - P3[1] * V2[5]) + (P1[1] *
      (P3[0] * V2[5] - P3[3] * V2[2]) + P1[3] * (P3[1] * V2[2] - P3[0] *
      V2[3]))))));
  T3[17] = denom * 2. * cI * (TMP23 * (OM3 * P3[3] * (TMP1 * (P2[0] * (P3[2] *
      V1[3] - P3[1] * V1[4]) + (P2[1] * (P3[0] * V1[4] - P3[2] * V1[2]) + P2[2]
      * (P3[1] * V1[2] - P3[0] * V1[3]))) + 0.333333333 * (P3[3] * TMP33)) +
      (P1[3] * (P2[0] * (P3[1] * V1[4] - P3[2] * V1[3]) + (P2[1] * (P3[2] *
      V1[2] - P3[0] * V1[4]) + P2[2] * (P3[0] * V1[3] - P3[1] * V1[2]))) +
      0.333333333 * (TMP33))) + TMP26 * (OM3 * P3[3] * (TMP2 * (P1[0] * (P3[2]
      * V2[3] - P3[1] * V2[4]) + (P1[1] * (P3[0] * V2[4] - P3[2] * V2[2]) +
      P1[2] * (P3[1] * V2[2] - P3[0] * V2[3]))) + 0.333333333 * (P3[3] *
      TMP34)) + (P2[3] * (P1[0] * (P3[1] * V2[4] - P3[2] * V2[3]) + (P1[1] *
      (P3[2] * V2[2] - P3[0] * V2[4]) + P1[2] * (P3[0] * V2[3] - P3[1] *
      V2[2]))) + 0.333333333 * (TMP34))));
}


void VVT11_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP36; 
  complex<double> TMP35; 
  TMP36 = (V1[2] * - 1. * (V2[3] * T3[3] + V2[4] * T3[4] + V2[5] * T3[5] -
      V2[2] * T3[2]) + (V1[3] * (V2[3] * T3[7] + V2[4] * T3[8] + V2[5] * T3[9]
      - V2[2] * T3[6]) + (V1[4] * (V2[3] * T3[11] + V2[4] * T3[12] + V2[5] *
      T3[13] - V2[2] * T3[10]) + V1[5] * (V2[3] * T3[15] + V2[4] * T3[16] +
      V2[5] * T3[17] - V2[2] * T3[14]))));
  TMP35 = (V1[2] * - 1. * (V2[3] * T3[6] + V2[4] * T3[10] + V2[5] * T3[14] -
      V2[2] * T3[2]) + (V1[3] * (V2[3] * T3[7] + V2[4] * T3[11] + V2[5] *
      T3[15] - V2[2] * T3[3]) + (V1[4] * (V2[3] * T3[8] + V2[4] * T3[12] +
      V2[5] * T3[16] - V2[2] * T3[4]) + V1[5] * (V2[3] * T3[9] + V2[4] * T3[13]
      + V2[5] * T3[17] - V2[2] * T3[5]))));
  vertex = COUP * - 1. * (+cI * (TMP35 + TMP36)); 
}


void FFT4_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP2; 
  complex<double> TMP1; 
  double P1[4]; 
  complex<double> TMP10; 
  double P2[4]; 
  double P3[4]; 
  double OM3; 
  complex<double> denom; 
  complex<double> TMP19; 
  P1[0] = F1[0].real(); 
  P1[1] = F1[1].real(); 
  P1[2] = F1[1].imag(); 
  P1[3] = F1[0].imag(); 
  P2[0] = F2[0].real(); 
  P2[1] = F2[1].real(); 
  P2[2] = F2[1].imag(); 
  P2[3] = F2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +F1[0] + F2[0]; 
  T3[1] = +F1[1] + F2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP19 = (F2[2] * F1[2] + F2[3] * F1[3] + F2[4] * F1[4] + F2[5] * F1[5]); 
  TMP1 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP10 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP2 = (P3[0] * P2[0] - P3[1] * P2[1] - P3[2] * P2[2] - P3[3] * P2[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * 2. * TMP19 * (OM3 * (P3[0] * (P3[0] * 0.333333333 * (+cI *
      (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (TMP1 * P2[0] + P1[0] *
      TMP2))) + 0.333333333 * cI * (TMP1 * TMP2)) + (-0.333333333 * cI *
      (TMP10) + cI * (P1[0] * P2[0])));
  T3[3] = denom * TMP19 * (OM3 * (P3[0] * (P3[1] * 0.666666667 * (+cI * (TMP10)
      + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (P1[1] * TMP2 + TMP1 * P2[1])))
      - P3[1] * (+cI * (TMP1 * P2[0] + P1[0] * TMP2))) + (+cI * (P1[1] * P2[0]
      + P1[0] * P2[1])));
  T3[4] = denom * TMP19 * (OM3 * (P3[0] * (P3[2] * 0.666666667 * (+cI * (TMP10)
      + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (P1[2] * TMP2 + TMP1 * P2[2])))
      - P3[2] * (+cI * (TMP1 * P2[0] + P1[0] * TMP2))) + (+cI * (P1[2] * P2[0]
      + P1[0] * P2[2])));
  T3[5] = denom * TMP19 * (OM3 * (P3[0] * (P3[3] * 0.666666667 * (+cI * (TMP10)
      + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (P1[3] * TMP2 + TMP1 * P2[3])))
      - P3[3] * (+cI * (TMP1 * P2[0] + P1[0] * TMP2))) + (+cI * (P1[3] * P2[0]
      + P1[0] * P2[3])));
  T3[6] = denom * TMP19 * (OM3 * (P3[0] * (P3[1] * 0.666666667 * (+cI * (TMP10)
      + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (TMP1 * P2[1] + P1[1] * TMP2)))
      - P3[1] * (+cI * (P1[0] * TMP2 + TMP1 * P2[0]))) + (+cI * (P1[0] * P2[1]
      + P1[1] * P2[0])));
  T3[7] = denom * 2. * TMP19 * (OM3 * (P3[1] * (P3[1] * 0.333333333 * (+cI *
      (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (TMP1 * P2[1] + P1[1] *
      TMP2))) - 0.333333333 * cI * (TMP1 * TMP2)) + (+cI * (P1[1] * P2[1]) +
      0.333333333 * cI * (TMP10)));
  T3[8] = denom * TMP19 * (OM3 * (P3[1] * (P3[2] * 0.666666667 * (+cI * (TMP10)
      + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (P1[2] * TMP2 + TMP1 * P2[2])))
      - P3[2] * (+cI * (TMP1 * P2[1] + P1[1] * TMP2))) + (+cI * (P1[2] * P2[1]
      + P1[1] * P2[2])));
  T3[9] = denom * TMP19 * (OM3 * (P3[1] * (P3[3] * 0.666666667 * (+cI * (TMP10)
      + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (P1[3] * TMP2 + TMP1 * P2[3])))
      - P3[3] * (+cI * (TMP1 * P2[1] + P1[1] * TMP2))) + (+cI * (P1[3] * P2[1]
      + P1[1] * P2[3])));
  T3[10] = denom * TMP19 * (OM3 * (P3[0] * (P3[2] * 0.666666667 * (+cI *
      (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (TMP1 * P2[2] + P1[2] *
      TMP2))) - P3[2] * (+cI * (P1[0] * TMP2 + TMP1 * P2[0]))) + (+cI * (P1[0]
      * P2[2] + P1[2] * P2[0])));
  T3[11] = denom * TMP19 * (OM3 * (P3[1] * (P3[2] * 0.666666667 * (+cI *
      (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (TMP1 * P2[2] + P1[2] *
      TMP2))) - P3[2] * (+cI * (P1[1] * TMP2 + TMP1 * P2[1]))) + (+cI * (P1[1]
      * P2[2] + P1[2] * P2[1])));
  T3[12] = denom * 2. * TMP19 * (OM3 * (P3[2] * (P3[2] * 0.333333333 * (+cI *
      (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (TMP1 * P2[2] + P1[2] *
      TMP2))) - 0.333333333 * cI * (TMP1 * TMP2)) + (+cI * (P1[2] * P2[2]) +
      0.333333333 * cI * (TMP10)));
  T3[13] = denom * TMP19 * (OM3 * (P3[2] * (P3[3] * 0.666666667 * (+cI *
      (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (P1[3] * TMP2 + TMP1 *
      P2[3]))) - P3[3] * (+cI * (TMP1 * P2[2] + P1[2] * TMP2))) + (+cI * (P1[3]
      * P2[2] + P1[2] * P2[3])));
  T3[14] = denom * TMP19 * (OM3 * (P3[0] * (P3[3] * 0.666666667 * (+cI *
      (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (TMP1 * P2[3] + P1[3] *
      TMP2))) - P3[3] * (+cI * (P1[0] * TMP2 + TMP1 * P2[0]))) + (+cI * (P1[0]
      * P2[3] + P1[3] * P2[0])));
  T3[15] = denom * TMP19 * (OM3 * (P3[1] * (P3[3] * 0.666666667 * (+cI *
      (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (TMP1 * P2[3] + P1[3] *
      TMP2))) - P3[3] * (+cI * (P1[1] * TMP2 + TMP1 * P2[1]))) + (+cI * (P1[1]
      * P2[3] + P1[3] * P2[1])));
  T3[16] = denom * TMP19 * (OM3 * (P3[2] * (P3[3] * 0.666666667 * (+cI *
      (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (TMP1 * P2[3] + P1[3] *
      TMP2))) - P3[3] * (+cI * (P1[2] * TMP2 + TMP1 * P2[2]))) + (+cI * (P1[2]
      * P2[3] + P1[3] * P2[2])));
  T3[17] = denom * 2. * TMP19 * (OM3 * (P3[3] * (P3[3] * 0.333333333 * (+cI *
      (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (TMP1 * P2[3] + P1[3] *
      TMP2))) - 0.333333333 * cI * (TMP1 * TMP2)) + (+cI * (P1[3] * P2[3]) +
      0.333333333 * cI * (TMP10)));
}


void VVT7_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> TMP22; 
  double P2[4]; 
  complex<double> TMP21; 
  complex<double> TMP4; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  TMP4 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP21 = (P1[0] * - 1. * (P2[1] * T3[6] + P2[2] * T3[10] + P2[3] * T3[14] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[11] + P2[3] *
      T3[15] - P2[0] * T3[3]) + (P1[2] * (P2[1] * T3[8] + P2[2] * T3[12] +
      P2[3] * T3[16] - P2[0] * T3[4]) + P1[3] * (P2[1] * T3[9] + P2[2] * T3[13]
      + P2[3] * T3[17] - P2[0] * T3[5]))));
  TMP22 = (P1[0] * - 1. * (P2[1] * T3[3] + P2[2] * T3[4] + P2[3] * T3[5] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[8] + P2[3] * T3[9]
      - P2[0] * T3[6]) + (P1[2] * (P2[1] * T3[11] + P2[2] * T3[12] + P2[3] *
      T3[13] - P2[0] * T3[10]) + P1[3] * (P2[1] * T3[15] + P2[2] * T3[16] +
      P2[3] * T3[17] - P2[0] * T3[14]))));
  vertex = COUP * - TMP4 * (+cI * (TMP21 + TMP22)); 
}


void FFV7_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  complex<double> TMP17; 
  double P3[4]; 
  double OM3; 
  complex<double> TMP18; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP17 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
      F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
  TMP18 = (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) +
      F1[5] * (F2[2] * (+cI * (P3[2]) - P3[1]) + F2[3] * (P3[0] + P3[3])));
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * - 2. * cI * (OM3 * - 0.500000000 * P3[0] * (TMP17 + 2. *
      (TMP18)) + (+0.500000000 * (F2[4] * F1[2] + F2[5] * F1[3]) + F2[2] *
      F1[4] + F2[3] * F1[5]));
  V3[3] = denom * - 2. * cI * (OM3 * - 0.500000000 * P3[1] * (TMP17 + 2. *
      (TMP18)) + (-0.500000000 * (F2[5] * F1[2] + F2[4] * F1[3]) + F2[3] *
      F1[4] + F2[2] * F1[5]));
  V3[4] = denom * 2. * cI * (OM3 * 0.500000000 * P3[2] * (TMP17 + 2. * (TMP18))
      + (+0.500000000 * cI * (F2[5] * F1[2]) - 0.500000000 * cI * (F2[4] *
      F1[3]) - cI * (F2[3] * F1[4]) + cI * (F2[2] * F1[5])));
  V3[5] = denom * 2. * cI * (OM3 * 0.500000000 * P3[3] * (TMP17 + 2. * (TMP18))
      + (+0.500000000 * (F2[4] * F1[2]) - 0.500000000 * (F2[5] * F1[3]) - F2[2]
      * F1[4] + F2[3] * F1[5]));
}


void VVT2_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> TMP22; 
  double P2[4]; 
  complex<double> TMP30; 
  complex<double> TMP21; 
  complex<double> TMP29; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  TMP30 = -1. * (P1[0] * (P2[1] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P2[2] *
      (V2[3] * V1[5] - V2[5] * V1[3]) + P2[3] * (V2[4] * V1[3] - V2[3] *
      V1[4]))) + (P1[1] * (P2[0] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P2[2] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P2[3] * (V2[2] * V1[4] - V2[4] *
      V1[2]))) + (P1[2] * (P2[0] * (V2[5] * V1[3] - V2[3] * V1[5]) + (P2[1] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P2[3] * (V2[3] * V1[2] - V2[2] *
      V1[3]))) + P1[3] * (P2[0] * (V2[3] * V1[4] - V2[4] * V1[3]) + (P2[1] *
      (V2[4] * V1[2] - V2[2] * V1[4]) + P2[2] * (V2[2] * V1[3] - V2[3] *
      V1[2]))))));
  TMP21 = (P1[0] * - 1. * (P2[1] * T3[6] + P2[2] * T3[10] + P2[3] * T3[14] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[11] + P2[3] *
      T3[15] - P2[0] * T3[3]) + (P1[2] * (P2[1] * T3[8] + P2[2] * T3[12] +
      P2[3] * T3[16] - P2[0] * T3[4]) + P1[3] * (P2[1] * T3[9] + P2[2] * T3[13]
      + P2[3] * T3[17] - P2[0] * T3[5]))));
  TMP22 = (P1[0] * - 1. * (P2[1] * T3[3] + P2[2] * T3[4] + P2[3] * T3[5] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[8] + P2[3] * T3[9]
      - P2[0] * T3[6]) + (P1[2] * (P2[1] * T3[11] + P2[2] * T3[12] + P2[3] *
      T3[13] - P2[0] * T3[10]) + P1[3] * (P2[1] * T3[15] + P2[2] * T3[16] +
      P2[3] * T3[17] - P2[0] * T3[14]))));
  TMP29 = -1. * (P1[0] * (P2[1] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P2[2] *
      (V2[5] * V1[3] - V2[3] * V1[5]) + P2[3] * (V2[3] * V1[4] - V2[4] *
      V1[3]))) + (P1[1] * (P2[0] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P2[2] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P2[3] * (V2[4] * V1[2] - V2[2] *
      V1[4]))) + (P1[2] * (P2[0] * (V2[3] * V1[5] - V2[5] * V1[3]) + (P2[1] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P2[3] * (V2[2] * V1[3] - V2[3] *
      V1[2]))) + P1[3] * (P2[0] * (V2[4] * V1[3] - V2[3] * V1[4]) + (P2[1] *
      (V2[2] * V1[4] - V2[4] * V1[2]) + P2[2] * (V2[3] * V1[2] - V2[2] *
      V1[3]))))));
  vertex = COUP * (TMP21 * (-cI * (TMP30) + cI * (TMP29)) + TMP22 * (-cI *
      (TMP30) + cI * (TMP29)));
}


void FFT3_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP2; 
  complex<double> TMP39; 
  complex<double> TMP1; 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  complex<double> TMP40; 
  complex<double> denom; 
  complex<double> TMP11; 
  double OM3; 
  P1[0] = F1[0].real(); 
  P1[1] = F1[1].real(); 
  P1[2] = F1[1].imag(); 
  P1[3] = F1[0].imag(); 
  P2[0] = F2[0].real(); 
  P2[1] = F2[1].real(); 
  P2[2] = F2[1].imag(); 
  P2[3] = F2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +F1[0] + F2[0]; 
  T3[1] = +F1[1] + F2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP1 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP39 = (F1[2] * (F2[4] * - 1. * (P1[0] + P1[3]) - F2[5] * (P1[1] + cI *
      (P1[2]))) + (F1[3] * (F2[4] * (+cI * (P1[2]) - P1[1]) + F2[5] * (P1[3] -
      P1[0])) + (F1[4] * (F2[2] * (P1[0] - P1[3]) - F2[3] * (P1[1] + cI *
      (P1[2]))) + F1[5] * (F2[2] * (+cI * (P1[2]) - P1[1]) + F2[3] * (P1[0] +
      P1[3])))));
  TMP11 = (F1[2] * (F2[4] * - 1. * (P3[0] + P3[3]) - F2[5] * (P3[1] + cI *
      (P3[2]))) + (F1[3] * (F2[4] * (+cI * (P3[2]) - P3[1]) + F2[5] * (P3[3] -
      P3[0])) + (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI *
      (P3[2]))) + F1[5] * (F2[2] * (+cI * (P3[2]) - P3[1]) + F2[3] * (P3[0] +
      P3[3])))));
  TMP40 = (F1[2] * (F2[4] * - 1. * (P2[0] + P2[3]) - F2[5] * (P2[1] + cI *
      (P2[2]))) + (F1[3] * (F2[4] * (+cI * (P2[2]) - P2[1]) + F2[5] * (P2[3] -
      P2[0])) + (F1[4] * (F2[2] * (P2[0] - P2[3]) - F2[3] * (P2[1] + cI *
      (P2[2]))) + F1[5] * (F2[2] * (+cI * (P2[2]) - P2[1]) + F2[3] * (P2[0] +
      P2[3])))));
  TMP2 = (P3[0] * P2[0] - P3[1] * P2[1] - P3[2] * P2[2] - P3[3] * P2[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * 2. * cI * (OM3 * (P3[0] * (TMP1 * (F2[4] * F1[2] + F2[5] *
      F1[3] + 0.666666667 * (P3[0] * OM3 * TMP11) - F2[2] * F1[4] - F2[3] *
      F1[5]) + (TMP2 * - 1. * (F2[4] * F1[2] + F2[5] * F1[3] + 0.666666667 *
      (P3[0] * OM3 * TMP11) - F2[2] * F1[4] - F2[3] * F1[5]) + (P3[0] *
      0.333333333 * (TMP39 - TMP40) + TMP11 * (P2[0] - P1[0])))) + 0.333333333
      * (TMP11 * (TMP1 - TMP2))) + (P1[0] * (F2[2] * F1[4] + F2[3] * F1[5] -
      F2[4] * F1[2] - F2[5] * F1[3]) + (P2[0] * (F2[4] * F1[2] + F2[5] * F1[3]
      - F2[2] * F1[4] - F2[3] * F1[5]) + (-0.333333333 * (TMP39) + 0.333333333
      * (TMP40)))));
  T3[6] = denom * cI * (OM3 * (P3[0] * (TMP1 * - 1. * (F2[5] * F1[2] + F2[4] *
      F1[3] + F2[3] * F1[4] + F2[2] * F1[5] - 1.333333333 * (P3[1] * OM3 *
      TMP11)) + (TMP2 * (F2[5] * F1[2] + F2[4] * F1[3] + F2[3] * F1[4] + F2[2]
      * F1[5] - 1.333333333 * (P3[1] * OM3 * TMP11)) + (P3[1] * 0.666666667 *
      (TMP39 - TMP40) + TMP11 * (P2[1] - P1[1])))) + P3[1] * (TMP1 * (F2[4] *
      F1[2] + F2[5] * F1[3] - F2[2] * F1[4] - F2[3] * F1[5]) + (TMP2 * (F2[2] *
      F1[4] + F2[3] * F1[5] - F2[4] * F1[2] - F2[5] * F1[3]) + TMP11 * (P2[0] -
      P1[0])))) + (F1[2] * (F2[4] * (P2[1] - P1[1]) + F2[5] * (P1[0] - P2[0]))
      + (F1[3] * (F2[4] * (P1[0] - P2[0]) + F2[5] * (P2[1] - P1[1])) + (F1[4] *
      (F2[2] * (P1[1] - P2[1]) + F2[3] * (P1[0] - P2[0])) + F1[5] * (F2[2] *
      (P1[0] - P2[0]) + F2[3] * (P1[1] - P2[1]))))));
  T3[10] = denom * cI * (OM3 * (P3[0] * (TMP1 * (-cI * (F2[5] * F1[2] + F2[3] *
      F1[4]) + cI * (F2[4] * F1[3] + F2[2] * F1[5]) + 1.333333333 * (P3[2] *
      OM3 * TMP11)) + (TMP2 * - 1. * (-cI * (F2[5] * F1[2] + F2[3] * F1[4]) +
      cI * (F2[4] * F1[3] + F2[2] * F1[5]) + 1.333333333 * (P3[2] * OM3 *
      TMP11)) + (P3[2] * 0.666666667 * (TMP39 - TMP40) + TMP11 * (P2[2] -
      P1[2])))) + P3[2] * (TMP1 * (F2[4] * F1[2] + F2[5] * F1[3] - F2[2] *
      F1[4] - F2[3] * F1[5]) + (TMP2 * (F2[2] * F1[4] + F2[3] * F1[5] - F2[4] *
      F1[2] - F2[5] * F1[3]) + TMP11 * (P2[0] - P1[0])))) + (F1[2] * (F2[4] *
      (P2[2] - P1[2]) + F2[5] * (-cI * (P2[0]) + cI * (P1[0]))) + (F1[3] *
      (F2[4] * (-cI * (P1[0]) + cI * (P2[0])) + F2[5] * (P2[2] - P1[2])) +
      (F1[4] * (F2[2] * (P1[2] - P2[2]) + F2[3] * (-cI * (P2[0]) + cI *
      (P1[0]))) + F1[5] * (F2[2] * (-cI * (P1[0]) + cI * (P2[0])) + F2[3] *
      (P1[2] - P2[2]))))));
  T3[14] = denom * cI * (OM3 * (P3[0] * (TMP1 * (F2[5] * F1[3] + F2[3] * F1[5]
      + 1.333333333 * (P3[3] * OM3 * TMP11) - F2[4] * F1[2] - F2[2] * F1[4]) +
      (TMP2 * - 1. * (F2[5] * F1[3] + F2[3] * F1[5] + 1.333333333 * (P3[3] *
      OM3 * TMP11) - F2[4] * F1[2] - F2[2] * F1[4]) + (P3[3] * 0.666666667 *
      (TMP39 - TMP40) + TMP11 * (P2[3] - P1[3])))) + P3[3] * (TMP1 * (F2[4] *
      F1[2] + F2[5] * F1[3] - F2[2] * F1[4] - F2[3] * F1[5]) + (TMP2 * (F2[2] *
      F1[4] + F2[3] * F1[5] - F2[4] * F1[2] - F2[5] * F1[3]) + TMP11 * (P2[0] -
      P1[0])))) + (F1[2] * F2[4] * (P1[0] + P2[3] - P1[3] - P2[0]) + (F1[3] *
      F2[5] * (P2[0] + P2[3] - P1[0] - P1[3]) + (F1[4] * F2[2] * (P1[0] + P1[3]
      - P2[0] - P2[3]) + F1[5] * F2[3] * (P1[3] + P2[0] - P1[0] - P2[3])))));
  T3[3] = denom * cI * (OM3 * (P3[0] * (TMP1 * - 1. * (F2[5] * F1[2] + F2[4] *
      F1[3] + F2[3] * F1[4] + F2[2] * F1[5] - 1.333333333 * (P3[1] * OM3 *
      TMP11)) + (TMP2 * (F2[5] * F1[2] + F2[4] * F1[3] + F2[3] * F1[4] + F2[2]
      * F1[5] - 1.333333333 * (P3[1] * OM3 * TMP11)) + (P3[1] * 0.666666667 *
      (TMP39 - TMP40) + TMP11 * (P2[1] - P1[1])))) + P3[1] * (TMP1 * (F2[4] *
      F1[2] + F2[5] * F1[3] - F2[2] * F1[4] - F2[3] * F1[5]) + (TMP2 * (F2[2] *
      F1[4] + F2[3] * F1[5] - F2[4] * F1[2] - F2[5] * F1[3]) + TMP11 * (P2[0] -
      P1[0])))) + (F1[2] * (F2[4] * (P2[1] - P1[1]) + F2[5] * (P1[0] - P2[0]))
      + (F1[3] * (F2[4] * (P1[0] - P2[0]) + F2[5] * (P2[1] - P1[1])) + (F1[4] *
      (F2[2] * (P1[1] - P2[1]) + F2[3] * (P1[0] - P2[0])) + F1[5] * (F2[2] *
      (P1[0] - P2[0]) + F2[3] * (P1[1] - P2[1]))))));
  T3[7] = denom * 2. * cI * (OM3 * (P3[1] * (TMP1 * - 1. * (F2[5] * F1[2] +
      F2[4] * F1[3] + F2[3] * F1[4] + F2[2] * F1[5] - 0.666666667 * (P3[1] *
      OM3 * TMP11)) + (TMP2 * (F2[5] * F1[2] + F2[4] * F1[3] + F2[3] * F1[4] +
      F2[2] * F1[5] - 0.666666667 * (P3[1] * OM3 * TMP11)) + (P3[1] *
      0.333333333 * (TMP39 - TMP40) + TMP11 * (P2[1] - P1[1])))) + 0.333333333
      * (TMP11 * (TMP2 - TMP1))) + (P1[1] * (F2[5] * F1[2] + F2[4] * F1[3] +
      F2[3] * F1[4] + F2[2] * F1[5]) + (P2[1] * - 1. * (F2[5] * F1[2] + F2[4] *
      F1[3] + F2[3] * F1[4] + F2[2] * F1[5]) + (-0.333333333 * (TMP40) +
      0.333333333 * (TMP39)))));
  T3[11] = denom * cI * (OM3 * (P3[1] * (TMP1 * (-cI * (F2[5] * F1[2] + F2[3] *
      F1[4]) + cI * (F2[4] * F1[3] + F2[2] * F1[5]) + 1.333333333 * (P3[2] *
      OM3 * TMP11)) + (TMP2 * - 1. * (-cI * (F2[5] * F1[2] + F2[3] * F1[4]) +
      cI * (F2[4] * F1[3] + F2[2] * F1[5]) + 1.333333333 * (P3[2] * OM3 *
      TMP11)) + (P3[2] * 0.666666667 * (TMP39 - TMP40) + TMP11 * (P2[2] -
      P1[2])))) + P3[2] * (TMP1 * - 1. * (F2[5] * F1[2] + F2[4] * F1[3] + F2[3]
      * F1[4] + F2[2] * F1[5]) + (TMP2 * (F2[5] * F1[2] + F2[4] * F1[3] + F2[3]
      * F1[4] + F2[2] * F1[5]) + TMP11 * (P2[1] - P1[1])))) + (F1[2] * F2[5] *
      (P1[2] - cI * (P2[1]) + cI * (P1[1]) - P2[2]) + (F1[3] * F2[4] * (P1[2] -
      cI * (P1[1]) + cI * (P2[1]) - P2[2]) + (F1[4] * F2[3] * (P1[2] - cI *
      (P2[1]) + cI * (P1[1]) - P2[2]) + F1[5] * F2[2] * (P1[2] - cI * (P1[1]) +
      cI * (P2[1]) - P2[2])))));
  T3[15] = denom * cI * (OM3 * (P3[1] * (TMP1 * (F2[5] * F1[3] + F2[3] * F1[5]
      + 1.333333333 * (P3[3] * OM3 * TMP11) - F2[4] * F1[2] - F2[2] * F1[4]) +
      (TMP2 * - 1. * (F2[5] * F1[3] + F2[3] * F1[5] + 1.333333333 * (P3[3] *
      OM3 * TMP11) - F2[4] * F1[2] - F2[2] * F1[4]) + (P3[3] * 0.666666667 *
      (TMP39 - TMP40) + TMP11 * (P2[3] - P1[3])))) + P3[3] * (TMP1 * - 1. *
      (F2[5] * F1[2] + F2[4] * F1[3] + F2[3] * F1[4] + F2[2] * F1[5]) + (TMP2 *
      (F2[5] * F1[2] + F2[4] * F1[3] + F2[3] * F1[4] + F2[2] * F1[5]) + TMP11 *
      (P2[1] - P1[1])))) + (F1[2] * (F2[4] * (P1[1] - P2[1]) + F2[5] * (P1[3] -
      P2[3])) + (F1[3] * (F2[4] * (P1[3] - P2[3]) + F2[5] * (P2[1] - P1[1])) +
      (F1[4] * (F2[2] * (P1[1] - P2[1]) + F2[3] * (P1[3] - P2[3])) + F1[5] *
      (F2[2] * (P1[3] - P2[3]) + F2[3] * (P2[1] - P1[1]))))));
  T3[4] = denom * cI * (OM3 * (P3[0] * (TMP1 * (-cI * (F2[5] * F1[2] + F2[3] *
      F1[4]) + cI * (F2[4] * F1[3] + F2[2] * F1[5]) + 1.333333333 * (P3[2] *
      OM3 * TMP11)) + (TMP2 * - 1. * (-cI * (F2[5] * F1[2] + F2[3] * F1[4]) +
      cI * (F2[4] * F1[3] + F2[2] * F1[5]) + 1.333333333 * (P3[2] * OM3 *
      TMP11)) + (P3[2] * 0.666666667 * (TMP39 - TMP40) + TMP11 * (P2[2] -
      P1[2])))) + P3[2] * (TMP1 * (F2[4] * F1[2] + F2[5] * F1[3] - F2[2] *
      F1[4] - F2[3] * F1[5]) + (TMP2 * (F2[2] * F1[4] + F2[3] * F1[5] - F2[4] *
      F1[2] - F2[5] * F1[3]) + TMP11 * (P2[0] - P1[0])))) + (F1[2] * (F2[4] *
      (P2[2] - P1[2]) + F2[5] * (-cI * (P2[0]) + cI * (P1[0]))) + (F1[3] *
      (F2[4] * (-cI * (P1[0]) + cI * (P2[0])) + F2[5] * (P2[2] - P1[2])) +
      (F1[4] * (F2[2] * (P1[2] - P2[2]) + F2[3] * (-cI * (P2[0]) + cI *
      (P1[0]))) + F1[5] * (F2[2] * (-cI * (P1[0]) + cI * (P2[0])) + F2[3] *
      (P1[2] - P2[2]))))));
  T3[8] = denom * cI * (OM3 * (P3[1] * (TMP1 * (-cI * (F2[5] * F1[2] + F2[3] *
      F1[4]) + cI * (F2[4] * F1[3] + F2[2] * F1[5]) + 1.333333333 * (P3[2] *
      OM3 * TMP11)) + (TMP2 * - 1. * (-cI * (F2[5] * F1[2] + F2[3] * F1[4]) +
      cI * (F2[4] * F1[3] + F2[2] * F1[5]) + 1.333333333 * (P3[2] * OM3 *
      TMP11)) + (P3[2] * 0.666666667 * (TMP39 - TMP40) + TMP11 * (P2[2] -
      P1[2])))) + P3[2] * (TMP1 * - 1. * (F2[5] * F1[2] + F2[4] * F1[3] + F2[3]
      * F1[4] + F2[2] * F1[5]) + (TMP2 * (F2[5] * F1[2] + F2[4] * F1[3] + F2[3]
      * F1[4] + F2[2] * F1[5]) + TMP11 * (P2[1] - P1[1])))) + (F1[2] * F2[5] *
      (P1[2] - cI * (P2[1]) + cI * (P1[1]) - P2[2]) + (F1[3] * F2[4] * (P1[2] -
      cI * (P1[1]) + cI * (P2[1]) - P2[2]) + (F1[4] * F2[3] * (P1[2] - cI *
      (P2[1]) + cI * (P1[1]) - P2[2]) + F1[5] * F2[2] * (P1[2] - cI * (P1[1]) +
      cI * (P2[1]) - P2[2])))));
  T3[12] = denom * 2. * cI * (OM3 * (P3[2] * (TMP1 * (-cI * (F2[5] * F1[2] +
      F2[3] * F1[4]) + cI * (F2[4] * F1[3] + F2[2] * F1[5]) + 0.666666667 *
      (P3[2] * OM3 * TMP11)) + (TMP2 * - 1. * (-cI * (F2[5] * F1[2] + F2[3] *
      F1[4]) + cI * (F2[4] * F1[3] + F2[2] * F1[5]) + 0.666666667 * (P3[2] *
      OM3 * TMP11)) + (P3[2] * 0.333333333 * (TMP39 - TMP40) + TMP11 * (P2[2] -
      P1[2])))) + 0.333333333 * (TMP11 * (TMP2 - TMP1))) + (P1[2] * (-cI *
      (F2[4] * F1[3] + F2[2] * F1[5]) + cI * (F2[5] * F1[2] + F2[3] * F1[4])) +
      (P2[2] * (-cI * (F2[5] * F1[2] + F2[3] * F1[4]) + cI * (F2[4] * F1[3] +
      F2[2] * F1[5])) + (-0.333333333 * (TMP40) + 0.333333333 * (TMP39)))));
  T3[16] = denom * cI * (OM3 * (P3[2] * (TMP1 * (F2[5] * F1[3] + F2[3] * F1[5]
      + 1.333333333 * (P3[3] * OM3 * TMP11) - F2[4] * F1[2] - F2[2] * F1[4]) +
      (TMP2 * - 1. * (F2[5] * F1[3] + F2[3] * F1[5] + 1.333333333 * (P3[3] *
      OM3 * TMP11) - F2[4] * F1[2] - F2[2] * F1[4]) + (P3[3] * 0.666666667 *
      (TMP39 - TMP40) + TMP11 * (P2[3] - P1[3])))) + P3[3] * (TMP1 * (-cI *
      (F2[5] * F1[2] + F2[3] * F1[4]) + cI * (F2[4] * F1[3] + F2[2] * F1[5])) +
      (TMP2 * (-cI * (F2[4] * F1[3] + F2[2] * F1[5]) + cI * (F2[5] * F1[2] +
      F2[3] * F1[4])) + TMP11 * (P2[2] - P1[2])))) + (F1[2] * (F2[4] * (P1[2] -
      P2[2]) + F2[5] * (-cI * (P2[3]) + cI * (P1[3]))) + (F1[3] * (F2[4] * (-cI
      * (P1[3]) + cI * (P2[3])) + F2[5] * (P2[2] - P1[2])) + (F1[4] * (F2[2] *
      (P1[2] - P2[2]) + F2[3] * (-cI * (P2[3]) + cI * (P1[3]))) + F1[5] *
      (F2[2] * (-cI * (P1[3]) + cI * (P2[3])) + F2[3] * (P2[2] - P1[2]))))));
  T3[5] = denom * cI * (OM3 * (P3[0] * (TMP1 * (F2[5] * F1[3] + F2[3] * F1[5] +
      1.333333333 * (P3[3] * OM3 * TMP11) - F2[4] * F1[2] - F2[2] * F1[4]) +
      (TMP2 * - 1. * (F2[5] * F1[3] + F2[3] * F1[5] + 1.333333333 * (P3[3] *
      OM3 * TMP11) - F2[4] * F1[2] - F2[2] * F1[4]) + (P3[3] * 0.666666667 *
      (TMP39 - TMP40) + TMP11 * (P2[3] - P1[3])))) + P3[3] * (TMP1 * (F2[4] *
      F1[2] + F2[5] * F1[3] - F2[2] * F1[4] - F2[3] * F1[5]) + (TMP2 * (F2[2] *
      F1[4] + F2[3] * F1[5] - F2[4] * F1[2] - F2[5] * F1[3]) + TMP11 * (P2[0] -
      P1[0])))) + (F1[2] * F2[4] * (P1[0] + P2[3] - P1[3] - P2[0]) + (F1[3] *
      F2[5] * (P2[3] + P2[0] - P1[3] - P1[0]) + (F1[4] * F2[2] * (P1[3] + P1[0]
      - P2[3] - P2[0]) + F1[5] * F2[3] * (P1[3] + P2[0] - P1[0] - P2[3])))));
  T3[9] = denom * cI * (OM3 * (P3[1] * (TMP1 * (F2[5] * F1[3] + F2[3] * F1[5] +
      1.333333333 * (P3[3] * OM3 * TMP11) - F2[4] * F1[2] - F2[2] * F1[4]) +
      (TMP2 * - 1. * (F2[5] * F1[3] + F2[3] * F1[5] + 1.333333333 * (P3[3] *
      OM3 * TMP11) - F2[4] * F1[2] - F2[2] * F1[4]) + (P3[3] * 0.666666667 *
      (TMP39 - TMP40) + TMP11 * (P2[3] - P1[3])))) + P3[3] * (TMP1 * - 1. *
      (F2[5] * F1[2] + F2[4] * F1[3] + F2[3] * F1[4] + F2[2] * F1[5]) + (TMP2 *
      (F2[5] * F1[2] + F2[4] * F1[3] + F2[3] * F1[4] + F2[2] * F1[5]) + TMP11 *
      (P2[1] - P1[1])))) + (F1[2] * (F2[4] * (P1[1] - P2[1]) + F2[5] * (P1[3] -
      P2[3])) + (F1[3] * (F2[4] * (P1[3] - P2[3]) + F2[5] * (P2[1] - P1[1])) +
      (F1[4] * (F2[2] * (P1[1] - P2[1]) + F2[3] * (P1[3] - P2[3])) + F1[5] *
      (F2[2] * (P1[3] - P2[3]) + F2[3] * (P2[1] - P1[1]))))));
  T3[13] = denom * cI * (OM3 * (P3[2] * (TMP1 * (F2[5] * F1[3] + F2[3] * F1[5]
      + 1.333333333 * (P3[3] * OM3 * TMP11) - F2[4] * F1[2] - F2[2] * F1[4]) +
      (TMP2 * - 1. * (F2[5] * F1[3] + F2[3] * F1[5] + 1.333333333 * (P3[3] *
      OM3 * TMP11) - F2[4] * F1[2] - F2[2] * F1[4]) + (P3[3] * 0.666666667 *
      (TMP39 - TMP40) + TMP11 * (P2[3] - P1[3])))) + P3[3] * (TMP1 * (-cI *
      (F2[5] * F1[2] + F2[3] * F1[4]) + cI * (F2[4] * F1[3] + F2[2] * F1[5])) +
      (TMP2 * (-cI * (F2[4] * F1[3] + F2[2] * F1[5]) + cI * (F2[5] * F1[2] +
      F2[3] * F1[4])) + TMP11 * (P2[2] - P1[2])))) + (F1[2] * (F2[4] * (P1[2] -
      P2[2]) + F2[5] * (-cI * (P2[3]) + cI * (P1[3]))) + (F1[3] * (F2[4] * (-cI
      * (P1[3]) + cI * (P2[3])) + F2[5] * (P2[2] - P1[2])) + (F1[4] * (F2[2] *
      (P1[2] - P2[2]) + F2[3] * (-cI * (P2[3]) + cI * (P1[3]))) + F1[5] *
      (F2[2] * (-cI * (P1[3]) + cI * (P2[3])) + F2[3] * (P2[2] - P1[2]))))));
  T3[17] = denom * 2. * cI * (OM3 * (P3[3] * (TMP1 * (F2[5] * F1[3] + F2[3] *
      F1[5] + 0.666666667 * (P3[3] * OM3 * TMP11) - F2[4] * F1[2] - F2[2] *
      F1[4]) + (TMP2 * - 1. * (F2[5] * F1[3] + F2[3] * F1[5] + 0.666666667 *
      (P3[3] * OM3 * TMP11) - F2[4] * F1[2] - F2[2] * F1[4]) + (P3[3] *
      0.333333333 * (TMP39 - TMP40) + TMP11 * (P2[3] - P1[3])))) + 0.333333333
      * (TMP11 * (TMP2 - TMP1))) + (P1[3] * (F2[4] * F1[2] + F2[2] * F1[4] -
      F2[5] * F1[3] - F2[3] * F1[5]) + (P2[3] * (F2[5] * F1[3] + F2[3] * F1[5]
      - F2[4] * F1[2] - F2[2] * F1[4]) + (-0.333333333 * (TMP40) + 0.333333333
      * (TMP39)))));
}


void VVT10_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP1; 
  complex<double> TMP10; 
  double P3[4]; 
  complex<double> TMP47; 
  complex<double> TMP9; 
  complex<double> TMP2; 
  double P2[4]; 
  double OM3; 
  complex<double> TMP49; 
  double P1[4]; 
  complex<double> TMP23; 
  complex<double> TMP7; 
  complex<double> denom; 
  complex<double> TMP48; 
  complex<double> TMP50; 
  complex<double> TMP26; 
  complex<double> TMP4; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +V1[0] + V2[0]; 
  T3[1] = +V1[1] + V2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP50 = (P2[0] * V2[2] - P2[1] * V2[3] - P2[2] * V2[4] - P2[3] * V2[5]); 
  TMP26 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP47 = (P1[0] * P1[0] - P1[1] * P1[1] - P1[2] * P1[2] - P1[3] * P1[3]); 
  TMP23 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP48 = (P2[0] * P2[0] - P2[1] * P2[1] - P2[2] * P2[2] - P2[3] * P2[3]); 
  TMP49 = (P1[0] * V1[2] - P1[1] * V1[3] - P1[2] * V1[4] - P1[3] * V1[5]); 
  TMP9 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP4 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP7 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP1 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP10 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP2 = (P3[0] * P2[0] - P3[1] * P2[1] - P3[2] * P2[2] - P3[3] * P2[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * 0.333333333 * (OM3 * (P3[0] * (P3[0] * (OM3 * (TMP1 * (TMP2 *
      (TMP4 * - 2. * (-2. * cI * (TMP10) + cI * (TMP1 + TMP2)) + (TMP7 * 2. *
      (-2. * cI * (TMP9) + cI * (TMP26)) + 2. * cI * (TMP9 * TMP23))) + (TMP1 *
      (TMP9 * 2. * (+cI * (TMP7 + TMP23)) - 2. * cI * (TMP4 * TMP10)) - 2. * cI
      * (TMP10 * TMP23 * TMP26))) + TMP2 * (TMP2 * (TMP7 * 2. * (+cI * (TMP9 +
      TMP26)) - 2. * cI * (TMP4 * TMP10)) - 2. * cI * (TMP10 * TMP23 * TMP26)))
      + (TMP10 * (TMP4 * - 1. * (-2. * cI * (TMP10) + cI * (TMP47 + TMP48)) +
      (-2. * cI * (TMP7 * TMP9) - cI * (TMP23 * TMP49 + TMP26 * TMP50))) +
      (TMP7 * (TMP48 * (+cI * (TMP9 + TMP26)) + (+cI * (TMP9 * TMP47 + TMP2 *
      TMP49))) + (TMP1 * (-cI * (TMP4 * TMP48) + cI * (TMP9 * TMP50)) + TMP47 *
      (-cI * (TMP2 * TMP4) + cI * (TMP9 * TMP23)))))) + (TMP1 * (TMP9 * (P1[0]
      * - 6. * (+cI * (TMP7 + TMP23)) + (P2[0] * 3. * (-cI * (TMP23) + 2. * cI
      * (TMP7)) - 3. * cI * (TMP2 * V2[2]))) + (TMP4 * (P1[0] * 6. * (+cI *
      (TMP10 + TMP2)) + 6. * (P2[0] * (-cI * (TMP10) + cI * (TMP2)))) + 3. *
      (V1[2] * (-cI * (TMP2 * TMP7) + cI * (TMP10 * TMP23))))) + (TMP2 * (TMP7
      * (P1[0] * 3. * (-cI * (TMP26) + 2. * cI * (TMP9)) - 6. * (P2[0] * (+cI *
      (TMP9 + TMP26)))) + TMP10 * (TMP4 * 6. * (-cI * (P1[0]) + cI * (P2[0])) +
      3. * cI * (V2[2] * TMP26))) + 3. * (TMP10 * TMP23 * TMP26 * (+cI * (P1[0]
      + P2[0])))))) + (TMP1 * (TMP2 * (TMP4 * - 1. * (-2. * cI * (TMP10) + cI *
      (TMP1 + TMP2)) + (TMP7 * (-2. * cI * (TMP9) + cI * (TMP26)) + cI * (TMP9
      * TMP23))) + (TMP1 * (TMP9 * (+cI * (TMP7 + TMP23)) - cI * (TMP4 *
      TMP10)) - cI * (TMP10 * TMP23 * TMP26))) + TMP2 * (TMP2 * (TMP7 * (+cI *
      (TMP9 + TMP26)) - cI * (TMP4 * TMP10)) - cI * (TMP10 * TMP23 * TMP26))))
      + (TMP10 * (TMP4 * (P1[0] * 3. * (-cI * (P1[0]) + 2. * cI * (P2[0])) +
      (+1.000000000 * cI * (TMP47 + TMP48) - 2.000000000 * cI * (TMP10) - 3. *
      cI * (P2[0] * P2[0]))) + (TMP23 * 1.000000000 * (+cI * (TMP49) -
      3.000000000 * cI * (P1[0] * V1[2])) + (+2.000000000 * cI * (TMP7 * TMP9)
      + TMP26 * 1.000000000 * (+cI * (TMP50) - 3.000000000 * cI * (P2[0] *
      V2[2]))))) + (TMP7 * (TMP9 * (P1[0] * 3. * (-2. * cI * (P2[0]) + cI *
      (P1[0])) + (-1.000000000 * cI * (TMP47 + TMP48) + 3. * cI * (P2[0] *
      P2[0]))) + (+1.000000000 * (TMP2 * (-cI * (TMP49) + 3.000000000 * cI *
      (P1[0] * V1[2]))) + TMP26 * 1.000000000 * (-cI * (TMP48) + 3.000000000 *
      cI * (P2[0] * P2[0])))) + (TMP1 * 1.000000000 * (TMP9 * (-cI * (TMP50) +
      3.000000000 * cI * (P2[0] * V2[2])) + 1.000000000 * (TMP4 * 1.000000000 *
      (+cI * (TMP48) - 3.000000000 * cI * (P2[0] * P2[0])))) + (+1.000000000 *
      (TMP23 * TMP9 * (-cI * (TMP47) + 3.000000000 * cI * (P1[0] * P1[0]))) +
      TMP2 * 1.000000000 * TMP4 * (+cI * (TMP47) - 3.000000000 * cI * (P1[0] *
      P1[0])))))));
  T3[6] = denom * (OM3 * (P3[0] * (P3[1] * (OM3 * (TMP1 * (TMP2 * (TMP4 * -
      0.666666667 * (-2. * cI * (TMP10) + cI * (TMP1 + TMP2)) + (TMP7 *
      1.333333333 * (+0.500000000 * cI * (TMP26) - cI * (TMP9)) + 0.666666667 *
      cI * (TMP9 * TMP23))) + (TMP1 * (TMP9 * 0.666666667 * (+cI * (TMP7 +
      TMP23)) - 0.666666667 * cI * (TMP4 * TMP10)) - 0.666666667 * cI * (TMP10
      * TMP23 * TMP26))) + TMP2 * (TMP2 * (TMP7 * 0.666666667 * (+cI * (TMP9 +
      TMP26)) - 0.666666667 * cI * (TMP4 * TMP10)) - 0.666666667 * cI * (TMP10
      * TMP23 * TMP26))) + (TMP10 * (TMP4 * - 0.333333333 * (-2. * cI * (TMP10)
      + cI * (TMP47 + TMP48)) + (-0.666666667 * cI * (TMP7 * TMP9) -
      0.333333333 * cI * (TMP23 * TMP49 + TMP26 * TMP50))) + (TMP7 * (TMP48 *
      0.333333333 * (+cI * (TMP9 + TMP26)) + (+0.333333333 * cI * (TMP9 * TMP47
      + TMP2 * TMP49))) + (TMP1 * 0.333333333 * (-cI * (TMP4 * TMP48) + cI *
      (TMP9 * TMP50)) + 0.333333333 * (TMP47 * (-cI * (TMP2 * TMP4) + cI *
      (TMP9 * TMP23))))))) + (TMP1 * (TMP9 * (P1[1] * - 1. * (+cI * (TMP7 +
      TMP23)) + (-0.500000000 * cI * (TMP2 * V2[3]) + P2[1] * 0.500000000 *
      (-cI * (TMP23) + 2. * cI * (TMP7)))) + (+0.500000000 * (V1[3] * (-cI *
      (TMP2 * TMP7) + cI * (TMP10 * TMP23))) + TMP4 * (P1[1] * (+cI * (TMP10 +
      TMP2)) + P2[1] * (-cI * (TMP10) + cI * (TMP2))))) + (+0.500000000 *
      (TMP10 * TMP23 * TMP26 * (+cI * (P1[1] + P2[1]))) + TMP2 * 0.500000000 *
      (TMP10 * (+2. * (TMP4 * (-cI * (P1[1]) + cI * (P2[1]))) + cI * (V2[3] *
      TMP26)) + 2. * (TMP7 * (P1[1] * 0.500000000 * (-cI * (TMP26) + 2. * cI *
      (TMP9)) - P2[1] * (+cI * (TMP9 + TMP26)))))))) + P3[1] * (TMP1 * (TMP9 *
      (P1[0] * - 1. * (+cI * (TMP7 + TMP23)) + (-0.500000000 * cI * (TMP2 *
      V2[2]) + P2[0] * 0.500000000 * (-cI * (TMP23) + 2. * cI * (TMP7)))) +
      (+0.500000000 * (V1[2] * (-cI * (TMP2 * TMP7) + cI * (TMP10 * TMP23))) +
      TMP4 * (P1[0] * (+cI * (TMP10 + TMP2)) + P2[0] * (-cI * (TMP10) + cI *
      (TMP2))))) + (+0.500000000 * (TMP10 * TMP23 * TMP26 * (+cI * (P1[0] +
      P2[0]))) + TMP2 * 0.500000000 * (TMP10 * (+2. * (TMP4 * (-cI * (P1[0]) +
      cI * (P2[0]))) + cI * (V2[2] * TMP26)) + 2. * (TMP7 * (P1[0] *
      0.500000000 * (-cI * (TMP26) + 2. * cI * (TMP9)) - P2[0] * (+cI * (TMP9 +
      TMP26)))))))) + (P1[0] * (P1[1] * (TMP4 * - 1. * (+cI * (TMP10 + TMP2)) +
      TMP9 * (+cI * (TMP7 + TMP23))) + (+0.500000000 * (V1[3] * (-cI * (TMP10 *
      TMP23) + cI * (TMP2 * TMP7))) + P2[1] * (-cI * (TMP7 * TMP9) + cI * (TMP4
      * TMP10)))) + (P2[0] * (P2[1] * (TMP4 * - 1. * (+cI * (TMP10 + TMP1)) +
      TMP7 * (+cI * (TMP9 + TMP26))) + (+0.500000000 * (V2[3] * (-cI * (TMP10 *
      TMP26) + cI * (TMP1 * TMP9))) + P1[1] * (-cI * (TMP7 * TMP9) + cI * (TMP4
      * TMP10)))) + (+0.500000000 * (P2[1] * V2[2] * (-cI * (TMP10 * TMP26) +
      cI * (TMP1 * TMP9))) + P1[1] * 0.500000000 * V1[2] * (-cI * (TMP10 *
      TMP23) + cI * (TMP2 * TMP7))))));
  T3[10] = denom * (OM3 * (P3[0] * (P3[2] * (OM3 * (TMP1 * (TMP2 * (TMP4 * -
      0.666666667 * (-2. * cI * (TMP10) + cI * (TMP1 + TMP2)) + (TMP7 *
      1.333333333 * (+0.500000000 * cI * (TMP26) - cI * (TMP9)) + 0.666666667 *
      cI * (TMP9 * TMP23))) + (TMP1 * (TMP9 * 0.666666667 * (+cI * (TMP7 +
      TMP23)) - 0.666666667 * cI * (TMP4 * TMP10)) - 0.666666667 * cI * (TMP10
      * TMP23 * TMP26))) + TMP2 * (TMP2 * (TMP7 * 0.666666667 * (+cI * (TMP9 +
      TMP26)) - 0.666666667 * cI * (TMP4 * TMP10)) - 0.666666667 * cI * (TMP10
      * TMP23 * TMP26))) + (TMP10 * (TMP4 * - 0.333333333 * (-2. * cI * (TMP10)
      + cI * (TMP47 + TMP48)) + (-0.666666667 * cI * (TMP7 * TMP9) -
      0.333333333 * cI * (TMP23 * TMP49 + TMP26 * TMP50))) + (TMP7 * (TMP48 *
      0.333333333 * (+cI * (TMP9 + TMP26)) + (+0.333333333 * cI * (TMP9 * TMP47
      + TMP2 * TMP49))) + (TMP1 * 0.333333333 * (-cI * (TMP4 * TMP48) + cI *
      (TMP9 * TMP50)) + 0.333333333 * (TMP47 * (-cI * (TMP2 * TMP4) + cI *
      (TMP9 * TMP23))))))) + (TMP1 * (TMP9 * (P1[2] * - 1. * (+cI * (TMP7 +
      TMP23)) + (-0.500000000 * cI * (TMP2 * V2[4]) + P2[2] * 0.500000000 *
      (-cI * (TMP23) + 2. * cI * (TMP7)))) + (+0.500000000 * (V1[4] * (-cI *
      (TMP2 * TMP7) + cI * (TMP10 * TMP23))) + TMP4 * (P1[2] * (+cI * (TMP10 +
      TMP2)) + P2[2] * (-cI * (TMP10) + cI * (TMP2))))) + (+0.500000000 *
      (TMP10 * TMP23 * TMP26 * (+cI * (P1[2] + P2[2]))) + TMP2 * 0.500000000 *
      (TMP10 * (+2. * (TMP4 * (-cI * (P1[2]) + cI * (P2[2]))) + cI * (V2[4] *
      TMP26)) + 2. * (TMP7 * (P1[2] * 0.500000000 * (-cI * (TMP26) + 2. * cI *
      (TMP9)) - P2[2] * (+cI * (TMP9 + TMP26)))))))) + P3[2] * (TMP1 * (TMP9 *
      (P1[0] * - 1. * (+cI * (TMP7 + TMP23)) + (-0.500000000 * cI * (TMP2 *
      V2[2]) + P2[0] * 0.500000000 * (-cI * (TMP23) + 2. * cI * (TMP7)))) +
      (+0.500000000 * (V1[2] * (-cI * (TMP2 * TMP7) + cI * (TMP10 * TMP23))) +
      TMP4 * (P1[0] * (+cI * (TMP10 + TMP2)) + P2[0] * (-cI * (TMP10) + cI *
      (TMP2))))) + (+0.500000000 * (TMP10 * TMP23 * TMP26 * (+cI * (P1[0] +
      P2[0]))) + TMP2 * 0.500000000 * (TMP10 * (+2. * (TMP4 * (-cI * (P1[0]) +
      cI * (P2[0]))) + cI * (V2[2] * TMP26)) + 2. * (TMP7 * (P1[0] *
      0.500000000 * (-cI * (TMP26) + 2. * cI * (TMP9)) - P2[0] * (+cI * (TMP9 +
      TMP26)))))))) + (P1[0] * (P1[2] * (TMP4 * - 1. * (+cI * (TMP10 + TMP2)) +
      TMP9 * (+cI * (TMP7 + TMP23))) + (+0.500000000 * (V1[4] * (-cI * (TMP10 *
      TMP23) + cI * (TMP2 * TMP7))) + P2[2] * (-cI * (TMP7 * TMP9) + cI * (TMP4
      * TMP10)))) + (P2[0] * (P2[2] * (TMP4 * - 1. * (+cI * (TMP10 + TMP1)) +
      TMP7 * (+cI * (TMP9 + TMP26))) + (+0.500000000 * (V2[4] * (-cI * (TMP10 *
      TMP26) + cI * (TMP1 * TMP9))) + P1[2] * (-cI * (TMP7 * TMP9) + cI * (TMP4
      * TMP10)))) + (+0.500000000 * (P2[2] * V2[2] * (-cI * (TMP10 * TMP26) +
      cI * (TMP1 * TMP9))) + P1[2] * 0.500000000 * V1[2] * (-cI * (TMP10 *
      TMP23) + cI * (TMP2 * TMP7))))));
  T3[14] = denom * (OM3 * (P3[0] * (P3[3] * (OM3 * (TMP1 * (TMP2 * (TMP4 * -
      0.666666667 * (-2. * cI * (TMP10) + cI * (TMP1 + TMP2)) + (TMP7 *
      1.333333333 * (+0.500000000 * cI * (TMP26) - cI * (TMP9)) + 0.666666667 *
      cI * (TMP9 * TMP23))) + (TMP1 * (TMP9 * 0.666666667 * (+cI * (TMP7 +
      TMP23)) - 0.666666667 * cI * (TMP4 * TMP10)) - 0.666666667 * cI * (TMP10
      * TMP23 * TMP26))) + TMP2 * (TMP2 * (TMP7 * 0.666666667 * (+cI * (TMP9 +
      TMP26)) - 0.666666667 * cI * (TMP4 * TMP10)) - 0.666666667 * cI * (TMP10
      * TMP23 * TMP26))) + (TMP10 * (TMP4 * - 0.333333333 * (-2. * cI * (TMP10)
      + cI * (TMP47 + TMP48)) + (-0.666666667 * cI * (TMP7 * TMP9) -
      0.333333333 * cI * (TMP23 * TMP49 + TMP26 * TMP50))) + (TMP7 * (TMP48 *
      0.333333333 * (+cI * (TMP9 + TMP26)) + (+0.333333333 * cI * (TMP9 * TMP47
      + TMP2 * TMP49))) + (TMP1 * 0.333333333 * (-cI * (TMP4 * TMP48) + cI *
      (TMP9 * TMP50)) + 0.333333333 * (TMP47 * (-cI * (TMP2 * TMP4) + cI *
      (TMP9 * TMP23))))))) + (TMP1 * (TMP9 * (P1[3] * - 1. * (+cI * (TMP7 +
      TMP23)) + (-0.500000000 * cI * (TMP2 * V2[5]) + P2[3] * 0.500000000 *
      (-cI * (TMP23) + 2. * cI * (TMP7)))) + (+0.500000000 * (V1[5] * (-cI *
      (TMP2 * TMP7) + cI * (TMP10 * TMP23))) + TMP4 * (P1[3] * (+cI * (TMP10 +
      TMP2)) + P2[3] * (-cI * (TMP10) + cI * (TMP2))))) + (+0.500000000 *
      (TMP10 * TMP23 * TMP26 * (+cI * (P1[3] + P2[3]))) + TMP2 * 0.500000000 *
      (TMP10 * (+2. * (TMP4 * (-cI * (P1[3]) + cI * (P2[3]))) + cI * (V2[5] *
      TMP26)) + 2. * (TMP7 * (P1[3] * 0.500000000 * (-cI * (TMP26) + 2. * cI *
      (TMP9)) - P2[3] * (+cI * (TMP9 + TMP26)))))))) + P3[3] * (TMP1 * (TMP9 *
      (P1[0] * - 1. * (+cI * (TMP7 + TMP23)) + (-0.500000000 * cI * (TMP2 *
      V2[2]) + P2[0] * 0.500000000 * (-cI * (TMP23) + 2. * cI * (TMP7)))) +
      (+0.500000000 * (V1[2] * (-cI * (TMP2 * TMP7) + cI * (TMP10 * TMP23))) +
      TMP4 * (P1[0] * (+cI * (TMP10 + TMP2)) + P2[0] * (-cI * (TMP10) + cI *
      (TMP2))))) + (+0.500000000 * (TMP10 * TMP23 * TMP26 * (+cI * (P1[0] +
      P2[0]))) + TMP2 * 0.500000000 * (TMP10 * (+2. * (TMP4 * (-cI * (P1[0]) +
      cI * (P2[0]))) + cI * (V2[2] * TMP26)) + 2. * (TMP7 * (P1[0] *
      0.500000000 * (-cI * (TMP26) + 2. * cI * (TMP9)) - P2[0] * (+cI * (TMP9 +
      TMP26)))))))) + (P1[0] * (P1[3] * (TMP4 * - 1. * (+cI * (TMP10 + TMP2)) +
      TMP9 * (+cI * (TMP7 + TMP23))) + (+0.500000000 * (V1[5] * (-cI * (TMP10 *
      TMP23) + cI * (TMP2 * TMP7))) + P2[3] * (-cI * (TMP7 * TMP9) + cI * (TMP4
      * TMP10)))) + (P2[0] * (P2[3] * (TMP4 * - 1. * (+cI * (TMP10 + TMP1)) +
      TMP7 * (+cI * (TMP9 + TMP26))) + (+0.500000000 * (V2[5] * (-cI * (TMP10 *
      TMP26) + cI * (TMP1 * TMP9))) + P1[3] * (-cI * (TMP7 * TMP9) + cI * (TMP4
      * TMP10)))) + (+0.500000000 * (P2[3] * V2[2] * (-cI * (TMP10 * TMP26) +
      cI * (TMP1 * TMP9))) + P1[3] * 0.500000000 * V1[2] * (-cI * (TMP10 *
      TMP23) + cI * (TMP2 * TMP7))))));
  T3[3] = denom * (OM3 * (P3[0] * (P3[1] * (OM3 * (TMP1 * (TMP2 * (TMP4 * -
      0.666666667 * (-2. * cI * (TMP10) + cI * (TMP1 + TMP2)) + (TMP7 *
      1.333333333 * (+0.500000000 * cI * (TMP26) - cI * (TMP9)) + 0.666666667 *
      cI * (TMP9 * TMP23))) + (TMP1 * (TMP9 * 0.666666667 * (+cI * (TMP7 +
      TMP23)) - 0.666666667 * cI * (TMP4 * TMP10)) - 0.666666667 * cI * (TMP10
      * TMP23 * TMP26))) + TMP2 * (TMP2 * (TMP7 * 0.666666667 * (+cI * (TMP9 +
      TMP26)) - 0.666666667 * cI * (TMP4 * TMP10)) - 0.666666667 * cI * (TMP10
      * TMP23 * TMP26))) + (TMP10 * (TMP4 * - 0.333333333 * (-2. * cI * (TMP10)
      + cI * (TMP47 + TMP48)) + (-0.666666667 * cI * (TMP7 * TMP9) -
      0.333333333 * cI * (TMP23 * TMP49 + TMP26 * TMP50))) + (TMP7 * (TMP48 *
      0.333333333 * (+cI * (TMP9 + TMP26)) + (+0.333333333 * cI * (TMP9 * TMP47
      + TMP2 * TMP49))) + (TMP1 * 0.333333333 * (-cI * (TMP4 * TMP48) + cI *
      (TMP9 * TMP50)) + 0.333333333 * (TMP47 * (-cI * (TMP2 * TMP4) + cI *
      (TMP9 * TMP23))))))) + (TMP1 * (TMP9 * (P1[1] * - 1. * (+cI * (TMP7 +
      TMP23)) + (-0.500000000 * cI * (TMP2 * V2[3]) + P2[1] * 0.500000000 *
      (-cI * (TMP23) + 2. * cI * (TMP7)))) + (+0.500000000 * (V1[3] * (-cI *
      (TMP2 * TMP7) + cI * (TMP10 * TMP23))) + TMP4 * (P1[1] * (+cI * (TMP10 +
      TMP2)) + P2[1] * (-cI * (TMP10) + cI * (TMP2))))) + (+0.500000000 *
      (TMP10 * TMP23 * TMP26 * (+cI * (P1[1] + P2[1]))) + TMP2 * 0.500000000 *
      (TMP10 * (+2. * (TMP4 * (-cI * (P1[1]) + cI * (P2[1]))) + cI * (V2[3] *
      TMP26)) + 2. * (TMP7 * (P1[1] * 0.500000000 * (-cI * (TMP26) + 2. * cI *
      (TMP9)) - P2[1] * (+cI * (TMP9 + TMP26)))))))) + P3[1] * (TMP1 * (TMP9 *
      (P1[0] * - 1. * (+cI * (TMP7 + TMP23)) + (-0.500000000 * cI * (TMP2 *
      V2[2]) + P2[0] * 0.500000000 * (-cI * (TMP23) + 2. * cI * (TMP7)))) +
      (+0.500000000 * (V1[2] * (-cI * (TMP2 * TMP7) + cI * (TMP10 * TMP23))) +
      TMP4 * (P1[0] * (+cI * (TMP10 + TMP2)) + P2[0] * (-cI * (TMP10) + cI *
      (TMP2))))) + (+0.500000000 * (TMP10 * TMP23 * TMP26 * (+cI * (P1[0] +
      P2[0]))) + TMP2 * 0.500000000 * (TMP10 * (+2. * (TMP4 * (-cI * (P1[0]) +
      cI * (P2[0]))) + cI * (V2[2] * TMP26)) + 2. * (TMP7 * (P1[0] *
      0.500000000 * (-cI * (TMP26) + 2. * cI * (TMP9)) - P2[0] * (+cI * (TMP9 +
      TMP26)))))))) + (P1[0] * (P1[1] * (TMP4 * - 1. * (+cI * (TMP10 + TMP2)) +
      TMP9 * (+cI * (TMP7 + TMP23))) + (+0.500000000 * (V1[3] * (-cI * (TMP10 *
      TMP23) + cI * (TMP2 * TMP7))) + P2[1] * (-cI * (TMP7 * TMP9) + cI * (TMP4
      * TMP10)))) + (P2[0] * (P2[1] * (TMP4 * - 1. * (+cI * (TMP10 + TMP1)) +
      TMP7 * (+cI * (TMP9 + TMP26))) + (+0.500000000 * (V2[3] * (-cI * (TMP10 *
      TMP26) + cI * (TMP1 * TMP9))) + P1[1] * (-cI * (TMP7 * TMP9) + cI * (TMP4
      * TMP10)))) + (+0.500000000 * (P2[1] * V2[2] * (-cI * (TMP10 * TMP26) +
      cI * (TMP1 * TMP9))) + P1[1] * 0.500000000 * V1[2] * (-cI * (TMP10 *
      TMP23) + cI * (TMP2 * TMP7))))));
  T3[7] = denom * 0.333333333 * (OM3 * (P3[1] * (P3[1] * (OM3 * (TMP1 * (TMP2 *
      (TMP4 * - 2. * (-2. * cI * (TMP10) + cI * (TMP1 + TMP2)) + (TMP7 * 2. *
      (-2. * cI * (TMP9) + cI * (TMP26)) + 2. * cI * (TMP9 * TMP23))) + (TMP1 *
      (TMP9 * 2. * (+cI * (TMP7 + TMP23)) - 2. * cI * (TMP4 * TMP10)) - 2. * cI
      * (TMP10 * TMP23 * TMP26))) + TMP2 * (TMP2 * (TMP7 * 2. * (+cI * (TMP9 +
      TMP26)) - 2. * cI * (TMP4 * TMP10)) - 2. * cI * (TMP10 * TMP23 * TMP26)))
      + (TMP10 * (TMP4 * - 1. * (-2. * cI * (TMP10) + cI * (TMP47 + TMP48)) +
      (-2. * cI * (TMP7 * TMP9) - cI * (TMP23 * TMP49 + TMP26 * TMP50))) +
      (TMP7 * (TMP48 * (+cI * (TMP9 + TMP26)) + (+cI * (TMP9 * TMP47 + TMP2 *
      TMP49))) + (TMP1 * (-cI * (TMP4 * TMP48) + cI * (TMP9 * TMP50)) + TMP47 *
      (-cI * (TMP2 * TMP4) + cI * (TMP9 * TMP23)))))) + (TMP1 * (TMP9 * (P1[1]
      * - 6. * (+cI * (TMP7 + TMP23)) + (P2[1] * 3. * (-cI * (TMP23) + 2. * cI
      * (TMP7)) - 3. * cI * (TMP2 * V2[3]))) + (TMP4 * (P1[1] * 6. * (+cI *
      (TMP10 + TMP2)) + 6. * (P2[1] * (-cI * (TMP10) + cI * (TMP2)))) + 3. *
      (V1[3] * (-cI * (TMP2 * TMP7) + cI * (TMP10 * TMP23))))) + (TMP2 * (TMP7
      * (P1[1] * 3. * (-cI * (TMP26) + 2. * cI * (TMP9)) - 6. * (P2[1] * (+cI *
      (TMP9 + TMP26)))) + TMP10 * (TMP4 * 6. * (-cI * (P1[1]) + cI * (P2[1])) +
      3. * cI * (V2[3] * TMP26))) + 3. * (TMP10 * TMP23 * TMP26 * (+cI * (P1[1]
      + P2[1])))))) + (TMP1 * (TMP2 * (TMP4 * (-2. * cI * (TMP10) + cI * (TMP1
      + TMP2)) + (TMP7 * (-cI * (TMP26) + 2. * cI * (TMP9)) - cI * (TMP9 *
      TMP23))) + (TMP1 * (TMP9 * - 1. * (+cI * (TMP7 + TMP23)) + cI * (TMP4 *
      TMP10)) + cI * (TMP10 * TMP23 * TMP26))) + TMP2 * (TMP2 * (TMP7 * - 1. *
      (+cI * (TMP9 + TMP26)) + cI * (TMP4 * TMP10)) + cI * (TMP10 * TMP23 *
      TMP26)))) + (TMP10 * (TMP4 * (P1[1] * 3. * (-cI * (P1[1]) + 2. * cI *
      (P2[1])) + (-1.000000000 * cI * (TMP47 + TMP48) + 2.000000000 * cI *
      (TMP10) - 3. * cI * (P2[1] * P2[1]))) + (TMP23 * - 1.000000000 * (+cI *
      (TMP49) + 3.000000000 * cI * (P1[1] * V1[3])) + (-2.000000000 * cI *
      (TMP7 * TMP9) + TMP26 * - 1.000000000 * (+cI * (TMP50) + 3.000000000 * cI
      * (P2[1] * V2[3]))))) + (TMP7 * (TMP9 * (P1[1] * 3. * (-2. * cI * (P2[1])
      + cI * (P1[1])) + (+1.000000000 * cI * (TMP47 + TMP48) + 3. * cI * (P2[1]
      * P2[1]))) + (+1.000000000 * (TMP2 * (+cI * (TMP49) + 3.000000000 * cI *
      (P1[1] * V1[3]))) + TMP26 * 1.000000000 * (+cI * (TMP48) + 3.000000000 *
      cI * (P2[1] * P2[1])))) + (TMP1 * 1.000000000 * (TMP9 * (+cI * (TMP50) +
      3.000000000 * cI * (P2[1] * V2[3])) + 1.000000000 * (TMP4 * - 1.000000000
      * (+cI * (TMP48) + 3.000000000 * cI * (P2[1] * P2[1])))) + (+1.000000000
      * (TMP23 * TMP9 * (+cI * (TMP47) + 3.000000000 * cI * (P1[1] * P1[1]))) +
      TMP2 * - 1.000000000 * TMP4 * (+cI * (TMP47) + 3.000000000 * cI * (P1[1]
      * P1[1])))))));
  T3[11] = denom * (OM3 * (P3[1] * (P3[2] * (OM3 * (TMP1 * (TMP2 * (TMP4 * -
      0.666666667 * (-2. * cI * (TMP10) + cI * (TMP1 + TMP2)) + (TMP7 *
      1.333333333 * (+0.500000000 * cI * (TMP26) - cI * (TMP9)) + 0.666666667 *
      cI * (TMP9 * TMP23))) + (TMP1 * (TMP9 * 0.666666667 * (+cI * (TMP7 +
      TMP23)) - 0.666666667 * cI * (TMP4 * TMP10)) - 0.666666667 * cI * (TMP10
      * TMP23 * TMP26))) + TMP2 * (TMP2 * (TMP7 * 0.666666667 * (+cI * (TMP9 +
      TMP26)) - 0.666666667 * cI * (TMP4 * TMP10)) - 0.666666667 * cI * (TMP10
      * TMP23 * TMP26))) + (TMP10 * (TMP4 * - 0.333333333 * (-2. * cI * (TMP10)
      + cI * (TMP47 + TMP48)) + (-0.666666667 * cI * (TMP7 * TMP9) -
      0.333333333 * cI * (TMP23 * TMP49 + TMP26 * TMP50))) + (TMP7 * (TMP48 *
      0.333333333 * (+cI * (TMP9 + TMP26)) + (+0.333333333 * cI * (TMP9 * TMP47
      + TMP2 * TMP49))) + (TMP1 * 0.333333333 * (-cI * (TMP4 * TMP48) + cI *
      (TMP9 * TMP50)) + 0.333333333 * (TMP47 * (-cI * (TMP2 * TMP4) + cI *
      (TMP9 * TMP23))))))) + (TMP1 * (TMP9 * (P1[2] * - 1. * (+cI * (TMP7 +
      TMP23)) + (-0.500000000 * cI * (TMP2 * V2[4]) + P2[2] * 0.500000000 *
      (-cI * (TMP23) + 2. * cI * (TMP7)))) + (+0.500000000 * (V1[4] * (-cI *
      (TMP2 * TMP7) + cI * (TMP10 * TMP23))) + TMP4 * (P1[2] * (+cI * (TMP10 +
      TMP2)) + P2[2] * (-cI * (TMP10) + cI * (TMP2))))) + (+0.500000000 *
      (TMP10 * TMP23 * TMP26 * (+cI * (P1[2] + P2[2]))) + TMP2 * 0.500000000 *
      (TMP10 * (+2. * (TMP4 * (-cI * (P1[2]) + cI * (P2[2]))) + cI * (V2[4] *
      TMP26)) + 2. * (TMP7 * (P1[2] * 0.500000000 * (-cI * (TMP26) + 2. * cI *
      (TMP9)) - P2[2] * (+cI * (TMP9 + TMP26)))))))) + P3[2] * (TMP1 * (TMP9 *
      (P1[1] * - 1. * (+cI * (TMP7 + TMP23)) + (-0.500000000 * cI * (TMP2 *
      V2[3]) + P2[1] * 0.500000000 * (-cI * (TMP23) + 2. * cI * (TMP7)))) +
      (+0.500000000 * (V1[3] * (-cI * (TMP2 * TMP7) + cI * (TMP10 * TMP23))) +
      TMP4 * (P1[1] * (+cI * (TMP10 + TMP2)) + P2[1] * (-cI * (TMP10) + cI *
      (TMP2))))) + (+0.500000000 * (TMP10 * TMP23 * TMP26 * (+cI * (P1[1] +
      P2[1]))) + TMP2 * 0.500000000 * (TMP10 * (+2. * (TMP4 * (-cI * (P1[1]) +
      cI * (P2[1]))) + cI * (V2[3] * TMP26)) + 2. * (TMP7 * (P1[1] *
      0.500000000 * (-cI * (TMP26) + 2. * cI * (TMP9)) - P2[1] * (+cI * (TMP9 +
      TMP26)))))))) + (P1[1] * (P1[2] * (TMP4 * - 1. * (+cI * (TMP10 + TMP2)) +
      TMP9 * (+cI * (TMP7 + TMP23))) + (+0.500000000 * (V1[4] * (-cI * (TMP10 *
      TMP23) + cI * (TMP2 * TMP7))) + P2[2] * (-cI * (TMP7 * TMP9) + cI * (TMP4
      * TMP10)))) + (P2[1] * (P2[2] * (TMP4 * - 1. * (+cI * (TMP10 + TMP1)) +
      TMP7 * (+cI * (TMP9 + TMP26))) + (+0.500000000 * (V2[4] * (-cI * (TMP10 *
      TMP26) + cI * (TMP1 * TMP9))) + P1[2] * (-cI * (TMP7 * TMP9) + cI * (TMP4
      * TMP10)))) + (+0.500000000 * (P2[2] * V2[3] * (-cI * (TMP10 * TMP26) +
      cI * (TMP1 * TMP9))) + P1[2] * 0.500000000 * V1[3] * (-cI * (TMP10 *
      TMP23) + cI * (TMP2 * TMP7))))));
  T3[15] = denom * (OM3 * (P3[1] * (P3[3] * (OM3 * (TMP1 * (TMP2 * (TMP4 * -
      0.666666667 * (-2. * cI * (TMP10) + cI * (TMP1 + TMP2)) + (TMP7 *
      1.333333333 * (+0.500000000 * cI * (TMP26) - cI * (TMP9)) + 0.666666667 *
      cI * (TMP9 * TMP23))) + (TMP1 * (TMP9 * 0.666666667 * (+cI * (TMP7 +
      TMP23)) - 0.666666667 * cI * (TMP4 * TMP10)) - 0.666666667 * cI * (TMP10
      * TMP23 * TMP26))) + TMP2 * (TMP2 * (TMP7 * 0.666666667 * (+cI * (TMP9 +
      TMP26)) - 0.666666667 * cI * (TMP4 * TMP10)) - 0.666666667 * cI * (TMP10
      * TMP23 * TMP26))) + (TMP10 * (TMP4 * - 0.333333333 * (-2. * cI * (TMP10)
      + cI * (TMP47 + TMP48)) + (-0.666666667 * cI * (TMP7 * TMP9) -
      0.333333333 * cI * (TMP23 * TMP49 + TMP26 * TMP50))) + (TMP7 * (TMP48 *
      0.333333333 * (+cI * (TMP9 + TMP26)) + (+0.333333333 * cI * (TMP9 * TMP47
      + TMP2 * TMP49))) + (TMP1 * 0.333333333 * (-cI * (TMP4 * TMP48) + cI *
      (TMP9 * TMP50)) + 0.333333333 * (TMP47 * (-cI * (TMP2 * TMP4) + cI *
      (TMP9 * TMP23))))))) + (TMP1 * (TMP9 * (P1[3] * - 1. * (+cI * (TMP7 +
      TMP23)) + (-0.500000000 * cI * (TMP2 * V2[5]) + P2[3] * 0.500000000 *
      (-cI * (TMP23) + 2. * cI * (TMP7)))) + (+0.500000000 * (V1[5] * (-cI *
      (TMP2 * TMP7) + cI * (TMP10 * TMP23))) + TMP4 * (P1[3] * (+cI * (TMP10 +
      TMP2)) + P2[3] * (-cI * (TMP10) + cI * (TMP2))))) + (+0.500000000 *
      (TMP10 * TMP23 * TMP26 * (+cI * (P1[3] + P2[3]))) + TMP2 * 0.500000000 *
      (TMP10 * (+2. * (TMP4 * (-cI * (P1[3]) + cI * (P2[3]))) + cI * (V2[5] *
      TMP26)) + 2. * (TMP7 * (P1[3] * 0.500000000 * (-cI * (TMP26) + 2. * cI *
      (TMP9)) - P2[3] * (+cI * (TMP9 + TMP26)))))))) + P3[3] * (TMP1 * (TMP9 *
      (P1[1] * - 1. * (+cI * (TMP7 + TMP23)) + (-0.500000000 * cI * (TMP2 *
      V2[3]) + P2[1] * 0.500000000 * (-cI * (TMP23) + 2. * cI * (TMP7)))) +
      (+0.500000000 * (V1[3] * (-cI * (TMP2 * TMP7) + cI * (TMP10 * TMP23))) +
      TMP4 * (P1[1] * (+cI * (TMP10 + TMP2)) + P2[1] * (-cI * (TMP10) + cI *
      (TMP2))))) + (+0.500000000 * (TMP10 * TMP23 * TMP26 * (+cI * (P1[1] +
      P2[1]))) + TMP2 * 0.500000000 * (TMP10 * (+2. * (TMP4 * (-cI * (P1[1]) +
      cI * (P2[1]))) + cI * (V2[3] * TMP26)) + 2. * (TMP7 * (P1[1] *
      0.500000000 * (-cI * (TMP26) + 2. * cI * (TMP9)) - P2[1] * (+cI * (TMP9 +
      TMP26)))))))) + (P1[1] * (P1[3] * (TMP4 * - 1. * (+cI * (TMP10 + TMP2)) +
      TMP9 * (+cI * (TMP7 + TMP23))) + (+0.500000000 * (V1[5] * (-cI * (TMP10 *
      TMP23) + cI * (TMP2 * TMP7))) + P2[3] * (-cI * (TMP7 * TMP9) + cI * (TMP4
      * TMP10)))) + (P2[1] * (P2[3] * (TMP4 * - 1. * (+cI * (TMP10 + TMP1)) +
      TMP7 * (+cI * (TMP9 + TMP26))) + (+0.500000000 * (V2[5] * (-cI * (TMP10 *
      TMP26) + cI * (TMP1 * TMP9))) + P1[3] * (-cI * (TMP7 * TMP9) + cI * (TMP4
      * TMP10)))) + (+0.500000000 * (P2[3] * V2[3] * (-cI * (TMP10 * TMP26) +
      cI * (TMP1 * TMP9))) + P1[3] * 0.500000000 * V1[3] * (-cI * (TMP10 *
      TMP23) + cI * (TMP2 * TMP7))))));
  T3[4] = denom * (OM3 * (P3[0] * (P3[2] * (OM3 * (TMP1 * (TMP2 * (TMP4 * -
      0.666666667 * (-2. * cI * (TMP10) + cI * (TMP1 + TMP2)) + (TMP7 *
      1.333333333 * (+0.500000000 * cI * (TMP26) - cI * (TMP9)) + 0.666666667 *
      cI * (TMP9 * TMP23))) + (TMP1 * (TMP9 * 0.666666667 * (+cI * (TMP7 +
      TMP23)) - 0.666666667 * cI * (TMP4 * TMP10)) - 0.666666667 * cI * (TMP10
      * TMP23 * TMP26))) + TMP2 * (TMP2 * (TMP7 * 0.666666667 * (+cI * (TMP9 +
      TMP26)) - 0.666666667 * cI * (TMP4 * TMP10)) - 0.666666667 * cI * (TMP10
      * TMP23 * TMP26))) + (TMP10 * (TMP4 * - 0.333333333 * (-2. * cI * (TMP10)
      + cI * (TMP47 + TMP48)) + (-0.666666667 * cI * (TMP7 * TMP9) -
      0.333333333 * cI * (TMP23 * TMP49 + TMP26 * TMP50))) + (TMP7 * (TMP48 *
      0.333333333 * (+cI * (TMP9 + TMP26)) + (+0.333333333 * cI * (TMP9 * TMP47
      + TMP2 * TMP49))) + (TMP1 * 0.333333333 * (-cI * (TMP4 * TMP48) + cI *
      (TMP9 * TMP50)) + 0.333333333 * (TMP47 * (-cI * (TMP2 * TMP4) + cI *
      (TMP9 * TMP23))))))) + (TMP1 * (TMP9 * (P1[2] * - 1. * (+cI * (TMP7 +
      TMP23)) + (-0.500000000 * cI * (TMP2 * V2[4]) + P2[2] * 0.500000000 *
      (-cI * (TMP23) + 2. * cI * (TMP7)))) + (+0.500000000 * (V1[4] * (-cI *
      (TMP2 * TMP7) + cI * (TMP10 * TMP23))) + TMP4 * (P1[2] * (+cI * (TMP10 +
      TMP2)) + P2[2] * (-cI * (TMP10) + cI * (TMP2))))) + (+0.500000000 *
      (TMP10 * TMP23 * TMP26 * (+cI * (P1[2] + P2[2]))) + TMP2 * 0.500000000 *
      (TMP10 * (+2. * (TMP4 * (-cI * (P1[2]) + cI * (P2[2]))) + cI * (V2[4] *
      TMP26)) + 2. * (TMP7 * (P1[2] * 0.500000000 * (-cI * (TMP26) + 2. * cI *
      (TMP9)) - P2[2] * (+cI * (TMP9 + TMP26)))))))) + P3[2] * (TMP1 * (TMP9 *
      (P1[0] * - 1. * (+cI * (TMP7 + TMP23)) + (-0.500000000 * cI * (TMP2 *
      V2[2]) + P2[0] * 0.500000000 * (-cI * (TMP23) + 2. * cI * (TMP7)))) +
      (+0.500000000 * (V1[2] * (-cI * (TMP2 * TMP7) + cI * (TMP10 * TMP23))) +
      TMP4 * (P1[0] * (+cI * (TMP10 + TMP2)) + P2[0] * (-cI * (TMP10) + cI *
      (TMP2))))) + (+0.500000000 * (TMP10 * TMP23 * TMP26 * (+cI * (P1[0] +
      P2[0]))) + TMP2 * 0.500000000 * (TMP10 * (+2. * (TMP4 * (-cI * (P1[0]) +
      cI * (P2[0]))) + cI * (V2[2] * TMP26)) + 2. * (TMP7 * (P1[0] *
      0.500000000 * (-cI * (TMP26) + 2. * cI * (TMP9)) - P2[0] * (+cI * (TMP9 +
      TMP26)))))))) + (P1[0] * (P1[2] * (TMP4 * - 1. * (+cI * (TMP10 + TMP2)) +
      TMP9 * (+cI * (TMP7 + TMP23))) + (+0.500000000 * (V1[4] * (-cI * (TMP10 *
      TMP23) + cI * (TMP2 * TMP7))) + P2[2] * (-cI * (TMP7 * TMP9) + cI * (TMP4
      * TMP10)))) + (P2[0] * (P2[2] * (TMP4 * - 1. * (+cI * (TMP10 + TMP1)) +
      TMP7 * (+cI * (TMP9 + TMP26))) + (+0.500000000 * (V2[4] * (-cI * (TMP10 *
      TMP26) + cI * (TMP1 * TMP9))) + P1[2] * (-cI * (TMP7 * TMP9) + cI * (TMP4
      * TMP10)))) + (+0.500000000 * (P2[2] * V2[2] * (-cI * (TMP10 * TMP26) +
      cI * (TMP1 * TMP9))) + P1[2] * 0.500000000 * V1[2] * (-cI * (TMP10 *
      TMP23) + cI * (TMP2 * TMP7))))));
  T3[8] = denom * (OM3 * (P3[1] * (P3[2] * (OM3 * (TMP1 * (TMP2 * (TMP4 * -
      0.666666667 * (-2. * cI * (TMP10) + cI * (TMP1 + TMP2)) + (TMP7 *
      1.333333333 * (+0.500000000 * cI * (TMP26) - cI * (TMP9)) + 0.666666667 *
      cI * (TMP9 * TMP23))) + (TMP1 * (TMP9 * 0.666666667 * (+cI * (TMP7 +
      TMP23)) - 0.666666667 * cI * (TMP4 * TMP10)) - 0.666666667 * cI * (TMP10
      * TMP23 * TMP26))) + TMP2 * (TMP2 * (TMP7 * 0.666666667 * (+cI * (TMP9 +
      TMP26)) - 0.666666667 * cI * (TMP4 * TMP10)) - 0.666666667 * cI * (TMP10
      * TMP23 * TMP26))) + (TMP10 * (TMP4 * - 0.333333333 * (-2. * cI * (TMP10)
      + cI * (TMP47 + TMP48)) + (-0.666666667 * cI * (TMP7 * TMP9) -
      0.333333333 * cI * (TMP23 * TMP49 + TMP26 * TMP50))) + (TMP7 * (TMP48 *
      0.333333333 * (+cI * (TMP9 + TMP26)) + (+0.333333333 * cI * (TMP9 * TMP47
      + TMP2 * TMP49))) + (TMP1 * 0.333333333 * (-cI * (TMP4 * TMP48) + cI *
      (TMP9 * TMP50)) + 0.333333333 * (TMP47 * (-cI * (TMP2 * TMP4) + cI *
      (TMP9 * TMP23))))))) + (TMP1 * (TMP9 * (P1[2] * - 1. * (+cI * (TMP7 +
      TMP23)) + (-0.500000000 * cI * (TMP2 * V2[4]) + P2[2] * 0.500000000 *
      (-cI * (TMP23) + 2. * cI * (TMP7)))) + (+0.500000000 * (V1[4] * (-cI *
      (TMP2 * TMP7) + cI * (TMP10 * TMP23))) + TMP4 * (P1[2] * (+cI * (TMP10 +
      TMP2)) + P2[2] * (-cI * (TMP10) + cI * (TMP2))))) + (+0.500000000 *
      (TMP10 * TMP23 * TMP26 * (+cI * (P1[2] + P2[2]))) + TMP2 * 0.500000000 *
      (TMP10 * (+2. * (TMP4 * (-cI * (P1[2]) + cI * (P2[2]))) + cI * (V2[4] *
      TMP26)) + 2. * (TMP7 * (P1[2] * 0.500000000 * (-cI * (TMP26) + 2. * cI *
      (TMP9)) - P2[2] * (+cI * (TMP9 + TMP26)))))))) + P3[2] * (TMP1 * (TMP9 *
      (P1[1] * - 1. * (+cI * (TMP7 + TMP23)) + (-0.500000000 * cI * (TMP2 *
      V2[3]) + P2[1] * 0.500000000 * (-cI * (TMP23) + 2. * cI * (TMP7)))) +
      (+0.500000000 * (V1[3] * (-cI * (TMP2 * TMP7) + cI * (TMP10 * TMP23))) +
      TMP4 * (P1[1] * (+cI * (TMP10 + TMP2)) + P2[1] * (-cI * (TMP10) + cI *
      (TMP2))))) + (+0.500000000 * (TMP10 * TMP23 * TMP26 * (+cI * (P1[1] +
      P2[1]))) + TMP2 * 0.500000000 * (TMP10 * (+2. * (TMP4 * (-cI * (P1[1]) +
      cI * (P2[1]))) + cI * (V2[3] * TMP26)) + 2. * (TMP7 * (P1[1] *
      0.500000000 * (-cI * (TMP26) + 2. * cI * (TMP9)) - P2[1] * (+cI * (TMP9 +
      TMP26)))))))) + (P1[1] * (P1[2] * (TMP4 * - 1. * (+cI * (TMP10 + TMP2)) +
      TMP9 * (+cI * (TMP7 + TMP23))) + (+0.500000000 * (V1[4] * (-cI * (TMP10 *
      TMP23) + cI * (TMP2 * TMP7))) + P2[2] * (-cI * (TMP7 * TMP9) + cI * (TMP4
      * TMP10)))) + (P2[1] * (P2[2] * (TMP4 * - 1. * (+cI * (TMP10 + TMP1)) +
      TMP7 * (+cI * (TMP9 + TMP26))) + (+0.500000000 * (V2[4] * (-cI * (TMP10 *
      TMP26) + cI * (TMP1 * TMP9))) + P1[2] * (-cI * (TMP7 * TMP9) + cI * (TMP4
      * TMP10)))) + (+0.500000000 * (P2[2] * V2[3] * (-cI * (TMP10 * TMP26) +
      cI * (TMP1 * TMP9))) + P1[2] * 0.500000000 * V1[3] * (-cI * (TMP10 *
      TMP23) + cI * (TMP2 * TMP7))))));
  T3[12] = denom * 0.333333333 * (OM3 * (P3[2] * (P3[2] * (OM3 * (TMP1 * (TMP2
      * (TMP4 * - 2. * (-2. * cI * (TMP10) + cI * (TMP1 + TMP2)) + (TMP7 * 2. *
      (-2. * cI * (TMP9) + cI * (TMP26)) + 2. * cI * (TMP9 * TMP23))) + (TMP1 *
      (TMP9 * 2. * (+cI * (TMP7 + TMP23)) - 2. * cI * (TMP4 * TMP10)) - 2. * cI
      * (TMP10 * TMP23 * TMP26))) + TMP2 * (TMP2 * (TMP7 * 2. * (+cI * (TMP9 +
      TMP26)) - 2. * cI * (TMP4 * TMP10)) - 2. * cI * (TMP10 * TMP23 * TMP26)))
      + (TMP10 * (TMP4 * - 1. * (-2. * cI * (TMP10) + cI * (TMP47 + TMP48)) +
      (-2. * cI * (TMP7 * TMP9) - cI * (TMP23 * TMP49 + TMP26 * TMP50))) +
      (TMP7 * (TMP48 * (+cI * (TMP9 + TMP26)) + (+cI * (TMP9 * TMP47 + TMP2 *
      TMP49))) + (TMP1 * (-cI * (TMP4 * TMP48) + cI * (TMP9 * TMP50)) + TMP47 *
      (-cI * (TMP2 * TMP4) + cI * (TMP9 * TMP23)))))) + (TMP1 * (TMP9 * (P1[2]
      * - 6. * (+cI * (TMP7 + TMP23)) + (P2[2] * 3. * (-cI * (TMP23) + 2. * cI
      * (TMP7)) - 3. * cI * (TMP2 * V2[4]))) + (TMP4 * (P1[2] * 6. * (+cI *
      (TMP10 + TMP2)) + 6. * (P2[2] * (-cI * (TMP10) + cI * (TMP2)))) + 3. *
      (V1[4] * (-cI * (TMP2 * TMP7) + cI * (TMP10 * TMP23))))) + (TMP2 * (TMP7
      * (P1[2] * 3. * (-cI * (TMP26) + 2. * cI * (TMP9)) - 6. * (P2[2] * (+cI *
      (TMP9 + TMP26)))) + TMP10 * (TMP4 * 6. * (-cI * (P1[2]) + cI * (P2[2])) +
      3. * cI * (V2[4] * TMP26))) + 3. * (TMP10 * TMP23 * TMP26 * (+cI * (P1[2]
      + P2[2])))))) + (TMP1 * (TMP2 * (TMP4 * (-2. * cI * (TMP10) + cI * (TMP1
      + TMP2)) + (TMP7 * (-cI * (TMP26) + 2. * cI * (TMP9)) - cI * (TMP9 *
      TMP23))) + (TMP1 * (TMP9 * - 1. * (+cI * (TMP7 + TMP23)) + cI * (TMP4 *
      TMP10)) + cI * (TMP10 * TMP23 * TMP26))) + TMP2 * (TMP2 * (TMP7 * - 1. *
      (+cI * (TMP9 + TMP26)) + cI * (TMP4 * TMP10)) + cI * (TMP10 * TMP23 *
      TMP26)))) + (TMP10 * (TMP4 * (P1[2] * 3. * (-cI * (P1[2]) + 2. * cI *
      (P2[2])) + (-1.000000000 * cI * (TMP47 + TMP48) + 2.000000000 * cI *
      (TMP10) - 3. * cI * (P2[2] * P2[2]))) + (TMP23 * - 1.000000000 * (+cI *
      (TMP49) + 3.000000000 * cI * (P1[2] * V1[4])) + (-2.000000000 * cI *
      (TMP7 * TMP9) + TMP26 * - 1.000000000 * (+cI * (TMP50) + 3.000000000 * cI
      * (P2[2] * V2[4]))))) + (TMP7 * (TMP9 * (P1[2] * 3. * (-2. * cI * (P2[2])
      + cI * (P1[2])) + (+1.000000000 * cI * (TMP47 + TMP48) + 3. * cI * (P2[2]
      * P2[2]))) + (+1.000000000 * (TMP2 * (+cI * (TMP49) + 3.000000000 * cI *
      (P1[2] * V1[4]))) + TMP26 * 1.000000000 * (+cI * (TMP48) + 3.000000000 *
      cI * (P2[2] * P2[2])))) + (TMP1 * 1.000000000 * (TMP9 * (+cI * (TMP50) +
      3.000000000 * cI * (P2[2] * V2[4])) + 1.000000000 * (TMP4 * - 1.000000000
      * (+cI * (TMP48) + 3.000000000 * cI * (P2[2] * P2[2])))) + (+1.000000000
      * (TMP23 * TMP9 * (+cI * (TMP47) + 3.000000000 * cI * (P1[2] * P1[2]))) +
      TMP2 * - 1.000000000 * TMP4 * (+cI * (TMP47) + 3.000000000 * cI * (P1[2]
      * P1[2])))))));
  T3[16] = denom * (OM3 * (P3[2] * (P3[3] * (OM3 * (TMP1 * (TMP2 * (TMP4 * -
      0.666666667 * (-2. * cI * (TMP10) + cI * (TMP1 + TMP2)) + (TMP7 *
      1.333333333 * (+0.500000000 * cI * (TMP26) - cI * (TMP9)) + 0.666666667 *
      cI * (TMP9 * TMP23))) + (TMP1 * (TMP9 * 0.666666667 * (+cI * (TMP7 +
      TMP23)) - 0.666666667 * cI * (TMP4 * TMP10)) - 0.666666667 * cI * (TMP10
      * TMP23 * TMP26))) + TMP2 * (TMP2 * (TMP7 * 0.666666667 * (+cI * (TMP9 +
      TMP26)) - 0.666666667 * cI * (TMP4 * TMP10)) - 0.666666667 * cI * (TMP10
      * TMP23 * TMP26))) + (TMP10 * (TMP4 * - 0.333333333 * (-2. * cI * (TMP10)
      + cI * (TMP47 + TMP48)) + (-0.666666667 * cI * (TMP7 * TMP9) -
      0.333333333 * cI * (TMP23 * TMP49 + TMP26 * TMP50))) + (TMP7 * (TMP48 *
      0.333333333 * (+cI * (TMP9 + TMP26)) + (+0.333333333 * cI * (TMP9 * TMP47
      + TMP2 * TMP49))) + (TMP1 * 0.333333333 * (-cI * (TMP4 * TMP48) + cI *
      (TMP9 * TMP50)) + 0.333333333 * (TMP47 * (-cI * (TMP2 * TMP4) + cI *
      (TMP9 * TMP23))))))) + (TMP1 * (TMP9 * (P1[3] * - 1. * (+cI * (TMP7 +
      TMP23)) + (-0.500000000 * cI * (TMP2 * V2[5]) + P2[3] * 0.500000000 *
      (-cI * (TMP23) + 2. * cI * (TMP7)))) + (+0.500000000 * (V1[5] * (-cI *
      (TMP2 * TMP7) + cI * (TMP10 * TMP23))) + TMP4 * (P1[3] * (+cI * (TMP10 +
      TMP2)) + P2[3] * (-cI * (TMP10) + cI * (TMP2))))) + (+0.500000000 *
      (TMP10 * TMP23 * TMP26 * (+cI * (P1[3] + P2[3]))) + TMP2 * 0.500000000 *
      (TMP10 * (+2. * (TMP4 * (-cI * (P1[3]) + cI * (P2[3]))) + cI * (V2[5] *
      TMP26)) + 2. * (TMP7 * (P1[3] * 0.500000000 * (-cI * (TMP26) + 2. * cI *
      (TMP9)) - P2[3] * (+cI * (TMP9 + TMP26)))))))) + P3[3] * (TMP1 * (TMP9 *
      (P1[2] * - 1. * (+cI * (TMP7 + TMP23)) + (-0.500000000 * cI * (TMP2 *
      V2[4]) + P2[2] * 0.500000000 * (-cI * (TMP23) + 2. * cI * (TMP7)))) +
      (+0.500000000 * (V1[4] * (-cI * (TMP2 * TMP7) + cI * (TMP10 * TMP23))) +
      TMP4 * (P1[2] * (+cI * (TMP10 + TMP2)) + P2[2] * (-cI * (TMP10) + cI *
      (TMP2))))) + (+0.500000000 * (TMP10 * TMP23 * TMP26 * (+cI * (P1[2] +
      P2[2]))) + TMP2 * 0.500000000 * (TMP10 * (+2. * (TMP4 * (-cI * (P1[2]) +
      cI * (P2[2]))) + cI * (V2[4] * TMP26)) + 2. * (TMP7 * (P1[2] *
      0.500000000 * (-cI * (TMP26) + 2. * cI * (TMP9)) - P2[2] * (+cI * (TMP9 +
      TMP26)))))))) + (P1[2] * (P1[3] * (TMP4 * - 1. * (+cI * (TMP10 + TMP2)) +
      TMP9 * (+cI * (TMP7 + TMP23))) + (+0.500000000 * (V1[5] * (-cI * (TMP10 *
      TMP23) + cI * (TMP2 * TMP7))) + P2[3] * (-cI * (TMP7 * TMP9) + cI * (TMP4
      * TMP10)))) + (P2[2] * (P2[3] * (TMP4 * - 1. * (+cI * (TMP10 + TMP1)) +
      TMP7 * (+cI * (TMP9 + TMP26))) + (+0.500000000 * (V2[5] * (-cI * (TMP10 *
      TMP26) + cI * (TMP1 * TMP9))) + P1[3] * (-cI * (TMP7 * TMP9) + cI * (TMP4
      * TMP10)))) + (+0.500000000 * (P2[3] * V2[4] * (-cI * (TMP10 * TMP26) +
      cI * (TMP1 * TMP9))) + P1[3] * 0.500000000 * V1[4] * (-cI * (TMP10 *
      TMP23) + cI * (TMP2 * TMP7))))));
  T3[5] = denom * (OM3 * (P3[0] * (P3[3] * (OM3 * (TMP1 * (TMP2 * (TMP4 * -
      0.666666667 * (-2. * cI * (TMP10) + cI * (TMP1 + TMP2)) + (TMP7 *
      1.333333333 * (+0.500000000 * cI * (TMP26) - cI * (TMP9)) + 0.666666667 *
      cI * (TMP9 * TMP23))) + (TMP1 * (TMP9 * 0.666666667 * (+cI * (TMP7 +
      TMP23)) - 0.666666667 * cI * (TMP4 * TMP10)) - 0.666666667 * cI * (TMP10
      * TMP23 * TMP26))) + TMP2 * (TMP2 * (TMP7 * 0.666666667 * (+cI * (TMP9 +
      TMP26)) - 0.666666667 * cI * (TMP4 * TMP10)) - 0.666666667 * cI * (TMP10
      * TMP23 * TMP26))) + (TMP10 * (TMP4 * - 0.333333333 * (-2. * cI * (TMP10)
      + cI * (TMP47 + TMP48)) + (-0.666666667 * cI * (TMP7 * TMP9) -
      0.333333333 * cI * (TMP23 * TMP49 + TMP26 * TMP50))) + (TMP7 * (TMP48 *
      0.333333333 * (+cI * (TMP9 + TMP26)) + (+0.333333333 * cI * (TMP9 * TMP47
      + TMP2 * TMP49))) + (TMP1 * 0.333333333 * (-cI * (TMP4 * TMP48) + cI *
      (TMP9 * TMP50)) + 0.333333333 * (TMP47 * (-cI * (TMP2 * TMP4) + cI *
      (TMP9 * TMP23))))))) + (TMP1 * (TMP9 * (P1[3] * - 1. * (+cI * (TMP7 +
      TMP23)) + (-0.500000000 * cI * (TMP2 * V2[5]) + P2[3] * 0.500000000 *
      (-cI * (TMP23) + 2. * cI * (TMP7)))) + (+0.500000000 * (V1[5] * (-cI *
      (TMP2 * TMP7) + cI * (TMP10 * TMP23))) + TMP4 * (P1[3] * (+cI * (TMP10 +
      TMP2)) + P2[3] * (-cI * (TMP10) + cI * (TMP2))))) + (+0.500000000 *
      (TMP10 * TMP23 * TMP26 * (+cI * (P1[3] + P2[3]))) + TMP2 * 0.500000000 *
      (TMP10 * (+2. * (TMP4 * (-cI * (P1[3]) + cI * (P2[3]))) + cI * (V2[5] *
      TMP26)) + 2. * (TMP7 * (P1[3] * 0.500000000 * (-cI * (TMP26) + 2. * cI *
      (TMP9)) - P2[3] * (+cI * (TMP9 + TMP26)))))))) + P3[3] * (TMP1 * (TMP9 *
      (P1[0] * - 1. * (+cI * (TMP7 + TMP23)) + (-0.500000000 * cI * (TMP2 *
      V2[2]) + P2[0] * 0.500000000 * (-cI * (TMP23) + 2. * cI * (TMP7)))) +
      (+0.500000000 * (V1[2] * (-cI * (TMP2 * TMP7) + cI * (TMP10 * TMP23))) +
      TMP4 * (P1[0] * (+cI * (TMP10 + TMP2)) + P2[0] * (-cI * (TMP10) + cI *
      (TMP2))))) + (+0.500000000 * (TMP10 * TMP23 * TMP26 * (+cI * (P1[0] +
      P2[0]))) + TMP2 * 0.500000000 * (TMP10 * (+2. * (TMP4 * (-cI * (P1[0]) +
      cI * (P2[0]))) + cI * (V2[2] * TMP26)) + 2. * (TMP7 * (P1[0] *
      0.500000000 * (-cI * (TMP26) + 2. * cI * (TMP9)) - P2[0] * (+cI * (TMP9 +
      TMP26)))))))) + (P1[0] * (P1[3] * (TMP4 * - 1. * (+cI * (TMP10 + TMP2)) +
      TMP9 * (+cI * (TMP7 + TMP23))) + (+0.500000000 * (V1[5] * (-cI * (TMP10 *
      TMP23) + cI * (TMP2 * TMP7))) + P2[3] * (-cI * (TMP7 * TMP9) + cI * (TMP4
      * TMP10)))) + (P2[0] * (P2[3] * (TMP4 * - 1. * (+cI * (TMP10 + TMP1)) +
      TMP7 * (+cI * (TMP9 + TMP26))) + (+0.500000000 * (V2[5] * (-cI * (TMP10 *
      TMP26) + cI * (TMP1 * TMP9))) + P1[3] * (-cI * (TMP7 * TMP9) + cI * (TMP4
      * TMP10)))) + (+0.500000000 * (P2[3] * V2[2] * (-cI * (TMP10 * TMP26) +
      cI * (TMP1 * TMP9))) + P1[3] * 0.500000000 * V1[2] * (-cI * (TMP10 *
      TMP23) + cI * (TMP2 * TMP7))))));
  T3[9] = denom * (OM3 * (P3[1] * (P3[3] * (OM3 * (TMP1 * (TMP2 * (TMP4 * -
      0.666666667 * (-2. * cI * (TMP10) + cI * (TMP1 + TMP2)) + (TMP7 *
      1.333333333 * (+0.500000000 * cI * (TMP26) - cI * (TMP9)) + 0.666666667 *
      cI * (TMP9 * TMP23))) + (TMP1 * (TMP9 * 0.666666667 * (+cI * (TMP7 +
      TMP23)) - 0.666666667 * cI * (TMP4 * TMP10)) - 0.666666667 * cI * (TMP10
      * TMP23 * TMP26))) + TMP2 * (TMP2 * (TMP7 * 0.666666667 * (+cI * (TMP9 +
      TMP26)) - 0.666666667 * cI * (TMP4 * TMP10)) - 0.666666667 * cI * (TMP10
      * TMP23 * TMP26))) + (TMP10 * (TMP4 * - 0.333333333 * (-2. * cI * (TMP10)
      + cI * (TMP47 + TMP48)) + (-0.666666667 * cI * (TMP7 * TMP9) -
      0.333333333 * cI * (TMP23 * TMP49 + TMP26 * TMP50))) + (TMP7 * (TMP48 *
      0.333333333 * (+cI * (TMP9 + TMP26)) + (+0.333333333 * cI * (TMP9 * TMP47
      + TMP2 * TMP49))) + (TMP1 * 0.333333333 * (-cI * (TMP4 * TMP48) + cI *
      (TMP9 * TMP50)) + 0.333333333 * (TMP47 * (-cI * (TMP2 * TMP4) + cI *
      (TMP9 * TMP23))))))) + (TMP1 * (TMP9 * (P1[3] * - 1. * (+cI * (TMP7 +
      TMP23)) + (-0.500000000 * cI * (TMP2 * V2[5]) + P2[3] * 0.500000000 *
      (-cI * (TMP23) + 2. * cI * (TMP7)))) + (+0.500000000 * (V1[5] * (-cI *
      (TMP2 * TMP7) + cI * (TMP10 * TMP23))) + TMP4 * (P1[3] * (+cI * (TMP10 +
      TMP2)) + P2[3] * (-cI * (TMP10) + cI * (TMP2))))) + (+0.500000000 *
      (TMP10 * TMP23 * TMP26 * (+cI * (P1[3] + P2[3]))) + TMP2 * 0.500000000 *
      (TMP10 * (+2. * (TMP4 * (-cI * (P1[3]) + cI * (P2[3]))) + cI * (V2[5] *
      TMP26)) + 2. * (TMP7 * (P1[3] * 0.500000000 * (-cI * (TMP26) + 2. * cI *
      (TMP9)) - P2[3] * (+cI * (TMP9 + TMP26)))))))) + P3[3] * (TMP1 * (TMP9 *
      (P1[1] * - 1. * (+cI * (TMP7 + TMP23)) + (-0.500000000 * cI * (TMP2 *
      V2[3]) + P2[1] * 0.500000000 * (-cI * (TMP23) + 2. * cI * (TMP7)))) +
      (+0.500000000 * (V1[3] * (-cI * (TMP2 * TMP7) + cI * (TMP10 * TMP23))) +
      TMP4 * (P1[1] * (+cI * (TMP10 + TMP2)) + P2[1] * (-cI * (TMP10) + cI *
      (TMP2))))) + (+0.500000000 * (TMP10 * TMP23 * TMP26 * (+cI * (P1[1] +
      P2[1]))) + TMP2 * 0.500000000 * (TMP10 * (+2. * (TMP4 * (-cI * (P1[1]) +
      cI * (P2[1]))) + cI * (V2[3] * TMP26)) + 2. * (TMP7 * (P1[1] *
      0.500000000 * (-cI * (TMP26) + 2. * cI * (TMP9)) - P2[1] * (+cI * (TMP9 +
      TMP26)))))))) + (P1[1] * (P1[3] * (TMP4 * - 1. * (+cI * (TMP10 + TMP2)) +
      TMP9 * (+cI * (TMP7 + TMP23))) + (+0.500000000 * (V1[5] * (-cI * (TMP10 *
      TMP23) + cI * (TMP2 * TMP7))) + P2[3] * (-cI * (TMP7 * TMP9) + cI * (TMP4
      * TMP10)))) + (P2[1] * (P2[3] * (TMP4 * - 1. * (+cI * (TMP10 + TMP1)) +
      TMP7 * (+cI * (TMP9 + TMP26))) + (+0.500000000 * (V2[5] * (-cI * (TMP10 *
      TMP26) + cI * (TMP1 * TMP9))) + P1[3] * (-cI * (TMP7 * TMP9) + cI * (TMP4
      * TMP10)))) + (+0.500000000 * (P2[3] * V2[3] * (-cI * (TMP10 * TMP26) +
      cI * (TMP1 * TMP9))) + P1[3] * 0.500000000 * V1[3] * (-cI * (TMP10 *
      TMP23) + cI * (TMP2 * TMP7))))));
  T3[13] = denom * (OM3 * (P3[2] * (P3[3] * (OM3 * (TMP1 * (TMP2 * (TMP4 * -
      0.666666667 * (-2. * cI * (TMP10) + cI * (TMP1 + TMP2)) + (TMP7 *
      1.333333333 * (+0.500000000 * cI * (TMP26) - cI * (TMP9)) + 0.666666667 *
      cI * (TMP9 * TMP23))) + (TMP1 * (TMP9 * 0.666666667 * (+cI * (TMP7 +
      TMP23)) - 0.666666667 * cI * (TMP4 * TMP10)) - 0.666666667 * cI * (TMP10
      * TMP23 * TMP26))) + TMP2 * (TMP2 * (TMP7 * 0.666666667 * (+cI * (TMP9 +
      TMP26)) - 0.666666667 * cI * (TMP4 * TMP10)) - 0.666666667 * cI * (TMP10
      * TMP23 * TMP26))) + (TMP10 * (TMP4 * - 0.333333333 * (-2. * cI * (TMP10)
      + cI * (TMP47 + TMP48)) + (-0.666666667 * cI * (TMP7 * TMP9) -
      0.333333333 * cI * (TMP23 * TMP49 + TMP26 * TMP50))) + (TMP7 * (TMP48 *
      0.333333333 * (+cI * (TMP9 + TMP26)) + (+0.333333333 * cI * (TMP9 * TMP47
      + TMP2 * TMP49))) + (TMP1 * 0.333333333 * (-cI * (TMP4 * TMP48) + cI *
      (TMP9 * TMP50)) + 0.333333333 * (TMP47 * (-cI * (TMP2 * TMP4) + cI *
      (TMP9 * TMP23))))))) + (TMP1 * (TMP9 * (P1[3] * - 1. * (+cI * (TMP7 +
      TMP23)) + (-0.500000000 * cI * (TMP2 * V2[5]) + P2[3] * 0.500000000 *
      (-cI * (TMP23) + 2. * cI * (TMP7)))) + (+0.500000000 * (V1[5] * (-cI *
      (TMP2 * TMP7) + cI * (TMP10 * TMP23))) + TMP4 * (P1[3] * (+cI * (TMP10 +
      TMP2)) + P2[3] * (-cI * (TMP10) + cI * (TMP2))))) + (+0.500000000 *
      (TMP10 * TMP23 * TMP26 * (+cI * (P1[3] + P2[3]))) + TMP2 * 0.500000000 *
      (TMP10 * (+2. * (TMP4 * (-cI * (P1[3]) + cI * (P2[3]))) + cI * (V2[5] *
      TMP26)) + 2. * (TMP7 * (P1[3] * 0.500000000 * (-cI * (TMP26) + 2. * cI *
      (TMP9)) - P2[3] * (+cI * (TMP9 + TMP26)))))))) + P3[3] * (TMP1 * (TMP9 *
      (P1[2] * - 1. * (+cI * (TMP7 + TMP23)) + (-0.500000000 * cI * (TMP2 *
      V2[4]) + P2[2] * 0.500000000 * (-cI * (TMP23) + 2. * cI * (TMP7)))) +
      (+0.500000000 * (V1[4] * (-cI * (TMP2 * TMP7) + cI * (TMP10 * TMP23))) +
      TMP4 * (P1[2] * (+cI * (TMP10 + TMP2)) + P2[2] * (-cI * (TMP10) + cI *
      (TMP2))))) + (+0.500000000 * (TMP10 * TMP23 * TMP26 * (+cI * (P1[2] +
      P2[2]))) + TMP2 * 0.500000000 * (TMP10 * (+2. * (TMP4 * (-cI * (P1[2]) +
      cI * (P2[2]))) + cI * (V2[4] * TMP26)) + 2. * (TMP7 * (P1[2] *
      0.500000000 * (-cI * (TMP26) + 2. * cI * (TMP9)) - P2[2] * (+cI * (TMP9 +
      TMP26)))))))) + (P1[2] * (P1[3] * (TMP4 * - 1. * (+cI * (TMP10 + TMP2)) +
      TMP9 * (+cI * (TMP7 + TMP23))) + (+0.500000000 * (V1[5] * (-cI * (TMP10 *
      TMP23) + cI * (TMP2 * TMP7))) + P2[3] * (-cI * (TMP7 * TMP9) + cI * (TMP4
      * TMP10)))) + (P2[2] * (P2[3] * (TMP4 * - 1. * (+cI * (TMP10 + TMP1)) +
      TMP7 * (+cI * (TMP9 + TMP26))) + (+0.500000000 * (V2[5] * (-cI * (TMP10 *
      TMP26) + cI * (TMP1 * TMP9))) + P1[3] * (-cI * (TMP7 * TMP9) + cI * (TMP4
      * TMP10)))) + (+0.500000000 * (P2[3] * V2[4] * (-cI * (TMP10 * TMP26) +
      cI * (TMP1 * TMP9))) + P1[3] * 0.500000000 * V1[4] * (-cI * (TMP10 *
      TMP23) + cI * (TMP2 * TMP7))))));
  T3[17] = denom * 0.333333333 * (OM3 * (P3[3] * (P3[3] * (OM3 * (TMP1 * (TMP2
      * (TMP4 * - 2. * (-2. * cI * (TMP10) + cI * (TMP1 + TMP2)) + (TMP7 * 2. *
      (-2. * cI * (TMP9) + cI * (TMP26)) + 2. * cI * (TMP9 * TMP23))) + (TMP1 *
      (TMP9 * 2. * (+cI * (TMP7 + TMP23)) - 2. * cI * (TMP4 * TMP10)) - 2. * cI
      * (TMP10 * TMP23 * TMP26))) + TMP2 * (TMP2 * (TMP7 * 2. * (+cI * (TMP9 +
      TMP26)) - 2. * cI * (TMP4 * TMP10)) - 2. * cI * (TMP10 * TMP23 * TMP26)))
      + (TMP10 * (TMP4 * - 1. * (-2. * cI * (TMP10) + cI * (TMP47 + TMP48)) +
      (-2. * cI * (TMP7 * TMP9) - cI * (TMP23 * TMP49 + TMP26 * TMP50))) +
      (TMP7 * (TMP48 * (+cI * (TMP9 + TMP26)) + (+cI * (TMP9 * TMP47 + TMP2 *
      TMP49))) + (TMP1 * (-cI * (TMP4 * TMP48) + cI * (TMP9 * TMP50)) + TMP47 *
      (-cI * (TMP2 * TMP4) + cI * (TMP9 * TMP23)))))) + (TMP1 * (TMP9 * (P1[3]
      * - 6. * (+cI * (TMP7 + TMP23)) + (P2[3] * 3. * (-cI * (TMP23) + 2. * cI
      * (TMP7)) - 3. * cI * (TMP2 * V2[5]))) + (TMP4 * (P1[3] * 6. * (+cI *
      (TMP10 + TMP2)) + 6. * (P2[3] * (-cI * (TMP10) + cI * (TMP2)))) + 3. *
      (V1[5] * (-cI * (TMP2 * TMP7) + cI * (TMP10 * TMP23))))) + (TMP2 * (TMP7
      * (P1[3] * 3. * (-cI * (TMP26) + 2. * cI * (TMP9)) - 6. * (P2[3] * (+cI *
      (TMP9 + TMP26)))) + TMP10 * (TMP4 * 6. * (-cI * (P1[3]) + cI * (P2[3])) +
      3. * cI * (V2[5] * TMP26))) + 3. * (TMP10 * TMP23 * TMP26 * (+cI * (P1[3]
      + P2[3])))))) + (TMP1 * (TMP2 * (TMP4 * (-2. * cI * (TMP10) + cI * (TMP1
      + TMP2)) + (TMP7 * (-cI * (TMP26) + 2. * cI * (TMP9)) - cI * (TMP9 *
      TMP23))) + (TMP1 * (TMP9 * - 1. * (+cI * (TMP7 + TMP23)) + cI * (TMP4 *
      TMP10)) + cI * (TMP10 * TMP23 * TMP26))) + TMP2 * (TMP2 * (TMP7 * - 1. *
      (+cI * (TMP9 + TMP26)) + cI * (TMP4 * TMP10)) + cI * (TMP10 * TMP23 *
      TMP26)))) + (TMP10 * (TMP4 * (P1[3] * 3. * (-cI * (P1[3]) + 2. * cI *
      (P2[3])) + (-1.000000000 * cI * (TMP47 + TMP48) + 2.000000000 * cI *
      (TMP10) - 3. * cI * (P2[3] * P2[3]))) + (TMP23 * - 1.000000000 * (+cI *
      (TMP49) + 3.000000000 * cI * (P1[3] * V1[5])) + (-2.000000000 * cI *
      (TMP7 * TMP9) + TMP26 * - 1.000000000 * (+cI * (TMP50) + 3.000000000 * cI
      * (P2[3] * V2[5]))))) + (TMP7 * (TMP9 * (P1[3] * 3. * (-2. * cI * (P2[3])
      + cI * (P1[3])) + (+1.000000000 * cI * (TMP47 + TMP48) + 3. * cI * (P2[3]
      * P2[3]))) + (+1.000000000 * (TMP2 * (+cI * (TMP49) + 3.000000000 * cI *
      (P1[3] * V1[5]))) + TMP26 * 1.000000000 * (+cI * (TMP48) + 3.000000000 *
      cI * (P2[3] * P2[3])))) + (TMP1 * 1.000000000 * (TMP9 * (+cI * (TMP50) +
      3.000000000 * cI * (P2[3] * V2[5])) + 1.000000000 * (TMP4 * - 1.000000000
      * (+cI * (TMP48) + 3.000000000 * cI * (P2[3] * P2[3])))) + (+1.000000000
      * (TMP23 * TMP9 * (+cI * (TMP47) + 3.000000000 * cI * (P1[3] * P1[3]))) +
      TMP2 * - 1.000000000 * TMP4 * (+cI * (TMP47) + 3.000000000 * cI * (P1[3]
      * P1[3])))))));
}


void VVT13_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP2; 
  complex<double> TMP1; 
  double P1[4]; 
  complex<double> TMP22; 
  complex<double> TMP36; 
  double P2[4]; 
  complex<double> TMP23; 
  double P3[4]; 
  complex<double> TMP21; 
  complex<double> TMP26; 
  complex<double> TMP28; 
  complex<double> TMP27; 
  complex<double> TMP24; 
  complex<double> TMP25; 
  complex<double> TMP35; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = T3[0].real(); 
  P3[1] = T3[1].real(); 
  P3[2] = T3[1].imag(); 
  P3[3] = T3[0].imag(); 
  TMP24 = (P2[0] * - 1. * (V1[3] * T3[6] + V1[4] * T3[10] + V1[5] * T3[14] -
      V1[2] * T3[2]) + (P2[1] * (V1[3] * T3[7] + V1[4] * T3[11] + V1[5] *
      T3[15] - V1[2] * T3[3]) + (P2[2] * (V1[3] * T3[8] + V1[4] * T3[12] +
      V1[5] * T3[16] - V1[2] * T3[4]) + P2[3] * (V1[3] * T3[9] + V1[4] * T3[13]
      + V1[5] * T3[17] - V1[2] * T3[5]))));
  TMP25 = (P2[0] * - 1. * (V1[3] * T3[3] + V1[4] * T3[4] + V1[5] * T3[5] -
      V1[2] * T3[2]) + (P2[1] * (V1[3] * T3[7] + V1[4] * T3[8] + V1[5] * T3[9]
      - V1[2] * T3[6]) + (P2[2] * (V1[3] * T3[11] + V1[4] * T3[12] + V1[5] *
      T3[13] - V1[2] * T3[10]) + P2[3] * (V1[3] * T3[15] + V1[4] * T3[16] +
      V1[5] * T3[17] - V1[2] * T3[14]))));
  TMP26 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP27 = (P1[0] * - 1. * (V2[3] * T3[6] + V2[4] * T3[10] + V2[5] * T3[14] -
      V2[2] * T3[2]) + (P1[1] * (V2[3] * T3[7] + V2[4] * T3[11] + V2[5] *
      T3[15] - V2[2] * T3[3]) + (P1[2] * (V2[3] * T3[8] + V2[4] * T3[12] +
      V2[5] * T3[16] - V2[2] * T3[4]) + P1[3] * (V2[3] * T3[9] + V2[4] * T3[13]
      + V2[5] * T3[17] - V2[2] * T3[5]))));
  TMP21 = (P1[0] * - 1. * (P2[1] * T3[6] + P2[2] * T3[10] + P2[3] * T3[14] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[11] + P2[3] *
      T3[15] - P2[0] * T3[3]) + (P1[2] * (P2[1] * T3[8] + P2[2] * T3[12] +
      P2[3] * T3[16] - P2[0] * T3[4]) + P1[3] * (P2[1] * T3[9] + P2[2] * T3[13]
      + P2[3] * T3[17] - P2[0] * T3[5]))));
  TMP22 = (P1[0] * - 1. * (P2[1] * T3[3] + P2[2] * T3[4] + P2[3] * T3[5] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[8] + P2[3] * T3[9]
      - P2[0] * T3[6]) + (P1[2] * (P2[1] * T3[11] + P2[2] * T3[12] + P2[3] *
      T3[13] - P2[0] * T3[10]) + P1[3] * (P2[1] * T3[15] + P2[2] * T3[16] +
      P2[3] * T3[17] - P2[0] * T3[14]))));
  TMP23 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP28 = (P1[0] * - 1. * (V2[3] * T3[3] + V2[4] * T3[4] + V2[5] * T3[5] -
      V2[2] * T3[2]) + (P1[1] * (V2[3] * T3[7] + V2[4] * T3[8] + V2[5] * T3[9]
      - V2[2] * T3[6]) + (P1[2] * (V2[3] * T3[11] + V2[4] * T3[12] + V2[5] *
      T3[13] - V2[2] * T3[10]) + P1[3] * (V2[3] * T3[15] + V2[4] * T3[16] +
      V2[5] * T3[17] - V2[2] * T3[14]))));
  TMP1 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP36 = (V1[2] * - 1. * (V2[3] * T3[3] + V2[4] * T3[4] + V2[5] * T3[5] -
      V2[2] * T3[2]) + (V1[3] * (V2[3] * T3[7] + V2[4] * T3[8] + V2[5] * T3[9]
      - V2[2] * T3[6]) + (V1[4] * (V2[3] * T3[11] + V2[4] * T3[12] + V2[5] *
      T3[13] - V2[2] * T3[10]) + V1[5] * (V2[3] * T3[15] + V2[4] * T3[16] +
      V2[5] * T3[17] - V2[2] * T3[14]))));
  TMP35 = (V1[2] * - 1. * (V2[3] * T3[6] + V2[4] * T3[10] + V2[5] * T3[14] -
      V2[2] * T3[2]) + (V1[3] * (V2[3] * T3[7] + V2[4] * T3[11] + V2[5] *
      T3[15] - V2[2] * T3[3]) + (V1[4] * (V2[3] * T3[8] + V2[4] * T3[12] +
      V2[5] * T3[16] - V2[2] * T3[4]) + V1[5] * (V2[3] * T3[9] + V2[4] * T3[13]
      + V2[5] * T3[17] - V2[2] * T3[5]))));
  TMP2 = (P3[0] * P2[0] - P3[1] * P2[1] - P3[2] * P2[2] - P3[3] * P2[3]); 
  vertex = COUP * (TMP1 * (TMP2 * - 1. * (+cI * (TMP35 + TMP36)) + TMP23 * (+cI
      * (TMP24 + TMP25))) + TMP26 * (TMP2 * (+cI * (TMP27 + TMP28)) - TMP23 *
      (+cI * (TMP21 + TMP22))));
}


void VVT9_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP2; 
  complex<double> TMP1; 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP23; 
  complex<double> TMP7; 
  double P3[4]; 
  double OM3; 
  complex<double> denom; 
  complex<double> TMP26; 
  complex<double> TMP9; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +V1[0] + V2[0]; 
  T3[1] = +V1[1] + V2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP26 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP23 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP9 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP7 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP1 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP2 = (P3[0] * P2[0] - P3[1] * P2[1] - P3[2] * P2[2] - P3[3] * P2[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * 2. * (OM3 * (TMP23 * (TMP26 * (P3[0] * (OM3 * 0.666666667 *
      P3[0] * (+cI * (TMP2 + TMP1)) + (-cI * (P2[0] + P1[0]))) + (+0.333333333
      * cI * (TMP2 + TMP1))) + P3[0] * (-cI * (TMP2 * V1[2]) + 0.333333333 * cI
      * (P3[0] * TMP9))) + P3[0] * TMP26 * (-cI * (TMP1 * V2[2]) + 0.333333333
      * cI * (P3[0] * TMP7))) + (TMP23 * (-0.333333333 * cI * (TMP9) + cI *
      (P2[0] * V1[2])) + TMP26 * (-0.333333333 * cI * (TMP7) + cI * (P1[0] *
      V2[2]))));
  T3[6] = denom * (OM3 * (TMP23 * (TMP26 * (P3[0] * (OM3 * 1.333333333 * P3[1]
      * (+cI * (TMP2 + TMP1)) + (-cI * (P2[1] + P1[1]))) - P3[1] * (+cI *
      (P2[0] + P1[0]))) + (P3[0] * (-cI * (TMP2 * V1[3]) + 0.666666667 * cI *
      (P3[1] * TMP9)) - cI * (P3[1] * TMP2 * V1[2]))) + TMP26 * (P3[0] * (-cI *
      (TMP1 * V2[3]) + 0.666666667 * cI * (P3[1] * TMP7)) - cI * (P3[1] * TMP1
      * V2[2]))) + (TMP23 * (+cI * (P2[0] * V1[3] + P2[1] * V1[2])) + TMP26 *
      (+cI * (P1[0] * V2[3] + P1[1] * V2[2]))));
  T3[10] = denom * (OM3 * (TMP23 * (TMP26 * (P3[0] * (OM3 * 1.333333333 * P3[2]
      * (+cI * (TMP2 + TMP1)) + (-cI * (P2[2] + P1[2]))) - P3[2] * (+cI *
      (P2[0] + P1[0]))) + (P3[0] * (-cI * (TMP2 * V1[4]) + 0.666666667 * cI *
      (P3[2] * TMP9)) - cI * (P3[2] * TMP2 * V1[2]))) + TMP26 * (P3[0] * (-cI *
      (TMP1 * V2[4]) + 0.666666667 * cI * (P3[2] * TMP7)) - cI * (P3[2] * TMP1
      * V2[2]))) + (TMP23 * (+cI * (P2[0] * V1[4] + P2[2] * V1[2])) + TMP26 *
      (+cI * (P1[0] * V2[4] + P1[2] * V2[2]))));
  T3[14] = denom * (OM3 * (TMP23 * (TMP26 * (P3[0] * (OM3 * 1.333333333 * P3[3]
      * (+cI * (TMP2 + TMP1)) + (-cI * (P2[3] + P1[3]))) - P3[3] * (+cI *
      (P2[0] + P1[0]))) + (P3[0] * (-cI * (TMP2 * V1[5]) + 0.666666667 * cI *
      (P3[3] * TMP9)) - cI * (P3[3] * TMP2 * V1[2]))) + TMP26 * (P3[0] * (-cI *
      (TMP1 * V2[5]) + 0.666666667 * cI * (P3[3] * TMP7)) - cI * (P3[3] * TMP1
      * V2[2]))) + (TMP23 * (+cI * (P2[0] * V1[5] + P2[3] * V1[2])) + TMP26 *
      (+cI * (P1[0] * V2[5] + P1[3] * V2[2]))));
  T3[3] = denom * (OM3 * (TMP23 * (TMP26 * (P3[0] * (OM3 * 1.333333333 * P3[1]
      * (+cI * (TMP2 + TMP1)) + (-cI * (P2[1] + P1[1]))) - P3[1] * (+cI *
      (P2[0] + P1[0]))) + (P3[0] * (-cI * (TMP2 * V1[3]) + 0.666666667 * cI *
      (P3[1] * TMP9)) - cI * (P3[1] * TMP2 * V1[2]))) + TMP26 * (P3[0] * (-cI *
      (TMP1 * V2[3]) + 0.666666667 * cI * (P3[1] * TMP7)) - cI * (P3[1] * TMP1
      * V2[2]))) + (TMP23 * (+cI * (P2[1] * V1[2] + P2[0] * V1[3])) + TMP26 *
      (+cI * (P1[1] * V2[2] + P1[0] * V2[3]))));
  T3[7] = denom * 2. * (OM3 * (TMP23 * (TMP26 * (P3[1] * (OM3 * 0.666666667 *
      P3[1] * (+cI * (TMP2 + TMP1)) + (-cI * (P2[1] + P1[1]))) + (-0.333333333
      * cI * (TMP2 + TMP1))) + P3[1] * (-cI * (TMP2 * V1[3]) + 0.333333333 * cI
      * (P3[1] * TMP9))) + P3[1] * TMP26 * (-cI * (TMP1 * V2[3]) + 0.333333333
      * cI * (P3[1] * TMP7))) + (TMP23 * (+cI * (P2[1] * V1[3]) + 0.333333333 *
      cI * (TMP9)) + TMP26 * (+cI * (P1[1] * V2[3]) + 0.333333333 * cI *
      (TMP7))));
  T3[11] = denom * (OM3 * (TMP23 * (TMP26 * (P3[1] * (OM3 * 1.333333333 * P3[2]
      * (+cI * (TMP2 + TMP1)) + (-cI * (P2[2] + P1[2]))) - P3[2] * (+cI *
      (P2[1] + P1[1]))) + (P3[1] * (-cI * (TMP2 * V1[4]) + 0.666666667 * cI *
      (P3[2] * TMP9)) - cI * (P3[2] * TMP2 * V1[3]))) + TMP26 * (P3[1] * (-cI *
      (TMP1 * V2[4]) + 0.666666667 * cI * (P3[2] * TMP7)) - cI * (P3[2] * TMP1
      * V2[3]))) + (TMP23 * (+cI * (P2[1] * V1[4] + P2[2] * V1[3])) + TMP26 *
      (+cI * (P1[1] * V2[4] + P1[2] * V2[3]))));
  T3[15] = denom * (OM3 * (TMP23 * (TMP26 * (P3[1] * (OM3 * 1.333333333 * P3[3]
      * (+cI * (TMP2 + TMP1)) + (-cI * (P2[3] + P1[3]))) - P3[3] * (+cI *
      (P2[1] + P1[1]))) + (P3[1] * (-cI * (TMP2 * V1[5]) + 0.666666667 * cI *
      (P3[3] * TMP9)) - cI * (P3[3] * TMP2 * V1[3]))) + TMP26 * (P3[1] * (-cI *
      (TMP1 * V2[5]) + 0.666666667 * cI * (P3[3] * TMP7)) - cI * (P3[3] * TMP1
      * V2[3]))) + (TMP23 * (+cI * (P2[1] * V1[5] + P2[3] * V1[3])) + TMP26 *
      (+cI * (P1[1] * V2[5] + P1[3] * V2[3]))));
  T3[4] = denom * (OM3 * (TMP23 * (TMP26 * (P3[0] * (OM3 * 1.333333333 * P3[2]
      * (+cI * (TMP2 + TMP1)) + (-cI * (P2[2] + P1[2]))) - P3[2] * (+cI *
      (P2[0] + P1[0]))) + (P3[0] * (-cI * (TMP2 * V1[4]) + 0.666666667 * cI *
      (P3[2] * TMP9)) - cI * (P3[2] * TMP2 * V1[2]))) + TMP26 * (P3[0] * (-cI *
      (TMP1 * V2[4]) + 0.666666667 * cI * (P3[2] * TMP7)) - cI * (P3[2] * TMP1
      * V2[2]))) + (TMP23 * (+cI * (P2[2] * V1[2] + P2[0] * V1[4])) + TMP26 *
      (+cI * (P1[2] * V2[2] + P1[0] * V2[4]))));
  T3[8] = denom * (OM3 * (TMP23 * (TMP26 * (P3[1] * (OM3 * 1.333333333 * P3[2]
      * (+cI * (TMP2 + TMP1)) + (-cI * (P2[2] + P1[2]))) - P3[2] * (+cI *
      (P2[1] + P1[1]))) + (P3[1] * (-cI * (TMP2 * V1[4]) + 0.666666667 * cI *
      (P3[2] * TMP9)) - cI * (P3[2] * TMP2 * V1[3]))) + TMP26 * (P3[1] * (-cI *
      (TMP1 * V2[4]) + 0.666666667 * cI * (P3[2] * TMP7)) - cI * (P3[2] * TMP1
      * V2[3]))) + (TMP23 * (+cI * (P2[2] * V1[3] + P2[1] * V1[4])) + TMP26 *
      (+cI * (P1[2] * V2[3] + P1[1] * V2[4]))));
  T3[12] = denom * 2. * (OM3 * (TMP23 * (TMP26 * (P3[2] * (OM3 * 0.666666667 *
      P3[2] * (+cI * (TMP2 + TMP1)) + (-cI * (P2[2] + P1[2]))) + (-0.333333333
      * cI * (TMP2 + TMP1))) + P3[2] * (-cI * (TMP2 * V1[4]) + 0.333333333 * cI
      * (P3[2] * TMP9))) + P3[2] * TMP26 * (-cI * (TMP1 * V2[4]) + 0.333333333
      * cI * (P3[2] * TMP7))) + (TMP23 * (+cI * (P2[2] * V1[4]) + 0.333333333 *
      cI * (TMP9)) + TMP26 * (+cI * (P1[2] * V2[4]) + 0.333333333 * cI *
      (TMP7))));
  T3[16] = denom * (OM3 * (TMP23 * (TMP26 * (P3[2] * (OM3 * 1.333333333 * P3[3]
      * (+cI * (TMP2 + TMP1)) + (-cI * (P2[3] + P1[3]))) - P3[3] * (+cI *
      (P2[2] + P1[2]))) + (P3[2] * (-cI * (TMP2 * V1[5]) + 0.666666667 * cI *
      (P3[3] * TMP9)) - cI * (P3[3] * TMP2 * V1[4]))) + TMP26 * (P3[2] * (-cI *
      (TMP1 * V2[5]) + 0.666666667 * cI * (P3[3] * TMP7)) - cI * (P3[3] * TMP1
      * V2[4]))) + (TMP23 * (+cI * (P2[2] * V1[5] + P2[3] * V1[4])) + TMP26 *
      (+cI * (P1[2] * V2[5] + P1[3] * V2[4]))));
  T3[5] = denom * (OM3 * (TMP23 * (TMP26 * (P3[0] * (OM3 * 1.333333333 * P3[3]
      * (+cI * (TMP2 + TMP1)) + (-cI * (P2[3] + P1[3]))) - P3[3] * (+cI *
      (P2[0] + P1[0]))) + (P3[0] * (-cI * (TMP2 * V1[5]) + 0.666666667 * cI *
      (P3[3] * TMP9)) - cI * (P3[3] * TMP2 * V1[2]))) + TMP26 * (P3[0] * (-cI *
      (TMP1 * V2[5]) + 0.666666667 * cI * (P3[3] * TMP7)) - cI * (P3[3] * TMP1
      * V2[2]))) + (TMP23 * (+cI * (P2[3] * V1[2] + P2[0] * V1[5])) + TMP26 *
      (+cI * (P1[3] * V2[2] + P1[0] * V2[5]))));
  T3[9] = denom * (OM3 * (TMP23 * (TMP26 * (P3[1] * (OM3 * 1.333333333 * P3[3]
      * (+cI * (TMP2 + TMP1)) + (-cI * (P2[3] + P1[3]))) - P3[3] * (+cI *
      (P2[1] + P1[1]))) + (P3[1] * (-cI * (TMP2 * V1[5]) + 0.666666667 * cI *
      (P3[3] * TMP9)) - cI * (P3[3] * TMP2 * V1[3]))) + TMP26 * (P3[1] * (-cI *
      (TMP1 * V2[5]) + 0.666666667 * cI * (P3[3] * TMP7)) - cI * (P3[3] * TMP1
      * V2[3]))) + (TMP23 * (+cI * (P2[3] * V1[3] + P2[1] * V1[5])) + TMP26 *
      (+cI * (P1[3] * V2[3] + P1[1] * V2[5]))));
  T3[13] = denom * (OM3 * (TMP23 * (TMP26 * (P3[2] * (OM3 * 1.333333333 * P3[3]
      * (+cI * (TMP2 + TMP1)) + (-cI * (P2[3] + P1[3]))) - P3[3] * (+cI *
      (P2[2] + P1[2]))) + (P3[2] * (-cI * (TMP2 * V1[5]) + 0.666666667 * cI *
      (P3[3] * TMP9)) - cI * (P3[3] * TMP2 * V1[4]))) + TMP26 * (P3[2] * (-cI *
      (TMP1 * V2[5]) + 0.666666667 * cI * (P3[3] * TMP7)) - cI * (P3[3] * TMP1
      * V2[4]))) + (TMP23 * (+cI * (P2[3] * V1[4] + P2[2] * V1[5])) + TMP26 *
      (+cI * (P1[3] * V2[4] + P1[2] * V2[5]))));
  T3[17] = denom * 2. * (OM3 * (TMP23 * (TMP26 * (P3[3] * (OM3 * 0.666666667 *
      P3[3] * (+cI * (TMP2 + TMP1)) + (-cI * (P2[3] + P1[3]))) + (-0.333333333
      * cI * (TMP2 + TMP1))) + P3[3] * (-cI * (TMP2 * V1[5]) + 0.333333333 * cI
      * (P3[3] * TMP9))) + P3[3] * TMP26 * (-cI * (TMP1 * V2[5]) + 0.333333333
      * cI * (P3[3] * TMP7))) + (TMP23 * (+cI * (P2[3] * V1[5]) + 0.333333333 *
      cI * (TMP9)) + TMP26 * (+cI * (P1[3] * V2[5]) + 0.333333333 * cI *
      (TMP7))));
}


void VVT4_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP2; 
  complex<double> TMP1; 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP23; 
  complex<double> TMP31; 
  double P3[4]; 
  double OM3; 
  complex<double> denom; 
  complex<double> TMP26; 
  complex<double> TMP32; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +V1[0] + V2[0]; 
  T3[1] = +V1[1] + V2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP26 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP23 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP32 = -1. * (P1[0] * (P2[1] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] *
      (P3[3] * V1[3] - P3[1] * V1[5]) + P2[3] * (P3[1] * V1[4] - P3[2] *
      V1[3]))) + (P1[1] * (P2[0] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] *
      (P3[0] * V1[5] - P3[3] * V1[2]) + P2[3] * (P3[2] * V1[2] - P3[0] *
      V1[4]))) + (P1[2] * (P2[0] * (P3[1] * V1[5] - P3[3] * V1[3]) + (P2[1] *
      (P3[3] * V1[2] - P3[0] * V1[5]) + P2[3] * (P3[0] * V1[3] - P3[1] *
      V1[2]))) + P1[3] * (P2[0] * (P3[2] * V1[3] - P3[1] * V1[4]) + (P2[1] *
      (P3[0] * V1[4] - P3[2] * V1[2]) + P2[2] * (P3[1] * V1[2] - P3[0] *
      V1[3]))))));
  TMP31 = -1. * (P1[0] * (P2[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P2[2] *
      (P3[1] * V2[5] - P3[3] * V2[3]) + P2[3] * (P3[2] * V2[3] - P3[1] *
      V2[4]))) + (P1[1] * (P2[0] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P2[2] *
      (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3] * (P3[0] * V2[4] - P3[2] *
      V2[2]))) + (P1[2] * (P2[0] * (P3[3] * V2[3] - P3[1] * V2[5]) + (P2[1] *
      (P3[0] * V2[5] - P3[3] * V2[2]) + P2[3] * (P3[1] * V2[2] - P3[0] *
      V2[3]))) + P1[3] * (P2[0] * (P3[1] * V2[4] - P3[2] * V2[3]) + (P2[1] *
      (P3[2] * V2[2] - P3[0] * V2[4]) + P2[2] * (P3[0] * V2[3] - P3[1] *
      V2[2]))))));
  TMP1 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP2 = (P3[0] * P2[0] - P3[1] * P2[1] - P3[2] * P2[2] - P3[3] * P2[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * - 2. * cI * (TMP23 * (OM3 * P3[0] * (TMP2 * (P1[1] * (P3[2] *
      V1[5] - P3[3] * V1[4]) + (P1[2] * (P3[3] * V1[3] - P3[1] * V1[5]) + P1[3]
      * (P3[1] * V1[4] - P3[2] * V1[3]))) - 0.333333333 * (P3[0] * TMP32)) +
      (P2[0] * (P1[1] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[1] *
      V1[5] - P3[3] * V1[3]) + P1[3] * (P3[2] * V1[3] - P3[1] * V1[4]))) +
      0.333333333 * (TMP32))) + TMP26 * (OM3 * P3[0] * (TMP1 * (P2[1] * (P3[2]
      * V2[5] - P3[3] * V2[4]) + (P2[2] * (P3[3] * V2[3] - P3[1] * V2[5]) +
      P2[3] * (P3[1] * V2[4] - P3[2] * V2[3]))) - 0.333333333 * (P3[0] *
      TMP31)) + (P1[0] * (P2[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P2[2] *
      (P3[1] * V2[5] - P3[3] * V2[3]) + P2[3] * (P3[2] * V2[3] - P3[1] *
      V2[4]))) + 0.333333333 * (TMP31))));
  T3[3] = denom * cI * (OM3 * (P3[0] * (TMP23 * (TMP2 * (P1[0] * (P3[3] * V1[4]
      - P3[2] * V1[5]) + (P1[2] * (P3[0] * V1[5] - P3[3] * V1[2]) + P1[3] *
      (P3[2] * V1[2] - P3[0] * V1[4]))) + 0.666666667 * (P3[1] * TMP32)) +
      TMP26 * (TMP1 * (P2[0] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P2[2] *
      (P3[0] * V2[5] - P3[3] * V2[2]) + P2[3] * (P3[2] * V2[2] - P3[0] *
      V2[4]))) + 0.666666667 * (P3[1] * TMP31))) + P3[1] * (TMP1 * TMP26 *
      (P2[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[1] * V2[5] -
      P3[3] * V2[3]) + P2[3] * (P3[2] * V2[3] - P3[1] * V2[4]))) + TMP2 * TMP23
      * (P1[1] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[1] * V1[5] -
      P3[3] * V1[3]) + P1[3] * (P3[2] * V1[3] - P3[1] * V1[4]))))) + (TMP23 *
      (P2[0] * (P1[0] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P1[2] * (P3[3] *
      V1[2] - P3[0] * V1[5]) + P1[3] * (P3[0] * V1[4] - P3[2] * V1[2]))) +
      P2[1] * (P1[1] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P1[2] * (P3[3] *
      V1[3] - P3[1] * V1[5]) + P1[3] * (P3[1] * V1[4] - P3[2] * V1[3])))) +
      TMP26 * (P1[0] * (P2[0] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P2[2] *
      (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3] * (P3[0] * V2[4] - P3[2] *
      V2[2]))) + P1[1] * (P2[1] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P2[2] *
      (P3[3] * V2[3] - P3[1] * V2[5]) + P2[3] * (P3[1] * V2[4] - P3[2] *
      V2[3]))))));
  T3[4] = denom * cI * (OM3 * (P3[0] * (TMP23 * (TMP2 * (P1[0] * (P3[1] * V1[5]
      - P3[3] * V1[3]) + (P1[1] * (P3[3] * V1[2] - P3[0] * V1[5]) + P1[3] *
      (P3[0] * V1[3] - P3[1] * V1[2]))) + 0.666666667 * (P3[2] * TMP32)) +
      TMP26 * (TMP1 * (P2[0] * (P3[1] * V2[5] - P3[3] * V2[3]) + (P2[1] *
      (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3] * (P3[0] * V2[3] - P3[1] *
      V2[2]))) + 0.666666667 * (P3[2] * TMP31))) + P3[2] * (TMP1 * TMP26 *
      (P2[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[1] * V2[5] -
      P3[3] * V2[3]) + P2[3] * (P3[2] * V2[3] - P3[1] * V2[4]))) + TMP2 * TMP23
      * (P1[1] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[1] * V1[5] -
      P3[3] * V1[3]) + P1[3] * (P3[2] * V1[3] - P3[1] * V1[4]))))) + (TMP23 *
      (P2[0] * (P1[0] * (P3[3] * V1[3] - P3[1] * V1[5]) + (P1[1] * (P3[0] *
      V1[5] - P3[3] * V1[2]) + P1[3] * (P3[1] * V1[2] - P3[0] * V1[3]))) +
      P2[2] * (P1[1] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P1[2] * (P3[3] *
      V1[3] - P3[1] * V1[5]) + P1[3] * (P3[1] * V1[4] - P3[2] * V1[3])))) +
      TMP26 * (P1[0] * (P2[0] * (P3[3] * V2[3] - P3[1] * V2[5]) + (P2[1] *
      (P3[0] * V2[5] - P3[3] * V2[2]) + P2[3] * (P3[1] * V2[2] - P3[0] *
      V2[3]))) + P1[2] * (P2[1] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P2[2] *
      (P3[3] * V2[3] - P3[1] * V2[5]) + P2[3] * (P3[1] * V2[4] - P3[2] *
      V2[3]))))));
  T3[5] = denom * cI * (OM3 * (P3[0] * (TMP23 * (TMP2 * (P1[0] * (P3[2] * V1[3]
      - P3[1] * V1[4]) + (P1[1] * (P3[0] * V1[4] - P3[2] * V1[2]) + P1[2] *
      (P3[1] * V1[2] - P3[0] * V1[3]))) + 0.666666667 * (P3[3] * TMP32)) +
      TMP26 * (TMP1 * (P2[0] * (P3[2] * V2[3] - P3[1] * V2[4]) + (P2[1] *
      (P3[0] * V2[4] - P3[2] * V2[2]) + P2[2] * (P3[1] * V2[2] - P3[0] *
      V2[3]))) + 0.666666667 * (P3[3] * TMP31))) + P3[3] * (TMP1 * TMP26 *
      (P2[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[1] * V2[5] -
      P3[3] * V2[3]) + P2[3] * (P3[2] * V2[3] - P3[1] * V2[4]))) + TMP2 * TMP23
      * (P1[1] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[1] * V1[5] -
      P3[3] * V1[3]) + P1[3] * (P3[2] * V1[3] - P3[1] * V1[4]))))) + (TMP23 *
      (P2[0] * (P1[0] * (P3[1] * V1[4] - P3[2] * V1[3]) + (P1[1] * (P3[2] *
      V1[2] - P3[0] * V1[4]) + P1[2] * (P3[0] * V1[3] - P3[1] * V1[2]))) +
      P2[3] * (P1[1] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P1[2] * (P3[3] *
      V1[3] - P3[1] * V1[5]) + P1[3] * (P3[1] * V1[4] - P3[2] * V1[3])))) +
      TMP26 * (P1[0] * (P2[0] * (P3[1] * V2[4] - P3[2] * V2[3]) + (P2[1] *
      (P3[2] * V2[2] - P3[0] * V2[4]) + P2[2] * (P3[0] * V2[3] - P3[1] *
      V2[2]))) + P1[3] * (P2[1] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P2[2] *
      (P3[3] * V2[3] - P3[1] * V2[5]) + P2[3] * (P3[1] * V2[4] - P3[2] *
      V2[3]))))));
  T3[6] = denom * cI * (OM3 * (P3[0] * (TMP23 * (TMP2 * (P1[0] * (P3[3] * V1[4]
      - P3[2] * V1[5]) + (P1[2] * (P3[0] * V1[5] - P3[3] * V1[2]) + P1[3] *
      (P3[2] * V1[2] - P3[0] * V1[4]))) + 0.666666667 * (P3[1] * TMP32)) +
      TMP26 * (TMP1 * (P2[0] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P2[2] *
      (P3[0] * V2[5] - P3[3] * V2[2]) + P2[3] * (P3[2] * V2[2] - P3[0] *
      V2[4]))) + 0.666666667 * (P3[1] * TMP31))) + P3[1] * (TMP1 * TMP26 *
      (P2[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[1] * V2[5] -
      P3[3] * V2[3]) + P2[3] * (P3[2] * V2[3] - P3[1] * V2[4]))) + TMP2 * TMP23
      * (P1[1] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[1] * V1[5] -
      P3[3] * V1[3]) + P1[3] * (P3[2] * V1[3] - P3[1] * V1[4]))))) + (TMP23 *
      (P2[0] * (P1[0] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P1[2] * (P3[3] *
      V1[2] - P3[0] * V1[5]) + P1[3] * (P3[0] * V1[4] - P3[2] * V1[2]))) +
      P2[1] * (P1[1] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P1[2] * (P3[3] *
      V1[3] - P3[1] * V1[5]) + P1[3] * (P3[1] * V1[4] - P3[2] * V1[3])))) +
      TMP26 * (P1[0] * (P2[0] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P2[2] *
      (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3] * (P3[0] * V2[4] - P3[2] *
      V2[2]))) + P1[1] * (P2[1] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P2[2] *
      (P3[3] * V2[3] - P3[1] * V2[5]) + P2[3] * (P3[1] * V2[4] - P3[2] *
      V2[3]))))));
  T3[7] = denom * 2. * cI * (TMP23 * (OM3 * P3[1] * (TMP2 * (P1[0] * (P3[3] *
      V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[0] * V1[5] - P3[3] * V1[2]) + P1[3]
      * (P3[2] * V1[2] - P3[0] * V1[4]))) + 0.333333333 * (P3[1] * TMP32)) +
      (P2[1] * (P1[0] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P1[2] * (P3[3] *
      V1[2] - P3[0] * V1[5]) + P1[3] * (P3[0] * V1[4] - P3[2] * V1[2]))) +
      0.333333333 * (TMP32))) + TMP26 * (OM3 * P3[1] * (TMP1 * (P2[0] * (P3[3]
      * V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[0] * V2[5] - P3[3] * V2[2]) +
      P2[3] * (P3[2] * V2[2] - P3[0] * V2[4]))) + 0.333333333 * (P3[1] *
      TMP31)) + (P1[1] * (P2[0] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P2[2] *
      (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3] * (P3[0] * V2[4] - P3[2] *
      V2[2]))) + 0.333333333 * (TMP31))));
  T3[8] = denom * cI * (OM3 * (P3[1] * (TMP23 * (TMP2 * (P1[0] * (P3[1] * V1[5]
      - P3[3] * V1[3]) + (P1[1] * (P3[3] * V1[2] - P3[0] * V1[5]) + P1[3] *
      (P3[0] * V1[3] - P3[1] * V1[2]))) + 0.666666667 * (P3[2] * TMP32)) +
      TMP26 * (TMP1 * (P2[0] * (P3[1] * V2[5] - P3[3] * V2[3]) + (P2[1] *
      (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3] * (P3[0] * V2[3] - P3[1] *
      V2[2]))) + 0.666666667 * (P3[2] * TMP31))) + P3[2] * (TMP1 * TMP26 *
      (P2[0] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[0] * V2[5] -
      P3[3] * V2[2]) + P2[3] * (P3[2] * V2[2] - P3[0] * V2[4]))) + TMP2 * TMP23
      * (P1[0] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[0] * V1[5] -
      P3[3] * V1[2]) + P1[3] * (P3[2] * V1[2] - P3[0] * V1[4]))))) + (TMP23 *
      (P2[1] * (P1[0] * (P3[3] * V1[3] - P3[1] * V1[5]) + (P1[1] * (P3[0] *
      V1[5] - P3[3] * V1[2]) + P1[3] * (P3[1] * V1[2] - P3[0] * V1[3]))) +
      P2[2] * (P1[0] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P1[2] * (P3[3] *
      V1[2] - P3[0] * V1[5]) + P1[3] * (P3[0] * V1[4] - P3[2] * V1[2])))) +
      TMP26 * (P1[1] * (P2[0] * (P3[3] * V2[3] - P3[1] * V2[5]) + (P2[1] *
      (P3[0] * V2[5] - P3[3] * V2[2]) + P2[3] * (P3[1] * V2[2] - P3[0] *
      V2[3]))) + P1[2] * (P2[0] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P2[2] *
      (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3] * (P3[0] * V2[4] - P3[2] *
      V2[2]))))));
  T3[9] = denom * cI * (OM3 * (P3[1] * (TMP23 * (TMP2 * (P1[0] * (P3[2] * V1[3]
      - P3[1] * V1[4]) + (P1[1] * (P3[0] * V1[4] - P3[2] * V1[2]) + P1[2] *
      (P3[1] * V1[2] - P3[0] * V1[3]))) + 0.666666667 * (P3[3] * TMP32)) +
      TMP26 * (TMP1 * (P2[0] * (P3[2] * V2[3] - P3[1] * V2[4]) + (P2[1] *
      (P3[0] * V2[4] - P3[2] * V2[2]) + P2[2] * (P3[1] * V2[2] - P3[0] *
      V2[3]))) + 0.666666667 * (P3[3] * TMP31))) + P3[3] * (TMP1 * TMP26 *
      (P2[0] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[0] * V2[5] -
      P3[3] * V2[2]) + P2[3] * (P3[2] * V2[2] - P3[0] * V2[4]))) + TMP2 * TMP23
      * (P1[0] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[0] * V1[5] -
      P3[3] * V1[2]) + P1[3] * (P3[2] * V1[2] - P3[0] * V1[4]))))) + (TMP23 *
      (P2[1] * (P1[0] * (P3[1] * V1[4] - P3[2] * V1[3]) + (P1[1] * (P3[2] *
      V1[2] - P3[0] * V1[4]) + P1[2] * (P3[0] * V1[3] - P3[1] * V1[2]))) +
      P2[3] * (P1[0] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P1[2] * (P3[3] *
      V1[2] - P3[0] * V1[5]) + P1[3] * (P3[0] * V1[4] - P3[2] * V1[2])))) +
      TMP26 * (P1[1] * (P2[0] * (P3[1] * V2[4] - P3[2] * V2[3]) + (P2[1] *
      (P3[2] * V2[2] - P3[0] * V2[4]) + P2[2] * (P3[0] * V2[3] - P3[1] *
      V2[2]))) + P1[3] * (P2[0] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P2[2] *
      (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3] * (P3[0] * V2[4] - P3[2] *
      V2[2]))))));
  T3[10] = denom * cI * (OM3 * (P3[0] * (TMP23 * (TMP2 * (P1[0] * (P3[1] *
      V1[5] - P3[3] * V1[3]) + (P1[1] * (P3[3] * V1[2] - P3[0] * V1[5]) + P1[3]
      * (P3[0] * V1[3] - P3[1] * V1[2]))) + 0.666666667 * (P3[2] * TMP32)) +
      TMP26 * (TMP1 * (P2[0] * (P3[1] * V2[5] - P3[3] * V2[3]) + (P2[1] *
      (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3] * (P3[0] * V2[3] - P3[1] *
      V2[2]))) + 0.666666667 * (P3[2] * TMP31))) + P3[2] * (TMP1 * TMP26 *
      (P2[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[1] * V2[5] -
      P3[3] * V2[3]) + P2[3] * (P3[2] * V2[3] - P3[1] * V2[4]))) + TMP2 * TMP23
      * (P1[1] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[1] * V1[5] -
      P3[3] * V1[3]) + P1[3] * (P3[2] * V1[3] - P3[1] * V1[4]))))) + (TMP23 *
      (P2[0] * (P1[0] * (P3[3] * V1[3] - P3[1] * V1[5]) + (P1[1] * (P3[0] *
      V1[5] - P3[3] * V1[2]) + P1[3] * (P3[1] * V1[2] - P3[0] * V1[3]))) +
      P2[2] * (P1[1] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P1[2] * (P3[3] *
      V1[3] - P3[1] * V1[5]) + P1[3] * (P3[1] * V1[4] - P3[2] * V1[3])))) +
      TMP26 * (P1[0] * (P2[0] * (P3[3] * V2[3] - P3[1] * V2[5]) + (P2[1] *
      (P3[0] * V2[5] - P3[3] * V2[2]) + P2[3] * (P3[1] * V2[2] - P3[0] *
      V2[3]))) + P1[2] * (P2[1] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P2[2] *
      (P3[3] * V2[3] - P3[1] * V2[5]) + P2[3] * (P3[1] * V2[4] - P3[2] *
      V2[3]))))));
  T3[11] = denom * cI * (OM3 * (P3[1] * (TMP23 * (TMP2 * (P1[0] * (P3[1] *
      V1[5] - P3[3] * V1[3]) + (P1[1] * (P3[3] * V1[2] - P3[0] * V1[5]) + P1[3]
      * (P3[0] * V1[3] - P3[1] * V1[2]))) + 0.666666667 * (P3[2] * TMP32)) +
      TMP26 * (TMP1 * (P2[0] * (P3[1] * V2[5] - P3[3] * V2[3]) + (P2[1] *
      (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3] * (P3[0] * V2[3] - P3[1] *
      V2[2]))) + 0.666666667 * (P3[2] * TMP31))) + P3[2] * (TMP1 * TMP26 *
      (P2[0] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[0] * V2[5] -
      P3[3] * V2[2]) + P2[3] * (P3[2] * V2[2] - P3[0] * V2[4]))) + TMP2 * TMP23
      * (P1[0] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[0] * V1[5] -
      P3[3] * V1[2]) + P1[3] * (P3[2] * V1[2] - P3[0] * V1[4]))))) + (TMP23 *
      (P2[1] * (P1[0] * (P3[3] * V1[3] - P3[1] * V1[5]) + (P1[1] * (P3[0] *
      V1[5] - P3[3] * V1[2]) + P1[3] * (P3[1] * V1[2] - P3[0] * V1[3]))) +
      P2[2] * (P1[0] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P1[2] * (P3[3] *
      V1[2] - P3[0] * V1[5]) + P1[3] * (P3[0] * V1[4] - P3[2] * V1[2])))) +
      TMP26 * (P1[1] * (P2[0] * (P3[3] * V2[3] - P3[1] * V2[5]) + (P2[1] *
      (P3[0] * V2[5] - P3[3] * V2[2]) + P2[3] * (P3[1] * V2[2] - P3[0] *
      V2[3]))) + P1[2] * (P2[0] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P2[2] *
      (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3] * (P3[0] * V2[4] - P3[2] *
      V2[2]))))));
  T3[12] = denom * 2. * cI * (TMP23 * (OM3 * P3[2] * (TMP2 * (P1[0] * (P3[1] *
      V1[5] - P3[3] * V1[3]) + (P1[1] * (P3[3] * V1[2] - P3[0] * V1[5]) + P1[3]
      * (P3[0] * V1[3] - P3[1] * V1[2]))) + 0.333333333 * (P3[2] * TMP32)) +
      (P2[2] * (P1[0] * (P3[3] * V1[3] - P3[1] * V1[5]) + (P1[1] * (P3[0] *
      V1[5] - P3[3] * V1[2]) + P1[3] * (P3[1] * V1[2] - P3[0] * V1[3]))) +
      0.333333333 * (TMP32))) + TMP26 * (OM3 * P3[2] * (TMP1 * (P2[0] * (P3[1]
      * V2[5] - P3[3] * V2[3]) + (P2[1] * (P3[3] * V2[2] - P3[0] * V2[5]) +
      P2[3] * (P3[0] * V2[3] - P3[1] * V2[2]))) + 0.333333333 * (P3[2] *
      TMP31)) + (P1[2] * (P2[0] * (P3[3] * V2[3] - P3[1] * V2[5]) + (P2[1] *
      (P3[0] * V2[5] - P3[3] * V2[2]) + P2[3] * (P3[1] * V2[2] - P3[0] *
      V2[3]))) + 0.333333333 * (TMP31))));
  T3[13] = denom * cI * (OM3 * (P3[2] * (TMP23 * (TMP2 * (P1[0] * (P3[2] *
      V1[3] - P3[1] * V1[4]) + (P1[1] * (P3[0] * V1[4] - P3[2] * V1[2]) + P1[2]
      * (P3[1] * V1[2] - P3[0] * V1[3]))) + 0.666666667 * (P3[3] * TMP32)) +
      TMP26 * (TMP1 * (P2[0] * (P3[2] * V2[3] - P3[1] * V2[4]) + (P2[1] *
      (P3[0] * V2[4] - P3[2] * V2[2]) + P2[2] * (P3[1] * V2[2] - P3[0] *
      V2[3]))) + 0.666666667 * (P3[3] * TMP31))) + P3[3] * (TMP1 * TMP26 *
      (P2[0] * (P3[1] * V2[5] - P3[3] * V2[3]) + (P2[1] * (P3[3] * V2[2] -
      P3[0] * V2[5]) + P2[3] * (P3[0] * V2[3] - P3[1] * V2[2]))) + TMP2 * TMP23
      * (P1[0] * (P3[1] * V1[5] - P3[3] * V1[3]) + (P1[1] * (P3[3] * V1[2] -
      P3[0] * V1[5]) + P1[3] * (P3[0] * V1[3] - P3[1] * V1[2]))))) + (TMP23 *
      (P2[2] * (P1[0] * (P3[1] * V1[4] - P3[2] * V1[3]) + (P1[1] * (P3[2] *
      V1[2] - P3[0] * V1[4]) + P1[2] * (P3[0] * V1[3] - P3[1] * V1[2]))) +
      P2[3] * (P1[0] * (P3[3] * V1[3] - P3[1] * V1[5]) + (P1[1] * (P3[0] *
      V1[5] - P3[3] * V1[2]) + P1[3] * (P3[1] * V1[2] - P3[0] * V1[3])))) +
      TMP26 * (P1[2] * (P2[0] * (P3[1] * V2[4] - P3[2] * V2[3]) + (P2[1] *
      (P3[2] * V2[2] - P3[0] * V2[4]) + P2[2] * (P3[0] * V2[3] - P3[1] *
      V2[2]))) + P1[3] * (P2[0] * (P3[3] * V2[3] - P3[1] * V2[5]) + (P2[1] *
      (P3[0] * V2[5] - P3[3] * V2[2]) + P2[3] * (P3[1] * V2[2] - P3[0] *
      V2[3]))))));
  T3[14] = denom * cI * (OM3 * (P3[0] * (TMP23 * (TMP2 * (P1[0] * (P3[2] *
      V1[3] - P3[1] * V1[4]) + (P1[1] * (P3[0] * V1[4] - P3[2] * V1[2]) + P1[2]
      * (P3[1] * V1[2] - P3[0] * V1[3]))) + 0.666666667 * (P3[3] * TMP32)) +
      TMP26 * (TMP1 * (P2[0] * (P3[2] * V2[3] - P3[1] * V2[4]) + (P2[1] *
      (P3[0] * V2[4] - P3[2] * V2[2]) + P2[2] * (P3[1] * V2[2] - P3[0] *
      V2[3]))) + 0.666666667 * (P3[3] * TMP31))) + P3[3] * (TMP1 * TMP26 *
      (P2[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[1] * V2[5] -
      P3[3] * V2[3]) + P2[3] * (P3[2] * V2[3] - P3[1] * V2[4]))) + TMP2 * TMP23
      * (P1[1] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[1] * V1[5] -
      P3[3] * V1[3]) + P1[3] * (P3[2] * V1[3] - P3[1] * V1[4]))))) + (TMP23 *
      (P2[0] * (P1[0] * (P3[1] * V1[4] - P3[2] * V1[3]) + (P1[1] * (P3[2] *
      V1[2] - P3[0] * V1[4]) + P1[2] * (P3[0] * V1[3] - P3[1] * V1[2]))) +
      P2[3] * (P1[1] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P1[2] * (P3[3] *
      V1[3] - P3[1] * V1[5]) + P1[3] * (P3[1] * V1[4] - P3[2] * V1[3])))) +
      TMP26 * (P1[0] * (P2[0] * (P3[1] * V2[4] - P3[2] * V2[3]) + (P2[1] *
      (P3[2] * V2[2] - P3[0] * V2[4]) + P2[2] * (P3[0] * V2[3] - P3[1] *
      V2[2]))) + P1[3] * (P2[1] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P2[2] *
      (P3[3] * V2[3] - P3[1] * V2[5]) + P2[3] * (P3[1] * V2[4] - P3[2] *
      V2[3]))))));
  T3[15] = denom * cI * (OM3 * (P3[1] * (TMP23 * (TMP2 * (P1[0] * (P3[2] *
      V1[3] - P3[1] * V1[4]) + (P1[1] * (P3[0] * V1[4] - P3[2] * V1[2]) + P1[2]
      * (P3[1] * V1[2] - P3[0] * V1[3]))) + 0.666666667 * (P3[3] * TMP32)) +
      TMP26 * (TMP1 * (P2[0] * (P3[2] * V2[3] - P3[1] * V2[4]) + (P2[1] *
      (P3[0] * V2[4] - P3[2] * V2[2]) + P2[2] * (P3[1] * V2[2] - P3[0] *
      V2[3]))) + 0.666666667 * (P3[3] * TMP31))) + P3[3] * (TMP1 * TMP26 *
      (P2[0] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[0] * V2[5] -
      P3[3] * V2[2]) + P2[3] * (P3[2] * V2[2] - P3[0] * V2[4]))) + TMP2 * TMP23
      * (P1[0] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[0] * V1[5] -
      P3[3] * V1[2]) + P1[3] * (P3[2] * V1[2] - P3[0] * V1[4]))))) + (TMP23 *
      (P2[1] * (P1[0] * (P3[1] * V1[4] - P3[2] * V1[3]) + (P1[1] * (P3[2] *
      V1[2] - P3[0] * V1[4]) + P1[2] * (P3[0] * V1[3] - P3[1] * V1[2]))) +
      P2[3] * (P1[0] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P1[2] * (P3[3] *
      V1[2] - P3[0] * V1[5]) + P1[3] * (P3[0] * V1[4] - P3[2] * V1[2])))) +
      TMP26 * (P1[1] * (P2[0] * (P3[1] * V2[4] - P3[2] * V2[3]) + (P2[1] *
      (P3[2] * V2[2] - P3[0] * V2[4]) + P2[2] * (P3[0] * V2[3] - P3[1] *
      V2[2]))) + P1[3] * (P2[0] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P2[2] *
      (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3] * (P3[0] * V2[4] - P3[2] *
      V2[2]))))));
  T3[16] = denom * cI * (OM3 * (P3[2] * (TMP23 * (TMP2 * (P1[0] * (P3[2] *
      V1[3] - P3[1] * V1[4]) + (P1[1] * (P3[0] * V1[4] - P3[2] * V1[2]) + P1[2]
      * (P3[1] * V1[2] - P3[0] * V1[3]))) + 0.666666667 * (P3[3] * TMP32)) +
      TMP26 * (TMP1 * (P2[0] * (P3[2] * V2[3] - P3[1] * V2[4]) + (P2[1] *
      (P3[0] * V2[4] - P3[2] * V2[2]) + P2[2] * (P3[1] * V2[2] - P3[0] *
      V2[3]))) + 0.666666667 * (P3[3] * TMP31))) + P3[3] * (TMP1 * TMP26 *
      (P2[0] * (P3[1] * V2[5] - P3[3] * V2[3]) + (P2[1] * (P3[3] * V2[2] -
      P3[0] * V2[5]) + P2[3] * (P3[0] * V2[3] - P3[1] * V2[2]))) + TMP2 * TMP23
      * (P1[0] * (P3[1] * V1[5] - P3[3] * V1[3]) + (P1[1] * (P3[3] * V1[2] -
      P3[0] * V1[5]) + P1[3] * (P3[0] * V1[3] - P3[1] * V1[2]))))) + (TMP23 *
      (P2[2] * (P1[0] * (P3[1] * V1[4] - P3[2] * V1[3]) + (P1[1] * (P3[2] *
      V1[2] - P3[0] * V1[4]) + P1[2] * (P3[0] * V1[3] - P3[1] * V1[2]))) +
      P2[3] * (P1[0] * (P3[3] * V1[3] - P3[1] * V1[5]) + (P1[1] * (P3[0] *
      V1[5] - P3[3] * V1[2]) + P1[3] * (P3[1] * V1[2] - P3[0] * V1[3])))) +
      TMP26 * (P1[2] * (P2[0] * (P3[1] * V2[4] - P3[2] * V2[3]) + (P2[1] *
      (P3[2] * V2[2] - P3[0] * V2[4]) + P2[2] * (P3[0] * V2[3] - P3[1] *
      V2[2]))) + P1[3] * (P2[0] * (P3[3] * V2[3] - P3[1] * V2[5]) + (P2[1] *
      (P3[0] * V2[5] - P3[3] * V2[2]) + P2[3] * (P3[1] * V2[2] - P3[0] *
      V2[3]))))));
  T3[17] = denom * 2. * cI * (TMP23 * (OM3 * P3[3] * (TMP2 * (P1[0] * (P3[2] *
      V1[3] - P3[1] * V1[4]) + (P1[1] * (P3[0] * V1[4] - P3[2] * V1[2]) + P1[2]
      * (P3[1] * V1[2] - P3[0] * V1[3]))) + 0.333333333 * (P3[3] * TMP32)) +
      (P2[3] * (P1[0] * (P3[1] * V1[4] - P3[2] * V1[3]) + (P1[1] * (P3[2] *
      V1[2] - P3[0] * V1[4]) + P1[2] * (P3[0] * V1[3] - P3[1] * V1[2]))) +
      0.333333333 * (TMP32))) + TMP26 * (OM3 * P3[3] * (TMP1 * (P2[0] * (P3[2]
      * V2[3] - P3[1] * V2[4]) + (P2[1] * (P3[0] * V2[4] - P3[2] * V2[2]) +
      P2[2] * (P3[1] * V2[2] - P3[0] * V2[3]))) + 0.333333333 * (P3[3] *
      TMP31)) + (P1[3] * (P2[0] * (P3[1] * V2[4] - P3[2] * V2[3]) + (P2[1] *
      (P3[2] * V2[2] - P3[0] * V2[4]) + P2[2] * (P3[0] * V2[3] - P3[1] *
      V2[2]))) + 0.333333333 * (TMP31))));
}


void FFT1_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP2; 
  complex<double> TMP1; 
  double P1[4]; 
  complex<double> TMP10; 
  double P2[4]; 
  double P3[4]; 
  double OM3; 
  complex<double> denom; 
  complex<double> TMP0; 
  P1[0] = F1[0].real(); 
  P1[1] = F1[1].real(); 
  P1[2] = F1[1].imag(); 
  P1[3] = F1[0].imag(); 
  P2[0] = F2[0].real(); 
  P2[1] = F2[1].real(); 
  P2[2] = F2[1].imag(); 
  P2[3] = F2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +F1[0] + F2[0]; 
  T3[1] = +F1[1] + F2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP10 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP1 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP0 = (F2[4] * F1[4] + F2[5] * F1[5] - F2[2] * F1[2] - F2[3] * F1[3]); 
  TMP2 = (P3[0] * P2[0] - P3[1] * P2[1] - P3[2] * P2[2] - P3[3] * P2[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * 2. * TMP0 * (OM3 * (P3[0] * (P3[0] * 0.333333333 * (+cI *
      (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (TMP1 * P2[0] + P1[0] *
      TMP2))) + 0.333333333 * cI * (TMP1 * TMP2)) + (-0.333333333 * cI *
      (TMP10) + cI * (P1[0] * P2[0])));
  T3[3] = denom * TMP0 * (OM3 * (P3[0] * (P3[1] * 0.666666667 * (+cI * (TMP10)
      + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (P1[1] * TMP2 + TMP1 * P2[1])))
      - P3[1] * (+cI * (TMP1 * P2[0] + P1[0] * TMP2))) + (+cI * (P1[1] * P2[0]
      + P1[0] * P2[1])));
  T3[4] = denom * TMP0 * (OM3 * (P3[0] * (P3[2] * 0.666666667 * (+cI * (TMP10)
      + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (P1[2] * TMP2 + TMP1 * P2[2])))
      - P3[2] * (+cI * (TMP1 * P2[0] + P1[0] * TMP2))) + (+cI * (P1[2] * P2[0]
      + P1[0] * P2[2])));
  T3[5] = denom * TMP0 * (OM3 * (P3[0] * (P3[3] * 0.666666667 * (+cI * (TMP10)
      + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (P1[3] * TMP2 + TMP1 * P2[3])))
      - P3[3] * (+cI * (TMP1 * P2[0] + P1[0] * TMP2))) + (+cI * (P1[3] * P2[0]
      + P1[0] * P2[3])));
  T3[6] = denom * TMP0 * (OM3 * (P3[0] * (P3[1] * 0.666666667 * (+cI * (TMP10)
      + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (TMP1 * P2[1] + P1[1] * TMP2)))
      - P3[1] * (+cI * (P1[0] * TMP2 + TMP1 * P2[0]))) + (+cI * (P1[0] * P2[1]
      + P1[1] * P2[0])));
  T3[7] = denom * 2. * TMP0 * (OM3 * (P3[1] * (P3[1] * 0.333333333 * (+cI *
      (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (TMP1 * P2[1] + P1[1] *
      TMP2))) - 0.333333333 * cI * (TMP1 * TMP2)) + (+cI * (P1[1] * P2[1]) +
      0.333333333 * cI * (TMP10)));
  T3[8] = denom * TMP0 * (OM3 * (P3[1] * (P3[2] * 0.666666667 * (+cI * (TMP10)
      + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (P1[2] * TMP2 + TMP1 * P2[2])))
      - P3[2] * (+cI * (TMP1 * P2[1] + P1[1] * TMP2))) + (+cI * (P1[2] * P2[1]
      + P1[1] * P2[2])));
  T3[9] = denom * TMP0 * (OM3 * (P3[1] * (P3[3] * 0.666666667 * (+cI * (TMP10)
      + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (P1[3] * TMP2 + TMP1 * P2[3])))
      - P3[3] * (+cI * (TMP1 * P2[1] + P1[1] * TMP2))) + (+cI * (P1[3] * P2[1]
      + P1[1] * P2[3])));
  T3[10] = denom * TMP0 * (OM3 * (P3[0] * (P3[2] * 0.666666667 * (+cI * (TMP10)
      + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (TMP1 * P2[2] + P1[2] * TMP2)))
      - P3[2] * (+cI * (P1[0] * TMP2 + TMP1 * P2[0]))) + (+cI * (P1[0] * P2[2]
      + P1[2] * P2[0])));
  T3[11] = denom * TMP0 * (OM3 * (P3[1] * (P3[2] * 0.666666667 * (+cI * (TMP10)
      + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (TMP1 * P2[2] + P1[2] * TMP2)))
      - P3[2] * (+cI * (P1[1] * TMP2 + TMP1 * P2[1]))) + (+cI * (P1[1] * P2[2]
      + P1[2] * P2[1])));
  T3[12] = denom * 2. * TMP0 * (OM3 * (P3[2] * (P3[2] * 0.333333333 * (+cI *
      (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (TMP1 * P2[2] + P1[2] *
      TMP2))) - 0.333333333 * cI * (TMP1 * TMP2)) + (+cI * (P1[2] * P2[2]) +
      0.333333333 * cI * (TMP10)));
  T3[13] = denom * TMP0 * (OM3 * (P3[2] * (P3[3] * 0.666666667 * (+cI * (TMP10)
      + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (P1[3] * TMP2 + TMP1 * P2[3])))
      - P3[3] * (+cI * (TMP1 * P2[2] + P1[2] * TMP2))) + (+cI * (P1[3] * P2[2]
      + P1[2] * P2[3])));
  T3[14] = denom * TMP0 * (OM3 * (P3[0] * (P3[3] * 0.666666667 * (+cI * (TMP10)
      + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (TMP1 * P2[3] + P1[3] * TMP2)))
      - P3[3] * (+cI * (P1[0] * TMP2 + TMP1 * P2[0]))) + (+cI * (P1[0] * P2[3]
      + P1[3] * P2[0])));
  T3[15] = denom * TMP0 * (OM3 * (P3[1] * (P3[3] * 0.666666667 * (+cI * (TMP10)
      + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (TMP1 * P2[3] + P1[3] * TMP2)))
      - P3[3] * (+cI * (P1[1] * TMP2 + TMP1 * P2[1]))) + (+cI * (P1[1] * P2[3]
      + P1[3] * P2[1])));
  T3[16] = denom * TMP0 * (OM3 * (P3[2] * (P3[3] * 0.666666667 * (+cI * (TMP10)
      + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (TMP1 * P2[3] + P1[3] * TMP2)))
      - P3[3] * (+cI * (P1[2] * TMP2 + TMP1 * P2[2]))) + (+cI * (P1[2] * P2[3]
      + P1[3] * P2[2])));
  T3[17] = denom * 2. * TMP0 * (OM3 * (P3[3] * (P3[3] * 0.333333333 * (+cI *
      (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (TMP1 * P2[3] + P1[3] *
      TMP2))) - 0.333333333 * cI * (TMP1 * TMP2)) + (+cI * (P1[3] * P2[3]) +
      0.333333333 * cI * (TMP10)));
}

void FFT1_2_3_4_3(complex<double> F1[], complex<double> F2[], complex<double>
    COUP1, complex<double> COUP2, complex<double> COUP3, complex<double> COUP4,
    double M3, double W3, complex<double> T3[])
{
//   complex<double> cI = complex<double> (0., 1.); 
  complex<double> Ttmp[18]; 
//   double P1[4]; 
//   double P2[4]; 
//   double P3[4]; 
//   double OM3; 
  complex<double> denom; 
  int i; 
  FFT1_3(F1, F2, COUP1, M3, W3, T3); 
  FFT2_3(F1, F2, COUP2, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  FFT3_3(F1, F2, COUP3, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  FFT4_3(F1, F2, COUP4, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
}

void VVT11_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP23; 
  double P3[4]; 
  double OM3; 
  complex<double> denom; 
  complex<double> TMP26; 
  complex<double> TMP4; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +V1[0] + V2[0]; 
  T3[1] = +V1[1] + V2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP4 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP26 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP23 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * 2. * (OM3 * (P3[0] * (P3[0] * 0.333333333 * (+cI * (TMP4) +
      2. * cI * (OM3 * TMP23 * TMP26)) + (-cI * (V2[2] * TMP26 + V1[2] *
      TMP23))) + 0.333333333 * cI * (TMP23 * TMP26)) + (-0.333333333 * cI *
      (TMP4) + cI * (V2[2] * V1[2])));
  T3[6] = denom * (OM3 * (P3[0] * (P3[1] * 0.666666667 * (+cI * (TMP4) + 2. *
      cI * (OM3 * TMP23 * TMP26)) + (-cI * (V2[3] * TMP26 + V1[3] * TMP23))) -
      P3[1] * (+cI * (V1[2] * TMP23 + V2[2] * TMP26))) + (+cI * (V2[3] * V1[2]
      + V2[2] * V1[3])));
  T3[10] = denom * (OM3 * (P3[0] * (P3[2] * 0.666666667 * (+cI * (TMP4) + 2. *
      cI * (OM3 * TMP23 * TMP26)) + (-cI * (V2[4] * TMP26 + V1[4] * TMP23))) -
      P3[2] * (+cI * (V1[2] * TMP23 + V2[2] * TMP26))) + (+cI * (V2[4] * V1[2]
      + V2[2] * V1[4])));
  T3[14] = denom * (OM3 * (P3[0] * (P3[3] * 0.666666667 * (+cI * (TMP4) + 2. *
      cI * (OM3 * TMP23 * TMP26)) + (-cI * (V2[5] * TMP26 + V1[5] * TMP23))) -
      P3[3] * (+cI * (V1[2] * TMP23 + V2[2] * TMP26))) + (+cI * (V2[5] * V1[2]
      + V2[2] * V1[5])));
  T3[3] = denom * (OM3 * (P3[0] * (P3[1] * 0.666666667 * (+cI * (TMP4) + 2. *
      cI * (OM3 * TMP23 * TMP26)) + (-cI * (V1[3] * TMP23 + V2[3] * TMP26))) -
      P3[1] * (+cI * (V2[2] * TMP26 + V1[2] * TMP23))) + (+cI * (V2[2] * V1[3]
      + V2[3] * V1[2])));
  T3[7] = denom * 2. * (OM3 * (P3[1] * (P3[1] * 0.333333333 * (+cI * (TMP4) +
      2. * cI * (OM3 * TMP23 * TMP26)) + (-cI * (V2[3] * TMP26 + V1[3] *
      TMP23))) - 0.333333333 * cI * (TMP23 * TMP26)) + (+cI * (V2[3] * V1[3]) +
      0.333333333 * cI * (TMP4)));
  T3[11] = denom * (OM3 * (P3[1] * (P3[2] * 0.666666667 * (+cI * (TMP4) + 2. *
      cI * (OM3 * TMP23 * TMP26)) + (-cI * (V2[4] * TMP26 + V1[4] * TMP23))) -
      P3[2] * (+cI * (V1[3] * TMP23 + V2[3] * TMP26))) + (+cI * (V2[4] * V1[3]
      + V2[3] * V1[4])));
  T3[15] = denom * (OM3 * (P3[1] * (P3[3] * 0.666666667 * (+cI * (TMP4) + 2. *
      cI * (OM3 * TMP23 * TMP26)) + (-cI * (V2[5] * TMP26 + V1[5] * TMP23))) -
      P3[3] * (+cI * (V1[3] * TMP23 + V2[3] * TMP26))) + (+cI * (V2[5] * V1[3]
      + V2[3] * V1[5])));
  T3[4] = denom * (OM3 * (P3[0] * (P3[2] * 0.666666667 * (+cI * (TMP4) + 2. *
      cI * (OM3 * TMP23 * TMP26)) + (-cI * (V1[4] * TMP23 + V2[4] * TMP26))) -
      P3[2] * (+cI * (V2[2] * TMP26 + V1[2] * TMP23))) + (+cI * (V2[2] * V1[4]
      + V2[4] * V1[2])));
  T3[8] = denom * (OM3 * (P3[1] * (P3[2] * 0.666666667 * (+cI * (TMP4) + 2. *
      cI * (OM3 * TMP23 * TMP26)) + (-cI * (V1[4] * TMP23 + V2[4] * TMP26))) -
      P3[2] * (+cI * (V2[3] * TMP26 + V1[3] * TMP23))) + (+cI * (V2[3] * V1[4]
      + V2[4] * V1[3])));
  T3[12] = denom * 2. * (OM3 * (P3[2] * (P3[2] * 0.333333333 * (+cI * (TMP4) +
      2. * cI * (OM3 * TMP23 * TMP26)) + (-cI * (V2[4] * TMP26 + V1[4] *
      TMP23))) - 0.333333333 * cI * (TMP23 * TMP26)) + (+cI * (V2[4] * V1[4]) +
      0.333333333 * cI * (TMP4)));
  T3[16] = denom * (OM3 * (P3[2] * (P3[3] * 0.666666667 * (+cI * (TMP4) + 2. *
      cI * (OM3 * TMP23 * TMP26)) + (-cI * (V2[5] * TMP26 + V1[5] * TMP23))) -
      P3[3] * (+cI * (V1[4] * TMP23 + V2[4] * TMP26))) + (+cI * (V2[5] * V1[4]
      + V2[4] * V1[5])));
  T3[5] = denom * (OM3 * (P3[0] * (P3[3] * 0.666666667 * (+cI * (TMP4) + 2. *
      cI * (OM3 * TMP23 * TMP26)) + (-cI * (V1[5] * TMP23 + V2[5] * TMP26))) -
      P3[3] * (+cI * (V2[2] * TMP26 + V1[2] * TMP23))) + (+cI * (V2[2] * V1[5]
      + V2[5] * V1[2])));
  T3[9] = denom * (OM3 * (P3[1] * (P3[3] * 0.666666667 * (+cI * (TMP4) + 2. *
      cI * (OM3 * TMP23 * TMP26)) + (-cI * (V1[5] * TMP23 + V2[5] * TMP26))) -
      P3[3] * (+cI * (V2[3] * TMP26 + V1[3] * TMP23))) + (+cI * (V2[3] * V1[5]
      + V2[5] * V1[3])));
  T3[13] = denom * (OM3 * (P3[2] * (P3[3] * 0.666666667 * (+cI * (TMP4) + 2. *
      cI * (OM3 * TMP23 * TMP26)) + (-cI * (V1[5] * TMP23 + V2[5] * TMP26))) -
      P3[3] * (+cI * (V2[4] * TMP26 + V1[4] * TMP23))) + (+cI * (V2[4] * V1[5]
      + V2[5] * V1[4])));
  T3[17] = denom * 2. * (OM3 * (P3[3] * (P3[3] * 0.333333333 * (+cI * (TMP4) +
      2. * cI * (OM3 * TMP23 * TMP26)) + (-cI * (V2[5] * TMP26 + V1[5] *
      TMP23))) - 0.333333333 * cI * (TMP23 * TMP26)) + (+cI * (V2[5] * V1[5]) +
      0.333333333 * cI * (TMP4)));
}


void VVT13_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP2; 
  complex<double> TMP1; 
  double P1[4]; 
  complex<double> TMP10; 
  double P2[4]; 
  complex<double> TMP23; 
  complex<double> TMP7; 
  double P3[4]; 
  double OM3; 
  complex<double> denom; 
  complex<double> TMP26; 
  complex<double> TMP4; 
  complex<double> TMP9; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +V1[0] + V2[0]; 
  T3[1] = +V1[1] + V2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP26 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP23 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP9 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP4 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP7 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP1 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP10 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP2 = (P3[0] * P2[0] - P3[1] * P2[1] - P3[2] * P2[2] - P3[3] * P2[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * 2. * (TMP1 * (TMP2 * (TMP4 * 0.333333333 * (+cI * (P3[0] *
      P3[0] * OM3) + - 1.000000000 * cI) + cI * (V2[2] * V1[2])) + TMP23 *
      (TMP9 * - 0.333333333 * (+cI * (P3[0] * P3[0] * OM3) + - 1.000000000 *
      cI) - cI * (P2[0] * V1[2]))) + TMP26 * (TMP2 * (TMP7 * - 0.333333333 *
      (+cI * (P3[0] * P3[0] * OM3) + - 1.000000000 * cI) - cI * (P1[0] *
      V2[2])) + TMP23 * (TMP10 * 0.333333333 * (+cI * (P3[0] * P3[0] * OM3) + -
      1.000000000 * cI) + cI * (P1[0] * P2[0]))));
  T3[3] = denom * (TMP1 * (TMP2 * (+cI * (V2[2] * V1[3] + V2[3] * V1[2]) +
      0.666666667 * cI * (P3[0] * P3[1] * OM3 * TMP4)) - TMP23 * (+cI * (P2[1]
      * V1[2] + P2[0] * V1[3]) + 0.666666667 * cI * (P3[0] * P3[1] * OM3 *
      TMP9))) + TMP26 * (TMP2 * - 1. * (+cI * (P1[1] * V2[2] + P1[0] * V2[3]) +
      0.666666667 * cI * (P3[0] * P3[1] * OM3 * TMP7)) + TMP23 * (+cI * (P1[1]
      * P2[0] + P1[0] * P2[1]) + 0.666666667 * cI * (P3[0] * P3[1] * OM3 *
      TMP10))));
  T3[4] = denom * (TMP1 * (TMP2 * (+cI * (V2[2] * V1[4] + V2[4] * V1[2]) +
      0.666666667 * cI * (P3[0] * P3[2] * OM3 * TMP4)) - TMP23 * (+cI * (P2[2]
      * V1[2] + P2[0] * V1[4]) + 0.666666667 * cI * (P3[0] * P3[2] * OM3 *
      TMP9))) + TMP26 * (TMP2 * - 1. * (+cI * (P1[2] * V2[2] + P1[0] * V2[4]) +
      0.666666667 * cI * (P3[0] * P3[2] * OM3 * TMP7)) + TMP23 * (+cI * (P1[2]
      * P2[0] + P1[0] * P2[2]) + 0.666666667 * cI * (P3[0] * P3[2] * OM3 *
      TMP10))));
  T3[5] = denom * (TMP1 * (TMP2 * (+cI * (V2[2] * V1[5] + V2[5] * V1[2]) +
      0.666666667 * cI * (P3[0] * P3[3] * OM3 * TMP4)) - TMP23 * (+cI * (P2[3]
      * V1[2] + P2[0] * V1[5]) + 0.666666667 * cI * (P3[0] * P3[3] * OM3 *
      TMP9))) + TMP26 * (TMP2 * - 1. * (+cI * (P1[3] * V2[2] + P1[0] * V2[5]) +
      0.666666667 * cI * (P3[0] * P3[3] * OM3 * TMP7)) + TMP23 * (+cI * (P1[3]
      * P2[0] + P1[0] * P2[3]) + 0.666666667 * cI * (P3[0] * P3[3] * OM3 *
      TMP10))));
  T3[6] = denom * (TMP1 * (TMP2 * (+cI * (V2[3] * V1[2] + V2[2] * V1[3]) +
      0.666666667 * cI * (P3[0] * P3[1] * OM3 * TMP4)) - TMP23 * (+cI * (P2[0]
      * V1[3] + P2[1] * V1[2]) + 0.666666667 * cI * (P3[0] * P3[1] * OM3 *
      TMP9))) + TMP26 * (TMP2 * - 1. * (+cI * (P1[0] * V2[3] + P1[1] * V2[2]) +
      0.666666667 * cI * (P3[0] * P3[1] * OM3 * TMP7)) + TMP23 * (+cI * (P1[0]
      * P2[1] + P1[1] * P2[0]) + 0.666666667 * cI * (P3[0] * P3[1] * OM3 *
      TMP10))));
  T3[7] = denom * 2. * (TMP1 * (TMP2 * (TMP4 * 0.333333333 * (+cI * (P3[1] *
      P3[1] * OM3) + 1.000000000 * cI) + cI * (V2[3] * V1[3])) + TMP23 * (TMP9
      * - 0.333333333 * (+cI * (P3[1] * P3[1] * OM3) + 1.000000000 * cI) - cI *
      (P2[1] * V1[3]))) + TMP26 * (TMP2 * (TMP7 * - 0.333333333 * (+cI * (P3[1]
      * P3[1] * OM3) + 1.000000000 * cI) - cI * (P1[1] * V2[3])) + TMP23 *
      (TMP10 * 0.333333333 * (+cI * (P3[1] * P3[1] * OM3) + 1.000000000 * cI) +
      cI * (P1[1] * P2[1]))));
  T3[8] = denom * (TMP1 * (TMP2 * (+cI * (V2[3] * V1[4] + V2[4] * V1[3]) +
      0.666666667 * cI * (P3[1] * P3[2] * OM3 * TMP4)) - TMP23 * (+cI * (P2[2]
      * V1[3] + P2[1] * V1[4]) + 0.666666667 * cI * (P3[1] * P3[2] * OM3 *
      TMP9))) + TMP26 * (TMP2 * - 1. * (+cI * (P1[2] * V2[3] + P1[1] * V2[4]) +
      0.666666667 * cI * (P3[1] * P3[2] * OM3 * TMP7)) + TMP23 * (+cI * (P1[2]
      * P2[1] + P1[1] * P2[2]) + 0.666666667 * cI * (P3[1] * P3[2] * OM3 *
      TMP10))));
  T3[9] = denom * (TMP1 * (TMP2 * (+cI * (V2[3] * V1[5] + V2[5] * V1[3]) +
      0.666666667 * cI * (P3[1] * P3[3] * OM3 * TMP4)) - TMP23 * (+cI * (P2[3]
      * V1[3] + P2[1] * V1[5]) + 0.666666667 * cI * (P3[1] * P3[3] * OM3 *
      TMP9))) + TMP26 * (TMP2 * - 1. * (+cI * (P1[3] * V2[3] + P1[1] * V2[5]) +
      0.666666667 * cI * (P3[1] * P3[3] * OM3 * TMP7)) + TMP23 * (+cI * (P1[3]
      * P2[1] + P1[1] * P2[3]) + 0.666666667 * cI * (P3[1] * P3[3] * OM3 *
      TMP10))));
  T3[10] = denom * (TMP1 * (TMP2 * (+cI * (V2[4] * V1[2] + V2[2] * V1[4]) +
      0.666666667 * cI * (P3[0] * P3[2] * OM3 * TMP4)) - TMP23 * (+cI * (P2[0]
      * V1[4] + P2[2] * V1[2]) + 0.666666667 * cI * (P3[0] * P3[2] * OM3 *
      TMP9))) + TMP26 * (TMP2 * - 1. * (+cI * (P1[0] * V2[4] + P1[2] * V2[2]) +
      0.666666667 * cI * (P3[0] * P3[2] * OM3 * TMP7)) + TMP23 * (+cI * (P1[0]
      * P2[2] + P1[2] * P2[0]) + 0.666666667 * cI * (P3[0] * P3[2] * OM3 *
      TMP10))));
  T3[11] = denom * (TMP1 * (TMP2 * (+cI * (V2[4] * V1[3] + V2[3] * V1[4]) +
      0.666666667 * cI * (P3[1] * P3[2] * OM3 * TMP4)) - TMP23 * (+cI * (P2[1]
      * V1[4] + P2[2] * V1[3]) + 0.666666667 * cI * (P3[1] * P3[2] * OM3 *
      TMP9))) + TMP26 * (TMP2 * - 1. * (+cI * (P1[1] * V2[4] + P1[2] * V2[3]) +
      0.666666667 * cI * (P3[1] * P3[2] * OM3 * TMP7)) + TMP23 * (+cI * (P1[1]
      * P2[2] + P1[2] * P2[1]) + 0.666666667 * cI * (P3[1] * P3[2] * OM3 *
      TMP10))));
  T3[12] = denom * 2. * (TMP1 * (TMP2 * (TMP4 * 0.333333333 * (+cI * (P3[2] *
      P3[2] * OM3) + 1.000000000 * cI) + cI * (V2[4] * V1[4])) + TMP23 * (TMP9
      * - 0.333333333 * (+cI * (P3[2] * P3[2] * OM3) + 1.000000000 * cI) - cI *
      (P2[2] * V1[4]))) + TMP26 * (TMP2 * (TMP7 * - 0.333333333 * (+cI * (P3[2]
      * P3[2] * OM3) + 1.000000000 * cI) - cI * (P1[2] * V2[4])) + TMP23 *
      (TMP10 * 0.333333333 * (+cI * (P3[2] * P3[2] * OM3) + 1.000000000 * cI) +
      cI * (P1[2] * P2[2]))));
  T3[13] = denom * (TMP1 * (TMP2 * (+cI * (V2[4] * V1[5] + V2[5] * V1[4]) +
      0.666666667 * cI * (P3[2] * P3[3] * OM3 * TMP4)) - TMP23 * (+cI * (P2[3]
      * V1[4] + P2[2] * V1[5]) + 0.666666667 * cI * (P3[2] * P3[3] * OM3 *
      TMP9))) + TMP26 * (TMP2 * - 1. * (+cI * (P1[3] * V2[4] + P1[2] * V2[5]) +
      0.666666667 * cI * (P3[2] * P3[3] * OM3 * TMP7)) + TMP23 * (+cI * (P1[3]
      * P2[2] + P1[2] * P2[3]) + 0.666666667 * cI * (P3[2] * P3[3] * OM3 *
      TMP10))));
  T3[14] = denom * (TMP1 * (TMP2 * (+cI * (V2[5] * V1[2] + V2[2] * V1[5]) +
      0.666666667 * cI * (P3[0] * P3[3] * OM3 * TMP4)) - TMP23 * (+cI * (P2[0]
      * V1[5] + P2[3] * V1[2]) + 0.666666667 * cI * (P3[0] * P3[3] * OM3 *
      TMP9))) + TMP26 * (TMP2 * - 1. * (+cI * (P1[0] * V2[5] + P1[3] * V2[2]) +
      0.666666667 * cI * (P3[0] * P3[3] * OM3 * TMP7)) + TMP23 * (+cI * (P1[0]
      * P2[3] + P1[3] * P2[0]) + 0.666666667 * cI * (P3[0] * P3[3] * OM3 *
      TMP10))));
  T3[15] = denom * (TMP1 * (TMP2 * (+cI * (V2[5] * V1[3] + V2[3] * V1[5]) +
      0.666666667 * cI * (P3[1] * P3[3] * OM3 * TMP4)) - TMP23 * (+cI * (P2[1]
      * V1[5] + P2[3] * V1[3]) + 0.666666667 * cI * (P3[1] * P3[3] * OM3 *
      TMP9))) + TMP26 * (TMP2 * - 1. * (+cI * (P1[1] * V2[5] + P1[3] * V2[3]) +
      0.666666667 * cI * (P3[1] * P3[3] * OM3 * TMP7)) + TMP23 * (+cI * (P1[1]
      * P2[3] + P1[3] * P2[1]) + 0.666666667 * cI * (P3[1] * P3[3] * OM3 *
      TMP10))));
  T3[16] = denom * (TMP1 * (TMP2 * (+cI * (V2[5] * V1[4] + V2[4] * V1[5]) +
      0.666666667 * cI * (P3[2] * P3[3] * OM3 * TMP4)) - TMP23 * (+cI * (P2[2]
      * V1[5] + P2[3] * V1[4]) + 0.666666667 * cI * (P3[2] * P3[3] * OM3 *
      TMP9))) + TMP26 * (TMP2 * - 1. * (+cI * (P1[2] * V2[5] + P1[3] * V2[4]) +
      0.666666667 * cI * (P3[2] * P3[3] * OM3 * TMP7)) + TMP23 * (+cI * (P1[2]
      * P2[3] + P1[3] * P2[2]) + 0.666666667 * cI * (P3[2] * P3[3] * OM3 *
      TMP10))));
  T3[17] = denom * 2. * (TMP1 * (TMP2 * (TMP4 * 0.333333333 * (+cI * (P3[3] *
      P3[3] * OM3) + 1.000000000 * cI) + cI * (V2[5] * V1[5])) + TMP23 * (TMP9
      * - 0.333333333 * (+cI * (P3[3] * P3[3] * OM3) + 1.000000000 * cI) - cI *
      (P2[3] * V1[5]))) + TMP26 * (TMP2 * (TMP7 * - 0.333333333 * (+cI * (P3[3]
      * P3[3] * OM3) + 1.000000000 * cI) - cI * (P1[3] * V2[5])) + TMP23 *
      (TMP10 * 0.333333333 * (+cI * (P3[3] * P3[3] * OM3) + 1.000000000 * cI) +
      cI * (P1[3] * P2[3]))));
}


void VVT3_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  complex<double> TMP54; 
  complex<double> TMP53; 
  complex<double> TMP52; 
  complex<double> TMP51; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = T3[0].real(); 
  P3[1] = T3[1].real(); 
  P3[2] = T3[1].imag(); 
  P3[3] = T3[0].imag(); 
  TMP51 = -1. * (P1[0] * (P3[0] * (T3[3] * (V2[5] * V1[4] - V2[4] * V1[5]) +
      (T3[4] * (V2[3] * V1[5] - V2[5] * V1[3]) + T3[5] * (V2[4] * V1[3] - V2[3]
      * V1[4]))) + (P3[1] * (T3[2] * (V2[4] * V1[5] - V2[5] * V1[4]) + (T3[4] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + T3[5] * (V2[2] * V1[4] - V2[4] *
      V1[2]))) + (P3[2] * (T3[2] * (V2[5] * V1[3] - V2[3] * V1[5]) + (T3[3] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + T3[5] * (V2[3] * V1[2] - V2[2] *
      V1[3]))) + P3[3] * (T3[2] * (V2[3] * V1[4] - V2[4] * V1[3]) + (T3[3] *
      (V2[4] * V1[2] - V2[2] * V1[4]) + T3[4] * (V2[2] * V1[3] - V2[3] *
      V1[2])))))) + (P1[1] * (P3[0] * (T3[7] * (V2[4] * V1[5] - V2[5] * V1[4])
      + (T3[8] * (V2[5] * V1[3] - V2[3] * V1[5]) + T3[9] * (V2[3] * V1[4] -
      V2[4] * V1[3]))) + (P3[1] * (T3[6] * (V2[5] * V1[4] - V2[4] * V1[5]) +
      (T3[8] * (V2[2] * V1[5] - V2[5] * V1[2]) + T3[9] * (V2[4] * V1[2] - V2[2]
      * V1[4]))) + (P3[2] * (T3[6] * (V2[3] * V1[5] - V2[5] * V1[3]) + (T3[7] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + T3[9] * (V2[2] * V1[3] - V2[3] *
      V1[2]))) + P3[3] * (T3[6] * (V2[4] * V1[3] - V2[3] * V1[4]) + (T3[7] *
      (V2[2] * V1[4] - V2[4] * V1[2]) + T3[8] * (V2[3] * V1[2] - V2[2] *
      V1[3])))))) + (P1[2] * (P3[0] * (T3[11] * (V2[4] * V1[5] - V2[5] * V1[4])
      + (T3[12] * (V2[5] * V1[3] - V2[3] * V1[5]) + T3[13] * (V2[3] * V1[4] -
      V2[4] * V1[3]))) + (P3[1] * (T3[12] * (V2[2] * V1[5] - V2[5] * V1[2]) +
      (T3[13] * (V2[4] * V1[2] - V2[2] * V1[4]) + T3[10] * (V2[5] * V1[4] -
      V2[4] * V1[5]))) + (P3[2] * (T3[11] * (V2[5] * V1[2] - V2[2] * V1[5]) +
      (T3[13] * (V2[2] * V1[3] - V2[3] * V1[2]) + T3[10] * (V2[3] * V1[5] -
      V2[5] * V1[3]))) + P3[3] * (T3[11] * (V2[2] * V1[4] - V2[4] * V1[2]) +
      (T3[12] * (V2[3] * V1[2] - V2[2] * V1[3]) + T3[10] * (V2[4] * V1[3] -
      V2[3] * V1[4])))))) + P1[3] * (P3[0] * (T3[15] * (V2[4] * V1[5] - V2[5] *
      V1[4]) + (T3[16] * (V2[5] * V1[3] - V2[3] * V1[5]) + T3[17] * (V2[3] *
      V1[4] - V2[4] * V1[3]))) + (P3[1] * (T3[14] * (V2[5] * V1[4] - V2[4] *
      V1[5]) + (T3[16] * (V2[2] * V1[5] - V2[5] * V1[2]) + T3[17] * (V2[4] *
      V1[2] - V2[2] * V1[4]))) + (P3[2] * (T3[14] * (V2[3] * V1[5] - V2[5] *
      V1[3]) + (T3[15] * (V2[5] * V1[2] - V2[2] * V1[5]) + T3[17] * (V2[2] *
      V1[3] - V2[3] * V1[2]))) + P3[3] * (T3[14] * (V2[4] * V1[3] - V2[3] *
      V1[4]) + (T3[15] * (V2[2] * V1[4] - V2[4] * V1[2]) + T3[16] * (V2[3] *
      V1[2] - V2[2] * V1[3])))))))));
  TMP53 = -1. * (P1[0] * (P3[0] * (T3[14] * (V2[4] * V1[3] - V2[3] * V1[4]) +
      (T3[6] * (V2[5] * V1[4] - V2[4] * V1[5]) + T3[10] * (V2[3] * V1[5] -
      V2[5] * V1[3]))) + (P3[1] * (T3[2] * (V2[4] * V1[5] - V2[5] * V1[4]) +
      (T3[14] * (V2[2] * V1[4] - V2[4] * V1[2]) + T3[10] * (V2[5] * V1[2] -
      V2[2] * V1[5]))) + (P3[2] * (T3[2] * (V2[5] * V1[3] - V2[3] * V1[5]) +
      (T3[14] * (V2[3] * V1[2] - V2[2] * V1[3]) + T3[6] * (V2[2] * V1[5] -
      V2[5] * V1[2]))) + P3[3] * (T3[2] * (V2[3] * V1[4] - V2[4] * V1[3]) +
      (T3[6] * (V2[4] * V1[2] - V2[2] * V1[4]) + T3[10] * (V2[2] * V1[3] -
      V2[3] * V1[2])))))) + (P1[1] * (P3[0] * (T3[11] * (V2[5] * V1[3] - V2[3]
      * V1[5]) + (T3[15] * (V2[3] * V1[4] - V2[4] * V1[3]) + T3[7] * (V2[4] *
      V1[5] - V2[5] * V1[4]))) + (P3[1] * (T3[11] * (V2[2] * V1[5] - V2[5] *
      V1[2]) + (T3[15] * (V2[4] * V1[2] - V2[2] * V1[4]) + T3[3] * (V2[5] *
      V1[4] - V2[4] * V1[5]))) + (P3[2] * (T3[15] * (V2[2] * V1[3] - V2[3] *
      V1[2]) + (T3[3] * (V2[3] * V1[5] - V2[5] * V1[3]) + T3[7] * (V2[5] *
      V1[2] - V2[2] * V1[5]))) + P3[3] * (T3[11] * (V2[3] * V1[2] - V2[2] *
      V1[3]) + (T3[3] * (V2[4] * V1[3] - V2[3] * V1[4]) + T3[7] * (V2[2] *
      V1[4] - V2[4] * V1[2])))))) + (P1[2] * (P3[0] * (T3[12] * (V2[5] * V1[3]
      - V2[3] * V1[5]) + (T3[16] * (V2[3] * V1[4] - V2[4] * V1[3]) + T3[8] *
      (V2[4] * V1[5] - V2[5] * V1[4]))) + (P3[1] * (T3[12] * (V2[2] * V1[5] -
      V2[5] * V1[2]) + (T3[16] * (V2[4] * V1[2] - V2[2] * V1[4]) + T3[4] *
      (V2[5] * V1[4] - V2[4] * V1[5]))) + (P3[2] * (T3[16] * (V2[2] * V1[3] -
      V2[3] * V1[2]) + (T3[4] * (V2[3] * V1[5] - V2[5] * V1[3]) + T3[8] *
      (V2[5] * V1[2] - V2[2] * V1[5]))) + P3[3] * (T3[12] * (V2[3] * V1[2] -
      V2[2] * V1[3]) + (T3[4] * (V2[4] * V1[3] - V2[3] * V1[4]) + T3[8] *
      (V2[2] * V1[4] - V2[4] * V1[2])))))) + P1[3] * (P3[0] * (T3[13] * (V2[5]
      * V1[3] - V2[3] * V1[5]) + (T3[17] * (V2[3] * V1[4] - V2[4] * V1[3]) +
      T3[9] * (V2[4] * V1[5] - V2[5] * V1[4]))) + (P3[1] * (T3[13] * (V2[2] *
      V1[5] - V2[5] * V1[2]) + (T3[17] * (V2[4] * V1[2] - V2[2] * V1[4]) +
      T3[5] * (V2[5] * V1[4] - V2[4] * V1[5]))) + (P3[2] * (T3[17] * (V2[2] *
      V1[3] - V2[3] * V1[2]) + (T3[5] * (V2[3] * V1[5] - V2[5] * V1[3]) + T3[9]
      * (V2[5] * V1[2] - V2[2] * V1[5]))) + P3[3] * (T3[13] * (V2[3] * V1[2] -
      V2[2] * V1[3]) + (T3[5] * (V2[4] * V1[3] - V2[3] * V1[4]) + T3[9] *
      (V2[2] * V1[4] - V2[4] * V1[2])))))))));
  TMP52 = -1. * (P2[0] * (P3[0] * (T3[3] * (V2[5] * V1[4] - V2[4] * V1[5]) +
      (T3[4] * (V2[3] * V1[5] - V2[5] * V1[3]) + T3[5] * (V2[4] * V1[3] - V2[3]
      * V1[4]))) + (P3[1] * (T3[2] * (V2[4] * V1[5] - V2[5] * V1[4]) + (T3[4] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + T3[5] * (V2[2] * V1[4] - V2[4] *
      V1[2]))) + (P3[2] * (T3[2] * (V2[5] * V1[3] - V2[3] * V1[5]) + (T3[3] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + T3[5] * (V2[3] * V1[2] - V2[2] *
      V1[3]))) + P3[3] * (T3[2] * (V2[3] * V1[4] - V2[4] * V1[3]) + (T3[3] *
      (V2[4] * V1[2] - V2[2] * V1[4]) + T3[4] * (V2[2] * V1[3] - V2[3] *
      V1[2])))))) + (P2[1] * (P3[0] * (T3[7] * (V2[4] * V1[5] - V2[5] * V1[4])
      + (T3[8] * (V2[5] * V1[3] - V2[3] * V1[5]) + T3[9] * (V2[3] * V1[4] -
      V2[4] * V1[3]))) + (P3[1] * (T3[6] * (V2[5] * V1[4] - V2[4] * V1[5]) +
      (T3[8] * (V2[2] * V1[5] - V2[5] * V1[2]) + T3[9] * (V2[4] * V1[2] - V2[2]
      * V1[4]))) + (P3[2] * (T3[6] * (V2[3] * V1[5] - V2[5] * V1[3]) + (T3[7] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + T3[9] * (V2[2] * V1[3] - V2[3] *
      V1[2]))) + P3[3] * (T3[6] * (V2[4] * V1[3] - V2[3] * V1[4]) + (T3[7] *
      (V2[2] * V1[4] - V2[4] * V1[2]) + T3[8] * (V2[3] * V1[2] - V2[2] *
      V1[3])))))) + (P2[2] * (P3[0] * (T3[11] * (V2[4] * V1[5] - V2[5] * V1[4])
      + (T3[12] * (V2[5] * V1[3] - V2[3] * V1[5]) + T3[13] * (V2[3] * V1[4] -
      V2[4] * V1[3]))) + (P3[1] * (T3[12] * (V2[2] * V1[5] - V2[5] * V1[2]) +
      (T3[13] * (V2[4] * V1[2] - V2[2] * V1[4]) + T3[10] * (V2[5] * V1[4] -
      V2[4] * V1[5]))) + (P3[2] * (T3[11] * (V2[5] * V1[2] - V2[2] * V1[5]) +
      (T3[13] * (V2[2] * V1[3] - V2[3] * V1[2]) + T3[10] * (V2[3] * V1[5] -
      V2[5] * V1[3]))) + P3[3] * (T3[11] * (V2[2] * V1[4] - V2[4] * V1[2]) +
      (T3[12] * (V2[3] * V1[2] - V2[2] * V1[3]) + T3[10] * (V2[4] * V1[3] -
      V2[3] * V1[4])))))) + P2[3] * (P3[0] * (T3[15] * (V2[4] * V1[5] - V2[5] *
      V1[4]) + (T3[16] * (V2[5] * V1[3] - V2[3] * V1[5]) + T3[17] * (V2[3] *
      V1[4] - V2[4] * V1[3]))) + (P3[1] * (T3[14] * (V2[5] * V1[4] - V2[4] *
      V1[5]) + (T3[16] * (V2[2] * V1[5] - V2[5] * V1[2]) + T3[17] * (V2[4] *
      V1[2] - V2[2] * V1[4]))) + (P3[2] * (T3[14] * (V2[3] * V1[5] - V2[5] *
      V1[3]) + (T3[15] * (V2[5] * V1[2] - V2[2] * V1[5]) + T3[17] * (V2[2] *
      V1[3] - V2[3] * V1[2]))) + P3[3] * (T3[14] * (V2[4] * V1[3] - V2[3] *
      V1[4]) + (T3[15] * (V2[2] * V1[4] - V2[4] * V1[2]) + T3[16] * (V2[3] *
      V1[2] - V2[2] * V1[3])))))))));
  TMP54 = -1. * (P2[0] * (P3[0] * (T3[14] * (V2[4] * V1[3] - V2[3] * V1[4]) +
      (T3[6] * (V2[5] * V1[4] - V2[4] * V1[5]) + T3[10] * (V2[3] * V1[5] -
      V2[5] * V1[3]))) + (P3[1] * (T3[2] * (V2[4] * V1[5] - V2[5] * V1[4]) +
      (T3[14] * (V2[2] * V1[4] - V2[4] * V1[2]) + T3[10] * (V2[5] * V1[2] -
      V2[2] * V1[5]))) + (P3[2] * (T3[2] * (V2[5] * V1[3] - V2[3] * V1[5]) +
      (T3[14] * (V2[3] * V1[2] - V2[2] * V1[3]) + T3[6] * (V2[2] * V1[5] -
      V2[5] * V1[2]))) + P3[3] * (T3[2] * (V2[3] * V1[4] - V2[4] * V1[3]) +
      (T3[6] * (V2[4] * V1[2] - V2[2] * V1[4]) + T3[10] * (V2[2] * V1[3] -
      V2[3] * V1[2])))))) + (P2[1] * (P3[0] * (T3[11] * (V2[5] * V1[3] - V2[3]
      * V1[5]) + (T3[15] * (V2[3] * V1[4] - V2[4] * V1[3]) + T3[7] * (V2[4] *
      V1[5] - V2[5] * V1[4]))) + (P3[1] * (T3[11] * (V2[2] * V1[5] - V2[5] *
      V1[2]) + (T3[15] * (V2[4] * V1[2] - V2[2] * V1[4]) + T3[3] * (V2[5] *
      V1[4] - V2[4] * V1[5]))) + (P3[2] * (T3[15] * (V2[2] * V1[3] - V2[3] *
      V1[2]) + (T3[3] * (V2[3] * V1[5] - V2[5] * V1[3]) + T3[7] * (V2[5] *
      V1[2] - V2[2] * V1[5]))) + P3[3] * (T3[11] * (V2[3] * V1[2] - V2[2] *
      V1[3]) + (T3[3] * (V2[4] * V1[3] - V2[3] * V1[4]) + T3[7] * (V2[2] *
      V1[4] - V2[4] * V1[2])))))) + (P2[2] * (P3[0] * (T3[12] * (V2[5] * V1[3]
      - V2[3] * V1[5]) + (T3[16] * (V2[3] * V1[4] - V2[4] * V1[3]) + T3[8] *
      (V2[4] * V1[5] - V2[5] * V1[4]))) + (P3[1] * (T3[12] * (V2[2] * V1[5] -
      V2[5] * V1[2]) + (T3[16] * (V2[4] * V1[2] - V2[2] * V1[4]) + T3[4] *
      (V2[5] * V1[4] - V2[4] * V1[5]))) + (P3[2] * (T3[16] * (V2[2] * V1[3] -
      V2[3] * V1[2]) + (T3[4] * (V2[3] * V1[5] - V2[5] * V1[3]) + T3[8] *
      (V2[5] * V1[2] - V2[2] * V1[5]))) + P3[3] * (T3[12] * (V2[3] * V1[2] -
      V2[2] * V1[3]) + (T3[4] * (V2[4] * V1[3] - V2[3] * V1[4]) + T3[8] *
      (V2[2] * V1[4] - V2[4] * V1[2])))))) + P2[3] * (P3[0] * (T3[13] * (V2[5]
      * V1[3] - V2[3] * V1[5]) + (T3[17] * (V2[3] * V1[4] - V2[4] * V1[3]) +
      T3[9] * (V2[4] * V1[5] - V2[5] * V1[4]))) + (P3[1] * (T3[13] * (V2[2] *
      V1[5] - V2[5] * V1[2]) + (T3[17] * (V2[4] * V1[2] - V2[2] * V1[4]) +
      T3[5] * (V2[5] * V1[4] - V2[4] * V1[5]))) + (P3[2] * (T3[17] * (V2[2] *
      V1[3] - V2[3] * V1[2]) + (T3[5] * (V2[3] * V1[5] - V2[5] * V1[3]) + T3[9]
      * (V2[5] * V1[2] - V2[2] * V1[5]))) + P3[3] * (T3[13] * (V2[3] * V1[2] -
      V2[2] * V1[3]) + (T3[5] * (V2[4] * V1[3] - V2[3] * V1[4]) + T3[9] *
      (V2[2] * V1[4] - V2[4] * V1[2])))))))));
  vertex = COUP * (-cI * (TMP51 + TMP53) + cI * (TMP52 + TMP54)); 
}


void VVT7_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP2; 
  complex<double> TMP1; 
  double P1[4]; 
  complex<double> TMP10; 
  double P2[4]; 
  double P3[4]; 
  complex<double> denom; 
  double OM3; 
  complex<double> TMP4; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +V1[0] + V2[0]; 
  T3[1] = +V1[1] + V2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP4 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP1 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP10 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP2 = (P3[0] * P2[0] - P3[1] * P2[1] - P3[2] * P2[2] - P3[3] * P2[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * 2. * TMP4 * (OM3 * (P3[0] * (P3[0] * 0.333333333 * (+cI *
      (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (TMP1 * P2[0] + P1[0] *
      TMP2))) + 0.333333333 * cI * (TMP1 * TMP2)) + (-0.333333333 * cI *
      (TMP10) + cI * (P1[0] * P2[0])));
  T3[3] = denom * TMP4 * (OM3 * (P3[0] * (P3[1] * 0.666666667 * (+cI * (TMP10)
      + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (P1[1] * TMP2 + TMP1 * P2[1])))
      - P3[1] * (+cI * (TMP1 * P2[0] + P1[0] * TMP2))) + (+cI * (P1[1] * P2[0]
      + P1[0] * P2[1])));
  T3[4] = denom * TMP4 * (OM3 * (P3[0] * (P3[2] * 0.666666667 * (+cI * (TMP10)
      + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (P1[2] * TMP2 + TMP1 * P2[2])))
      - P3[2] * (+cI * (TMP1 * P2[0] + P1[0] * TMP2))) + (+cI * (P1[2] * P2[0]
      + P1[0] * P2[2])));
  T3[5] = denom * TMP4 * (OM3 * (P3[0] * (P3[3] * 0.666666667 * (+cI * (TMP10)
      + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (P1[3] * TMP2 + TMP1 * P2[3])))
      - P3[3] * (+cI * (TMP1 * P2[0] + P1[0] * TMP2))) + (+cI * (P1[3] * P2[0]
      + P1[0] * P2[3])));
  T3[6] = denom * TMP4 * (OM3 * (P3[0] * (P3[1] * 0.666666667 * (+cI * (TMP10)
      + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (TMP1 * P2[1] + P1[1] * TMP2)))
      - P3[1] * (+cI * (P1[0] * TMP2 + TMP1 * P2[0]))) + (+cI * (P1[0] * P2[1]
      + P1[1] * P2[0])));
  T3[7] = denom * 2. * TMP4 * (OM3 * (P3[1] * (P3[1] * 0.333333333 * (+cI *
      (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (TMP1 * P2[1] + P1[1] *
      TMP2))) - 0.333333333 * cI * (TMP1 * TMP2)) + (+cI * (P1[1] * P2[1]) +
      0.333333333 * cI * (TMP10)));
  T3[8] = denom * TMP4 * (OM3 * (P3[1] * (P3[2] * 0.666666667 * (+cI * (TMP10)
      + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (P1[2] * TMP2 + TMP1 * P2[2])))
      - P3[2] * (+cI * (TMP1 * P2[1] + P1[1] * TMP2))) + (+cI * (P1[2] * P2[1]
      + P1[1] * P2[2])));
  T3[9] = denom * TMP4 * (OM3 * (P3[1] * (P3[3] * 0.666666667 * (+cI * (TMP10)
      + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (P1[3] * TMP2 + TMP1 * P2[3])))
      - P3[3] * (+cI * (TMP1 * P2[1] + P1[1] * TMP2))) + (+cI * (P1[3] * P2[1]
      + P1[1] * P2[3])));
  T3[10] = denom * TMP4 * (OM3 * (P3[0] * (P3[2] * 0.666666667 * (+cI * (TMP10)
      + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (TMP1 * P2[2] + P1[2] * TMP2)))
      - P3[2] * (+cI * (P1[0] * TMP2 + TMP1 * P2[0]))) + (+cI * (P1[0] * P2[2]
      + P1[2] * P2[0])));
  T3[11] = denom * TMP4 * (OM3 * (P3[1] * (P3[2] * 0.666666667 * (+cI * (TMP10)
      + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (TMP1 * P2[2] + P1[2] * TMP2)))
      - P3[2] * (+cI * (P1[1] * TMP2 + TMP1 * P2[1]))) + (+cI * (P1[1] * P2[2]
      + P1[2] * P2[1])));
  T3[12] = denom * 2. * TMP4 * (OM3 * (P3[2] * (P3[2] * 0.333333333 * (+cI *
      (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (TMP1 * P2[2] + P1[2] *
      TMP2))) - 0.333333333 * cI * (TMP1 * TMP2)) + (+cI * (P1[2] * P2[2]) +
      0.333333333 * cI * (TMP10)));
  T3[13] = denom * TMP4 * (OM3 * (P3[2] * (P3[3] * 0.666666667 * (+cI * (TMP10)
      + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (P1[3] * TMP2 + TMP1 * P2[3])))
      - P3[3] * (+cI * (TMP1 * P2[2] + P1[2] * TMP2))) + (+cI * (P1[3] * P2[2]
      + P1[2] * P2[3])));
  T3[14] = denom * TMP4 * (OM3 * (P3[0] * (P3[3] * 0.666666667 * (+cI * (TMP10)
      + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (TMP1 * P2[3] + P1[3] * TMP2)))
      - P3[3] * (+cI * (P1[0] * TMP2 + TMP1 * P2[0]))) + (+cI * (P1[0] * P2[3]
      + P1[3] * P2[0])));
  T3[15] = denom * TMP4 * (OM3 * (P3[1] * (P3[3] * 0.666666667 * (+cI * (TMP10)
      + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (TMP1 * P2[3] + P1[3] * TMP2)))
      - P3[3] * (+cI * (P1[1] * TMP2 + TMP1 * P2[1]))) + (+cI * (P1[1] * P2[3]
      + P1[3] * P2[1])));
  T3[16] = denom * TMP4 * (OM3 * (P3[2] * (P3[3] * 0.666666667 * (+cI * (TMP10)
      + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (TMP1 * P2[3] + P1[3] * TMP2)))
      - P3[3] * (+cI * (P1[2] * TMP2 + TMP1 * P2[2]))) + (+cI * (P1[2] * P2[3]
      + P1[3] * P2[2])));
  T3[17] = denom * 2. * TMP4 * (OM3 * (P3[3] * (P3[3] * 0.333333333 * (+cI *
      (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (-cI * (TMP1 * P2[3] + P1[3] *
      TMP2))) - 0.333333333 * cI * (TMP1 * TMP2)) + (+cI * (P1[3] * P2[3]) +
      0.333333333 * cI * (TMP10)));
}


void VVT10_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP1; 
  complex<double> TMP10; 
  double P3[4]; 
  complex<double> TMP21; 
  complex<double> TMP43; 
  complex<double> TMP9; 
  complex<double> TMP2; 
  double P2[4]; 
  complex<double> TMP46; 
  complex<double> TMP45; 
  complex<double> TMP42; 
  double P1[4]; 
  complex<double> TMP23; 
  complex<double> TMP7; 
  complex<double> TMP41; 
  complex<double> TMP22; 
  complex<double> TMP44; 
  complex<double> TMP26; 
  complex<double> TMP4; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = T3[0].real(); 
  P3[1] = T3[1].real(); 
  P3[2] = T3[1].imag(); 
  P3[3] = T3[0].imag(); 
  TMP42 = (P2[0] * (P2[1] * - 1. * (T3[3] + T3[6]) + (P2[2] * - 1. * (T3[4] +
      T3[10]) + (P2[3] * - 1. * (T3[5] + T3[14]) + P2[0] * T3[2]))) + (P2[1] *
      (P2[2] * (T3[8] + T3[11]) + (P2[3] * (T3[9] + T3[15]) + P2[1] * T3[7])) +
      (P2[2] * (P2[3] * (T3[13] + T3[16]) + P2[2] * T3[12]) + P2[3] * P2[3] *
      T3[17])));
  TMP43 = (P1[0] * - 1. * (V1[3] * T3[6] + V1[4] * T3[10] + V1[5] * T3[14] -
      V1[2] * T3[2]) + (P1[1] * (V1[3] * T3[7] + V1[4] * T3[11] + V1[5] *
      T3[15] - V1[2] * T3[3]) + (P1[2] * (V1[3] * T3[8] + V1[4] * T3[12] +
      V1[5] * T3[16] - V1[2] * T3[4]) + P1[3] * (V1[3] * T3[9] + V1[4] * T3[13]
      + V1[5] * T3[17] - V1[2] * T3[5]))));
  TMP26 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP41 = (P1[0] * (P1[1] * - 1. * (T3[6] + T3[3]) + (P1[2] * - 1. * (T3[10] +
      T3[4]) + (P1[3] * - 1. * (T3[14] + T3[5]) + P1[0] * T3[2]))) + (P1[1] *
      (P1[2] * (T3[11] + T3[8]) + (P1[3] * (T3[15] + T3[9]) + P1[1] * T3[7])) +
      (P1[2] * (P1[3] * (T3[16] + T3[13]) + P1[2] * T3[12]) + P1[3] * P1[3] *
      T3[17])));
  TMP46 = (P2[0] * - 1. * (V2[3] * T3[3] + V2[4] * T3[4] + V2[5] * T3[5] -
      V2[2] * T3[2]) + (P2[1] * (V2[3] * T3[7] + V2[4] * T3[8] + V2[5] * T3[9]
      - V2[2] * T3[6]) + (P2[2] * (V2[3] * T3[11] + V2[4] * T3[12] + V2[5] *
      T3[13] - V2[2] * T3[10]) + P2[3] * (V2[3] * T3[15] + V2[4] * T3[16] +
      V2[5] * T3[17] - V2[2] * T3[14]))));
  TMP21 = (P1[0] * - 1. * (P2[1] * T3[6] + P2[2] * T3[10] + P2[3] * T3[14] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[11] + P2[3] *
      T3[15] - P2[0] * T3[3]) + (P1[2] * (P2[1] * T3[8] + P2[2] * T3[12] +
      P2[3] * T3[16] - P2[0] * T3[4]) + P1[3] * (P2[1] * T3[9] + P2[2] * T3[13]
      + P2[3] * T3[17] - P2[0] * T3[5]))));
  TMP22 = (P1[0] * - 1. * (P2[1] * T3[3] + P2[2] * T3[4] + P2[3] * T3[5] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[8] + P2[3] * T3[9]
      - P2[0] * T3[6]) + (P1[2] * (P2[1] * T3[11] + P2[2] * T3[12] + P2[3] *
      T3[13] - P2[0] * T3[10]) + P1[3] * (P2[1] * T3[15] + P2[2] * T3[16] +
      P2[3] * T3[17] - P2[0] * T3[14]))));
  TMP23 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP45 = (P2[0] * - 1. * (V2[3] * T3[6] + V2[4] * T3[10] + V2[5] * T3[14] -
      V2[2] * T3[2]) + (P2[1] * (V2[3] * T3[7] + V2[4] * T3[11] + V2[5] *
      T3[15] - V2[2] * T3[3]) + (P2[2] * (V2[3] * T3[8] + V2[4] * T3[12] +
      V2[5] * T3[16] - V2[2] * T3[4]) + P2[3] * (V2[3] * T3[9] + V2[4] * T3[13]
      + V2[5] * T3[17] - V2[2] * T3[5]))));
  TMP44 = (P1[0] * - 1. * (V1[3] * T3[3] + V1[4] * T3[4] + V1[5] * T3[5] -
      V1[2] * T3[2]) + (P1[1] * (V1[3] * T3[7] + V1[4] * T3[8] + V1[5] * T3[9]
      - V1[2] * T3[6]) + (P1[2] * (V1[3] * T3[11] + V1[4] * T3[12] + V1[5] *
      T3[13] - V1[2] * T3[10]) + P1[3] * (V1[3] * T3[15] + V1[4] * T3[16] +
      V1[5] * T3[17] - V1[2] * T3[14]))));
  TMP9 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP4 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP7 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP1 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP10 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP2 = (P3[0] * P2[0] - P3[1] * P2[1] - P3[2] * P2[2] - P3[3] * P2[3]); 
  vertex = COUP * (TMP10 * (TMP4 * (-cI * (TMP21 + TMP22) + cI * (TMP41 +
      TMP42)) + (+0.500000000 * (TMP26 * (+cI * (TMP45 + TMP46))) + TMP23 *
      0.500000000 * (+cI * (TMP43 + TMP44)))) + (TMP7 * (TMP9 * (-cI * (TMP41 +
      TMP42) + cI * (TMP21 + TMP22)) + (TMP2 * - 0.500000000 * (+cI * (TMP43 +
      TMP44)) - cI * (TMP26 * TMP42))) + (TMP1 * (TMP9 * - 0.500000000 * (+cI *
      (TMP45 + TMP46)) + cI * (TMP4 * TMP42)) + TMP41 * (-cI * (TMP9 * TMP23) +
      cI * (TMP2 * TMP4)))));
}

void VVT10_11_12_13_2_3_6_7_8_9_0(complex<double> V1[], complex<double> V2[],
    complex<double> T3[], complex<double> COUP1, complex<double> COUP2,
    complex<double> COUP3, complex<double> COUP4, complex<double> COUP5,
    complex<double> COUP6, complex<double> COUP7, complex<double> COUP8,
    complex<double> COUP9, complex<double> COUP10, complex<double> & vertex)
{
//   complex<double> cI = complex<double> (0., 1.); 
//   double P3[4]; 
  complex<double> tmp; 
//   double P2[4]; 
//   double P1[4]; 
  VVT10_0(V1, V2, T3, COUP1, vertex); 
  VVT11_0(V1, V2, T3, COUP2, tmp); 
  vertex = vertex + tmp; 
  VVT12_0(V1, V2, T3, COUP3, tmp); 
  vertex = vertex + tmp; 
  VVT13_0(V1, V2, T3, COUP4, tmp); 
  vertex = vertex + tmp; 
  VVT2_0(V1, V2, T3, COUP5, tmp); 
  vertex = vertex + tmp; 
  VVT3_0(V1, V2, T3, COUP6, tmp); 
  vertex = vertex + tmp; 
  VVT6_0(V1, V2, T3, COUP7, tmp); 
  vertex = vertex + tmp; 
  VVT7_0(V1, V2, T3, COUP8, tmp); 
  vertex = vertex + tmp; 
  VVT8_0(V1, V2, T3, COUP9, tmp); 
  vertex = vertex + tmp; 
  VVT9_0(V1, V2, T3, COUP10, tmp); 
  vertex = vertex + tmp; 
}

void VVT8_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP2; 
  complex<double> TMP1; 
  double P1[4]; 
  complex<double> TMP10; 
  double P2[4]; 
  complex<double> TMP7; 
  double P3[4]; 
  complex<double> denom; 
  double OM3; 
  complex<double> TMP4; 
  complex<double> TMP9; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +V1[0] + V2[0]; 
  T3[1] = +V1[1] + V2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP9 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP4 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP7 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP1 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP10 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP2 = (P3[0] * P2[0] - P3[1] * P2[1] - P3[2] * P2[2] - P3[3] * P2[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * 2. * (OM3 * (P3[0] * (TMP10 * (TMP4 * (P3[0] * - 0.333333333
      * (+cI * (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (+cI * (TMP1 * P2[0]
      + P1[0] * TMP2))) + 0.333333333 * cI * (P3[0] * TMP7 * TMP9)) + TMP7 *
      TMP9 * (TMP1 * (-cI * (P2[0]) + 0.666666667 * cI * (P3[0] * OM3 * TMP2))
      - cI * (P1[0] * TMP2))) + 0.333333333 * (TMP1 * TMP2 * (-cI * (TMP4 *
      TMP10) + cI * (TMP7 * TMP9)))) + (TMP10 * (TMP4 * (-cI * (P1[0] * P2[0])
      + 0.333333333 * cI * (TMP10)) - 0.333333333 * cI * (TMP7 * TMP9)) + cI *
      (P1[0] * P2[0] * TMP7 * TMP9)));
  T3[3] = denom * (OM3 * (P3[0] * (TMP10 * (TMP4 * (P3[1] * - 0.666666667 *
      (+cI * (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (+cI * (P1[1] * TMP2 +
      TMP1 * P2[1]))) + 0.666666667 * cI * (P3[1] * TMP7 * TMP9)) + TMP7 * TMP9
      * (TMP1 * (-cI * (P2[1]) + 1.333333333 * cI * (P3[1] * OM3 * TMP2)) - cI
      * (P1[1] * TMP2))) + P3[1] * (P1[0] * TMP2 * (-cI * (TMP7 * TMP9) + cI *
      (TMP4 * TMP10)) + P2[0] * TMP1 * (-cI * (TMP7 * TMP9) + cI * (TMP4 *
      TMP10)))) + (P1[0] * P2[1] * (-cI * (TMP4 * TMP10) + cI * (TMP7 * TMP9))
      + P1[1] * P2[0] * (-cI * (TMP4 * TMP10) + cI * (TMP7 * TMP9))));
  T3[4] = denom * (OM3 * (P3[0] * (TMP10 * (TMP4 * (P3[2] * - 0.666666667 *
      (+cI * (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (+cI * (P1[2] * TMP2 +
      TMP1 * P2[2]))) + 0.666666667 * cI * (P3[2] * TMP7 * TMP9)) + TMP7 * TMP9
      * (TMP1 * (-cI * (P2[2]) + 1.333333333 * cI * (P3[2] * OM3 * TMP2)) - cI
      * (P1[2] * TMP2))) + P3[2] * (P1[0] * TMP2 * (-cI * (TMP7 * TMP9) + cI *
      (TMP4 * TMP10)) + P2[0] * TMP1 * (-cI * (TMP7 * TMP9) + cI * (TMP4 *
      TMP10)))) + (P1[0] * P2[2] * (-cI * (TMP4 * TMP10) + cI * (TMP7 * TMP9))
      + P1[2] * P2[0] * (-cI * (TMP4 * TMP10) + cI * (TMP7 * TMP9))));
  T3[5] = denom * (OM3 * (P3[0] * (TMP10 * (TMP4 * (P3[3] * - 0.666666667 *
      (+cI * (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (+cI * (P1[3] * TMP2 +
      TMP1 * P2[3]))) + 0.666666667 * cI * (P3[3] * TMP7 * TMP9)) + TMP7 * TMP9
      * (TMP1 * (-cI * (P2[3]) + 1.333333333 * cI * (P3[3] * OM3 * TMP2)) - cI
      * (P1[3] * TMP2))) + P3[3] * (P1[0] * TMP2 * (-cI * (TMP7 * TMP9) + cI *
      (TMP4 * TMP10)) + P2[0] * TMP1 * (-cI * (TMP7 * TMP9) + cI * (TMP4 *
      TMP10)))) + (P1[0] * P2[3] * (-cI * (TMP4 * TMP10) + cI * (TMP7 * TMP9))
      + P1[3] * P2[0] * (-cI * (TMP4 * TMP10) + cI * (TMP7 * TMP9))));
  T3[6] = denom * (OM3 * (P3[0] * (TMP10 * (TMP4 * (P3[1] * - 0.666666667 *
      (+cI * (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (+cI * (TMP1 * P2[1] +
      P1[1] * TMP2))) + 0.666666667 * cI * (P3[1] * TMP7 * TMP9)) + TMP7 * TMP9
      * (TMP1 * (-cI * (P2[1]) + 1.333333333 * cI * (P3[1] * OM3 * TMP2)) - cI
      * (P1[1] * TMP2))) + P3[1] * (P1[0] * TMP2 * (-cI * (TMP7 * TMP9) + cI *
      (TMP4 * TMP10)) + P2[0] * TMP1 * (-cI * (TMP7 * TMP9) + cI * (TMP4 *
      TMP10)))) + (P1[0] * P2[1] * (-cI * (TMP4 * TMP10) + cI * (TMP7 * TMP9))
      + P1[1] * P2[0] * (-cI * (TMP4 * TMP10) + cI * (TMP7 * TMP9))));
  T3[7] = denom * 2. * (OM3 * (P3[1] * (TMP10 * (TMP4 * (P3[1] * - 0.333333333
      * (+cI * (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (+cI * (TMP1 * P2[1]
      + P1[1] * TMP2))) + 0.333333333 * cI * (P3[1] * TMP7 * TMP9)) + TMP7 *
      TMP9 * (TMP1 * (-cI * (P2[1]) + 0.666666667 * cI * (P3[1] * OM3 * TMP2))
      - cI * (P1[1] * TMP2))) + 0.333333333 * (TMP1 * TMP2 * (-cI * (TMP7 *
      TMP9) + cI * (TMP4 * TMP10)))) + (TMP10 * (TMP4 * - 1. * (+cI * (P1[1] *
      P2[1]) + 0.333333333 * cI * (TMP10)) + 0.333333333 * cI * (TMP7 * TMP9))
      + cI * (P1[1] * P2[1] * TMP7 * TMP9)));
  T3[8] = denom * (OM3 * (P3[1] * (TMP10 * (TMP4 * (P3[2] * - 0.666666667 *
      (+cI * (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (+cI * (P1[2] * TMP2 +
      TMP1 * P2[2]))) + 0.666666667 * cI * (P3[2] * TMP7 * TMP9)) + TMP7 * TMP9
      * (TMP1 * (-cI * (P2[2]) + 1.333333333 * cI * (P3[2] * OM3 * TMP2)) - cI
      * (P1[2] * TMP2))) + P3[2] * (P1[1] * TMP2 * (-cI * (TMP7 * TMP9) + cI *
      (TMP4 * TMP10)) + P2[1] * TMP1 * (-cI * (TMP7 * TMP9) + cI * (TMP4 *
      TMP10)))) + (P1[1] * P2[2] * (-cI * (TMP4 * TMP10) + cI * (TMP7 * TMP9))
      + P1[2] * P2[1] * (-cI * (TMP4 * TMP10) + cI * (TMP7 * TMP9))));
  T3[9] = denom * (OM3 * (P3[1] * (TMP10 * (TMP4 * (P3[3] * - 0.666666667 *
      (+cI * (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (+cI * (P1[3] * TMP2 +
      TMP1 * P2[3]))) + 0.666666667 * cI * (P3[3] * TMP7 * TMP9)) + TMP7 * TMP9
      * (TMP1 * (-cI * (P2[3]) + 1.333333333 * cI * (P3[3] * OM3 * TMP2)) - cI
      * (P1[3] * TMP2))) + P3[3] * (P1[1] * TMP2 * (-cI * (TMP7 * TMP9) + cI *
      (TMP4 * TMP10)) + P2[1] * TMP1 * (-cI * (TMP7 * TMP9) + cI * (TMP4 *
      TMP10)))) + (P1[1] * P2[3] * (-cI * (TMP4 * TMP10) + cI * (TMP7 * TMP9))
      + P1[3] * P2[1] * (-cI * (TMP4 * TMP10) + cI * (TMP7 * TMP9))));
  T3[10] = denom * (OM3 * (P3[0] * (TMP10 * (TMP4 * (P3[2] * - 0.666666667 *
      (+cI * (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (+cI * (TMP1 * P2[2] +
      P1[2] * TMP2))) + 0.666666667 * cI * (P3[2] * TMP7 * TMP9)) + TMP7 * TMP9
      * (TMP1 * (-cI * (P2[2]) + 1.333333333 * cI * (P3[2] * OM3 * TMP2)) - cI
      * (P1[2] * TMP2))) + P3[2] * (P1[0] * TMP2 * (-cI * (TMP7 * TMP9) + cI *
      (TMP4 * TMP10)) + P2[0] * TMP1 * (-cI * (TMP7 * TMP9) + cI * (TMP4 *
      TMP10)))) + (P1[0] * P2[2] * (-cI * (TMP4 * TMP10) + cI * (TMP7 * TMP9))
      + P1[2] * P2[0] * (-cI * (TMP4 * TMP10) + cI * (TMP7 * TMP9))));
  T3[11] = denom * (OM3 * (P3[1] * (TMP10 * (TMP4 * (P3[2] * - 0.666666667 *
      (+cI * (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (+cI * (TMP1 * P2[2] +
      P1[2] * TMP2))) + 0.666666667 * cI * (P3[2] * TMP7 * TMP9)) + TMP7 * TMP9
      * (TMP1 * (-cI * (P2[2]) + 1.333333333 * cI * (P3[2] * OM3 * TMP2)) - cI
      * (P1[2] * TMP2))) + P3[2] * (P1[1] * TMP2 * (-cI * (TMP7 * TMP9) + cI *
      (TMP4 * TMP10)) + P2[1] * TMP1 * (-cI * (TMP7 * TMP9) + cI * (TMP4 *
      TMP10)))) + (P1[1] * P2[2] * (-cI * (TMP4 * TMP10) + cI * (TMP7 * TMP9))
      + P1[2] * P2[1] * (-cI * (TMP4 * TMP10) + cI * (TMP7 * TMP9))));
  T3[12] = denom * 2. * (OM3 * (P3[2] * (TMP10 * (TMP4 * (P3[2] * - 0.333333333
      * (+cI * (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (+cI * (TMP1 * P2[2]
      + P1[2] * TMP2))) + 0.333333333 * cI * (P3[2] * TMP7 * TMP9)) + TMP7 *
      TMP9 * (TMP1 * (-cI * (P2[2]) + 0.666666667 * cI * (P3[2] * OM3 * TMP2))
      - cI * (P1[2] * TMP2))) + 0.333333333 * (TMP1 * TMP2 * (-cI * (TMP7 *
      TMP9) + cI * (TMP4 * TMP10)))) + (TMP10 * (TMP4 * - 1. * (+cI * (P1[2] *
      P2[2]) + 0.333333333 * cI * (TMP10)) + 0.333333333 * cI * (TMP7 * TMP9))
      + cI * (P1[2] * P2[2] * TMP7 * TMP9)));
  T3[13] = denom * (OM3 * (P3[2] * (TMP10 * (TMP4 * (P3[3] * - 0.666666667 *
      (+cI * (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (+cI * (P1[3] * TMP2 +
      TMP1 * P2[3]))) + 0.666666667 * cI * (P3[3] * TMP7 * TMP9)) + TMP7 * TMP9
      * (TMP1 * (-cI * (P2[3]) + 1.333333333 * cI * (P3[3] * OM3 * TMP2)) - cI
      * (P1[3] * TMP2))) + P3[3] * (P1[2] * TMP2 * (-cI * (TMP7 * TMP9) + cI *
      (TMP4 * TMP10)) + P2[2] * TMP1 * (-cI * (TMP7 * TMP9) + cI * (TMP4 *
      TMP10)))) + (P1[2] * P2[3] * (-cI * (TMP4 * TMP10) + cI * (TMP7 * TMP9))
      + P1[3] * P2[2] * (-cI * (TMP4 * TMP10) + cI * (TMP7 * TMP9))));
  T3[14] = denom * (OM3 * (P3[0] * (TMP10 * (TMP4 * (P3[3] * - 0.666666667 *
      (+cI * (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (+cI * (TMP1 * P2[3] +
      P1[3] * TMP2))) + 0.666666667 * cI * (P3[3] * TMP7 * TMP9)) + TMP7 * TMP9
      * (TMP1 * (-cI * (P2[3]) + 1.333333333 * cI * (P3[3] * OM3 * TMP2)) - cI
      * (P1[3] * TMP2))) + P3[3] * (P1[0] * TMP2 * (-cI * (TMP7 * TMP9) + cI *
      (TMP4 * TMP10)) + P2[0] * TMP1 * (-cI * (TMP7 * TMP9) + cI * (TMP4 *
      TMP10)))) + (P1[0] * P2[3] * (-cI * (TMP4 * TMP10) + cI * (TMP7 * TMP9))
      + P1[3] * P2[0] * (-cI * (TMP4 * TMP10) + cI * (TMP7 * TMP9))));
  T3[15] = denom * (OM3 * (P3[1] * (TMP10 * (TMP4 * (P3[3] * - 0.666666667 *
      (+cI * (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (+cI * (TMP1 * P2[3] +
      P1[3] * TMP2))) + 0.666666667 * cI * (P3[3] * TMP7 * TMP9)) + TMP7 * TMP9
      * (TMP1 * (-cI * (P2[3]) + 1.333333333 * cI * (P3[3] * OM3 * TMP2)) - cI
      * (P1[3] * TMP2))) + P3[3] * (P1[1] * TMP2 * (-cI * (TMP7 * TMP9) + cI *
      (TMP4 * TMP10)) + P2[1] * TMP1 * (-cI * (TMP7 * TMP9) + cI * (TMP4 *
      TMP10)))) + (P1[1] * P2[3] * (-cI * (TMP4 * TMP10) + cI * (TMP7 * TMP9))
      + P1[3] * P2[1] * (-cI * (TMP4 * TMP10) + cI * (TMP7 * TMP9))));
  T3[16] = denom * (OM3 * (P3[2] * (TMP10 * (TMP4 * (P3[3] * - 0.666666667 *
      (+cI * (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (+cI * (TMP1 * P2[3] +
      P1[3] * TMP2))) + 0.666666667 * cI * (P3[3] * TMP7 * TMP9)) + TMP7 * TMP9
      * (TMP1 * (-cI * (P2[3]) + 1.333333333 * cI * (P3[3] * OM3 * TMP2)) - cI
      * (P1[3] * TMP2))) + P3[3] * (P1[2] * TMP2 * (-cI * (TMP7 * TMP9) + cI *
      (TMP4 * TMP10)) + P2[2] * TMP1 * (-cI * (TMP7 * TMP9) + cI * (TMP4 *
      TMP10)))) + (P1[2] * P2[3] * (-cI * (TMP4 * TMP10) + cI * (TMP7 * TMP9))
      + P1[3] * P2[2] * (-cI * (TMP4 * TMP10) + cI * (TMP7 * TMP9))));
  T3[17] = denom * 2. * (OM3 * (P3[3] * (TMP10 * (TMP4 * (P3[3] * - 0.333333333
      * (+cI * (TMP10) + 2. * cI * (OM3 * TMP1 * TMP2)) + (+cI * (TMP1 * P2[3]
      + P1[3] * TMP2))) + 0.333333333 * cI * (P3[3] * TMP7 * TMP9)) + TMP7 *
      TMP9 * (TMP1 * (-cI * (P2[3]) + 0.666666667 * cI * (P3[3] * OM3 * TMP2))
      - cI * (P1[3] * TMP2))) + 0.333333333 * (TMP1 * TMP2 * (-cI * (TMP7 *
      TMP9) + cI * (TMP4 * TMP10)))) + (TMP10 * (TMP4 * - 1. * (+cI * (P1[3] *
      P2[3]) + 0.333333333 * cI * (TMP10)) + 0.333333333 * cI * (TMP7 * TMP9))
      + cI * (P1[3] * P2[3] * TMP7 * TMP9)));
}


void FFV5_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  complex<double> TMP17; 
  double P3[4]; 
  double OM3; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP17 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
      F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * - cI * (F2[4] * F1[2] + F2[5] * F1[3] - P3[0] * OM3 * TMP17); 
  V3[3] = denom * - cI * (-F2[5] * F1[2] - F2[4] * F1[3] - P3[1] * OM3 *
      TMP17);
  V3[4] = denom * - cI * (-cI * (F2[5] * F1[2]) + cI * (F2[4] * F1[3]) - P3[2]
      * OM3 * TMP17);
  V3[5] = denom * - cI * (F2[5] * F1[3] - F2[4] * F1[2] - P3[3] * OM3 * TMP17); 
}

void FFV5_7_3(complex<double> F1[], complex<double> F2[], complex<double>
    COUP1, complex<double> COUP2, double M3, double W3, complex<double> V3[])
{
//   complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
//   double P3[4]; 
//   double OM3; 
  int i; 
  complex<double> Vtmp[6]; 
  FFV5_3(F1, F2, COUP1, M3, W3, V3); 
  FFV7_3(F1, F2, COUP2, M3, W3, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V3[i] = V3[i] + Vtmp[i]; 
    i++; 
  }
}

void FFT2_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP2; 
  complex<double> TMP37; 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  complex<double> TMP20; 
  double OM3; 
  complex<double> denom; 
  complex<double> TMP1; 
  complex<double> TMP38; 
  P1[0] = F1[0].real(); 
  P1[1] = F1[1].real(); 
  P1[2] = F1[1].imag(); 
  P1[3] = F1[0].imag(); 
  P2[0] = F2[0].real(); 
  P2[1] = F2[1].real(); 
  P2[2] = F2[1].imag(); 
  P2[3] = F2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +F1[0] + F2[0]; 
  T3[1] = +F1[1] + F2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP37 = (F1[2] * (F2[4] * (P1[0] + P1[3]) + F2[5] * (P1[1] + cI * (P1[2]))) +
      (F1[3] * (F2[4] * (P1[1] - cI * (P1[2])) + F2[5] * (P1[0] - P1[3])) +
      (F1[4] * (F2[2] * (P1[0] - P1[3]) - F2[3] * (P1[1] + cI * (P1[2]))) +
      F1[5] * (F2[2] * (+cI * (P1[2]) - P1[1]) + F2[3] * (P1[0] + P1[3])))));
  TMP20 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
      (F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])) +
      (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) +
      F1[5] * (F2[2] * (+cI * (P3[2]) - P3[1]) + F2[3] * (P3[0] + P3[3])))));
  TMP38 = (F1[2] * (F2[4] * (P2[0] + P2[3]) + F2[5] * (P2[1] + cI * (P2[2]))) +
      (F1[3] * (F2[4] * (P2[1] - cI * (P2[2])) + F2[5] * (P2[0] - P2[3])) +
      (F1[4] * (F2[2] * (P2[0] - P2[3]) - F2[3] * (P2[1] + cI * (P2[2]))) +
      F1[5] * (F2[2] * (+cI * (P2[2]) - P2[1]) + F2[3] * (P2[0] + P2[3])))));
  TMP1 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP2 = (P3[0] * P2[0] - P3[1] * P2[1] - P3[2] * P2[2] - P3[3] * P2[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * 2. * cI * (OM3 * (P3[0] * (TMP1 * - 1. * (F2[4] * F1[2] +
      F2[5] * F1[3] + F2[2] * F1[4] + F2[3] * F1[5] - 0.666666667 * (P3[0] *
      OM3 * TMP20)) + (TMP2 * (F2[4] * F1[2] + F2[5] * F1[3] + F2[2] * F1[4] +
      F2[3] * F1[5] - 0.666666667 * (P3[0] * OM3 * TMP20)) + (P3[0] *
      0.333333333 * (TMP37 - TMP38) + TMP20 * (P2[0] - P1[0])))) + 0.333333333
      * (TMP20 * (TMP1 - TMP2))) + (P1[0] * (F2[4] * F1[2] + F2[5] * F1[3] +
      F2[2] * F1[4] + F2[3] * F1[5]) + (P2[0] * - 1. * (F2[4] * F1[2] + F2[5] *
      F1[3] + F2[2] * F1[4] + F2[3] * F1[5]) + (-0.333333333 * (TMP37) +
      0.333333333 * (TMP38)))));
  T3[6] = denom * cI * (OM3 * (P3[0] * (TMP1 * (F2[5] * F1[2] + F2[4] * F1[3] +
      1.333333333 * (P3[1] * OM3 * TMP20) - F2[3] * F1[4] - F2[2] * F1[5]) +
      (TMP2 * - 1. * (F2[5] * F1[2] + F2[4] * F1[3] + 1.333333333 * (P3[1] *
      OM3 * TMP20) - F2[3] * F1[4] - F2[2] * F1[5]) + (P3[1] * 0.666666667 *
      (TMP37 - TMP38) + TMP20 * (P2[1] - P1[1])))) + P3[1] * (TMP1 * - 1. *
      (F2[4] * F1[2] + F2[5] * F1[3] + F2[2] * F1[4] + F2[3] * F1[5]) + (TMP2 *
      (F2[4] * F1[2] + F2[5] * F1[3] + F2[2] * F1[4] + F2[3] * F1[5]) + TMP20 *
      (P2[0] - P1[0])))) + (F1[2] * (F2[4] * (P1[1] - P2[1]) + F2[5] * (P2[0] -
      P1[0])) + (F1[3] * (F2[4] * (P2[0] - P1[0]) + F2[5] * (P1[1] - P2[1])) +
      (F1[4] * (F2[2] * (P1[1] - P2[1]) + F2[3] * (P1[0] - P2[0])) + F1[5] *
      (F2[2] * (P1[0] - P2[0]) + F2[3] * (P1[1] - P2[1]))))));
  T3[10] = denom * cI * (OM3 * (P3[0] * (TMP1 * (-cI * (F2[4] * F1[3] + F2[3] *
      F1[4]) + cI * (F2[5] * F1[2] + F2[2] * F1[5]) + 1.333333333 * (P3[2] *
      OM3 * TMP20)) + (TMP2 * - 1. * (-cI * (F2[4] * F1[3] + F2[3] * F1[4]) +
      cI * (F2[5] * F1[2] + F2[2] * F1[5]) + 1.333333333 * (P3[2] * OM3 *
      TMP20)) + (P3[2] * 0.666666667 * (TMP37 - TMP38) + TMP20 * (P2[2] -
      P1[2])))) + P3[2] * (TMP1 * - 1. * (F2[4] * F1[2] + F2[5] * F1[3] + F2[2]
      * F1[4] + F2[3] * F1[5]) + (TMP2 * (F2[4] * F1[2] + F2[5] * F1[3] + F2[2]
      * F1[4] + F2[3] * F1[5]) + TMP20 * (P2[0] - P1[0])))) + (F1[2] * (F2[4] *
      (P1[2] - P2[2]) + F2[5] * (-cI * (P1[0]) + cI * (P2[0]))) + (F1[3] *
      (F2[4] * (-cI * (P2[0]) + cI * (P1[0])) + F2[5] * (P1[2] - P2[2])) +
      (F1[4] * (F2[2] * (P1[2] - P2[2]) + F2[3] * (-cI * (P2[0]) + cI *
      (P1[0]))) + F1[5] * (F2[2] * (-cI * (P1[0]) + cI * (P2[0])) + F2[3] *
      (P1[2] - P2[2]))))));
  T3[14] = denom * cI * (OM3 * (P3[0] * (TMP1 * (F2[4] * F1[2] + F2[3] * F1[5]
      + 1.333333333 * (P3[3] * OM3 * TMP20) - F2[5] * F1[3] - F2[2] * F1[4]) +
      (TMP2 * - 1. * (F2[4] * F1[2] + F2[3] * F1[5] + 1.333333333 * (P3[3] *
      OM3 * TMP20) - F2[5] * F1[3] - F2[2] * F1[4]) + (P3[3] * 0.666666667 *
      (TMP37 - TMP38) + TMP20 * (P2[3] - P1[3])))) + P3[3] * (TMP1 * - 1. *
      (F2[4] * F1[2] + F2[5] * F1[3] + F2[2] * F1[4] + F2[3] * F1[5]) + (TMP2 *
      (F2[4] * F1[2] + F2[5] * F1[3] + F2[2] * F1[4] + F2[3] * F1[5]) + TMP20 *
      (P2[0] - P1[0])))) + (F1[2] * F2[4] * (P1[3] + P2[0] - P1[0] - P2[3]) +
      (F1[3] * F2[5] * (P1[0] + P1[3] - P2[0] - P2[3]) + (F1[4] * F2[2] *
      (P1[0] + P1[3] - P2[0] - P2[3]) + F1[5] * F2[3] * (P1[3] + P2[0] - P1[0]
      - P2[3])))));
  T3[3] = denom * cI * (OM3 * (P3[0] * (TMP1 * (F2[5] * F1[2] + F2[4] * F1[3] +
      1.333333333 * (P3[1] * OM3 * TMP20) - F2[3] * F1[4] - F2[2] * F1[5]) +
      (TMP2 * - 1. * (F2[5] * F1[2] + F2[4] * F1[3] + 1.333333333 * (P3[1] *
      OM3 * TMP20) - F2[3] * F1[4] - F2[2] * F1[5]) + (P3[1] * 0.666666667 *
      (TMP37 - TMP38) + TMP20 * (P2[1] - P1[1])))) + P3[1] * (TMP1 * - 1. *
      (F2[4] * F1[2] + F2[5] * F1[3] + F2[2] * F1[4] + F2[3] * F1[5]) + (TMP2 *
      (F2[4] * F1[2] + F2[5] * F1[3] + F2[2] * F1[4] + F2[3] * F1[5]) + TMP20 *
      (P2[0] - P1[0])))) + (F1[2] * (F2[4] * (P1[1] - P2[1]) + F2[5] * (P2[0] -
      P1[0])) + (F1[3] * (F2[4] * (P2[0] - P1[0]) + F2[5] * (P1[1] - P2[1])) +
      (F1[4] * (F2[2] * (P1[1] - P2[1]) + F2[3] * (P1[0] - P2[0])) + F1[5] *
      (F2[2] * (P1[0] - P2[0]) + F2[3] * (P1[1] - P2[1]))))));
  T3[7] = denom * 2. * cI * (OM3 * (P3[1] * (TMP1 * (F2[5] * F1[2] + F2[4] *
      F1[3] + 0.666666667 * (P3[1] * OM3 * TMP20) - F2[3] * F1[4] - F2[2] *
      F1[5]) + (TMP2 * - 1. * (F2[5] * F1[2] + F2[4] * F1[3] + 0.666666667 *
      (P3[1] * OM3 * TMP20) - F2[3] * F1[4] - F2[2] * F1[5]) + (P3[1] *
      0.333333333 * (TMP37 - TMP38) + TMP20 * (P2[1] - P1[1])))) + 0.333333333
      * (TMP20 * (TMP2 - TMP1))) + (P1[1] * (F2[3] * F1[4] + F2[2] * F1[5] -
      F2[5] * F1[2] - F2[4] * F1[3]) + (P2[1] * (F2[5] * F1[2] + F2[4] * F1[3]
      - F2[3] * F1[4] - F2[2] * F1[5]) + (-0.333333333 * (TMP38) + 0.333333333
      * (TMP37)))));
  T3[11] = denom * cI * (OM3 * (P3[1] * (TMP1 * (-cI * (F2[4] * F1[3] + F2[3] *
      F1[4]) + cI * (F2[5] * F1[2] + F2[2] * F1[5]) + 1.333333333 * (P3[2] *
      OM3 * TMP20)) + (TMP2 * - 1. * (-cI * (F2[4] * F1[3] + F2[3] * F1[4]) +
      cI * (F2[5] * F1[2] + F2[2] * F1[5]) + 1.333333333 * (P3[2] * OM3 *
      TMP20)) + (P3[2] * 0.666666667 * (TMP37 - TMP38) + TMP20 * (P2[2] -
      P1[2])))) + P3[2] * (TMP1 * (F2[5] * F1[2] + F2[4] * F1[3] - F2[3] *
      F1[4] - F2[2] * F1[5]) + (TMP2 * (F2[3] * F1[4] + F2[2] * F1[5] - F2[5] *
      F1[2] - F2[4] * F1[3]) + TMP20 * (P2[1] - P1[1])))) + (F1[2] * F2[5] *
      (P2[2] - cI * (P1[1]) + cI * (P2[1]) - P1[2]) + (F1[3] * F2[4] * (P2[2] -
      cI * (P2[1]) + cI * (P1[1]) - P1[2]) + (F1[4] * F2[3] * (P1[2] - cI *
      (P2[1]) + cI * (P1[1]) - P2[2]) + F1[5] * F2[2] * (P1[2] - cI * (P1[1]) +
      cI * (P2[1]) - P2[2])))));
  T3[15] = denom * cI * (OM3 * (P3[1] * (TMP1 * (F2[4] * F1[2] + F2[3] * F1[5]
      + 1.333333333 * (P3[3] * OM3 * TMP20) - F2[5] * F1[3] - F2[2] * F1[4]) +
      (TMP2 * - 1. * (F2[4] * F1[2] + F2[3] * F1[5] + 1.333333333 * (P3[3] *
      OM3 * TMP20) - F2[5] * F1[3] - F2[2] * F1[4]) + (P3[3] * 0.666666667 *
      (TMP37 - TMP38) + TMP20 * (P2[3] - P1[3])))) + P3[3] * (TMP1 * (F2[5] *
      F1[2] + F2[4] * F1[3] - F2[3] * F1[4] - F2[2] * F1[5]) + (TMP2 * (F2[3] *
      F1[4] + F2[2] * F1[5] - F2[5] * F1[2] - F2[4] * F1[3]) + TMP20 * (P2[1] -
      P1[1])))) + (F1[2] * (F2[4] * (P2[1] - P1[1]) + F2[5] * (P2[3] - P1[3]))
      + (F1[3] * (F2[4] * (P2[3] - P1[3]) + F2[5] * (P1[1] - P2[1])) + (F1[4] *
      (F2[2] * (P1[1] - P2[1]) + F2[3] * (P1[3] - P2[3])) + F1[5] * (F2[2] *
      (P1[3] - P2[3]) + F2[3] * (P2[1] - P1[1]))))));
  T3[4] = denom * cI * (OM3 * (P3[0] * (TMP1 * (-cI * (F2[4] * F1[3] + F2[3] *
      F1[4]) + cI * (F2[5] * F1[2] + F2[2] * F1[5]) + 1.333333333 * (P3[2] *
      OM3 * TMP20)) + (TMP2 * - 1. * (-cI * (F2[4] * F1[3] + F2[3] * F1[4]) +
      cI * (F2[5] * F1[2] + F2[2] * F1[5]) + 1.333333333 * (P3[2] * OM3 *
      TMP20)) + (P3[2] * 0.666666667 * (TMP37 - TMP38) + TMP20 * (P2[2] -
      P1[2])))) + P3[2] * (TMP1 * - 1. * (F2[4] * F1[2] + F2[5] * F1[3] + F2[2]
      * F1[4] + F2[3] * F1[5]) + (TMP2 * (F2[4] * F1[2] + F2[5] * F1[3] + F2[2]
      * F1[4] + F2[3] * F1[5]) + TMP20 * (P2[0] - P1[0])))) + (F1[2] * (F2[4] *
      (P1[2] - P2[2]) + F2[5] * (-cI * (P1[0]) + cI * (P2[0]))) + (F1[3] *
      (F2[4] * (-cI * (P2[0]) + cI * (P1[0])) + F2[5] * (P1[2] - P2[2])) +
      (F1[4] * (F2[2] * (P1[2] - P2[2]) + F2[3] * (-cI * (P2[0]) + cI *
      (P1[0]))) + F1[5] * (F2[2] * (-cI * (P1[0]) + cI * (P2[0])) + F2[3] *
      (P1[2] - P2[2]))))));
  T3[8] = denom * cI * (OM3 * (P3[1] * (TMP1 * (-cI * (F2[4] * F1[3] + F2[3] *
      F1[4]) + cI * (F2[5] * F1[2] + F2[2] * F1[5]) + 1.333333333 * (P3[2] *
      OM3 * TMP20)) + (TMP2 * - 1. * (-cI * (F2[4] * F1[3] + F2[3] * F1[4]) +
      cI * (F2[5] * F1[2] + F2[2] * F1[5]) + 1.333333333 * (P3[2] * OM3 *
      TMP20)) + (P3[2] * 0.666666667 * (TMP37 - TMP38) + TMP20 * (P2[2] -
      P1[2])))) + P3[2] * (TMP1 * (F2[5] * F1[2] + F2[4] * F1[3] - F2[3] *
      F1[4] - F2[2] * F1[5]) + (TMP2 * (F2[3] * F1[4] + F2[2] * F1[5] - F2[5] *
      F1[2] - F2[4] * F1[3]) + TMP20 * (P2[1] - P1[1])))) + (F1[2] * F2[5] *
      (P2[2] - cI * (P1[1]) + cI * (P2[1]) - P1[2]) + (F1[3] * F2[4] * (P2[2] -
      cI * (P2[1]) + cI * (P1[1]) - P1[2]) + (F1[4] * F2[3] * (P1[2] - cI *
      (P2[1]) + cI * (P1[1]) - P2[2]) + F1[5] * F2[2] * (P1[2] - cI * (P1[1]) +
      cI * (P2[1]) - P2[2])))));
  T3[12] = denom * 2. * cI * (OM3 * (P3[2] * (TMP1 * (-cI * (F2[4] * F1[3] +
      F2[3] * F1[4]) + cI * (F2[5] * F1[2] + F2[2] * F1[5]) + 0.666666667 *
      (P3[2] * OM3 * TMP20)) + (TMP2 * - 1. * (-cI * (F2[4] * F1[3] + F2[3] *
      F1[4]) + cI * (F2[5] * F1[2] + F2[2] * F1[5]) + 0.666666667 * (P3[2] *
      OM3 * TMP20)) + (P3[2] * 0.333333333 * (TMP37 - TMP38) + TMP20 * (P2[2] -
      P1[2])))) + 0.333333333 * (TMP20 * (TMP2 - TMP1))) + (P1[2] * (-cI *
      (F2[5] * F1[2] + F2[2] * F1[5]) + cI * (F2[4] * F1[3] + F2[3] * F1[4])) +
      (P2[2] * (-cI * (F2[4] * F1[3] + F2[3] * F1[4]) + cI * (F2[5] * F1[2] +
      F2[2] * F1[5])) + (-0.333333333 * (TMP38) + 0.333333333 * (TMP37)))));
  T3[16] = denom * cI * (OM3 * (P3[2] * (TMP1 * (F2[4] * F1[2] + F2[3] * F1[5]
      + 1.333333333 * (P3[3] * OM3 * TMP20) - F2[5] * F1[3] - F2[2] * F1[4]) +
      (TMP2 * - 1. * (F2[4] * F1[2] + F2[3] * F1[5] + 1.333333333 * (P3[3] *
      OM3 * TMP20) - F2[5] * F1[3] - F2[2] * F1[4]) + (P3[3] * 0.666666667 *
      (TMP37 - TMP38) + TMP20 * (P2[3] - P1[3])))) + P3[3] * (TMP1 * (-cI *
      (F2[4] * F1[3] + F2[3] * F1[4]) + cI * (F2[5] * F1[2] + F2[2] * F1[5])) +
      (TMP2 * (-cI * (F2[5] * F1[2] + F2[2] * F1[5]) + cI * (F2[4] * F1[3] +
      F2[3] * F1[4])) + TMP20 * (P2[2] - P1[2])))) + (F1[2] * (F2[4] * (P2[2] -
      P1[2]) + F2[5] * (-cI * (P1[3]) + cI * (P2[3]))) + (F1[3] * (F2[4] * (-cI
      * (P2[3]) + cI * (P1[3])) + F2[5] * (P1[2] - P2[2])) + (F1[4] * (F2[2] *
      (P1[2] - P2[2]) + F2[3] * (-cI * (P2[3]) + cI * (P1[3]))) + F1[5] *
      (F2[2] * (-cI * (P1[3]) + cI * (P2[3])) + F2[3] * (P2[2] - P1[2]))))));
  T3[5] = denom * cI * (OM3 * (P3[0] * (TMP1 * (F2[4] * F1[2] + F2[3] * F1[5] +
      1.333333333 * (P3[3] * OM3 * TMP20) - F2[5] * F1[3] - F2[2] * F1[4]) +
      (TMP2 * - 1. * (F2[4] * F1[2] + F2[3] * F1[5] + 1.333333333 * (P3[3] *
      OM3 * TMP20) - F2[5] * F1[3] - F2[2] * F1[4]) + (P3[3] * 0.666666667 *
      (TMP37 - TMP38) + TMP20 * (P2[3] - P1[3])))) + P3[3] * (TMP1 * - 1. *
      (F2[4] * F1[2] + F2[5] * F1[3] + F2[2] * F1[4] + F2[3] * F1[5]) + (TMP2 *
      (F2[4] * F1[2] + F2[5] * F1[3] + F2[2] * F1[4] + F2[3] * F1[5]) + TMP20 *
      (P2[0] - P1[0])))) + (F1[2] * F2[4] * (P1[3] + P2[0] - P1[0] - P2[3]) +
      (F1[3] * F2[5] * (P1[3] + P1[0] - P2[3] - P2[0]) + (F1[4] * F2[2] *
      (P1[3] + P1[0] - P2[3] - P2[0]) + F1[5] * F2[3] * (P1[3] + P2[0] - P1[0]
      - P2[3])))));
  T3[9] = denom * cI * (OM3 * (P3[1] * (TMP1 * (F2[4] * F1[2] + F2[3] * F1[5] +
      1.333333333 * (P3[3] * OM3 * TMP20) - F2[5] * F1[3] - F2[2] * F1[4]) +
      (TMP2 * - 1. * (F2[4] * F1[2] + F2[3] * F1[5] + 1.333333333 * (P3[3] *
      OM3 * TMP20) - F2[5] * F1[3] - F2[2] * F1[4]) + (P3[3] * 0.666666667 *
      (TMP37 - TMP38) + TMP20 * (P2[3] - P1[3])))) + P3[3] * (TMP1 * (F2[5] *
      F1[2] + F2[4] * F1[3] - F2[3] * F1[4] - F2[2] * F1[5]) + (TMP2 * (F2[3] *
      F1[4] + F2[2] * F1[5] - F2[5] * F1[2] - F2[4] * F1[3]) + TMP20 * (P2[1] -
      P1[1])))) + (F1[2] * (F2[4] * (P2[1] - P1[1]) + F2[5] * (P2[3] - P1[3]))
      + (F1[3] * (F2[4] * (P2[3] - P1[3]) + F2[5] * (P1[1] - P2[1])) + (F1[4] *
      (F2[2] * (P1[1] - P2[1]) + F2[3] * (P1[3] - P2[3])) + F1[5] * (F2[2] *
      (P1[3] - P2[3]) + F2[3] * (P2[1] - P1[1]))))));
  T3[13] = denom * cI * (OM3 * (P3[2] * (TMP1 * (F2[4] * F1[2] + F2[3] * F1[5]
      + 1.333333333 * (P3[3] * OM3 * TMP20) - F2[5] * F1[3] - F2[2] * F1[4]) +
      (TMP2 * - 1. * (F2[4] * F1[2] + F2[3] * F1[5] + 1.333333333 * (P3[3] *
      OM3 * TMP20) - F2[5] * F1[3] - F2[2] * F1[4]) + (P3[3] * 0.666666667 *
      (TMP37 - TMP38) + TMP20 * (P2[3] - P1[3])))) + P3[3] * (TMP1 * (-cI *
      (F2[4] * F1[3] + F2[3] * F1[4]) + cI * (F2[5] * F1[2] + F2[2] * F1[5])) +
      (TMP2 * (-cI * (F2[5] * F1[2] + F2[2] * F1[5]) + cI * (F2[4] * F1[3] +
      F2[3] * F1[4])) + TMP20 * (P2[2] - P1[2])))) + (F1[2] * (F2[4] * (P2[2] -
      P1[2]) + F2[5] * (-cI * (P1[3]) + cI * (P2[3]))) + (F1[3] * (F2[4] * (-cI
      * (P2[3]) + cI * (P1[3])) + F2[5] * (P1[2] - P2[2])) + (F1[4] * (F2[2] *
      (P1[2] - P2[2]) + F2[3] * (-cI * (P2[3]) + cI * (P1[3]))) + F1[5] *
      (F2[2] * (-cI * (P1[3]) + cI * (P2[3])) + F2[3] * (P2[2] - P1[2]))))));
  T3[17] = denom * 2. * cI * (OM3 * (P3[3] * (TMP1 * (F2[4] * F1[2] + F2[3] *
      F1[5] + 0.666666667 * (P3[3] * OM3 * TMP20) - F2[5] * F1[3] - F2[2] *
      F1[4]) + (TMP2 * - 1. * (F2[4] * F1[2] + F2[3] * F1[5] + 0.666666667 *
      (P3[3] * OM3 * TMP20) - F2[5] * F1[3] - F2[2] * F1[4]) + (P3[3] *
      0.333333333 * (TMP37 - TMP38) + TMP20 * (P2[3] - P1[3])))) + 0.333333333
      * (TMP20 * (TMP2 - TMP1))) + (P1[3] * (F2[5] * F1[3] + F2[2] * F1[4] -
      F2[4] * F1[2] - F2[3] * F1[5]) + (P2[3] * (F2[4] * F1[2] + F2[3] * F1[5]
      - F2[5] * F1[3] - F2[2] * F1[4]) + (-0.333333333 * (TMP38) + 0.333333333
      * (TMP37)))));
}


void VVT6_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP59; 
  complex<double> TMP58; 
  double P1[4]; 
  complex<double> TMP57; 
  double P2[4]; 
  complex<double> TMP23; 
  double P3[4]; 
  complex<double> TMP63; 
  complex<double> TMP60; 
  complex<double> TMP26; 
  complex<double> TMP61; 
  complex<double> TMP64; 
  complex<double> TMP62; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = T3[0].real(); 
  P3[1] = T3[1].real(); 
  P3[2] = T3[1].imag(); 
  P3[3] = T3[0].imag(); 
  TMP26 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP57 = -1. * (P1[0] * (P2[0] * (P3[1] * (V2[4] * T3[5] - V2[5] * T3[4]) +
      (P3[2] * (V2[5] * T3[3] - V2[3] * T3[5]) + P3[3] * (V2[3] * T3[4] - V2[4]
      * T3[3]))) + (P2[1] * (P3[0] * (V2[5] * T3[4] - V2[4] * T3[5]) + (P3[2] *
      (V2[2] * T3[5] - V2[5] * T3[2]) + P3[3] * (V2[4] * T3[2] - V2[2] *
      T3[4]))) + (P2[2] * (P3[0] * (V2[3] * T3[5] - V2[5] * T3[3]) + (P3[1] *
      (V2[5] * T3[2] - V2[2] * T3[5]) + P3[3] * (V2[2] * T3[3] - V2[3] *
      T3[2]))) + P2[3] * (P3[0] * (V2[4] * T3[3] - V2[3] * T3[4]) + (P3[1] *
      (V2[2] * T3[4] - V2[4] * T3[2]) + P3[2] * (V2[3] * T3[2] - V2[2] *
      T3[3])))))) + (P1[1] * (P2[0] * (P3[1] * (V2[5] * T3[8] - V2[4] * T3[9])
      + (P3[2] * (V2[3] * T3[9] - V2[5] * T3[7]) + P3[3] * (V2[4] * T3[7] -
      V2[3] * T3[8]))) + (P2[1] * (P3[0] * (V2[4] * T3[9] - V2[5] * T3[8]) +
      (P3[2] * (V2[5] * T3[6] - V2[2] * T3[9]) + P3[3] * (V2[2] * T3[8] - V2[4]
      * T3[6]))) + (P2[2] * (P3[0] * (V2[5] * T3[7] - V2[3] * T3[9]) + (P3[1] *
      (V2[2] * T3[9] - V2[5] * T3[6]) + P3[3] * (V2[3] * T3[6] - V2[2] *
      T3[7]))) + P2[3] * (P3[0] * (V2[3] * T3[8] - V2[4] * T3[7]) + (P3[1] *
      (V2[4] * T3[6] - V2[2] * T3[8]) + P3[2] * (V2[2] * T3[7] - V2[3] *
      T3[6])))))) + (P1[2] * (P2[0] * (P3[1] * (V2[5] * T3[12] - V2[4] *
      T3[13]) + (P3[2] * (V2[3] * T3[13] - V2[5] * T3[11]) + P3[3] * (V2[4] *
      T3[11] - V2[3] * T3[12]))) + (P2[1] * (P3[0] * (V2[4] * T3[13] - V2[5] *
      T3[12]) + (P3[2] * (V2[5] * T3[10] - V2[2] * T3[13]) + P3[3] * (V2[2] *
      T3[12] - V2[4] * T3[10]))) + (P2[2] * (P3[0] * (V2[5] * T3[11] - V2[3] *
      T3[13]) + (P3[1] * (V2[2] * T3[13] - V2[5] * T3[10]) + P3[3] * (V2[3] *
      T3[10] - V2[2] * T3[11]))) + P2[3] * (P3[0] * (V2[3] * T3[12] - V2[4] *
      T3[11]) + (P3[1] * (V2[4] * T3[10] - V2[2] * T3[12]) + P3[2] * (V2[2] *
      T3[11] - V2[3] * T3[10])))))) + P1[3] * (P2[0] * (P3[1] * (V2[5] * T3[16]
      - V2[4] * T3[17]) + (P3[2] * (V2[3] * T3[17] - V2[5] * T3[15]) + P3[3] *
      (V2[4] * T3[15] - V2[3] * T3[16]))) + (P2[1] * (P3[0] * (V2[4] * T3[17] -
      V2[5] * T3[16]) + (P3[2] * (V2[5] * T3[14] - V2[2] * T3[17]) + P3[3] *
      (V2[2] * T3[16] - V2[4] * T3[14]))) + (P2[2] * (P3[0] * (V2[5] * T3[15] -
      V2[3] * T3[17]) + (P3[1] * (V2[2] * T3[17] - V2[5] * T3[14]) + P3[3] *
      (V2[3] * T3[14] - V2[2] * T3[15]))) + P2[3] * (P3[0] * (V2[3] * T3[16] -
      V2[4] * T3[15]) + (P3[1] * (V2[4] * T3[14] - V2[2] * T3[16]) + P3[2] *
      (V2[2] * T3[15] - V2[3] * T3[14])))))))));
  TMP23 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP59 = -1. * (P2[0] * (P1[0] * (P3[1] * (V2[4] * T3[5] - V2[5] * T3[4]) +
      (P3[2] * (V2[5] * T3[3] - V2[3] * T3[5]) + P3[3] * (V2[3] * T3[4] - V2[4]
      * T3[3]))) + (P1[1] * (P3[0] * (V2[5] * T3[4] - V2[4] * T3[5]) + (P3[2] *
      (V2[2] * T3[5] - V2[5] * T3[2]) + P3[3] * (V2[4] * T3[2] - V2[2] *
      T3[4]))) + (P1[2] * (P3[0] * (V2[3] * T3[5] - V2[5] * T3[3]) + (P3[1] *
      (V2[5] * T3[2] - V2[2] * T3[5]) + P3[3] * (V2[2] * T3[3] - V2[3] *
      T3[2]))) + P1[3] * (P3[0] * (V2[4] * T3[3] - V2[3] * T3[4]) + (P3[1] *
      (V2[2] * T3[4] - V2[4] * T3[2]) + P3[2] * (V2[3] * T3[2] - V2[2] *
      T3[3])))))) + (P2[1] * (P1[0] * (P3[1] * (V2[5] * T3[8] - V2[4] * T3[9])
      + (P3[2] * (V2[3] * T3[9] - V2[5] * T3[7]) + P3[3] * (V2[4] * T3[7] -
      V2[3] * T3[8]))) + (P1[1] * (P3[0] * (V2[4] * T3[9] - V2[5] * T3[8]) +
      (P3[2] * (V2[5] * T3[6] - V2[2] * T3[9]) + P3[3] * (V2[2] * T3[8] - V2[4]
      * T3[6]))) + (P1[2] * (P3[0] * (V2[5] * T3[7] - V2[3] * T3[9]) + (P3[1] *
      (V2[2] * T3[9] - V2[5] * T3[6]) + P3[3] * (V2[3] * T3[6] - V2[2] *
      T3[7]))) + P1[3] * (P3[0] * (V2[3] * T3[8] - V2[4] * T3[7]) + (P3[1] *
      (V2[4] * T3[6] - V2[2] * T3[8]) + P3[2] * (V2[2] * T3[7] - V2[3] *
      T3[6])))))) + (P2[2] * (P1[0] * (P3[1] * (V2[5] * T3[12] - V2[4] *
      T3[13]) + (P3[2] * (V2[3] * T3[13] - V2[5] * T3[11]) + P3[3] * (V2[4] *
      T3[11] - V2[3] * T3[12]))) + (P1[1] * (P3[0] * (V2[4] * T3[13] - V2[5] *
      T3[12]) + (P3[2] * (V2[5] * T3[10] - V2[2] * T3[13]) + P3[3] * (V2[2] *
      T3[12] - V2[4] * T3[10]))) + (P1[2] * (P3[0] * (V2[5] * T3[11] - V2[3] *
      T3[13]) + (P3[1] * (V2[2] * T3[13] - V2[5] * T3[10]) + P3[3] * (V2[3] *
      T3[10] - V2[2] * T3[11]))) + P1[3] * (P3[0] * (V2[3] * T3[12] - V2[4] *
      T3[11]) + (P3[1] * (V2[4] * T3[10] - V2[2] * T3[12]) + P3[2] * (V2[2] *
      T3[11] - V2[3] * T3[10])))))) + P2[3] * (P1[0] * (P3[1] * (V2[5] * T3[16]
      - V2[4] * T3[17]) + (P3[2] * (V2[3] * T3[17] - V2[5] * T3[15]) + P3[3] *
      (V2[4] * T3[15] - V2[3] * T3[16]))) + (P1[1] * (P3[0] * (V2[4] * T3[17] -
      V2[5] * T3[16]) + (P3[2] * (V2[5] * T3[14] - V2[2] * T3[17]) + P3[3] *
      (V2[2] * T3[16] - V2[4] * T3[14]))) + (P1[2] * (P3[0] * (V2[5] * T3[15] -
      V2[3] * T3[17]) + (P3[1] * (V2[2] * T3[17] - V2[5] * T3[14]) + P3[3] *
      (V2[3] * T3[14] - V2[2] * T3[15]))) + P1[3] * (P3[0] * (V2[3] * T3[16] -
      V2[4] * T3[15]) + (P3[1] * (V2[4] * T3[14] - V2[2] * T3[16]) + P3[2] *
      (V2[2] * T3[15] - V2[3] * T3[14])))))))));
  TMP58 = -1. * (P1[0] * (P2[0] * (P3[1] * (V1[4] * T3[5] - V1[5] * T3[4]) +
      (P3[2] * (V1[5] * T3[3] - V1[3] * T3[5]) + P3[3] * (V1[3] * T3[4] - V1[4]
      * T3[3]))) + (P2[1] * (P3[0] * (V1[5] * T3[4] - V1[4] * T3[5]) + (P3[2] *
      (V1[2] * T3[5] - V1[5] * T3[2]) + P3[3] * (V1[4] * T3[2] - V1[2] *
      T3[4]))) + (P2[2] * (P3[0] * (V1[3] * T3[5] - V1[5] * T3[3]) + (P3[1] *
      (V1[5] * T3[2] - V1[2] * T3[5]) + P3[3] * (V1[2] * T3[3] - V1[3] *
      T3[2]))) + P2[3] * (P3[0] * (V1[4] * T3[3] - V1[3] * T3[4]) + (P3[1] *
      (V1[2] * T3[4] - V1[4] * T3[2]) + P3[2] * (V1[3] * T3[2] - V1[2] *
      T3[3])))))) + (P1[1] * (P2[0] * (P3[1] * (V1[5] * T3[8] - V1[4] * T3[9])
      + (P3[2] * (V1[3] * T3[9] - V1[5] * T3[7]) + P3[3] * (V1[4] * T3[7] -
      V1[3] * T3[8]))) + (P2[1] * (P3[0] * (V1[4] * T3[9] - V1[5] * T3[8]) +
      (P3[2] * (V1[5] * T3[6] - V1[2] * T3[9]) + P3[3] * (V1[2] * T3[8] - V1[4]
      * T3[6]))) + (P2[2] * (P3[0] * (V1[5] * T3[7] - V1[3] * T3[9]) + (P3[1] *
      (V1[2] * T3[9] - V1[5] * T3[6]) + P3[3] * (V1[3] * T3[6] - V1[2] *
      T3[7]))) + P2[3] * (P3[0] * (V1[3] * T3[8] - V1[4] * T3[7]) + (P3[1] *
      (V1[4] * T3[6] - V1[2] * T3[8]) + P3[2] * (V1[2] * T3[7] - V1[3] *
      T3[6])))))) + (P1[2] * (P2[0] * (P3[1] * (V1[5] * T3[12] - V1[4] *
      T3[13]) + (P3[2] * (V1[3] * T3[13] - V1[5] * T3[11]) + P3[3] * (V1[4] *
      T3[11] - V1[3] * T3[12]))) + (P2[1] * (P3[0] * (V1[4] * T3[13] - V1[5] *
      T3[12]) + (P3[2] * (V1[5] * T3[10] - V1[2] * T3[13]) + P3[3] * (V1[2] *
      T3[12] - V1[4] * T3[10]))) + (P2[2] * (P3[0] * (V1[5] * T3[11] - V1[3] *
      T3[13]) + (P3[1] * (V1[2] * T3[13] - V1[5] * T3[10]) + P3[3] * (V1[3] *
      T3[10] - V1[2] * T3[11]))) + P2[3] * (P3[0] * (V1[3] * T3[12] - V1[4] *
      T3[11]) + (P3[1] * (V1[4] * T3[10] - V1[2] * T3[12]) + P3[2] * (V1[2] *
      T3[11] - V1[3] * T3[10])))))) + P1[3] * (P2[0] * (P3[1] * (V1[5] * T3[16]
      - V1[4] * T3[17]) + (P3[2] * (V1[3] * T3[17] - V1[5] * T3[15]) + P3[3] *
      (V1[4] * T3[15] - V1[3] * T3[16]))) + (P2[1] * (P3[0] * (V1[4] * T3[17] -
      V1[5] * T3[16]) + (P3[2] * (V1[5] * T3[14] - V1[2] * T3[17]) + P3[3] *
      (V1[2] * T3[16] - V1[4] * T3[14]))) + (P2[2] * (P3[0] * (V1[5] * T3[15] -
      V1[3] * T3[17]) + (P3[1] * (V1[2] * T3[17] - V1[5] * T3[14]) + P3[3] *
      (V1[3] * T3[14] - V1[2] * T3[15]))) + P2[3] * (P3[0] * (V1[3] * T3[16] -
      V1[4] * T3[15]) + (P3[1] * (V1[4] * T3[14] - V1[2] * T3[16]) + P3[2] *
      (V1[2] * T3[15] - V1[3] * T3[14])))))))));
  TMP62 = -1. * (P1[0] * (P2[0] * (P3[1] * (V1[4] * T3[14] - V1[5] * T3[10]) +
      (P3[2] * (V1[5] * T3[6] - V1[3] * T3[14]) + P3[3] * (V1[3] * T3[10] -
      V1[4] * T3[6]))) + (P2[1] * (P3[0] * (V1[5] * T3[10] - V1[4] * T3[14]) +
      (P3[2] * (V1[2] * T3[14] - V1[5] * T3[2]) + P3[3] * (V1[4] * T3[2] -
      V1[2] * T3[10]))) + (P2[2] * (P3[0] * (V1[3] * T3[14] - V1[5] * T3[6]) +
      (P3[1] * (V1[5] * T3[2] - V1[2] * T3[14]) + P3[3] * (V1[2] * T3[6] -
      V1[3] * T3[2]))) + P2[3] * (P3[0] * (V1[4] * T3[6] - V1[3] * T3[10]) +
      (P3[1] * (V1[2] * T3[10] - V1[4] * T3[2]) + P3[2] * (V1[3] * T3[2] -
      V1[2] * T3[6])))))) + (P1[1] * (P2[0] * (P3[1] * (V1[5] * T3[11] - V1[4]
      * T3[15]) + (P3[2] * (V1[3] * T3[15] - V1[5] * T3[7]) + P3[3] * (V1[4] *
      T3[7] - V1[3] * T3[11]))) + (P2[1] * (P3[0] * (V1[4] * T3[15] - V1[5] *
      T3[11]) + (P3[2] * (V1[5] * T3[3] - V1[2] * T3[15]) + P3[3] * (V1[2] *
      T3[11] - V1[4] * T3[3]))) + (P2[2] * (P3[0] * (V1[5] * T3[7] - V1[3] *
      T3[15]) + (P3[1] * (V1[2] * T3[15] - V1[5] * T3[3]) + P3[3] * (V1[3] *
      T3[3] - V1[2] * T3[7]))) + P2[3] * (P3[0] * (V1[3] * T3[11] - V1[4] *
      T3[7]) + (P3[1] * (V1[4] * T3[3] - V1[2] * T3[11]) + P3[2] * (V1[2] *
      T3[7] - V1[3] * T3[3])))))) + (P1[2] * (P2[0] * (P3[1] * (V1[5] * T3[12]
      - V1[4] * T3[16]) + (P3[2] * (V1[3] * T3[16] - V1[5] * T3[8]) + P3[3] *
      (V1[4] * T3[8] - V1[3] * T3[12]))) + (P2[1] * (P3[0] * (V1[4] * T3[16] -
      V1[5] * T3[12]) + (P3[2] * (V1[5] * T3[4] - V1[2] * T3[16]) + P3[3] *
      (V1[2] * T3[12] - V1[4] * T3[4]))) + (P2[2] * (P3[0] * (V1[5] * T3[8] -
      V1[3] * T3[16]) + (P3[1] * (V1[2] * T3[16] - V1[5] * T3[4]) + P3[3] *
      (V1[3] * T3[4] - V1[2] * T3[8]))) + P2[3] * (P3[0] * (V1[3] * T3[12] -
      V1[4] * T3[8]) + (P3[1] * (V1[4] * T3[4] - V1[2] * T3[12]) + P3[2] *
      (V1[2] * T3[8] - V1[3] * T3[4])))))) + P1[3] * (P2[0] * (P3[1] * (V1[5] *
      T3[13] - V1[4] * T3[17]) + (P3[2] * (V1[3] * T3[17] - V1[5] * T3[9]) +
      P3[3] * (V1[4] * T3[9] - V1[3] * T3[13]))) + (P2[1] * (P3[0] * (V1[4] *
      T3[17] - V1[5] * T3[13]) + (P3[2] * (V1[5] * T3[5] - V1[2] * T3[17]) +
      P3[3] * (V1[2] * T3[13] - V1[4] * T3[5]))) + (P2[2] * (P3[0] * (V1[5] *
      T3[9] - V1[3] * T3[17]) + (P3[1] * (V1[2] * T3[17] - V1[5] * T3[5]) +
      P3[3] * (V1[3] * T3[5] - V1[2] * T3[9]))) + P2[3] * (P3[0] * (V1[3] *
      T3[13] - V1[4] * T3[9]) + (P3[1] * (V1[4] * T3[5] - V1[2] * T3[13]) +
      P3[2] * (V1[2] * T3[9] - V1[3] * T3[5])))))))));
  TMP63 = -1. * (P2[0] * (P1[0] * (P3[1] * (V2[4] * T3[14] - V2[5] * T3[10]) +
      (P3[2] * (V2[5] * T3[6] - V2[3] * T3[14]) + P3[3] * (V2[3] * T3[10] -
      V2[4] * T3[6]))) + (P1[1] * (P3[0] * (V2[5] * T3[10] - V2[4] * T3[14]) +
      (P3[2] * (V2[2] * T3[14] - V2[5] * T3[2]) + P3[3] * (V2[4] * T3[2] -
      V2[2] * T3[10]))) + (P1[2] * (P3[0] * (V2[3] * T3[14] - V2[5] * T3[6]) +
      (P3[1] * (V2[5] * T3[2] - V2[2] * T3[14]) + P3[3] * (V2[2] * T3[6] -
      V2[3] * T3[2]))) + P1[3] * (P3[0] * (V2[4] * T3[6] - V2[3] * T3[10]) +
      (P3[1] * (V2[2] * T3[10] - V2[4] * T3[2]) + P3[2] * (V2[3] * T3[2] -
      V2[2] * T3[6])))))) + (P2[1] * (P1[0] * (P3[1] * (V2[5] * T3[11] - V2[4]
      * T3[15]) + (P3[2] * (V2[3] * T3[15] - V2[5] * T3[7]) + P3[3] * (V2[4] *
      T3[7] - V2[3] * T3[11]))) + (P1[1] * (P3[0] * (V2[4] * T3[15] - V2[5] *
      T3[11]) + (P3[2] * (V2[5] * T3[3] - V2[2] * T3[15]) + P3[3] * (V2[2] *
      T3[11] - V2[4] * T3[3]))) + (P1[2] * (P3[0] * (V2[5] * T3[7] - V2[3] *
      T3[15]) + (P3[1] * (V2[2] * T3[15] - V2[5] * T3[3]) + P3[3] * (V2[3] *
      T3[3] - V2[2] * T3[7]))) + P1[3] * (P3[0] * (V2[3] * T3[11] - V2[4] *
      T3[7]) + (P3[1] * (V2[4] * T3[3] - V2[2] * T3[11]) + P3[2] * (V2[2] *
      T3[7] - V2[3] * T3[3])))))) + (P2[2] * (P1[0] * (P3[1] * (V2[5] * T3[12]
      - V2[4] * T3[16]) + (P3[2] * (V2[3] * T3[16] - V2[5] * T3[8]) + P3[3] *
      (V2[4] * T3[8] - V2[3] * T3[12]))) + (P1[1] * (P3[0] * (V2[4] * T3[16] -
      V2[5] * T3[12]) + (P3[2] * (V2[5] * T3[4] - V2[2] * T3[16]) + P3[3] *
      (V2[2] * T3[12] - V2[4] * T3[4]))) + (P1[2] * (P3[0] * (V2[5] * T3[8] -
      V2[3] * T3[16]) + (P3[1] * (V2[2] * T3[16] - V2[5] * T3[4]) + P3[3] *
      (V2[3] * T3[4] - V2[2] * T3[8]))) + P1[3] * (P3[0] * (V2[3] * T3[12] -
      V2[4] * T3[8]) + (P3[1] * (V2[4] * T3[4] - V2[2] * T3[12]) + P3[2] *
      (V2[2] * T3[8] - V2[3] * T3[4])))))) + P2[3] * (P1[0] * (P3[1] * (V2[5] *
      T3[13] - V2[4] * T3[17]) + (P3[2] * (V2[3] * T3[17] - V2[5] * T3[9]) +
      P3[3] * (V2[4] * T3[9] - V2[3] * T3[13]))) + (P1[1] * (P3[0] * (V2[4] *
      T3[17] - V2[5] * T3[13]) + (P3[2] * (V2[5] * T3[5] - V2[2] * T3[17]) +
      P3[3] * (V2[2] * T3[13] - V2[4] * T3[5]))) + (P1[2] * (P3[0] * (V2[5] *
      T3[9] - V2[3] * T3[17]) + (P3[1] * (V2[2] * T3[17] - V2[5] * T3[5]) +
      P3[3] * (V2[3] * T3[5] - V2[2] * T3[9]))) + P1[3] * (P3[0] * (V2[3] *
      T3[13] - V2[4] * T3[9]) + (P3[1] * (V2[4] * T3[5] - V2[2] * T3[13]) +
      P3[2] * (V2[2] * T3[9] - V2[3] * T3[5])))))))));
  TMP64 = -1. * (P2[0] * (P1[0] * (P3[1] * (V1[4] * T3[14] - V1[5] * T3[10]) +
      (P3[2] * (V1[5] * T3[6] - V1[3] * T3[14]) + P3[3] * (V1[3] * T3[10] -
      V1[4] * T3[6]))) + (P1[1] * (P3[0] * (V1[5] * T3[10] - V1[4] * T3[14]) +
      (P3[2] * (V1[2] * T3[14] - V1[5] * T3[2]) + P3[3] * (V1[4] * T3[2] -
      V1[2] * T3[10]))) + (P1[2] * (P3[0] * (V1[3] * T3[14] - V1[5] * T3[6]) +
      (P3[1] * (V1[5] * T3[2] - V1[2] * T3[14]) + P3[3] * (V1[2] * T3[6] -
      V1[3] * T3[2]))) + P1[3] * (P3[0] * (V1[4] * T3[6] - V1[3] * T3[10]) +
      (P3[1] * (V1[2] * T3[10] - V1[4] * T3[2]) + P3[2] * (V1[3] * T3[2] -
      V1[2] * T3[6])))))) + (P2[1] * (P1[0] * (P3[1] * (V1[5] * T3[11] - V1[4]
      * T3[15]) + (P3[2] * (V1[3] * T3[15] - V1[5] * T3[7]) + P3[3] * (V1[4] *
      T3[7] - V1[3] * T3[11]))) + (P1[1] * (P3[0] * (V1[4] * T3[15] - V1[5] *
      T3[11]) + (P3[2] * (V1[5] * T3[3] - V1[2] * T3[15]) + P3[3] * (V1[2] *
      T3[11] - V1[4] * T3[3]))) + (P1[2] * (P3[0] * (V1[5] * T3[7] - V1[3] *
      T3[15]) + (P3[1] * (V1[2] * T3[15] - V1[5] * T3[3]) + P3[3] * (V1[3] *
      T3[3] - V1[2] * T3[7]))) + P1[3] * (P3[0] * (V1[3] * T3[11] - V1[4] *
      T3[7]) + (P3[1] * (V1[4] * T3[3] - V1[2] * T3[11]) + P3[2] * (V1[2] *
      T3[7] - V1[3] * T3[3])))))) + (P2[2] * (P1[0] * (P3[1] * (V1[5] * T3[12]
      - V1[4] * T3[16]) + (P3[2] * (V1[3] * T3[16] - V1[5] * T3[8]) + P3[3] *
      (V1[4] * T3[8] - V1[3] * T3[12]))) + (P1[1] * (P3[0] * (V1[4] * T3[16] -
      V1[5] * T3[12]) + (P3[2] * (V1[5] * T3[4] - V1[2] * T3[16]) + P3[3] *
      (V1[2] * T3[12] - V1[4] * T3[4]))) + (P1[2] * (P3[0] * (V1[5] * T3[8] -
      V1[3] * T3[16]) + (P3[1] * (V1[2] * T3[16] - V1[5] * T3[4]) + P3[3] *
      (V1[3] * T3[4] - V1[2] * T3[8]))) + P1[3] * (P3[0] * (V1[3] * T3[12] -
      V1[4] * T3[8]) + (P3[1] * (V1[4] * T3[4] - V1[2] * T3[12]) + P3[2] *
      (V1[2] * T3[8] - V1[3] * T3[4])))))) + P2[3] * (P1[0] * (P3[1] * (V1[5] *
      T3[13] - V1[4] * T3[17]) + (P3[2] * (V1[3] * T3[17] - V1[5] * T3[9]) +
      P3[3] * (V1[4] * T3[9] - V1[3] * T3[13]))) + (P1[1] * (P3[0] * (V1[4] *
      T3[17] - V1[5] * T3[13]) + (P3[2] * (V1[5] * T3[5] - V1[2] * T3[17]) +
      P3[3] * (V1[2] * T3[13] - V1[4] * T3[5]))) + (P1[2] * (P3[0] * (V1[5] *
      T3[9] - V1[3] * T3[17]) + (P3[1] * (V1[2] * T3[17] - V1[5] * T3[5]) +
      P3[3] * (V1[3] * T3[5] - V1[2] * T3[9]))) + P1[3] * (P3[0] * (V1[3] *
      T3[13] - V1[4] * T3[9]) + (P3[1] * (V1[4] * T3[5] - V1[2] * T3[13]) +
      P3[2] * (V1[2] * T3[9] - V1[3] * T3[5])))))))));
  TMP60 = -1. * (P2[0] * (P1[0] * (P3[1] * (V1[4] * T3[5] - V1[5] * T3[4]) +
      (P3[2] * (V1[5] * T3[3] - V1[3] * T3[5]) + P3[3] * (V1[3] * T3[4] - V1[4]
      * T3[3]))) + (P1[1] * (P3[0] * (V1[5] * T3[4] - V1[4] * T3[5]) + (P3[2] *
      (V1[2] * T3[5] - V1[5] * T3[2]) + P3[3] * (V1[4] * T3[2] - V1[2] *
      T3[4]))) + (P1[2] * (P3[0] * (V1[3] * T3[5] - V1[5] * T3[3]) + (P3[1] *
      (V1[5] * T3[2] - V1[2] * T3[5]) + P3[3] * (V1[2] * T3[3] - V1[3] *
      T3[2]))) + P1[3] * (P3[0] * (V1[4] * T3[3] - V1[3] * T3[4]) + (P3[1] *
      (V1[2] * T3[4] - V1[4] * T3[2]) + P3[2] * (V1[3] * T3[2] - V1[2] *
      T3[3])))))) + (P2[1] * (P1[0] * (P3[1] * (V1[5] * T3[8] - V1[4] * T3[9])
      + (P3[2] * (V1[3] * T3[9] - V1[5] * T3[7]) + P3[3] * (V1[4] * T3[7] -
      V1[3] * T3[8]))) + (P1[1] * (P3[0] * (V1[4] * T3[9] - V1[5] * T3[8]) +
      (P3[2] * (V1[5] * T3[6] - V1[2] * T3[9]) + P3[3] * (V1[2] * T3[8] - V1[4]
      * T3[6]))) + (P1[2] * (P3[0] * (V1[5] * T3[7] - V1[3] * T3[9]) + (P3[1] *
      (V1[2] * T3[9] - V1[5] * T3[6]) + P3[3] * (V1[3] * T3[6] - V1[2] *
      T3[7]))) + P1[3] * (P3[0] * (V1[3] * T3[8] - V1[4] * T3[7]) + (P3[1] *
      (V1[4] * T3[6] - V1[2] * T3[8]) + P3[2] * (V1[2] * T3[7] - V1[3] *
      T3[6])))))) + (P2[2] * (P1[0] * (P3[1] * (V1[5] * T3[12] - V1[4] *
      T3[13]) + (P3[2] * (V1[3] * T3[13] - V1[5] * T3[11]) + P3[3] * (V1[4] *
      T3[11] - V1[3] * T3[12]))) + (P1[1] * (P3[0] * (V1[4] * T3[13] - V1[5] *
      T3[12]) + (P3[2] * (V1[5] * T3[10] - V1[2] * T3[13]) + P3[3] * (V1[2] *
      T3[12] - V1[4] * T3[10]))) + (P1[2] * (P3[0] * (V1[5] * T3[11] - V1[3] *
      T3[13]) + (P3[1] * (V1[2] * T3[13] - V1[5] * T3[10]) + P3[3] * (V1[3] *
      T3[10] - V1[2] * T3[11]))) + P1[3] * (P3[0] * (V1[3] * T3[12] - V1[4] *
      T3[11]) + (P3[1] * (V1[4] * T3[10] - V1[2] * T3[12]) + P3[2] * (V1[2] *
      T3[11] - V1[3] * T3[10])))))) + P2[3] * (P1[0] * (P3[1] * (V1[5] * T3[16]
      - V1[4] * T3[17]) + (P3[2] * (V1[3] * T3[17] - V1[5] * T3[15]) + P3[3] *
      (V1[4] * T3[15] - V1[3] * T3[16]))) + (P1[1] * (P3[0] * (V1[4] * T3[17] -
      V1[5] * T3[16]) + (P3[2] * (V1[5] * T3[14] - V1[2] * T3[17]) + P3[3] *
      (V1[2] * T3[16] - V1[4] * T3[14]))) + (P1[2] * (P3[0] * (V1[5] * T3[15] -
      V1[3] * T3[17]) + (P3[1] * (V1[2] * T3[17] - V1[5] * T3[14]) + P3[3] *
      (V1[3] * T3[14] - V1[2] * T3[15]))) + P1[3] * (P3[0] * (V1[3] * T3[16] -
      V1[4] * T3[15]) + (P3[1] * (V1[4] * T3[14] - V1[2] * T3[16]) + P3[2] *
      (V1[2] * T3[15] - V1[3] * T3[14])))))))));
  TMP61 = -1. * (P1[0] * (P2[0] * (P3[1] * (V2[4] * T3[14] - V2[5] * T3[10]) +
      (P3[2] * (V2[5] * T3[6] - V2[3] * T3[14]) + P3[3] * (V2[3] * T3[10] -
      V2[4] * T3[6]))) + (P2[1] * (P3[0] * (V2[5] * T3[10] - V2[4] * T3[14]) +
      (P3[2] * (V2[2] * T3[14] - V2[5] * T3[2]) + P3[3] * (V2[4] * T3[2] -
      V2[2] * T3[10]))) + (P2[2] * (P3[0] * (V2[3] * T3[14] - V2[5] * T3[6]) +
      (P3[1] * (V2[5] * T3[2] - V2[2] * T3[14]) + P3[3] * (V2[2] * T3[6] -
      V2[3] * T3[2]))) + P2[3] * (P3[0] * (V2[4] * T3[6] - V2[3] * T3[10]) +
      (P3[1] * (V2[2] * T3[10] - V2[4] * T3[2]) + P3[2] * (V2[3] * T3[2] -
      V2[2] * T3[6])))))) + (P1[1] * (P2[0] * (P3[1] * (V2[5] * T3[11] - V2[4]
      * T3[15]) + (P3[2] * (V2[3] * T3[15] - V2[5] * T3[7]) + P3[3] * (V2[4] *
      T3[7] - V2[3] * T3[11]))) + (P2[1] * (P3[0] * (V2[4] * T3[15] - V2[5] *
      T3[11]) + (P3[2] * (V2[5] * T3[3] - V2[2] * T3[15]) + P3[3] * (V2[2] *
      T3[11] - V2[4] * T3[3]))) + (P2[2] * (P3[0] * (V2[5] * T3[7] - V2[3] *
      T3[15]) + (P3[1] * (V2[2] * T3[15] - V2[5] * T3[3]) + P3[3] * (V2[3] *
      T3[3] - V2[2] * T3[7]))) + P2[3] * (P3[0] * (V2[3] * T3[11] - V2[4] *
      T3[7]) + (P3[1] * (V2[4] * T3[3] - V2[2] * T3[11]) + P3[2] * (V2[2] *
      T3[7] - V2[3] * T3[3])))))) + (P1[2] * (P2[0] * (P3[1] * (V2[5] * T3[12]
      - V2[4] * T3[16]) + (P3[2] * (V2[3] * T3[16] - V2[5] * T3[8]) + P3[3] *
      (V2[4] * T3[8] - V2[3] * T3[12]))) + (P2[1] * (P3[0] * (V2[4] * T3[16] -
      V2[5] * T3[12]) + (P3[2] * (V2[5] * T3[4] - V2[2] * T3[16]) + P3[3] *
      (V2[2] * T3[12] - V2[4] * T3[4]))) + (P2[2] * (P3[0] * (V2[5] * T3[8] -
      V2[3] * T3[16]) + (P3[1] * (V2[2] * T3[16] - V2[5] * T3[4]) + P3[3] *
      (V2[3] * T3[4] - V2[2] * T3[8]))) + P2[3] * (P3[0] * (V2[3] * T3[12] -
      V2[4] * T3[8]) + (P3[1] * (V2[4] * T3[4] - V2[2] * T3[12]) + P3[2] *
      (V2[2] * T3[8] - V2[3] * T3[4])))))) + P1[3] * (P2[0] * (P3[1] * (V2[5] *
      T3[13] - V2[4] * T3[17]) + (P3[2] * (V2[3] * T3[17] - V2[5] * T3[9]) +
      P3[3] * (V2[4] * T3[9] - V2[3] * T3[13]))) + (P2[1] * (P3[0] * (V2[4] *
      T3[17] - V2[5] * T3[13]) + (P3[2] * (V2[5] * T3[5] - V2[2] * T3[17]) +
      P3[3] * (V2[2] * T3[13] - V2[4] * T3[5]))) + (P2[2] * (P3[0] * (V2[5] *
      T3[9] - V2[3] * T3[17]) + (P3[1] * (V2[2] * T3[17] - V2[5] * T3[5]) +
      P3[3] * (V2[3] * T3[5] - V2[2] * T3[9]))) + P2[3] * (P3[0] * (V2[3] *
      T3[13] - V2[4] * T3[9]) + (P3[1] * (V2[4] * T3[5] - V2[2] * T3[13]) +
      P3[2] * (V2[2] * T3[9] - V2[3] * T3[5])))))))));
  vertex = COUP * - 1. * (TMP23 * (+cI * (TMP58 + TMP60 + TMP62 + TMP64)) +
      TMP26 * (+cI * (TMP57 + TMP59 + TMP61 + TMP63)));
}


void VVT8_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> TMP22; 
  complex<double> TMP10; 
  double P2[4]; 
  complex<double> TMP7; 
  complex<double> TMP21; 
  complex<double> TMP4; 
  complex<double> TMP9; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  TMP21 = (P1[0] * - 1. * (P2[1] * T3[6] + P2[2] * T3[10] + P2[3] * T3[14] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[11] + P2[3] *
      T3[15] - P2[0] * T3[3]) + (P1[2] * (P2[1] * T3[8] + P2[2] * T3[12] +
      P2[3] * T3[16] - P2[0] * T3[4]) + P1[3] * (P2[1] * T3[9] + P2[2] * T3[13]
      + P2[3] * T3[17] - P2[0] * T3[5]))));
  TMP22 = (P1[0] * - 1. * (P2[1] * T3[3] + P2[2] * T3[4] + P2[3] * T3[5] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[8] + P2[3] * T3[9]
      - P2[0] * T3[6]) + (P1[2] * (P2[1] * T3[11] + P2[2] * T3[12] + P2[3] *
      T3[13] - P2[0] * T3[10]) + P1[3] * (P2[1] * T3[15] + P2[2] * T3[16] +
      P2[3] * T3[17] - P2[0] * T3[14]))));
  TMP9 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP4 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP7 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP10 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  vertex = COUP * (TMP10 * TMP4 * (+cI * (TMP21 + TMP22)) - TMP7 * TMP9 * (+cI
      * (TMP21 + TMP22)));
}


void VVT1_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP2; 
  complex<double> TMP1; 
  double P1[4]; 
  complex<double> TMP10; 
  double P2[4]; 
  double P3[4]; 
  complex<double> TMP30; 
  complex<double> denom; 
  double OM3; 
  complex<double> TMP29; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +V1[0] + V2[0]; 
  T3[1] = +V1[1] + V2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP10 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP30 = -1. * (P1[0] * (P2[1] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P2[2] *
      (V2[3] * V1[5] - V2[5] * V1[3]) + P2[3] * (V2[4] * V1[3] - V2[3] *
      V1[4]))) + (P1[1] * (P2[0] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P2[2] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P2[3] * (V2[2] * V1[4] - V2[4] *
      V1[2]))) + (P1[2] * (P2[0] * (V2[5] * V1[3] - V2[3] * V1[5]) + (P2[1] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P2[3] * (V2[3] * V1[2] - V2[2] *
      V1[3]))) + P1[3] * (P2[0] * (V2[3] * V1[4] - V2[4] * V1[3]) + (P2[1] *
      (V2[4] * V1[2] - V2[2] * V1[4]) + P2[2] * (V2[2] * V1[3] - V2[3] *
      V1[2]))))));
  TMP1 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP29 = -1. * (P1[0] * (P2[1] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P2[2] *
      (V2[5] * V1[3] - V2[3] * V1[5]) + P2[3] * (V2[3] * V1[4] - V2[4] *
      V1[3]))) + (P1[1] * (P2[0] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P2[2] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P2[3] * (V2[4] * V1[2] - V2[2] *
      V1[4]))) + (P1[2] * (P2[0] * (V2[3] * V1[5] - V2[5] * V1[3]) + (P2[1] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P2[3] * (V2[2] * V1[3] - V2[3] *
      V1[2]))) + P1[3] * (P2[0] * (V2[4] * V1[3] - V2[3] * V1[4]) + (P2[1] *
      (V2[2] * V1[4] - V2[4] * V1[2]) + P2[2] * (V2[3] * V1[2] - V2[2] *
      V1[3]))))));
  TMP2 = (P3[0] * P2[0] - P3[1] * P2[1] - P3[2] * P2[2] - P3[3] * P2[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * 8. * (OM3 * (P3[0] * (P3[0] * (OM3 * 0.666666667 * TMP1 *
      TMP2 * (-cI * (TMP29) + cI * (TMP30)) + 0.333333333 * (TMP10 * (-cI *
      (TMP29) + cI * (TMP30)))) + (P1[0] * TMP2 * (-cI * (TMP30) + cI *
      (TMP29)) + P2[0] * TMP1 * (-cI * (TMP30) + cI * (TMP29)))) + 0.333333333
      * (TMP1 * TMP2 * (-cI * (TMP29) + cI * (TMP30)))) + (P1[0] * P2[0] * (-cI
      * (TMP29) + cI * (TMP30)) + 0.333333333 * (TMP10 * (-cI * (TMP30) + cI *
      (TMP29)))));
  T3[3] = denom * 4. * (OM3 * (P3[0] * (P3[1] * (OM3 * 1.333333333 * TMP1 *
      TMP2 * (-cI * (TMP29) + cI * (TMP30)) + 0.666666667 * (TMP10 * (-cI *
      (TMP29) + cI * (TMP30)))) + (P1[1] * TMP2 * (-cI * (TMP30) + cI *
      (TMP29)) + P2[1] * TMP1 * (-cI * (TMP30) + cI * (TMP29)))) + P3[1] *
      (P1[0] * TMP2 * (-cI * (TMP30) + cI * (TMP29)) + P2[0] * TMP1 * (-cI *
      (TMP30) + cI * (TMP29)))) + (P1[0] * P2[1] * (-cI * (TMP29) + cI *
      (TMP30)) + P1[1] * P2[0] * (-cI * (TMP29) + cI * (TMP30))));
  T3[4] = denom * 4. * (OM3 * (P3[0] * (P3[2] * (OM3 * 1.333333333 * TMP1 *
      TMP2 * (-cI * (TMP29) + cI * (TMP30)) + 0.666666667 * (TMP10 * (-cI *
      (TMP29) + cI * (TMP30)))) + (P1[2] * TMP2 * (-cI * (TMP30) + cI *
      (TMP29)) + P2[2] * TMP1 * (-cI * (TMP30) + cI * (TMP29)))) + P3[2] *
      (P1[0] * TMP2 * (-cI * (TMP30) + cI * (TMP29)) + P2[0] * TMP1 * (-cI *
      (TMP30) + cI * (TMP29)))) + (P1[0] * P2[2] * (-cI * (TMP29) + cI *
      (TMP30)) + P1[2] * P2[0] * (-cI * (TMP29) + cI * (TMP30))));
  T3[5] = denom * 4. * (OM3 * (P3[0] * (P3[3] * (OM3 * 1.333333333 * TMP1 *
      TMP2 * (-cI * (TMP29) + cI * (TMP30)) + 0.666666667 * (TMP10 * (-cI *
      (TMP29) + cI * (TMP30)))) + (P1[3] * TMP2 * (-cI * (TMP30) + cI *
      (TMP29)) + P2[3] * TMP1 * (-cI * (TMP30) + cI * (TMP29)))) + P3[3] *
      (P1[0] * TMP2 * (-cI * (TMP30) + cI * (TMP29)) + P2[0] * TMP1 * (-cI *
      (TMP30) + cI * (TMP29)))) + (P1[0] * P2[3] * (-cI * (TMP29) + cI *
      (TMP30)) + P1[3] * P2[0] * (-cI * (TMP29) + cI * (TMP30))));
  T3[6] = denom * 4. * (OM3 * (P3[0] * (P3[1] * (OM3 * 1.333333333 * TMP1 *
      TMP2 * (-cI * (TMP29) + cI * (TMP30)) + 0.666666667 * (TMP10 * (-cI *
      (TMP29) + cI * (TMP30)))) + (P1[1] * TMP2 * (-cI * (TMP30) + cI *
      (TMP29)) + P2[1] * TMP1 * (-cI * (TMP30) + cI * (TMP29)))) + P3[1] *
      (P1[0] * TMP2 * (-cI * (TMP30) + cI * (TMP29)) + P2[0] * TMP1 * (-cI *
      (TMP30) + cI * (TMP29)))) + (P1[0] * P2[1] * (-cI * (TMP29) + cI *
      (TMP30)) + P1[1] * P2[0] * (-cI * (TMP29) + cI * (TMP30))));
  T3[7] = denom * 8. * (OM3 * (P3[1] * (P3[1] * (OM3 * 0.666666667 * TMP1 *
      TMP2 * (-cI * (TMP29) + cI * (TMP30)) + 0.333333333 * (TMP10 * (-cI *
      (TMP29) + cI * (TMP30)))) + (P1[1] * TMP2 * (-cI * (TMP30) + cI *
      (TMP29)) + P2[1] * TMP1 * (-cI * (TMP30) + cI * (TMP29)))) + 0.333333333
      * (TMP1 * TMP2 * (-cI * (TMP30) + cI * (TMP29)))) + (P1[1] * P2[1] * (-cI
      * (TMP29) + cI * (TMP30)) + 0.333333333 * (TMP10 * (-cI * (TMP29) + cI *
      (TMP30)))));
  T3[8] = denom * 4. * (OM3 * (P3[1] * (P3[2] * (OM3 * 1.333333333 * TMP1 *
      TMP2 * (-cI * (TMP29) + cI * (TMP30)) + 0.666666667 * (TMP10 * (-cI *
      (TMP29) + cI * (TMP30)))) + (P1[2] * TMP2 * (-cI * (TMP30) + cI *
      (TMP29)) + P2[2] * TMP1 * (-cI * (TMP30) + cI * (TMP29)))) + P3[2] *
      (P1[1] * TMP2 * (-cI * (TMP30) + cI * (TMP29)) + P2[1] * TMP1 * (-cI *
      (TMP30) + cI * (TMP29)))) + (P1[1] * P2[2] * (-cI * (TMP29) + cI *
      (TMP30)) + P1[2] * P2[1] * (-cI * (TMP29) + cI * (TMP30))));
  T3[9] = denom * 4. * (OM3 * (P3[1] * (P3[3] * (OM3 * 1.333333333 * TMP1 *
      TMP2 * (-cI * (TMP29) + cI * (TMP30)) + 0.666666667 * (TMP10 * (-cI *
      (TMP29) + cI * (TMP30)))) + (P1[3] * TMP2 * (-cI * (TMP30) + cI *
      (TMP29)) + P2[3] * TMP1 * (-cI * (TMP30) + cI * (TMP29)))) + P3[3] *
      (P1[1] * TMP2 * (-cI * (TMP30) + cI * (TMP29)) + P2[1] * TMP1 * (-cI *
      (TMP30) + cI * (TMP29)))) + (P1[1] * P2[3] * (-cI * (TMP29) + cI *
      (TMP30)) + P1[3] * P2[1] * (-cI * (TMP29) + cI * (TMP30))));
  T3[10] = denom * 4. * (OM3 * (P3[0] * (P3[2] * (OM3 * 1.333333333 * TMP1 *
      TMP2 * (-cI * (TMP29) + cI * (TMP30)) + 0.666666667 * (TMP10 * (-cI *
      (TMP29) + cI * (TMP30)))) + (P1[2] * TMP2 * (-cI * (TMP30) + cI *
      (TMP29)) + P2[2] * TMP1 * (-cI * (TMP30) + cI * (TMP29)))) + P3[2] *
      (P1[0] * TMP2 * (-cI * (TMP30) + cI * (TMP29)) + P2[0] * TMP1 * (-cI *
      (TMP30) + cI * (TMP29)))) + (P1[0] * P2[2] * (-cI * (TMP29) + cI *
      (TMP30)) + P1[2] * P2[0] * (-cI * (TMP29) + cI * (TMP30))));
  T3[11] = denom * 4. * (OM3 * (P3[1] * (P3[2] * (OM3 * 1.333333333 * TMP1 *
      TMP2 * (-cI * (TMP29) + cI * (TMP30)) + 0.666666667 * (TMP10 * (-cI *
      (TMP29) + cI * (TMP30)))) + (P1[2] * TMP2 * (-cI * (TMP30) + cI *
      (TMP29)) + P2[2] * TMP1 * (-cI * (TMP30) + cI * (TMP29)))) + P3[2] *
      (P1[1] * TMP2 * (-cI * (TMP30) + cI * (TMP29)) + P2[1] * TMP1 * (-cI *
      (TMP30) + cI * (TMP29)))) + (P1[1] * P2[2] * (-cI * (TMP29) + cI *
      (TMP30)) + P1[2] * P2[1] * (-cI * (TMP29) + cI * (TMP30))));
  T3[12] = denom * 8. * (OM3 * (P3[2] * (P3[2] * (OM3 * 0.666666667 * TMP1 *
      TMP2 * (-cI * (TMP29) + cI * (TMP30)) + 0.333333333 * (TMP10 * (-cI *
      (TMP29) + cI * (TMP30)))) + (P1[2] * TMP2 * (-cI * (TMP30) + cI *
      (TMP29)) + P2[2] * TMP1 * (-cI * (TMP30) + cI * (TMP29)))) + 0.333333333
      * (TMP1 * TMP2 * (-cI * (TMP30) + cI * (TMP29)))) + (P1[2] * P2[2] * (-cI
      * (TMP29) + cI * (TMP30)) + 0.333333333 * (TMP10 * (-cI * (TMP29) + cI *
      (TMP30)))));
  T3[13] = denom * 4. * (OM3 * (P3[2] * (P3[3] * (OM3 * 1.333333333 * TMP1 *
      TMP2 * (-cI * (TMP29) + cI * (TMP30)) + 0.666666667 * (TMP10 * (-cI *
      (TMP29) + cI * (TMP30)))) + (P1[3] * TMP2 * (-cI * (TMP30) + cI *
      (TMP29)) + P2[3] * TMP1 * (-cI * (TMP30) + cI * (TMP29)))) + P3[3] *
      (P1[2] * TMP2 * (-cI * (TMP30) + cI * (TMP29)) + P2[2] * TMP1 * (-cI *
      (TMP30) + cI * (TMP29)))) + (P1[2] * P2[3] * (-cI * (TMP29) + cI *
      (TMP30)) + P1[3] * P2[2] * (-cI * (TMP29) + cI * (TMP30))));
  T3[14] = denom * 4. * (OM3 * (P3[0] * (P3[3] * (OM3 * 1.333333333 * TMP1 *
      TMP2 * (-cI * (TMP29) + cI * (TMP30)) + 0.666666667 * (TMP10 * (-cI *
      (TMP29) + cI * (TMP30)))) + (P1[3] * TMP2 * (-cI * (TMP30) + cI *
      (TMP29)) + P2[3] * TMP1 * (-cI * (TMP30) + cI * (TMP29)))) + P3[3] *
      (P1[0] * TMP2 * (-cI * (TMP30) + cI * (TMP29)) + P2[0] * TMP1 * (-cI *
      (TMP30) + cI * (TMP29)))) + (P1[0] * P2[3] * (-cI * (TMP29) + cI *
      (TMP30)) + P1[3] * P2[0] * (-cI * (TMP29) + cI * (TMP30))));
  T3[15] = denom * 4. * (OM3 * (P3[1] * (P3[3] * (OM3 * 1.333333333 * TMP1 *
      TMP2 * (-cI * (TMP29) + cI * (TMP30)) + 0.666666667 * (TMP10 * (-cI *
      (TMP29) + cI * (TMP30)))) + (P1[3] * TMP2 * (-cI * (TMP30) + cI *
      (TMP29)) + P2[3] * TMP1 * (-cI * (TMP30) + cI * (TMP29)))) + P3[3] *
      (P1[1] * TMP2 * (-cI * (TMP30) + cI * (TMP29)) + P2[1] * TMP1 * (-cI *
      (TMP30) + cI * (TMP29)))) + (P1[1] * P2[3] * (-cI * (TMP29) + cI *
      (TMP30)) + P1[3] * P2[1] * (-cI * (TMP29) + cI * (TMP30))));
  T3[16] = denom * 4. * (OM3 * (P3[2] * (P3[3] * (OM3 * 1.333333333 * TMP1 *
      TMP2 * (-cI * (TMP29) + cI * (TMP30)) + 0.666666667 * (TMP10 * (-cI *
      (TMP29) + cI * (TMP30)))) + (P1[3] * TMP2 * (-cI * (TMP30) + cI *
      (TMP29)) + P2[3] * TMP1 * (-cI * (TMP30) + cI * (TMP29)))) + P3[3] *
      (P1[2] * TMP2 * (-cI * (TMP30) + cI * (TMP29)) + P2[2] * TMP1 * (-cI *
      (TMP30) + cI * (TMP29)))) + (P1[2] * P2[3] * (-cI * (TMP29) + cI *
      (TMP30)) + P1[3] * P2[2] * (-cI * (TMP29) + cI * (TMP30))));
  T3[17] = denom * 8. * (OM3 * (P3[3] * (P3[3] * (OM3 * 0.666666667 * TMP1 *
      TMP2 * (-cI * (TMP29) + cI * (TMP30)) + 0.333333333 * (TMP10 * (-cI *
      (TMP29) + cI * (TMP30)))) + (P1[3] * TMP2 * (-cI * (TMP30) + cI *
      (TMP29)) + P2[3] * TMP1 * (-cI * (TMP30) + cI * (TMP29)))) + 0.333333333
      * (TMP1 * TMP2 * (-cI * (TMP30) + cI * (TMP29)))) + (P1[3] * P2[3] * (-cI
      * (TMP29) + cI * (TMP30)) + 0.333333333 * (TMP10 * (-cI * (TMP29) + cI *
      (TMP30)))));
}

void VVT1_10_11_12_13_3_5_7_8_9_3(complex<double> V1[], complex<double> V2[],
    complex<double> COUP1, complex<double> COUP2, complex<double> COUP3,
    complex<double> COUP4, complex<double> COUP5, complex<double> COUP6,
    complex<double> COUP7, complex<double> COUP8, complex<double> COUP9,
    complex<double> COUP10, double M3, double W3, complex<double> T3[])
{
//   complex<double> cI = complex<double> (0., 1.); 
//   double P3[4]; 
//   double P2[4]; 
//   double OM3; 
//   double P1[4]; 
  complex<double> Ttmp[18]; 
  complex<double> denom; 
  int i; 
  VVT1_3(V1, V2, COUP1, M3, W3, T3); 
  VVT10_3(V1, V2, COUP2, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT11_3(V1, V2, COUP3, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT12_3(V1, V2, COUP4, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT13_3(V1, V2, COUP5, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT3_3(V1, V2, COUP6, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT5_3(V1, V2, COUP7, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT7_3(V1, V2, COUP8, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT8_3(V1, V2, COUP9, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT9_3(V1, V2, COUP10, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
}

}  // end namespace $(namespace)s_HEF_MEK

