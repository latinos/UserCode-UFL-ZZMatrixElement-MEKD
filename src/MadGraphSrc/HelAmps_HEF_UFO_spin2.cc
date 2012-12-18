//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph 5 v. 1.5.3, 2012-11-01
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#include "HelAmps_HEF_UFO_spin2.h"

namespace MG5_HEF_UFO_spin2
{


double Sgn(double a, double b)
{
  return (b < 0)? - abs(a):abs(a); 
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

void VVT16_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP33; 
  double P1[4]; 
  complex<double> TMP22; 
  complex<double> TMP0; 
  double P2[4]; 
  complex<double> TMP23; 
  double P3[4]; 
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
  TMP33 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP0 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP22 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  TMP23 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * TMP0 * (OM3 * (P3[0] * (P3[0] * 0.666666667 * (+cI * (TMP33)
      + 2. * cI * (OM3 * TMP22 * TMP23)) + (-2. * cI * (P2[0] * TMP23 + P1[0] *
      TMP22))) + 0.666666667 * cI * (TMP22 * TMP23)) + (-0.666666667 * cI *
      (TMP33) + 2. * cI * (P1[0] * P2[0])));
  T3[3] = denom * TMP0 * (OM3 * (P3[0] * (P3[1] * 0.666666667 * (+cI * (TMP33)
      + 2. * cI * (OM3 * TMP22 * TMP23)) + (-cI * (P1[1] * TMP22 + P2[1] *
      TMP23))) - P3[1] * (+cI * (P2[0] * TMP23 + P1[0] * TMP22))) + (+cI *
      (P1[1] * P2[0] + P1[0] * P2[1])));
  T3[4] = denom * TMP0 * (OM3 * (P3[0] * (P3[2] * 0.666666667 * (+cI * (TMP33)
      + 2. * cI * (OM3 * TMP22 * TMP23)) + (-cI * (P1[2] * TMP22 + P2[2] *
      TMP23))) - P3[2] * (+cI * (P2[0] * TMP23 + P1[0] * TMP22))) + (+cI *
      (P1[2] * P2[0] + P1[0] * P2[2])));
  T3[5] = denom * TMP0 * (OM3 * (P3[0] * (P3[3] * 0.666666667 * (+cI * (TMP33)
      + 2. * cI * (OM3 * TMP22 * TMP23)) + (-cI * (P1[3] * TMP22 + P2[3] *
      TMP23))) - P3[3] * (+cI * (P2[0] * TMP23 + P1[0] * TMP22))) + (+cI *
      (P1[3] * P2[0] + P1[0] * P2[3])));
  T3[6] = denom * TMP0 * (OM3 * (P3[0] * (P3[1] * 0.666666667 * (+cI * (TMP33)
      + 2. * cI * (OM3 * TMP22 * TMP23)) + (-cI * (P2[1] * TMP23 + P1[1] *
      TMP22))) - P3[1] * (+cI * (P1[0] * TMP22 + P2[0] * TMP23))) + (+cI *
      (P1[0] * P2[1] + P1[1] * P2[0])));
  T3[7] = denom * TMP0 * (OM3 * (P3[1] * (P3[1] * 0.666666667 * (+cI * (TMP33)
      + 2. * cI * (OM3 * TMP22 * TMP23)) + (-2. * cI * (P2[1] * TMP23 + P1[1] *
      TMP22))) - 0.666666667 * cI * (TMP22 * TMP23)) + (+2. * cI * (P1[1] *
      P2[1]) + 0.666666667 * cI * (TMP33)));
  T3[8] = denom * TMP0 * (OM3 * (P3[1] * (P3[2] * 0.666666667 * (+cI * (TMP33)
      + 2. * cI * (OM3 * TMP22 * TMP23)) + (-cI * (P1[2] * TMP22 + P2[2] *
      TMP23))) - P3[2] * (+cI * (P2[1] * TMP23 + P1[1] * TMP22))) + (+cI *
      (P1[2] * P2[1] + P1[1] * P2[2])));
  T3[9] = denom * TMP0 * (OM3 * (P3[1] * (P3[3] * 0.666666667 * (+cI * (TMP33)
      + 2. * cI * (OM3 * TMP22 * TMP23)) + (-cI * (P1[3] * TMP22 + P2[3] *
      TMP23))) - P3[3] * (+cI * (P2[1] * TMP23 + P1[1] * TMP22))) + (+cI *
      (P1[3] * P2[1] + P1[1] * P2[3])));
  T3[10] = denom * TMP0 * (OM3 * (P3[0] * (P3[2] * 0.666666667 * (+cI * (TMP33)
      + 2. * cI * (OM3 * TMP22 * TMP23)) + (-cI * (P2[2] * TMP23 + P1[2] *
      TMP22))) - P3[2] * (+cI * (P1[0] * TMP22 + P2[0] * TMP23))) + (+cI *
      (P1[0] * P2[2] + P1[2] * P2[0])));
  T3[11] = denom * TMP0 * (OM3 * (P3[1] * (P3[2] * 0.666666667 * (+cI * (TMP33)
      + 2. * cI * (OM3 * TMP22 * TMP23)) + (-cI * (P2[2] * TMP23 + P1[2] *
      TMP22))) - P3[2] * (+cI * (P1[1] * TMP22 + P2[1] * TMP23))) + (+cI *
      (P1[1] * P2[2] + P1[2] * P2[1])));
  T3[12] = denom * TMP0 * (OM3 * (P3[2] * (P3[2] * 0.666666667 * (+cI * (TMP33)
      + 2. * cI * (OM3 * TMP22 * TMP23)) + (-2. * cI * (P2[2] * TMP23 + P1[2] *
      TMP22))) - 0.666666667 * cI * (TMP22 * TMP23)) + (+2. * cI * (P1[2] *
      P2[2]) + 0.666666667 * cI * (TMP33)));
  T3[13] = denom * TMP0 * (OM3 * (P3[2] * (P3[3] * 0.666666667 * (+cI * (TMP33)
      + 2. * cI * (OM3 * TMP22 * TMP23)) + (-cI * (P1[3] * TMP22 + P2[3] *
      TMP23))) - P3[3] * (+cI * (P2[2] * TMP23 + P1[2] * TMP22))) + (+cI *
      (P1[3] * P2[2] + P1[2] * P2[3])));
  T3[14] = denom * TMP0 * (OM3 * (P3[0] * (P3[3] * 0.666666667 * (+cI * (TMP33)
      + 2. * cI * (OM3 * TMP22 * TMP23)) + (-cI * (P2[3] * TMP23 + P1[3] *
      TMP22))) - P3[3] * (+cI * (P1[0] * TMP22 + P2[0] * TMP23))) + (+cI *
      (P1[0] * P2[3] + P1[3] * P2[0])));
  T3[15] = denom * TMP0 * (OM3 * (P3[1] * (P3[3] * 0.666666667 * (+cI * (TMP33)
      + 2. * cI * (OM3 * TMP22 * TMP23)) + (-cI * (P2[3] * TMP23 + P1[3] *
      TMP22))) - P3[3] * (+cI * (P1[1] * TMP22 + P2[1] * TMP23))) + (+cI *
      (P1[1] * P2[3] + P1[3] * P2[1])));
  T3[16] = denom * TMP0 * (OM3 * (P3[2] * (P3[3] * 0.666666667 * (+cI * (TMP33)
      + 2. * cI * (OM3 * TMP22 * TMP23)) + (-cI * (P2[3] * TMP23 + P1[3] *
      TMP22))) - P3[3] * (+cI * (P1[2] * TMP22 + P2[2] * TMP23))) + (+cI *
      (P1[2] * P2[3] + P1[3] * P2[2])));
  T3[17] = denom * TMP0 * (OM3 * (P3[3] * (P3[3] * 0.666666667 * (+cI * (TMP33)
      + 2. * cI * (OM3 * TMP22 * TMP23)) + (-2. * cI * (P2[3] * TMP23 + P1[3] *
      TMP22))) - 0.666666667 * cI * (TMP22 * TMP23)) + (+2. * cI * (P1[3] *
      P2[3]) + 0.666666667 * cI * (TMP33)));
}


void FFV42_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  complex<double> TMP10; 
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
  TMP10 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
      F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * - cI * (F2[4] * F1[2] + F2[5] * F1[3] - P3[0] * OM3 * TMP10); 
  V3[3] = denom * - cI * (-F2[5] * F1[2] - F2[4] * F1[3] - P3[1] * OM3 *
      TMP10);
  V3[4] = denom * - cI * (-cI * (F2[5] * F1[2]) + cI * (F2[4] * F1[3]) - P3[2]
      * OM3 * TMP10);
  V3[5] = denom * - cI * (F2[5] * F1[3] - F2[4] * F1[2] - P3[3] * OM3 * TMP10); 
}

void FFV42_44_3(complex<double> F1[], complex<double> F2[], complex<double>
    COUP1, complex<double> COUP2, double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  double P3[4]; 
  double OM3; 
  int i; 
  complex<double> Vtmp[6]; 
  FFV42_3(F1, F2, COUP1, M3, W3, V3); 
  FFV44_3(F1, F2, COUP2, M3, W3, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V3[i] = V3[i] + Vtmp[i]; 
    i++; 
  }
}

void VVT19_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  double P3[4]; 
  complex<double> TMP5; 
  complex<double> TMP28; 
  complex<double> TMP25; 
  complex<double> TMP9; 
  double P2[4]; 
  complex<double> TMP6; 
  complex<double> TMP24; 
  complex<double> TMP3; 
  double P1[4]; 
  complex<double> TMP7; 
  complex<double> TMP30; 
  complex<double> TMP27; 
  complex<double> TMP0; 
  complex<double> TMP31; 
  complex<double> TMP26; 
  complex<double> TMP4; 
  complex<double> TMP29; 
  complex<double> TMP8; 
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
  TMP24 = (P1[0] * (P1[1] * - 1. * (T3[6] + T3[3]) + (P1[2] * - 1. * (T3[10] +
      T3[4]) + (P1[3] * - 1. * (T3[14] + T3[5]) + P1[0] * T3[2]))) + (P1[1] *
      (P1[2] * (T3[11] + T3[8]) + (P1[3] * (T3[15] + T3[9]) + P1[1] * T3[7])) +
      (P1[2] * (P1[3] * (T3[16] + T3[13]) + P1[2] * T3[12]) + P1[3] * P1[3] *
      T3[17])));
  TMP25 = (P1[0] * - 1. * (P2[1] * T3[6] + P2[2] * T3[10] + P2[3] * T3[14] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[11] + P2[3] *
      T3[15] - P2[0] * T3[3]) + (P1[2] * (P2[1] * T3[8] + P2[2] * T3[12] +
      P2[3] * T3[16] - P2[0] * T3[4]) + P1[3] * (P2[1] * T3[9] + P2[2] * T3[13]
      + P2[3] * T3[17] - P2[0] * T3[5]))));
  TMP26 = (P1[0] * - 1. * (P2[1] * T3[3] + P2[2] * T3[4] + P2[3] * T3[5] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[8] + P2[3] * T3[9]
      - P2[0] * T3[6]) + (P1[2] * (P2[1] * T3[11] + P2[2] * T3[12] + P2[3] *
      T3[13] - P2[0] * T3[10]) + P1[3] * (P2[1] * T3[15] + P2[2] * T3[16] +
      P2[3] * T3[17] - P2[0] * T3[14]))));
  TMP27 = (P2[0] * (P2[1] * - 1. * (T3[6] + T3[3]) + (P2[2] * - 1. * (T3[10] +
      T3[4]) + (P2[3] * - 1. * (T3[14] + T3[5]) + P2[0] * T3[2]))) + (P2[1] *
      (P2[2] * (T3[11] + T3[8]) + (P2[3] * (T3[15] + T3[9]) + P2[1] * T3[7])) +
      (P2[2] * (P2[3] * (T3[16] + T3[13]) + P2[2] * T3[12]) + P2[3] * P2[3] *
      T3[17])));
  TMP28 = (P1[0] * - 1. * (V1[3] * T3[6] + V1[4] * T3[10] + V1[5] * T3[14] -
      V1[2] * T3[2]) + (P1[1] * (V1[3] * T3[7] + V1[4] * T3[11] + V1[5] *
      T3[15] - V1[2] * T3[3]) + (P1[2] * (V1[3] * T3[8] + V1[4] * T3[12] +
      V1[5] * T3[16] - V1[2] * T3[4]) + P1[3] * (V1[3] * T3[9] + V1[4] * T3[13]
      + V1[5] * T3[17] - V1[2] * T3[5]))));
  TMP29 = (P1[0] * - 1. * (V1[3] * T3[3] + V1[4] * T3[4] + V1[5] * T3[5] -
      V1[2] * T3[2]) + (P1[1] * (V1[3] * T3[7] + V1[4] * T3[8] + V1[5] * T3[9]
      - V1[2] * T3[6]) + (P1[2] * (V1[3] * T3[11] + V1[4] * T3[12] + V1[5] *
      T3[13] - V1[2] * T3[10]) + P1[3] * (V1[3] * T3[15] + V1[4] * T3[16] +
      V1[5] * T3[17] - V1[2] * T3[14]))));
  TMP0 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP9 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP8 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP5 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  TMP4 = (V1[2] * P3[0] - V1[3] * P3[1] - V1[4] * P3[2] - V1[5] * P3[3]); 
  TMP31 = (P2[0] * - 1. * (V2[3] * T3[3] + V2[4] * T3[4] + V2[5] * T3[5] -
      V2[2] * T3[2]) + (P2[1] * (V2[3] * T3[7] + V2[4] * T3[8] + V2[5] * T3[9]
      - V2[2] * T3[6]) + (P2[2] * (V2[3] * T3[11] + V2[4] * T3[12] + V2[5] *
      T3[13] - V2[2] * T3[10]) + P2[3] * (V2[3] * T3[15] + V2[4] * T3[16] +
      V2[5] * T3[17] - V2[2] * T3[14]))));
  TMP30 = (P2[0] * - 1. * (V2[3] * T3[6] + V2[4] * T3[10] + V2[5] * T3[14] -
      V2[2] * T3[2]) + (P2[1] * (V2[3] * T3[7] + V2[4] * T3[11] + V2[5] *
      T3[15] - V2[2] * T3[3]) + (P2[2] * (V2[3] * T3[8] + V2[4] * T3[12] +
      V2[5] * T3[16] - V2[2] * T3[4]) + P2[3] * (V2[3] * T3[9] + V2[4] * T3[13]
      + V2[5] * T3[17] - V2[2] * T3[5]))));
  TMP7 = (V1[2] * P2[0] - V1[3] * P2[1] - V1[4] * P2[2] - V1[5] * P2[3]); 
  TMP3 = (V2[2] * P1[0] - V2[3] * P1[1] - V2[4] * P1[2] - V2[5] * P1[3]); 
  TMP6 = (V2[2] * P3[0] - V2[3] * P3[1] - V2[4] * P3[2] - V2[5] * P3[3]); 
  vertex = COUP * (TMP9 * (TMP0 * (-cI * (TMP25 + TMP26) + cI * (TMP24 +
      TMP27)) + (+0.500000000 * (TMP6 * (+cI * (TMP28 + TMP29))) + TMP4 *
      0.500000000 * (+cI * (TMP30 + TMP31)))) + (TMP3 * (TMP7 * (-cI * (TMP24 +
      TMP27) + cI * (TMP25 + TMP26)) + (TMP5 * - 0.500000000 * (+cI * (TMP28 +
      TMP29)) - cI * (TMP4 * TMP27))) + (TMP7 * (TMP8 * - 0.500000000 * (+cI *
      (TMP30 + TMP31)) - cI * (TMP6 * TMP24)) + TMP0 * (+cI * (TMP5 * TMP24 +
      TMP8 * TMP27)))));
}


void VVT18_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  complex<double> TMP20; 
  complex<double> TMP6; 
  complex<double> TMP21; 
  complex<double> TMP4; 
  complex<double> TMP19; 
  complex<double> TMP18; 
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
  TMP20 = (P1[0] * - 1. * (V2[3] * T3[6] + V2[4] * T3[10] + V2[5] * T3[14] -
      V2[2] * T3[2]) + (P1[1] * (V2[3] * T3[7] + V2[4] * T3[11] + V2[5] *
      T3[15] - V2[2] * T3[3]) + (P1[2] * (V2[3] * T3[8] + V2[4] * T3[12] +
      V2[5] * T3[16] - V2[2] * T3[4]) + P1[3] * (V2[3] * T3[9] + V2[4] * T3[13]
      + V2[5] * T3[17] - V2[2] * T3[5]))));
  TMP21 = (P1[0] * - 1. * (V2[3] * T3[3] + V2[4] * T3[4] + V2[5] * T3[5] -
      V2[2] * T3[2]) + (P1[1] * (V2[3] * T3[7] + V2[4] * T3[8] + V2[5] * T3[9]
      - V2[2] * T3[6]) + (P1[2] * (V2[3] * T3[11] + V2[4] * T3[12] + V2[5] *
      T3[13] - V2[2] * T3[10]) + P1[3] * (V2[3] * T3[15] + V2[4] * T3[16] +
      V2[5] * T3[17] - V2[2] * T3[14]))));
  TMP19 = (P2[0] * - 1. * (V1[3] * T3[3] + V1[4] * T3[4] + V1[5] * T3[5] -
      V1[2] * T3[2]) + (P2[1] * (V1[3] * T3[7] + V1[4] * T3[8] + V1[5] * T3[9]
      - V1[2] * T3[6]) + (P2[2] * (V1[3] * T3[11] + V1[4] * T3[12] + V1[5] *
      T3[13] - V1[2] * T3[10]) + P2[3] * (V1[3] * T3[15] + V1[4] * T3[16] +
      V1[5] * T3[17] - V1[2] * T3[14]))));
  TMP18 = (P2[0] * - 1. * (V1[3] * T3[6] + V1[4] * T3[10] + V1[5] * T3[14] -
      V1[2] * T3[2]) + (P2[1] * (V1[3] * T3[7] + V1[4] * T3[11] + V1[5] *
      T3[15] - V1[2] * T3[3]) + (P2[2] * (V1[3] * T3[8] + V1[4] * T3[12] +
      V1[5] * T3[16] - V1[2] * T3[4]) + P2[3] * (V1[3] * T3[9] + V1[4] * T3[13]
      + V1[5] * T3[17] - V1[2] * T3[5]))));
  TMP4 = (V1[2] * P3[0] - V1[3] * P3[1] - V1[4] * P3[2] - V1[5] * P3[3]); 
  TMP6 = (V2[2] * P3[0] - V2[3] * P3[1] - V2[4] * P3[2] - V2[5] * P3[3]); 
  vertex = COUP * (TMP4 * - 1. * (+cI * (TMP20 + TMP21)) - TMP6 * (+cI * (TMP18
      + TMP19)));
}


void VVT12_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP45; 
  double P3[4]; 
  complex<double> TMP46; 
  complex<double> TMP47; 
  complex<double> TMP48; 
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
  TMP48 = (P2[0] * (P3[0] * (T3[14] * (V2[3] * V1[4] - V2[4] * V1[3]) + (T3[6]
      * (V2[4] * V1[5] - V2[5] * V1[4]) + T3[10] * (V2[5] * V1[3] - V2[3] *
      V1[5]))) + (P3[1] * (T3[2] * (V2[5] * V1[4] - V2[4] * V1[5]) + (T3[14] *
      (V2[4] * V1[2] - V2[2] * V1[4]) + T3[10] * (V2[2] * V1[5] - V2[5] *
      V1[2]))) + (P3[2] * (T3[2] * (V2[3] * V1[5] - V2[5] * V1[3]) + (T3[14] *
      (V2[2] * V1[3] - V2[3] * V1[2]) + T3[6] * (V2[5] * V1[2] - V2[2] *
      V1[5]))) + P3[3] * (T3[2] * (V2[4] * V1[3] - V2[3] * V1[4]) + (T3[6] *
      (V2[2] * V1[4] - V2[4] * V1[2]) + T3[10] * (V2[3] * V1[2] - V2[2] *
      V1[3])))))) + (P2[1] * (P3[0] * (T3[11] * (V2[3] * V1[5] - V2[5] * V1[3])
      + (T3[15] * (V2[4] * V1[3] - V2[3] * V1[4]) + T3[7] * (V2[5] * V1[4] -
      V2[4] * V1[5]))) + (P3[1] * (T3[11] * (V2[5] * V1[2] - V2[2] * V1[5]) +
      (T3[15] * (V2[2] * V1[4] - V2[4] * V1[2]) + T3[3] * (V2[4] * V1[5] -
      V2[5] * V1[4]))) + (P3[2] * (T3[15] * (V2[3] * V1[2] - V2[2] * V1[3]) +
      (T3[3] * (V2[5] * V1[3] - V2[3] * V1[5]) + T3[7] * (V2[2] * V1[5] - V2[5]
      * V1[2]))) + P3[3] * (T3[11] * (V2[2] * V1[3] - V2[3] * V1[2]) + (T3[3] *
      (V2[3] * V1[4] - V2[4] * V1[3]) + T3[7] * (V2[4] * V1[2] - V2[2] *
      V1[4])))))) + (P2[2] * (P3[0] * (T3[12] * (V2[3] * V1[5] - V2[5] * V1[3])
      + (T3[16] * (V2[4] * V1[3] - V2[3] * V1[4]) + T3[8] * (V2[5] * V1[4] -
      V2[4] * V1[5]))) + (P3[1] * (T3[12] * (V2[5] * V1[2] - V2[2] * V1[5]) +
      (T3[16] * (V2[2] * V1[4] - V2[4] * V1[2]) + T3[4] * (V2[4] * V1[5] -
      V2[5] * V1[4]))) + (P3[2] * (T3[16] * (V2[3] * V1[2] - V2[2] * V1[3]) +
      (T3[4] * (V2[5] * V1[3] - V2[3] * V1[5]) + T3[8] * (V2[2] * V1[5] - V2[5]
      * V1[2]))) + P3[3] * (T3[12] * (V2[2] * V1[3] - V2[3] * V1[2]) + (T3[4] *
      (V2[3] * V1[4] - V2[4] * V1[3]) + T3[8] * (V2[4] * V1[2] - V2[2] *
      V1[4])))))) + P2[3] * (P3[0] * (T3[13] * (V2[3] * V1[5] - V2[5] * V1[3])
      + (T3[17] * (V2[4] * V1[3] - V2[3] * V1[4]) + T3[9] * (V2[5] * V1[4] -
      V2[4] * V1[5]))) + (P3[1] * (T3[13] * (V2[5] * V1[2] - V2[2] * V1[5]) +
      (T3[17] * (V2[2] * V1[4] - V2[4] * V1[2]) + T3[5] * (V2[4] * V1[5] -
      V2[5] * V1[4]))) + (P3[2] * (T3[17] * (V2[3] * V1[2] - V2[2] * V1[3]) +
      (T3[5] * (V2[5] * V1[3] - V2[3] * V1[5]) + T3[9] * (V2[2] * V1[5] - V2[5]
      * V1[2]))) + P3[3] * (T3[13] * (V2[2] * V1[3] - V2[3] * V1[2]) + (T3[5] *
      (V2[3] * V1[4] - V2[4] * V1[3]) + T3[9] * (V2[4] * V1[2] - V2[2] *
      V1[4])))))))));
  TMP46 = (P2[0] * (P3[0] * (T3[3] * (V2[4] * V1[5] - V2[5] * V1[4]) + (T3[4] *
      (V2[5] * V1[3] - V2[3] * V1[5]) + T3[5] * (V2[3] * V1[4] - V2[4] *
      V1[3]))) + (P3[1] * (T3[2] * (V2[5] * V1[4] - V2[4] * V1[5]) + (T3[4] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + T3[5] * (V2[4] * V1[2] - V2[2] *
      V1[4]))) + (P3[2] * (T3[2] * (V2[3] * V1[5] - V2[5] * V1[3]) + (T3[3] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + T3[5] * (V2[2] * V1[3] - V2[3] *
      V1[2]))) + P3[3] * (T3[2] * (V2[4] * V1[3] - V2[3] * V1[4]) + (T3[3] *
      (V2[2] * V1[4] - V2[4] * V1[2]) + T3[4] * (V2[3] * V1[2] - V2[2] *
      V1[3])))))) + (P2[1] * (P3[0] * (T3[7] * (V2[5] * V1[4] - V2[4] * V1[5])
      + (T3[8] * (V2[3] * V1[5] - V2[5] * V1[3]) + T3[9] * (V2[4] * V1[3] -
      V2[3] * V1[4]))) + (P3[1] * (T3[6] * (V2[4] * V1[5] - V2[5] * V1[4]) +
      (T3[8] * (V2[5] * V1[2] - V2[2] * V1[5]) + T3[9] * (V2[2] * V1[4] - V2[4]
      * V1[2]))) + (P3[2] * (T3[6] * (V2[5] * V1[3] - V2[3] * V1[5]) + (T3[7] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + T3[9] * (V2[3] * V1[2] - V2[2] *
      V1[3]))) + P3[3] * (T3[6] * (V2[3] * V1[4] - V2[4] * V1[3]) + (T3[7] *
      (V2[4] * V1[2] - V2[2] * V1[4]) + T3[8] * (V2[2] * V1[3] - V2[3] *
      V1[2])))))) + (P2[2] * (P3[0] * (T3[11] * (V2[5] * V1[4] - V2[4] * V1[5])
      + (T3[12] * (V2[3] * V1[5] - V2[5] * V1[3]) + T3[13] * (V2[4] * V1[3] -
      V2[3] * V1[4]))) + (P3[1] * (T3[12] * (V2[5] * V1[2] - V2[2] * V1[5]) +
      (T3[13] * (V2[2] * V1[4] - V2[4] * V1[2]) + T3[10] * (V2[4] * V1[5] -
      V2[5] * V1[4]))) + (P3[2] * (T3[11] * (V2[2] * V1[5] - V2[5] * V1[2]) +
      (T3[13] * (V2[3] * V1[2] - V2[2] * V1[3]) + T3[10] * (V2[5] * V1[3] -
      V2[3] * V1[5]))) + P3[3] * (T3[11] * (V2[4] * V1[2] - V2[2] * V1[4]) +
      (T3[12] * (V2[2] * V1[3] - V2[3] * V1[2]) + T3[10] * (V2[3] * V1[4] -
      V2[4] * V1[3])))))) + P2[3] * (P3[0] * (T3[15] * (V2[5] * V1[4] - V2[4] *
      V1[5]) + (T3[16] * (V2[3] * V1[5] - V2[5] * V1[3]) + T3[17] * (V2[4] *
      V1[3] - V2[3] * V1[4]))) + (P3[1] * (T3[14] * (V2[4] * V1[5] - V2[5] *
      V1[4]) + (T3[16] * (V2[5] * V1[2] - V2[2] * V1[5]) + T3[17] * (V2[2] *
      V1[4] - V2[4] * V1[2]))) + (P3[2] * (T3[14] * (V2[5] * V1[3] - V2[3] *
      V1[5]) + (T3[15] * (V2[2] * V1[5] - V2[5] * V1[2]) + T3[17] * (V2[3] *
      V1[2] - V2[2] * V1[3]))) + P3[3] * (T3[14] * (V2[3] * V1[4] - V2[4] *
      V1[3]) + (T3[15] * (V2[4] * V1[2] - V2[2] * V1[4]) + T3[16] * (V2[2] *
      V1[3] - V2[3] * V1[2])))))))));
  TMP47 = (P1[0] * (P3[0] * (T3[14] * (V2[3] * V1[4] - V2[4] * V1[3]) + (T3[6]
      * (V2[4] * V1[5] - V2[5] * V1[4]) + T3[10] * (V2[5] * V1[3] - V2[3] *
      V1[5]))) + (P3[1] * (T3[2] * (V2[5] * V1[4] - V2[4] * V1[5]) + (T3[14] *
      (V2[4] * V1[2] - V2[2] * V1[4]) + T3[10] * (V2[2] * V1[5] - V2[5] *
      V1[2]))) + (P3[2] * (T3[2] * (V2[3] * V1[5] - V2[5] * V1[3]) + (T3[14] *
      (V2[2] * V1[3] - V2[3] * V1[2]) + T3[6] * (V2[5] * V1[2] - V2[2] *
      V1[5]))) + P3[3] * (T3[2] * (V2[4] * V1[3] - V2[3] * V1[4]) + (T3[6] *
      (V2[2] * V1[4] - V2[4] * V1[2]) + T3[10] * (V2[3] * V1[2] - V2[2] *
      V1[3])))))) + (P1[1] * (P3[0] * (T3[11] * (V2[3] * V1[5] - V2[5] * V1[3])
      + (T3[15] * (V2[4] * V1[3] - V2[3] * V1[4]) + T3[7] * (V2[5] * V1[4] -
      V2[4] * V1[5]))) + (P3[1] * (T3[11] * (V2[5] * V1[2] - V2[2] * V1[5]) +
      (T3[15] * (V2[2] * V1[4] - V2[4] * V1[2]) + T3[3] * (V2[4] * V1[5] -
      V2[5] * V1[4]))) + (P3[2] * (T3[15] * (V2[3] * V1[2] - V2[2] * V1[3]) +
      (T3[3] * (V2[5] * V1[3] - V2[3] * V1[5]) + T3[7] * (V2[2] * V1[5] - V2[5]
      * V1[2]))) + P3[3] * (T3[11] * (V2[2] * V1[3] - V2[3] * V1[2]) + (T3[3] *
      (V2[3] * V1[4] - V2[4] * V1[3]) + T3[7] * (V2[4] * V1[2] - V2[2] *
      V1[4])))))) + (P1[2] * (P3[0] * (T3[12] * (V2[3] * V1[5] - V2[5] * V1[3])
      + (T3[16] * (V2[4] * V1[3] - V2[3] * V1[4]) + T3[8] * (V2[5] * V1[4] -
      V2[4] * V1[5]))) + (P3[1] * (T3[12] * (V2[5] * V1[2] - V2[2] * V1[5]) +
      (T3[16] * (V2[2] * V1[4] - V2[4] * V1[2]) + T3[4] * (V2[4] * V1[5] -
      V2[5] * V1[4]))) + (P3[2] * (T3[16] * (V2[3] * V1[2] - V2[2] * V1[3]) +
      (T3[4] * (V2[5] * V1[3] - V2[3] * V1[5]) + T3[8] * (V2[2] * V1[5] - V2[5]
      * V1[2]))) + P3[3] * (T3[12] * (V2[2] * V1[3] - V2[3] * V1[2]) + (T3[4] *
      (V2[3] * V1[4] - V2[4] * V1[3]) + T3[8] * (V2[4] * V1[2] - V2[2] *
      V1[4])))))) + P1[3] * (P3[0] * (T3[13] * (V2[3] * V1[5] - V2[5] * V1[3])
      + (T3[17] * (V2[4] * V1[3] - V2[3] * V1[4]) + T3[9] * (V2[5] * V1[4] -
      V2[4] * V1[5]))) + (P3[1] * (T3[13] * (V2[5] * V1[2] - V2[2] * V1[5]) +
      (T3[17] * (V2[2] * V1[4] - V2[4] * V1[2]) + T3[5] * (V2[4] * V1[5] -
      V2[5] * V1[4]))) + (P3[2] * (T3[17] * (V2[3] * V1[2] - V2[2] * V1[3]) +
      (T3[5] * (V2[5] * V1[3] - V2[3] * V1[5]) + T3[9] * (V2[2] * V1[5] - V2[5]
      * V1[2]))) + P3[3] * (T3[13] * (V2[2] * V1[3] - V2[3] * V1[2]) + (T3[5] *
      (V2[3] * V1[4] - V2[4] * V1[3]) + T3[9] * (V2[4] * V1[2] - V2[2] *
      V1[4])))))))));
  TMP45 = (P1[0] * (P3[0] * (T3[3] * (V2[4] * V1[5] - V2[5] * V1[4]) + (T3[4] *
      (V2[5] * V1[3] - V2[3] * V1[5]) + T3[5] * (V2[3] * V1[4] - V2[4] *
      V1[3]))) + (P3[1] * (T3[2] * (V2[5] * V1[4] - V2[4] * V1[5]) + (T3[4] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + T3[5] * (V2[4] * V1[2] - V2[2] *
      V1[4]))) + (P3[2] * (T3[2] * (V2[3] * V1[5] - V2[5] * V1[3]) + (T3[3] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + T3[5] * (V2[2] * V1[3] - V2[3] *
      V1[2]))) + P3[3] * (T3[2] * (V2[4] * V1[3] - V2[3] * V1[4]) + (T3[3] *
      (V2[2] * V1[4] - V2[4] * V1[2]) + T3[4] * (V2[3] * V1[2] - V2[2] *
      V1[3])))))) + (P1[1] * (P3[0] * (T3[7] * (V2[5] * V1[4] - V2[4] * V1[5])
      + (T3[8] * (V2[3] * V1[5] - V2[5] * V1[3]) + T3[9] * (V2[4] * V1[3] -
      V2[3] * V1[4]))) + (P3[1] * (T3[6] * (V2[4] * V1[5] - V2[5] * V1[4]) +
      (T3[8] * (V2[5] * V1[2] - V2[2] * V1[5]) + T3[9] * (V2[2] * V1[4] - V2[4]
      * V1[2]))) + (P3[2] * (T3[6] * (V2[5] * V1[3] - V2[3] * V1[5]) + (T3[7] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + T3[9] * (V2[3] * V1[2] - V2[2] *
      V1[3]))) + P3[3] * (T3[6] * (V2[3] * V1[4] - V2[4] * V1[3]) + (T3[7] *
      (V2[4] * V1[2] - V2[2] * V1[4]) + T3[8] * (V2[2] * V1[3] - V2[3] *
      V1[2])))))) + (P1[2] * (P3[0] * (T3[11] * (V2[5] * V1[4] - V2[4] * V1[5])
      + (T3[12] * (V2[3] * V1[5] - V2[5] * V1[3]) + T3[13] * (V2[4] * V1[3] -
      V2[3] * V1[4]))) + (P3[1] * (T3[12] * (V2[5] * V1[2] - V2[2] * V1[5]) +
      (T3[13] * (V2[2] * V1[4] - V2[4] * V1[2]) + T3[10] * (V2[4] * V1[5] -
      V2[5] * V1[4]))) + (P3[2] * (T3[11] * (V2[2] * V1[5] - V2[5] * V1[2]) +
      (T3[13] * (V2[3] * V1[2] - V2[2] * V1[3]) + T3[10] * (V2[5] * V1[3] -
      V2[3] * V1[5]))) + P3[3] * (T3[11] * (V2[4] * V1[2] - V2[2] * V1[4]) +
      (T3[12] * (V2[2] * V1[3] - V2[3] * V1[2]) + T3[10] * (V2[3] * V1[4] -
      V2[4] * V1[3])))))) + P1[3] * (P3[0] * (T3[15] * (V2[5] * V1[4] - V2[4] *
      V1[5]) + (T3[16] * (V2[3] * V1[5] - V2[5] * V1[3]) + T3[17] * (V2[4] *
      V1[3] - V2[3] * V1[4]))) + (P3[1] * (T3[14] * (V2[4] * V1[5] - V2[5] *
      V1[4]) + (T3[16] * (V2[5] * V1[2] - V2[2] * V1[5]) + T3[17] * (V2[2] *
      V1[4] - V2[4] * V1[2]))) + (P3[2] * (T3[14] * (V2[5] * V1[3] - V2[3] *
      V1[5]) + (T3[15] * (V2[2] * V1[5] - V2[5] * V1[2]) + T3[17] * (V2[3] *
      V1[2] - V2[2] * V1[3]))) + P3[3] * (T3[14] * (V2[3] * V1[4] - V2[4] *
      V1[3]) + (T3[15] * (V2[4] * V1[2] - V2[2] * V1[4]) + T3[16] * (V2[2] *
      V1[3] - V2[3] * V1[2])))))))));
  vertex = COUP * (-cI * (TMP45 + TMP47) + cI * (TMP46 + TMP48)); 
}


void VVT11_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP2; 
  complex<double> TMP1; 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP26; 
  complex<double> TMP25; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  TMP25 = (P1[0] * - 1. * (P2[1] * T3[6] + P2[2] * T3[10] + P2[3] * T3[14] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[11] + P2[3] *
      T3[15] - P2[0] * T3[3]) + (P1[2] * (P2[1] * T3[8] + P2[2] * T3[12] +
      P2[3] * T3[16] - P2[0] * T3[4]) + P1[3] * (P2[1] * T3[9] + P2[2] * T3[13]
      + P2[3] * T3[17] - P2[0] * T3[5]))));
  TMP26 = (P1[0] * - 1. * (P2[1] * T3[3] + P2[2] * T3[4] + P2[3] * T3[5] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[8] + P2[3] * T3[9]
      - P2[0] * T3[6]) + (P1[2] * (P2[1] * T3[11] + P2[2] * T3[12] + P2[3] *
      T3[13] - P2[0] * T3[10]) + P1[3] * (P2[1] * T3[15] + P2[2] * T3[16] +
      P2[3] * T3[17] - P2[0] * T3[14]))));
  TMP1 = -1. * (P1[0] * (P2[1] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P2[2] *
      (V2[5] * V1[3] - V2[3] * V1[5]) + P2[3] * (V2[3] * V1[4] - V2[4] *
      V1[3]))) + (P1[1] * (P2[0] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P2[2] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P2[3] * (V2[4] * V1[2] - V2[2] *
      V1[4]))) + (P1[2] * (P2[0] * (V2[3] * V1[5] - V2[5] * V1[3]) + (P2[1] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P2[3] * (V2[2] * V1[3] - V2[3] *
      V1[2]))) + P1[3] * (P2[0] * (V2[4] * V1[3] - V2[3] * V1[4]) + (P2[1] *
      (V2[2] * V1[4] - V2[4] * V1[2]) + P2[2] * (V2[3] * V1[2] - V2[2] *
      V1[3]))))));
  TMP2 = -1. * (P1[0] * (P2[1] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P2[2] *
      (V2[3] * V1[5] - V2[5] * V1[3]) + P2[3] * (V2[4] * V1[3] - V2[3] *
      V1[4]))) + (P1[1] * (P2[0] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P2[2] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P2[3] * (V2[2] * V1[4] - V2[4] *
      V1[2]))) + (P1[2] * (P2[0] * (V2[5] * V1[3] - V2[3] * V1[5]) + (P2[1] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P2[3] * (V2[3] * V1[2] - V2[2] *
      V1[3]))) + P1[3] * (P2[0] * (V2[3] * V1[4] - V2[4] * V1[3]) + (P2[1] *
      (V2[4] * V1[2] - V2[2] * V1[4]) + P2[2] * (V2[2] * V1[3] - V2[3] *
      V1[2]))))));
  vertex = COUP * (TMP1 * (+cI * (TMP25 + TMP26)) - TMP2 * (+cI * (TMP25 +
      TMP26)));
}

void VVT11_12_15_16_17_18_19_20_21_22_0(complex<double> V1[], complex<double>
    V2[], complex<double> T3[], complex<double> COUP1, complex<double> COUP2,
    complex<double> COUP3, complex<double> COUP4, complex<double> COUP5,
    complex<double> COUP6, complex<double> COUP7, complex<double> COUP8,
    complex<double> COUP9, complex<double> COUP10, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  complex<double> tmp; 
  VVT11_0(V1, V2, T3, COUP1, vertex); 
  VVT12_0(V1, V2, T3, COUP2, tmp); 
  vertex = vertex + tmp; 
  VVT15_0(V1, V2, T3, COUP3, tmp); 
  vertex = vertex + tmp; 
  VVT16_0(V1, V2, T3, COUP4, tmp); 
  vertex = vertex + tmp; 
  VVT17_0(V1, V2, T3, COUP5, tmp); 
  vertex = vertex + tmp; 
  VVT18_0(V1, V2, T3, COUP6, tmp); 
  vertex = vertex + tmp; 
  VVT19_0(V1, V2, T3, COUP7, tmp); 
  vertex = vertex + tmp; 
  VVT20_0(V1, V2, T3, COUP8, tmp); 
  vertex = vertex + tmp; 
  VVT21_0(V1, V2, T3, COUP9, tmp); 
  vertex = vertex + tmp; 
  VVT22_0(V1, V2, T3, COUP10, tmp); 
  vertex = vertex + tmp; 
}

void VVT22_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP17; 
  double P3[4]; 
  complex<double> TMP20; 
  complex<double> TMP6; 
  complex<double> TMP21; 
  complex<double> TMP5; 
  complex<double> TMP26; 
  complex<double> TMP4; 
  complex<double> TMP18; 
  complex<double> TMP16; 
  complex<double> TMP25; 
  complex<double> TMP19; 
  complex<double> TMP8; 
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
  TMP25 = (P1[0] * - 1. * (P2[1] * T3[6] + P2[2] * T3[10] + P2[3] * T3[14] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[11] + P2[3] *
      T3[15] - P2[0] * T3[3]) + (P1[2] * (P2[1] * T3[8] + P2[2] * T3[12] +
      P2[3] * T3[16] - P2[0] * T3[4]) + P1[3] * (P2[1] * T3[9] + P2[2] * T3[13]
      + P2[3] * T3[17] - P2[0] * T3[5]))));
  TMP26 = (P1[0] * - 1. * (P2[1] * T3[3] + P2[2] * T3[4] + P2[3] * T3[5] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[8] + P2[3] * T3[9]
      - P2[0] * T3[6]) + (P1[2] * (P2[1] * T3[11] + P2[2] * T3[12] + P2[3] *
      T3[13] - P2[0] * T3[10]) + P1[3] * (P2[1] * T3[15] + P2[2] * T3[16] +
      P2[3] * T3[17] - P2[0] * T3[14]))));
  TMP6 = (V2[2] * P3[0] - V2[3] * P3[1] - V2[4] * P3[2] - V2[5] * P3[3]); 
  TMP20 = (P1[0] * - 1. * (V2[3] * T3[6] + V2[4] * T3[10] + V2[5] * T3[14] -
      V2[2] * T3[2]) + (P1[1] * (V2[3] * T3[7] + V2[4] * T3[11] + V2[5] *
      T3[15] - V2[2] * T3[3]) + (P1[2] * (V2[3] * T3[8] + V2[4] * T3[12] +
      V2[5] * T3[16] - V2[2] * T3[4]) + P1[3] * (V2[3] * T3[9] + V2[4] * T3[13]
      + V2[5] * T3[17] - V2[2] * T3[5]))));
  TMP21 = (P1[0] * - 1. * (V2[3] * T3[3] + V2[4] * T3[4] + V2[5] * T3[5] -
      V2[2] * T3[2]) + (P1[1] * (V2[3] * T3[7] + V2[4] * T3[8] + V2[5] * T3[9]
      - V2[2] * T3[6]) + (P1[2] * (V2[3] * T3[11] + V2[4] * T3[12] + V2[5] *
      T3[13] - V2[2] * T3[10]) + P1[3] * (V2[3] * T3[15] + V2[4] * T3[16] +
      V2[5] * T3[17] - V2[2] * T3[14]))));
  TMP19 = (P2[0] * - 1. * (V1[3] * T3[3] + V1[4] * T3[4] + V1[5] * T3[5] -
      V1[2] * T3[2]) + (P2[1] * (V1[3] * T3[7] + V1[4] * T3[8] + V1[5] * T3[9]
      - V1[2] * T3[6]) + (P2[2] * (V1[3] * T3[11] + V1[4] * T3[12] + V1[5] *
      T3[13] - V1[2] * T3[10]) + P2[3] * (V1[3] * T3[15] + V1[4] * T3[16] +
      V1[5] * T3[17] - V1[2] * T3[14]))));
  TMP18 = (P2[0] * - 1. * (V1[3] * T3[6] + V1[4] * T3[10] + V1[5] * T3[14] -
      V1[2] * T3[2]) + (P2[1] * (V1[3] * T3[7] + V1[4] * T3[11] + V1[5] *
      T3[15] - V1[2] * T3[3]) + (P2[2] * (V1[3] * T3[8] + V1[4] * T3[12] +
      V1[5] * T3[16] - V1[2] * T3[4]) + P2[3] * (V1[3] * T3[9] + V1[4] * T3[13]
      + V1[5] * T3[17] - V1[2] * T3[5]))));
  TMP5 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  TMP4 = (V1[2] * P3[0] - V1[3] * P3[1] - V1[4] * P3[2] - V1[5] * P3[3]); 
  TMP17 = (V1[2] * - 1. * (V2[3] * T3[3] + V2[4] * T3[4] + V2[5] * T3[5] -
      V2[2] * T3[2]) + (V1[3] * (V2[3] * T3[7] + V2[4] * T3[8] + V2[5] * T3[9]
      - V2[2] * T3[6]) + (V1[4] * (V2[3] * T3[11] + V2[4] * T3[12] + V2[5] *
      T3[13] - V2[2] * T3[10]) + V1[5] * (V2[3] * T3[15] + V2[4] * T3[16] +
      V2[5] * T3[17] - V2[2] * T3[14]))));
  TMP16 = (V1[2] * - 1. * (V2[3] * T3[6] + V2[4] * T3[10] + V2[5] * T3[14] -
      V2[2] * T3[2]) + (V1[3] * (V2[3] * T3[7] + V2[4] * T3[11] + V2[5] *
      T3[15] - V2[2] * T3[3]) + (V1[4] * (V2[3] * T3[8] + V2[4] * T3[12] +
      V2[5] * T3[16] - V2[2] * T3[4]) + V1[5] * (V2[3] * T3[9] + V2[4] * T3[13]
      + V2[5] * T3[17] - V2[2] * T3[5]))));
  TMP8 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  vertex = COUP * (TMP4 * (TMP5 * (+cI * (TMP20 + TMP21)) - TMP6 * (+cI *
      (TMP25 + TMP26))) + TMP8 * (TMP5 * - 1. * (+cI * (TMP16 + TMP17)) + TMP6
      * (+cI * (TMP18 + TMP19))));
}


void VVT20_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP0; 
  double P3[4]; 
  complex<double> TMP6; 
  complex<double> denom; 
  double OM3; 
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
  TMP4 = (V1[2] * P3[0] - V1[3] * P3[1] - V1[4] * P3[2] - V1[5] * P3[3]); 
  TMP6 = (V2[2] * P3[0] - V2[3] * P3[1] - V2[4] * P3[2] - V2[5] * P3[3]); 
  TMP0 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * (OM3 * (P3[0] * (P3[0] * 0.666666667 * (+cI * (TMP0) + 2. *
      cI * (TMP4 * TMP6 * OM3)) + (-2. * cI * (V2[2] * TMP4 + V1[2] * TMP6))) +
      0.666666667 * cI * (TMP4 * TMP6)) + (-0.666666667 * cI * (TMP0) + 2. * cI
      * (V2[2] * V1[2])));
  T3[6] = denom * (OM3 * (P3[0] * (P3[1] * 0.666666667 * (+cI * (TMP0) + 2. *
      cI * (TMP4 * TMP6 * OM3)) + (-cI * (V2[3] * TMP4 + V1[3] * TMP6))) -
      P3[1] * (+cI * (V1[2] * TMP6 + V2[2] * TMP4))) + (+cI * (V2[3] * V1[2] +
      V2[2] * V1[3])));
  T3[10] = denom * (OM3 * (P3[0] * (P3[2] * 0.666666667 * (+cI * (TMP0) + 2. *
      cI * (TMP4 * TMP6 * OM3)) + (-cI * (V2[4] * TMP4 + V1[4] * TMP6))) -
      P3[2] * (+cI * (V1[2] * TMP6 + V2[2] * TMP4))) + (+cI * (V2[4] * V1[2] +
      V2[2] * V1[4])));
  T3[14] = denom * (OM3 * (P3[0] * (P3[3] * 0.666666667 * (+cI * (TMP0) + 2. *
      cI * (TMP4 * TMP6 * OM3)) + (-cI * (V2[5] * TMP4 + V1[5] * TMP6))) -
      P3[3] * (+cI * (V1[2] * TMP6 + V2[2] * TMP4))) + (+cI * (V2[5] * V1[2] +
      V2[2] * V1[5])));
  T3[3] = denom * (OM3 * (P3[0] * (P3[1] * 0.666666667 * (+cI * (TMP0) + 2. *
      cI * (TMP4 * TMP6 * OM3)) + (-cI * (V1[3] * TMP6 + V2[3] * TMP4))) -
      P3[1] * (+cI * (V2[2] * TMP4 + V1[2] * TMP6))) + (+cI * (V2[2] * V1[3] +
      V2[3] * V1[2])));
  T3[7] = denom * (OM3 * (P3[1] * (P3[1] * 0.666666667 * (+cI * (TMP0) + 2. *
      cI * (TMP4 * TMP6 * OM3)) + (-2. * cI * (V2[3] * TMP4 + V1[3] * TMP6))) -
      0.666666667 * cI * (TMP4 * TMP6)) + (+2. * cI * (V2[3] * V1[3]) +
      0.666666667 * cI * (TMP0)));
  T3[11] = denom * (OM3 * (P3[1] * (P3[2] * 0.666666667 * (+cI * (TMP0) + 2. *
      cI * (TMP4 * TMP6 * OM3)) + (-cI * (V2[4] * TMP4 + V1[4] * TMP6))) -
      P3[2] * (+cI * (V1[3] * TMP6 + V2[3] * TMP4))) + (+cI * (V2[4] * V1[3] +
      V2[3] * V1[4])));
  T3[15] = denom * (OM3 * (P3[1] * (P3[3] * 0.666666667 * (+cI * (TMP0) + 2. *
      cI * (TMP4 * TMP6 * OM3)) + (-cI * (V2[5] * TMP4 + V1[5] * TMP6))) -
      P3[3] * (+cI * (V1[3] * TMP6 + V2[3] * TMP4))) + (+cI * (V2[5] * V1[3] +
      V2[3] * V1[5])));
  T3[4] = denom * (OM3 * (P3[0] * (P3[2] * 0.666666667 * (+cI * (TMP0) + 2. *
      cI * (TMP4 * TMP6 * OM3)) + (-cI * (V1[4] * TMP6 + V2[4] * TMP4))) -
      P3[2] * (+cI * (V2[2] * TMP4 + V1[2] * TMP6))) + (+cI * (V2[2] * V1[4] +
      V2[4] * V1[2])));
  T3[8] = denom * (OM3 * (P3[1] * (P3[2] * 0.666666667 * (+cI * (TMP0) + 2. *
      cI * (TMP4 * TMP6 * OM3)) + (-cI * (V1[4] * TMP6 + V2[4] * TMP4))) -
      P3[2] * (+cI * (V2[3] * TMP4 + V1[3] * TMP6))) + (+cI * (V2[3] * V1[4] +
      V2[4] * V1[3])));
  T3[12] = denom * (OM3 * (P3[2] * (P3[2] * 0.666666667 * (+cI * (TMP0) + 2. *
      cI * (TMP4 * TMP6 * OM3)) + (-2. * cI * (V2[4] * TMP4 + V1[4] * TMP6))) -
      0.666666667 * cI * (TMP4 * TMP6)) + (+2. * cI * (V2[4] * V1[4]) +
      0.666666667 * cI * (TMP0)));
  T3[16] = denom * (OM3 * (P3[2] * (P3[3] * 0.666666667 * (+cI * (TMP0) + 2. *
      cI * (TMP4 * TMP6 * OM3)) + (-cI * (V2[5] * TMP4 + V1[5] * TMP6))) -
      P3[3] * (+cI * (V1[4] * TMP6 + V2[4] * TMP4))) + (+cI * (V2[5] * V1[4] +
      V2[4] * V1[5])));
  T3[5] = denom * (OM3 * (P3[0] * (P3[3] * 0.666666667 * (+cI * (TMP0) + 2. *
      cI * (TMP4 * TMP6 * OM3)) + (-cI * (V1[5] * TMP6 + V2[5] * TMP4))) -
      P3[3] * (+cI * (V2[2] * TMP4 + V1[2] * TMP6))) + (+cI * (V2[2] * V1[5] +
      V2[5] * V1[2])));
  T3[9] = denom * (OM3 * (P3[1] * (P3[3] * 0.666666667 * (+cI * (TMP0) + 2. *
      cI * (TMP4 * TMP6 * OM3)) + (-cI * (V1[5] * TMP6 + V2[5] * TMP4))) -
      P3[3] * (+cI * (V2[3] * TMP4 + V1[3] * TMP6))) + (+cI * (V2[3] * V1[5] +
      V2[5] * V1[3])));
  T3[13] = denom * (OM3 * (P3[2] * (P3[3] * 0.666666667 * (+cI * (TMP0) + 2. *
      cI * (TMP4 * TMP6 * OM3)) + (-cI * (V1[5] * TMP6 + V2[5] * TMP4))) -
      P3[3] * (+cI * (V2[4] * TMP4 + V1[4] * TMP6))) + (+cI * (V2[4] * V1[5] +
      V2[5] * V1[4])));
  T3[17] = denom * (OM3 * (P3[3] * (P3[3] * 0.666666667 * (+cI * (TMP0) + 2. *
      cI * (TMP4 * TMP6 * OM3)) + (-2. * cI * (V2[5] * TMP4 + V1[5] * TMP6))) -
      0.666666667 * cI * (TMP4 * TMP6)) + (+2. * cI * (V2[5] * V1[5]) +
      0.666666667 * cI * (TMP0)));
}


void VVT18_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> TMP22; 
  double P2[4]; 
  complex<double> TMP23; 
  complex<double> TMP7; 
  double P3[4]; 
  complex<double> TMP6; 
  complex<double> denom; 
  double OM3; 
  complex<double> TMP4; 
  complex<double> TMP3; 
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
  TMP22 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  TMP23 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP4 = (V1[2] * P3[0] - V1[3] * P3[1] - V1[4] * P3[2] - V1[5] * P3[3]); 
  TMP7 = (V1[2] * P2[0] - V1[3] * P2[1] - V1[4] * P2[2] - V1[5] * P2[3]); 
  TMP6 = (V2[2] * P3[0] - V2[3] * P3[1] - V2[4] * P3[2] - V2[5] * P3[3]); 
  TMP3 = (V2[2] * P1[0] - V2[3] * P1[1] - V2[4] * P1[2] - V2[5] * P1[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * (OM3 * (TMP4 * (TMP6 * (P3[0] * (OM3 * 1.333333333 * P3[0] *
      (+cI * (TMP22 + TMP23)) + (-2. * cI * (P2[0] + P1[0]))) + (+0.666666667 *
      cI * (TMP22 + TMP23))) + 2. * (P3[0] * (-cI * (V2[2] * TMP23) +
      0.333333333 * cI * (P3[0] * TMP3)))) + 2. * (P3[0] * TMP6 * (-cI * (V1[2]
      * TMP22) + 0.333333333 * cI * (P3[0] * TMP7)))) + (TMP4 * 2. *
      (-0.333333333 * cI * (TMP3) + cI * (V2[2] * P1[0])) + 2. * (TMP6 *
      (-0.333333333 * cI * (TMP7) + cI * (V1[2] * P2[0])))));
  T3[3] = denom * (OM3 * (TMP4 * (TMP6 * (P3[0] * (OM3 * 1.333333333 * P3[1] *
      (+cI * (TMP22 + TMP23)) + (-cI * (P2[1] + P1[1]))) - P3[1] * (+cI *
      (P2[0] + P1[0]))) + (P3[0] * (-cI * (V2[3] * TMP23) + 0.666666667 * cI *
      (P3[1] * TMP3)) - cI * (V2[2] * P3[1] * TMP23))) + TMP6 * (P3[0] * (-cI *
      (V1[3] * TMP22) + 0.666666667 * cI * (P3[1] * TMP7)) - cI * (V1[2] *
      P3[1] * TMP22))) + (TMP4 * (+cI * (V2[2] * P1[1] + V2[3] * P1[0])) + TMP6
      * (+cI * (V1[2] * P2[1] + V1[3] * P2[0]))));
  T3[4] = denom * (OM3 * (TMP4 * (TMP6 * (P3[0] * (OM3 * 1.333333333 * P3[2] *
      (+cI * (TMP22 + TMP23)) + (-cI * (P2[2] + P1[2]))) - P3[2] * (+cI *
      (P2[0] + P1[0]))) + (P3[0] * (-cI * (V2[4] * TMP23) + 0.666666667 * cI *
      (P3[2] * TMP3)) - cI * (V2[2] * P3[2] * TMP23))) + TMP6 * (P3[0] * (-cI *
      (V1[4] * TMP22) + 0.666666667 * cI * (P3[2] * TMP7)) - cI * (V1[2] *
      P3[2] * TMP22))) + (TMP4 * (+cI * (V2[2] * P1[2] + V2[4] * P1[0])) + TMP6
      * (+cI * (V1[2] * P2[2] + V1[4] * P2[0]))));
  T3[5] = denom * (OM3 * (TMP4 * (TMP6 * (P3[0] * (OM3 * 1.333333333 * P3[3] *
      (+cI * (TMP22 + TMP23)) + (-cI * (P2[3] + P1[3]))) - P3[3] * (+cI *
      (P2[0] + P1[0]))) + (P3[0] * (-cI * (V2[5] * TMP23) + 0.666666667 * cI *
      (P3[3] * TMP3)) - cI * (V2[2] * P3[3] * TMP23))) + TMP6 * (P3[0] * (-cI *
      (V1[5] * TMP22) + 0.666666667 * cI * (P3[3] * TMP7)) - cI * (V1[2] *
      P3[3] * TMP22))) + (TMP4 * (+cI * (V2[2] * P1[3] + V2[5] * P1[0])) + TMP6
      * (+cI * (V1[2] * P2[3] + V1[5] * P2[0]))));
  T3[6] = denom * (OM3 * (TMP4 * (TMP6 * (P3[0] * (OM3 * 1.333333333 * P3[1] *
      (+cI * (TMP22 + TMP23)) + (-cI * (P2[1] + P1[1]))) - P3[1] * (+cI *
      (P2[0] + P1[0]))) + (P3[0] * (-cI * (V2[3] * TMP23) + 0.666666667 * cI *
      (P3[1] * TMP3)) - cI * (V2[2] * P3[1] * TMP23))) + TMP6 * (P3[0] * (-cI *
      (V1[3] * TMP22) + 0.666666667 * cI * (P3[1] * TMP7)) - cI * (V1[2] *
      P3[1] * TMP22))) + (TMP4 * (+cI * (V2[3] * P1[0] + V2[2] * P1[1])) + TMP6
      * (+cI * (V1[3] * P2[0] + V1[2] * P2[1]))));
  T3[7] = denom * (OM3 * (TMP4 * (TMP6 * (P3[1] * (OM3 * 1.333333333 * P3[1] *
      (+cI * (TMP22 + TMP23)) + (-2. * cI * (P2[1] + P1[1]))) + (-0.666666667 *
      cI * (TMP22 + TMP23))) + 2. * (P3[1] * (-cI * (V2[3] * TMP23) +
      0.333333333 * cI * (P3[1] * TMP3)))) + 2. * (P3[1] * TMP6 * (-cI * (V1[3]
      * TMP22) + 0.333333333 * cI * (P3[1] * TMP7)))) + (TMP4 * 2. * (+cI *
      (V2[3] * P1[1]) + 0.333333333 * cI * (TMP3)) + 2. * (TMP6 * (+cI * (V1[3]
      * P2[1]) + 0.333333333 * cI * (TMP7)))));
  T3[8] = denom * (OM3 * (TMP4 * (TMP6 * (P3[1] * (OM3 * 1.333333333 * P3[2] *
      (+cI * (TMP22 + TMP23)) + (-cI * (P2[2] + P1[2]))) - P3[2] * (+cI *
      (P2[1] + P1[1]))) + (P3[1] * (-cI * (V2[4] * TMP23) + 0.666666667 * cI *
      (P3[2] * TMP3)) - cI * (V2[3] * P3[2] * TMP23))) + TMP6 * (P3[1] * (-cI *
      (V1[4] * TMP22) + 0.666666667 * cI * (P3[2] * TMP7)) - cI * (V1[3] *
      P3[2] * TMP22))) + (TMP4 * (+cI * (V2[3] * P1[2] + V2[4] * P1[1])) + TMP6
      * (+cI * (V1[3] * P2[2] + V1[4] * P2[1]))));
  T3[9] = denom * (OM3 * (TMP4 * (TMP6 * (P3[1] * (OM3 * 1.333333333 * P3[3] *
      (+cI * (TMP22 + TMP23)) + (-cI * (P2[3] + P1[3]))) - P3[3] * (+cI *
      (P2[1] + P1[1]))) + (P3[1] * (-cI * (V2[5] * TMP23) + 0.666666667 * cI *
      (P3[3] * TMP3)) - cI * (V2[3] * P3[3] * TMP23))) + TMP6 * (P3[1] * (-cI *
      (V1[5] * TMP22) + 0.666666667 * cI * (P3[3] * TMP7)) - cI * (V1[3] *
      P3[3] * TMP22))) + (TMP4 * (+cI * (V2[3] * P1[3] + V2[5] * P1[1])) + TMP6
      * (+cI * (V1[3] * P2[3] + V1[5] * P2[1]))));
  T3[10] = denom * (OM3 * (TMP4 * (TMP6 * (P3[0] * (OM3 * 1.333333333 * P3[2] *
      (+cI * (TMP22 + TMP23)) + (-cI * (P2[2] + P1[2]))) - P3[2] * (+cI *
      (P2[0] + P1[0]))) + (P3[0] * (-cI * (V2[4] * TMP23) + 0.666666667 * cI *
      (P3[2] * TMP3)) - cI * (V2[2] * P3[2] * TMP23))) + TMP6 * (P3[0] * (-cI *
      (V1[4] * TMP22) + 0.666666667 * cI * (P3[2] * TMP7)) - cI * (V1[2] *
      P3[2] * TMP22))) + (TMP4 * (+cI * (V2[4] * P1[0] + V2[2] * P1[2])) + TMP6
      * (+cI * (V1[4] * P2[0] + V1[2] * P2[2]))));
  T3[11] = denom * (OM3 * (TMP4 * (TMP6 * (P3[1] * (OM3 * 1.333333333 * P3[2] *
      (+cI * (TMP22 + TMP23)) + (-cI * (P2[2] + P1[2]))) - P3[2] * (+cI *
      (P2[1] + P1[1]))) + (P3[1] * (-cI * (V2[4] * TMP23) + 0.666666667 * cI *
      (P3[2] * TMP3)) - cI * (V2[3] * P3[2] * TMP23))) + TMP6 * (P3[1] * (-cI *
      (V1[4] * TMP22) + 0.666666667 * cI * (P3[2] * TMP7)) - cI * (V1[3] *
      P3[2] * TMP22))) + (TMP4 * (+cI * (V2[4] * P1[1] + V2[3] * P1[2])) + TMP6
      * (+cI * (V1[4] * P2[1] + V1[3] * P2[2]))));
  T3[12] = denom * (OM3 * (TMP4 * (TMP6 * (P3[2] * (OM3 * 1.333333333 * P3[2] *
      (+cI * (TMP22 + TMP23)) + (-2. * cI * (P2[2] + P1[2]))) + (-0.666666667 *
      cI * (TMP22 + TMP23))) + 2. * (P3[2] * (-cI * (V2[4] * TMP23) +
      0.333333333 * cI * (P3[2] * TMP3)))) + 2. * (P3[2] * TMP6 * (-cI * (V1[4]
      * TMP22) + 0.333333333 * cI * (P3[2] * TMP7)))) + (TMP4 * 2. * (+cI *
      (V2[4] * P1[2]) + 0.333333333 * cI * (TMP3)) + 2. * (TMP6 * (+cI * (V1[4]
      * P2[2]) + 0.333333333 * cI * (TMP7)))));
  T3[13] = denom * (OM3 * (TMP4 * (TMP6 * (P3[2] * (OM3 * 1.333333333 * P3[3] *
      (+cI * (TMP22 + TMP23)) + (-cI * (P2[3] + P1[3]))) - P3[3] * (+cI *
      (P2[2] + P1[2]))) + (P3[2] * (-cI * (V2[5] * TMP23) + 0.666666667 * cI *
      (P3[3] * TMP3)) - cI * (V2[4] * P3[3] * TMP23))) + TMP6 * (P3[2] * (-cI *
      (V1[5] * TMP22) + 0.666666667 * cI * (P3[3] * TMP7)) - cI * (V1[4] *
      P3[3] * TMP22))) + (TMP4 * (+cI * (V2[4] * P1[3] + V2[5] * P1[2])) + TMP6
      * (+cI * (V1[4] * P2[3] + V1[5] * P2[2]))));
  T3[14] = denom * (OM3 * (TMP4 * (TMP6 * (P3[0] * (OM3 * 1.333333333 * P3[3] *
      (+cI * (TMP22 + TMP23)) + (-cI * (P2[3] + P1[3]))) - P3[3] * (+cI *
      (P2[0] + P1[0]))) + (P3[0] * (-cI * (V2[5] * TMP23) + 0.666666667 * cI *
      (P3[3] * TMP3)) - cI * (V2[2] * P3[3] * TMP23))) + TMP6 * (P3[0] * (-cI *
      (V1[5] * TMP22) + 0.666666667 * cI * (P3[3] * TMP7)) - cI * (V1[2] *
      P3[3] * TMP22))) + (TMP4 * (+cI * (V2[5] * P1[0] + V2[2] * P1[3])) + TMP6
      * (+cI * (V1[5] * P2[0] + V1[2] * P2[3]))));
  T3[15] = denom * (OM3 * (TMP4 * (TMP6 * (P3[1] * (OM3 * 1.333333333 * P3[3] *
      (+cI * (TMP22 + TMP23)) + (-cI * (P2[3] + P1[3]))) - P3[3] * (+cI *
      (P2[1] + P1[1]))) + (P3[1] * (-cI * (V2[5] * TMP23) + 0.666666667 * cI *
      (P3[3] * TMP3)) - cI * (V2[3] * P3[3] * TMP23))) + TMP6 * (P3[1] * (-cI *
      (V1[5] * TMP22) + 0.666666667 * cI * (P3[3] * TMP7)) - cI * (V1[3] *
      P3[3] * TMP22))) + (TMP4 * (+cI * (V2[5] * P1[1] + V2[3] * P1[3])) + TMP6
      * (+cI * (V1[5] * P2[1] + V1[3] * P2[3]))));
  T3[16] = denom * (OM3 * (TMP4 * (TMP6 * (P3[2] * (OM3 * 1.333333333 * P3[3] *
      (+cI * (TMP22 + TMP23)) + (-cI * (P2[3] + P1[3]))) - P3[3] * (+cI *
      (P2[2] + P1[2]))) + (P3[2] * (-cI * (V2[5] * TMP23) + 0.666666667 * cI *
      (P3[3] * TMP3)) - cI * (V2[4] * P3[3] * TMP23))) + TMP6 * (P3[2] * (-cI *
      (V1[5] * TMP22) + 0.666666667 * cI * (P3[3] * TMP7)) - cI * (V1[4] *
      P3[3] * TMP22))) + (TMP4 * (+cI * (V2[5] * P1[2] + V2[4] * P1[3])) + TMP6
      * (+cI * (V1[5] * P2[2] + V1[4] * P2[3]))));
  T3[17] = denom * (OM3 * (TMP4 * (TMP6 * (P3[3] * (OM3 * 1.333333333 * P3[3] *
      (+cI * (TMP22 + TMP23)) + (-2. * cI * (P2[3] + P1[3]))) + (-0.666666667 *
      cI * (TMP22 + TMP23))) + 2. * (P3[3] * (-cI * (V2[5] * TMP23) +
      0.333333333 * cI * (P3[3] * TMP3)))) + 2. * (P3[3] * TMP6 * (-cI * (V1[5]
      * TMP22) + 0.333333333 * cI * (P3[3] * TMP7)))) + (TMP4 * 2. * (+cI *
      (V2[5] * P1[3]) + 0.333333333 * cI * (TMP3)) + 2. * (TMP6 * (+cI * (V1[5]
      * P2[3]) + 0.333333333 * cI * (TMP7)))));
}


void VVT21_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP33; 
  double P1[4]; 
  complex<double> TMP22; 
  complex<double> TMP0; 
  double P2[4]; 
  complex<double> TMP23; 
  complex<double> TMP7; 
  double P3[4]; 
  complex<double> TMP6; 
  complex<double> denom; 
  double OM3; 
  complex<double> TMP4; 
  complex<double> TMP9; 
  complex<double> TMP3; 
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
  TMP22 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  TMP23 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP9 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP33 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP4 = (V1[2] * P3[0] - V1[3] * P3[1] - V1[4] * P3[2] - V1[5] * P3[3]); 
  TMP7 = (V1[2] * P2[0] - V1[3] * P2[1] - V1[4] * P2[2] - V1[5] * P2[3]); 
  TMP6 = (V2[2] * P3[0] - V2[3] * P3[1] - V2[4] * P3[2] - V2[5] * P3[3]); 
  TMP0 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP3 = (V2[2] * P1[0] - V2[3] * P1[1] - V2[4] * P1[2] - V2[5] * P1[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * (OM3 * (P3[0] * (P3[0] * (OM3 * (TMP22 * 1.333333333 * (-cI *
      (TMP3 * TMP4) + cI * (TMP0 * TMP23)) + 1.333333333 * (TMP6 * (-cI * (TMP7
      * TMP23) + cI * (TMP4 * TMP9)))) + (TMP0 * 0.666666667 * (+cI * (TMP33 +
      TMP9)) - 1.333333333 * cI * (TMP3 * TMP7))) + (P1[0] * 2. * (-cI * (TMP0
      * TMP22) + cI * (TMP6 * TMP7)) + (P2[0] * 2. * (-cI * (TMP0 * TMP23) + cI
      * (TMP3 * TMP4)) + (TMP9 * - 2. * (+cI * (V2[2] * TMP4 + V1[2] * TMP6)) +
      (+2. * cI * (V1[2] * TMP3 * TMP22 + V2[2] * TMP7 * TMP23)))))) + (TMP22 *
      0.666666667 * (-cI * (TMP3 * TMP4) + cI * (TMP0 * TMP23)) + 0.666666667 *
      (TMP6 * (-cI * (TMP7 * TMP23) + cI * (TMP4 * TMP9))))) + (TMP0 * -
      0.666666667 * (+cI * (TMP33 + TMP9) - 3.000000000 * cI * (P1[0] * P2[0]))
      + (TMP3 * 2. * (-cI * (V1[2] * P2[0]) + 0.666666667 * cI * (TMP7)) + 2. *
      (V2[2] * (-cI * (P1[0] * TMP7) + cI * (V1[2] * TMP9))))));
  T3[3] = denom * (OM3 * (P3[0] * (P3[1] * (OM3 * (TMP22 * 1.333333333 * (-cI *
      (TMP3 * TMP4) + cI * (TMP0 * TMP23)) + 1.333333333 * (TMP6 * (-cI * (TMP7
      * TMP23) + cI * (TMP4 * TMP9)))) + (TMP0 * 0.666666667 * (+cI * (TMP33 +
      TMP9)) - 1.333333333 * cI * (TMP3 * TMP7))) + (P1[1] * (-cI * (TMP0 *
      TMP22) + cI * (TMP6 * TMP7)) + (P2[1] * (-cI * (TMP0 * TMP23) + cI *
      (TMP3 * TMP4)) + (TMP9 * - 1. * (+cI * (V1[3] * TMP6 + V2[3] * TMP4)) +
      (+cI * (V1[3] * TMP3 * TMP22 + V2[3] * TMP7 * TMP23)))))) + P3[1] *
      (P1[0] * (-cI * (TMP0 * TMP22) + cI * (TMP6 * TMP7)) + (P2[0] * (-cI *
      (TMP0 * TMP23) + cI * (TMP3 * TMP4)) + (TMP9 * - 1. * (+cI * (V2[2] *
      TMP4 + V1[2] * TMP6)) + (+cI * (V1[2] * TMP3 * TMP22 + V2[2] * TMP7 *
      TMP23)))))) + (P1[0] * (-cI * (V2[3] * TMP7) + cI * (TMP0 * P2[1])) +
      (P1[1] * (-cI * (V2[2] * TMP7) + cI * (TMP0 * P2[0])) + (TMP3 * - 1. *
      (+cI * (V1[2] * P2[1] + V1[3] * P2[0])) + TMP9 * (+cI * (V2[2] * V1[3] +
      V2[3] * V1[2]))))));
  T3[4] = denom * (OM3 * (P3[0] * (P3[2] * (OM3 * (TMP22 * 1.333333333 * (-cI *
      (TMP3 * TMP4) + cI * (TMP0 * TMP23)) + 1.333333333 * (TMP6 * (-cI * (TMP7
      * TMP23) + cI * (TMP4 * TMP9)))) + (TMP0 * 0.666666667 * (+cI * (TMP33 +
      TMP9)) - 1.333333333 * cI * (TMP3 * TMP7))) + (P1[2] * (-cI * (TMP0 *
      TMP22) + cI * (TMP6 * TMP7)) + (P2[2] * (-cI * (TMP0 * TMP23) + cI *
      (TMP3 * TMP4)) + (TMP9 * - 1. * (+cI * (V1[4] * TMP6 + V2[4] * TMP4)) +
      (+cI * (V1[4] * TMP3 * TMP22 + V2[4] * TMP7 * TMP23)))))) + P3[2] *
      (P1[0] * (-cI * (TMP0 * TMP22) + cI * (TMP6 * TMP7)) + (P2[0] * (-cI *
      (TMP0 * TMP23) + cI * (TMP3 * TMP4)) + (TMP9 * - 1. * (+cI * (V2[2] *
      TMP4 + V1[2] * TMP6)) + (+cI * (V1[2] * TMP3 * TMP22 + V2[2] * TMP7 *
      TMP23)))))) + (P1[0] * (-cI * (V2[4] * TMP7) + cI * (TMP0 * P2[2])) +
      (P1[2] * (-cI * (V2[2] * TMP7) + cI * (TMP0 * P2[0])) + (TMP3 * - 1. *
      (+cI * (V1[2] * P2[2] + V1[4] * P2[0])) + TMP9 * (+cI * (V2[2] * V1[4] +
      V2[4] * V1[2]))))));
  T3[5] = denom * (OM3 * (P3[0] * (P3[3] * (OM3 * (TMP22 * 1.333333333 * (-cI *
      (TMP3 * TMP4) + cI * (TMP0 * TMP23)) + 1.333333333 * (TMP6 * (-cI * (TMP7
      * TMP23) + cI * (TMP4 * TMP9)))) + (TMP0 * 0.666666667 * (+cI * (TMP33 +
      TMP9)) - 1.333333333 * cI * (TMP3 * TMP7))) + (P1[3] * (-cI * (TMP0 *
      TMP22) + cI * (TMP6 * TMP7)) + (P2[3] * (-cI * (TMP0 * TMP23) + cI *
      (TMP3 * TMP4)) + (TMP9 * - 1. * (+cI * (V1[5] * TMP6 + V2[5] * TMP4)) +
      (+cI * (V1[5] * TMP3 * TMP22 + V2[5] * TMP7 * TMP23)))))) + P3[3] *
      (P1[0] * (-cI * (TMP0 * TMP22) + cI * (TMP6 * TMP7)) + (P2[0] * (-cI *
      (TMP0 * TMP23) + cI * (TMP3 * TMP4)) + (TMP9 * - 1. * (+cI * (V2[2] *
      TMP4 + V1[2] * TMP6)) + (+cI * (V1[2] * TMP3 * TMP22 + V2[2] * TMP7 *
      TMP23)))))) + (P1[0] * (-cI * (V2[5] * TMP7) + cI * (TMP0 * P2[3])) +
      (P1[3] * (-cI * (V2[2] * TMP7) + cI * (TMP0 * P2[0])) + (TMP3 * - 1. *
      (+cI * (V1[2] * P2[3] + V1[5] * P2[0])) + TMP9 * (+cI * (V2[2] * V1[5] +
      V2[5] * V1[2]))))));
  T3[6] = denom * (OM3 * (P3[0] * (P3[1] * (OM3 * (TMP22 * 1.333333333 * (-cI *
      (TMP3 * TMP4) + cI * (TMP0 * TMP23)) + 1.333333333 * (TMP6 * (-cI * (TMP7
      * TMP23) + cI * (TMP4 * TMP9)))) + (TMP0 * 0.666666667 * (+cI * (TMP33 +
      TMP9)) - 1.333333333 * cI * (TMP3 * TMP7))) + (P1[1] * (-cI * (TMP0 *
      TMP22) + cI * (TMP6 * TMP7)) + (P2[1] * (-cI * (TMP0 * TMP23) + cI *
      (TMP3 * TMP4)) + (TMP9 * - 1. * (+cI * (V2[3] * TMP4 + V1[3] * TMP6)) +
      (+cI * (V1[3] * TMP3 * TMP22 + V2[3] * TMP7 * TMP23)))))) + P3[1] *
      (P1[0] * (-cI * (TMP0 * TMP22) + cI * (TMP6 * TMP7)) + (P2[0] * (-cI *
      (TMP0 * TMP23) + cI * (TMP3 * TMP4)) + (TMP9 * - 1. * (+cI * (V1[2] *
      TMP6 + V2[2] * TMP4)) + (+cI * (V1[2] * TMP3 * TMP22 + V2[2] * TMP7 *
      TMP23)))))) + (P1[0] * (-cI * (V2[3] * TMP7) + cI * (TMP0 * P2[1])) +
      (P1[1] * (-cI * (V2[2] * TMP7) + cI * (TMP0 * P2[0])) + (TMP3 * - 1. *
      (+cI * (V1[3] * P2[0] + V1[2] * P2[1])) + TMP9 * (+cI * (V2[3] * V1[2] +
      V2[2] * V1[3]))))));
  T3[7] = denom * (OM3 * (P3[1] * (P3[1] * (OM3 * (TMP22 * 1.333333333 * (-cI *
      (TMP3 * TMP4) + cI * (TMP0 * TMP23)) + 1.333333333 * (TMP6 * (-cI * (TMP7
      * TMP23) + cI * (TMP4 * TMP9)))) + (TMP0 * 0.666666667 * (+cI * (TMP33 +
      TMP9)) - 1.333333333 * cI * (TMP3 * TMP7))) + (P1[1] * 2. * (-cI * (TMP0
      * TMP22) + cI * (TMP6 * TMP7)) + (P2[1] * 2. * (-cI * (TMP0 * TMP23) + cI
      * (TMP3 * TMP4)) + (TMP9 * - 2. * (+cI * (V2[3] * TMP4 + V1[3] * TMP6)) +
      (+2. * cI * (V1[3] * TMP3 * TMP22 + V2[3] * TMP7 * TMP23)))))) + (TMP22 *
      0.666666667 * (-cI * (TMP0 * TMP23) + cI * (TMP3 * TMP4)) + 0.666666667 *
      (TMP6 * (-cI * (TMP4 * TMP9) + cI * (TMP7 * TMP23))))) + (TMP0 *
      0.666666667 * (+cI * (TMP33 + TMP9) + 3.000000000 * cI * (P1[1] * P2[1]))
      + (TMP3 * - 2. * (+cI * (V1[3] * P2[1]) + 0.666666667 * cI * (TMP7)) + 2.
      * (V2[3] * (-cI * (P1[1] * TMP7) + cI * (V1[3] * TMP9))))));
  T3[8] = denom * (OM3 * (P3[1] * (P3[2] * (OM3 * (TMP22 * 1.333333333 * (-cI *
      (TMP3 * TMP4) + cI * (TMP0 * TMP23)) + 1.333333333 * (TMP6 * (-cI * (TMP7
      * TMP23) + cI * (TMP4 * TMP9)))) + (TMP0 * 0.666666667 * (+cI * (TMP33 +
      TMP9)) - 1.333333333 * cI * (TMP3 * TMP7))) + (P1[2] * (-cI * (TMP0 *
      TMP22) + cI * (TMP6 * TMP7)) + (P2[2] * (-cI * (TMP0 * TMP23) + cI *
      (TMP3 * TMP4)) + (TMP9 * - 1. * (+cI * (V1[4] * TMP6 + V2[4] * TMP4)) +
      (+cI * (V1[4] * TMP3 * TMP22 + V2[4] * TMP7 * TMP23)))))) + P3[2] *
      (P1[1] * (-cI * (TMP0 * TMP22) + cI * (TMP6 * TMP7)) + (P2[1] * (-cI *
      (TMP0 * TMP23) + cI * (TMP3 * TMP4)) + (TMP9 * - 1. * (+cI * (V2[3] *
      TMP4 + V1[3] * TMP6)) + (+cI * (V1[3] * TMP3 * TMP22 + V2[3] * TMP7 *
      TMP23)))))) + (P1[1] * (-cI * (V2[4] * TMP7) + cI * (TMP0 * P2[2])) +
      (P1[2] * (-cI * (V2[3] * TMP7) + cI * (TMP0 * P2[1])) + (TMP3 * - 1. *
      (+cI * (V1[3] * P2[2] + V1[4] * P2[1])) + TMP9 * (+cI * (V2[3] * V1[4] +
      V2[4] * V1[3]))))));
  T3[9] = denom * (OM3 * (P3[1] * (P3[3] * (OM3 * (TMP22 * 1.333333333 * (-cI *
      (TMP3 * TMP4) + cI * (TMP0 * TMP23)) + 1.333333333 * (TMP6 * (-cI * (TMP7
      * TMP23) + cI * (TMP4 * TMP9)))) + (TMP0 * 0.666666667 * (+cI * (TMP33 +
      TMP9)) - 1.333333333 * cI * (TMP3 * TMP7))) + (P1[3] * (-cI * (TMP0 *
      TMP22) + cI * (TMP6 * TMP7)) + (P2[3] * (-cI * (TMP0 * TMP23) + cI *
      (TMP3 * TMP4)) + (TMP9 * - 1. * (+cI * (V1[5] * TMP6 + V2[5] * TMP4)) +
      (+cI * (V1[5] * TMP3 * TMP22 + V2[5] * TMP7 * TMP23)))))) + P3[3] *
      (P1[1] * (-cI * (TMP0 * TMP22) + cI * (TMP6 * TMP7)) + (P2[1] * (-cI *
      (TMP0 * TMP23) + cI * (TMP3 * TMP4)) + (TMP9 * - 1. * (+cI * (V2[3] *
      TMP4 + V1[3] * TMP6)) + (+cI * (V1[3] * TMP3 * TMP22 + V2[3] * TMP7 *
      TMP23)))))) + (P1[1] * (-cI * (V2[5] * TMP7) + cI * (TMP0 * P2[3])) +
      (P1[3] * (-cI * (V2[3] * TMP7) + cI * (TMP0 * P2[1])) + (TMP3 * - 1. *
      (+cI * (V1[3] * P2[3] + V1[5] * P2[1])) + TMP9 * (+cI * (V2[3] * V1[5] +
      V2[5] * V1[3]))))));
  T3[10] = denom * (OM3 * (P3[0] * (P3[2] * (OM3 * (TMP22 * 1.333333333 * (-cI
      * (TMP3 * TMP4) + cI * (TMP0 * TMP23)) + 1.333333333 * (TMP6 * (-cI *
      (TMP7 * TMP23) + cI * (TMP4 * TMP9)))) + (TMP0 * 0.666666667 * (+cI *
      (TMP33 + TMP9)) - 1.333333333 * cI * (TMP3 * TMP7))) + (P1[2] * (-cI *
      (TMP0 * TMP22) + cI * (TMP6 * TMP7)) + (P2[2] * (-cI * (TMP0 * TMP23) +
      cI * (TMP3 * TMP4)) + (TMP9 * - 1. * (+cI * (V2[4] * TMP4 + V1[4] *
      TMP6)) + (+cI * (V1[4] * TMP3 * TMP22 + V2[4] * TMP7 * TMP23)))))) +
      P3[2] * (P1[0] * (-cI * (TMP0 * TMP22) + cI * (TMP6 * TMP7)) + (P2[0] *
      (-cI * (TMP0 * TMP23) + cI * (TMP3 * TMP4)) + (TMP9 * - 1. * (+cI *
      (V1[2] * TMP6 + V2[2] * TMP4)) + (+cI * (V1[2] * TMP3 * TMP22 + V2[2] *
      TMP7 * TMP23)))))) + (P1[0] * (-cI * (V2[4] * TMP7) + cI * (TMP0 *
      P2[2])) + (P1[2] * (-cI * (V2[2] * TMP7) + cI * (TMP0 * P2[0])) + (TMP3 *
      - 1. * (+cI * (V1[4] * P2[0] + V1[2] * P2[2])) + TMP9 * (+cI * (V2[4] *
      V1[2] + V2[2] * V1[4]))))));
  T3[11] = denom * (OM3 * (P3[1] * (P3[2] * (OM3 * (TMP22 * 1.333333333 * (-cI
      * (TMP3 * TMP4) + cI * (TMP0 * TMP23)) + 1.333333333 * (TMP6 * (-cI *
      (TMP7 * TMP23) + cI * (TMP4 * TMP9)))) + (TMP0 * 0.666666667 * (+cI *
      (TMP33 + TMP9)) - 1.333333333 * cI * (TMP3 * TMP7))) + (P1[2] * (-cI *
      (TMP0 * TMP22) + cI * (TMP6 * TMP7)) + (P2[2] * (-cI * (TMP0 * TMP23) +
      cI * (TMP3 * TMP4)) + (TMP9 * - 1. * (+cI * (V2[4] * TMP4 + V1[4] *
      TMP6)) + (+cI * (V1[4] * TMP3 * TMP22 + V2[4] * TMP7 * TMP23)))))) +
      P3[2] * (P1[1] * (-cI * (TMP0 * TMP22) + cI * (TMP6 * TMP7)) + (P2[1] *
      (-cI * (TMP0 * TMP23) + cI * (TMP3 * TMP4)) + (TMP9 * - 1. * (+cI *
      (V1[3] * TMP6 + V2[3] * TMP4)) + (+cI * (V1[3] * TMP3 * TMP22 + V2[3] *
      TMP7 * TMP23)))))) + (P1[1] * (-cI * (V2[4] * TMP7) + cI * (TMP0 *
      P2[2])) + (P1[2] * (-cI * (V2[3] * TMP7) + cI * (TMP0 * P2[1])) + (TMP3 *
      - 1. * (+cI * (V1[4] * P2[1] + V1[3] * P2[2])) + TMP9 * (+cI * (V2[4] *
      V1[3] + V2[3] * V1[4]))))));
  T3[12] = denom * (OM3 * (P3[2] * (P3[2] * (OM3 * (TMP22 * 1.333333333 * (-cI
      * (TMP3 * TMP4) + cI * (TMP0 * TMP23)) + 1.333333333 * (TMP6 * (-cI *
      (TMP7 * TMP23) + cI * (TMP4 * TMP9)))) + (TMP0 * 0.666666667 * (+cI *
      (TMP33 + TMP9)) - 1.333333333 * cI * (TMP3 * TMP7))) + (P1[2] * 2. * (-cI
      * (TMP0 * TMP22) + cI * (TMP6 * TMP7)) + (P2[2] * 2. * (-cI * (TMP0 *
      TMP23) + cI * (TMP3 * TMP4)) + (TMP9 * - 2. * (+cI * (V2[4] * TMP4 +
      V1[4] * TMP6)) + (+2. * cI * (V1[4] * TMP3 * TMP22 + V2[4] * TMP7 *
      TMP23)))))) + (TMP22 * 0.666666667 * (-cI * (TMP0 * TMP23) + cI * (TMP3 *
      TMP4)) + 0.666666667 * (TMP6 * (-cI * (TMP4 * TMP9) + cI * (TMP7 *
      TMP23))))) + (TMP0 * 0.666666667 * (+cI * (TMP33 + TMP9) + 3.000000000 *
      cI * (P1[2] * P2[2])) + (TMP3 * - 2. * (+cI * (V1[4] * P2[2]) +
      0.666666667 * cI * (TMP7)) + 2. * (V2[4] * (-cI * (P1[2] * TMP7) + cI *
      (V1[4] * TMP9))))));
  T3[13] = denom * (OM3 * (P3[2] * (P3[3] * (OM3 * (TMP22 * 1.333333333 * (-cI
      * (TMP3 * TMP4) + cI * (TMP0 * TMP23)) + 1.333333333 * (TMP6 * (-cI *
      (TMP7 * TMP23) + cI * (TMP4 * TMP9)))) + (TMP0 * 0.666666667 * (+cI *
      (TMP33 + TMP9)) - 1.333333333 * cI * (TMP3 * TMP7))) + (P1[3] * (-cI *
      (TMP0 * TMP22) + cI * (TMP6 * TMP7)) + (P2[3] * (-cI * (TMP0 * TMP23) +
      cI * (TMP3 * TMP4)) + (TMP9 * - 1. * (+cI * (V1[5] * TMP6 + V2[5] *
      TMP4)) + (+cI * (V1[5] * TMP3 * TMP22 + V2[5] * TMP7 * TMP23)))))) +
      P3[3] * (P1[2] * (-cI * (TMP0 * TMP22) + cI * (TMP6 * TMP7)) + (P2[2] *
      (-cI * (TMP0 * TMP23) + cI * (TMP3 * TMP4)) + (TMP9 * - 1. * (+cI *
      (V2[4] * TMP4 + V1[4] * TMP6)) + (+cI * (V1[4] * TMP3 * TMP22 + V2[4] *
      TMP7 * TMP23)))))) + (P1[2] * (-cI * (V2[5] * TMP7) + cI * (TMP0 *
      P2[3])) + (P1[3] * (-cI * (V2[4] * TMP7) + cI * (TMP0 * P2[2])) + (TMP3 *
      - 1. * (+cI * (V1[4] * P2[3] + V1[5] * P2[2])) + TMP9 * (+cI * (V2[4] *
      V1[5] + V2[5] * V1[4]))))));
  T3[14] = denom * (OM3 * (P3[0] * (P3[3] * (OM3 * (TMP22 * 1.333333333 * (-cI
      * (TMP3 * TMP4) + cI * (TMP0 * TMP23)) + 1.333333333 * (TMP6 * (-cI *
      (TMP7 * TMP23) + cI * (TMP4 * TMP9)))) + (TMP0 * 0.666666667 * (+cI *
      (TMP33 + TMP9)) - 1.333333333 * cI * (TMP3 * TMP7))) + (P1[3] * (-cI *
      (TMP0 * TMP22) + cI * (TMP6 * TMP7)) + (P2[3] * (-cI * (TMP0 * TMP23) +
      cI * (TMP3 * TMP4)) + (TMP9 * - 1. * (+cI * (V2[5] * TMP4 + V1[5] *
      TMP6)) + (+cI * (V1[5] * TMP3 * TMP22 + V2[5] * TMP7 * TMP23)))))) +
      P3[3] * (P1[0] * (-cI * (TMP0 * TMP22) + cI * (TMP6 * TMP7)) + (P2[0] *
      (-cI * (TMP0 * TMP23) + cI * (TMP3 * TMP4)) + (TMP9 * - 1. * (+cI *
      (V1[2] * TMP6 + V2[2] * TMP4)) + (+cI * (V1[2] * TMP3 * TMP22 + V2[2] *
      TMP7 * TMP23)))))) + (P1[0] * (-cI * (V2[5] * TMP7) + cI * (TMP0 *
      P2[3])) + (P1[3] * (-cI * (V2[2] * TMP7) + cI * (TMP0 * P2[0])) + (TMP3 *
      - 1. * (+cI * (V1[5] * P2[0] + V1[2] * P2[3])) + TMP9 * (+cI * (V2[5] *
      V1[2] + V2[2] * V1[5]))))));
  T3[15] = denom * (OM3 * (P3[1] * (P3[3] * (OM3 * (TMP22 * 1.333333333 * (-cI
      * (TMP3 * TMP4) + cI * (TMP0 * TMP23)) + 1.333333333 * (TMP6 * (-cI *
      (TMP7 * TMP23) + cI * (TMP4 * TMP9)))) + (TMP0 * 0.666666667 * (+cI *
      (TMP33 + TMP9)) - 1.333333333 * cI * (TMP3 * TMP7))) + (P1[3] * (-cI *
      (TMP0 * TMP22) + cI * (TMP6 * TMP7)) + (P2[3] * (-cI * (TMP0 * TMP23) +
      cI * (TMP3 * TMP4)) + (TMP9 * - 1. * (+cI * (V2[5] * TMP4 + V1[5] *
      TMP6)) + (+cI * (V1[5] * TMP3 * TMP22 + V2[5] * TMP7 * TMP23)))))) +
      P3[3] * (P1[1] * (-cI * (TMP0 * TMP22) + cI * (TMP6 * TMP7)) + (P2[1] *
      (-cI * (TMP0 * TMP23) + cI * (TMP3 * TMP4)) + (TMP9 * - 1. * (+cI *
      (V1[3] * TMP6 + V2[3] * TMP4)) + (+cI * (V1[3] * TMP3 * TMP22 + V2[3] *
      TMP7 * TMP23)))))) + (P1[1] * (-cI * (V2[5] * TMP7) + cI * (TMP0 *
      P2[3])) + (P1[3] * (-cI * (V2[3] * TMP7) + cI * (TMP0 * P2[1])) + (TMP3 *
      - 1. * (+cI * (V1[5] * P2[1] + V1[3] * P2[3])) + TMP9 * (+cI * (V2[5] *
      V1[3] + V2[3] * V1[5]))))));
  T3[16] = denom * (OM3 * (P3[2] * (P3[3] * (OM3 * (TMP22 * 1.333333333 * (-cI
      * (TMP3 * TMP4) + cI * (TMP0 * TMP23)) + 1.333333333 * (TMP6 * (-cI *
      (TMP7 * TMP23) + cI * (TMP4 * TMP9)))) + (TMP0 * 0.666666667 * (+cI *
      (TMP33 + TMP9)) - 1.333333333 * cI * (TMP3 * TMP7))) + (P1[3] * (-cI *
      (TMP0 * TMP22) + cI * (TMP6 * TMP7)) + (P2[3] * (-cI * (TMP0 * TMP23) +
      cI * (TMP3 * TMP4)) + (TMP9 * - 1. * (+cI * (V2[5] * TMP4 + V1[5] *
      TMP6)) + (+cI * (V1[5] * TMP3 * TMP22 + V2[5] * TMP7 * TMP23)))))) +
      P3[3] * (P1[2] * (-cI * (TMP0 * TMP22) + cI * (TMP6 * TMP7)) + (P2[2] *
      (-cI * (TMP0 * TMP23) + cI * (TMP3 * TMP4)) + (TMP9 * - 1. * (+cI *
      (V1[4] * TMP6 + V2[4] * TMP4)) + (+cI * (V1[4] * TMP3 * TMP22 + V2[4] *
      TMP7 * TMP23)))))) + (P1[2] * (-cI * (V2[5] * TMP7) + cI * (TMP0 *
      P2[3])) + (P1[3] * (-cI * (V2[4] * TMP7) + cI * (TMP0 * P2[2])) + (TMP3 *
      - 1. * (+cI * (V1[5] * P2[2] + V1[4] * P2[3])) + TMP9 * (+cI * (V2[5] *
      V1[4] + V2[4] * V1[5]))))));
  T3[17] = denom * (OM3 * (P3[3] * (P3[3] * (OM3 * (TMP22 * 1.333333333 * (-cI
      * (TMP3 * TMP4) + cI * (TMP0 * TMP23)) + 1.333333333 * (TMP6 * (-cI *
      (TMP7 * TMP23) + cI * (TMP4 * TMP9)))) + (TMP0 * 0.666666667 * (+cI *
      (TMP33 + TMP9)) - 1.333333333 * cI * (TMP3 * TMP7))) + (P1[3] * 2. * (-cI
      * (TMP0 * TMP22) + cI * (TMP6 * TMP7)) + (P2[3] * 2. * (-cI * (TMP0 *
      TMP23) + cI * (TMP3 * TMP4)) + (TMP9 * - 2. * (+cI * (V2[5] * TMP4 +
      V1[5] * TMP6)) + (+2. * cI * (V1[5] * TMP3 * TMP22 + V2[5] * TMP7 *
      TMP23)))))) + (TMP22 * 0.666666667 * (-cI * (TMP0 * TMP23) + cI * (TMP3 *
      TMP4)) + 0.666666667 * (TMP6 * (-cI * (TMP4 * TMP9) + cI * (TMP7 *
      TMP23))))) + (TMP0 * 0.666666667 * (+cI * (TMP33 + TMP9) + 3.000000000 *
      cI * (P1[3] * P2[3])) + (TMP3 * - 2. * (+cI * (V1[5] * P2[3]) +
      0.666666667 * cI * (TMP7)) + 2. * (V2[5] * (-cI * (P1[3] * TMP7) + cI *
      (V1[5] * TMP9))))));
}


void VVT10_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP2; 
  complex<double> TMP1; 
  double P1[4]; 
  complex<double> TMP22; 
  double P2[4]; 
  complex<double> TMP23; 
  double P3[4]; 
  complex<double> denom; 
  double OM3; 
  complex<double> TMP33; 
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
  TMP33 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP23 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP1 = -1. * (P1[0] * (P2[1] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P2[2] *
      (V2[5] * V1[3] - V2[3] * V1[5]) + P2[3] * (V2[3] * V1[4] - V2[4] *
      V1[3]))) + (P1[1] * (P2[0] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P2[2] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P2[3] * (V2[4] * V1[2] - V2[2] *
      V1[4]))) + (P1[2] * (P2[0] * (V2[3] * V1[5] - V2[5] * V1[3]) + (P2[1] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P2[3] * (V2[2] * V1[3] - V2[3] *
      V1[2]))) + P1[3] * (P2[0] * (V2[4] * V1[3] - V2[3] * V1[4]) + (P2[1] *
      (V2[2] * V1[4] - V2[4] * V1[2]) + P2[2] * (V2[3] * V1[2] - V2[2] *
      V1[3]))))));
  TMP22 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  TMP2 = -1. * (P1[0] * (P2[1] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P2[2] *
      (V2[3] * V1[5] - V2[5] * V1[3]) + P2[3] * (V2[4] * V1[3] - V2[3] *
      V1[4]))) + (P1[1] * (P2[0] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P2[2] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P2[3] * (V2[2] * V1[4] - V2[4] *
      V1[2]))) + (P1[2] * (P2[0] * (V2[5] * V1[3] - V2[3] * V1[5]) + (P2[1] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P2[3] * (V2[3] * V1[2] - V2[2] *
      V1[3]))) + P1[3] * (P2[0] * (V2[3] * V1[4] - V2[4] * V1[3]) + (P2[1] *
      (V2[4] * V1[2] - V2[2] * V1[4]) + P2[2] * (V2[2] * V1[3] - V2[3] *
      V1[2]))))));
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * (OM3 * (P3[0] * (P3[0] * (OM3 * 5.333333333 * TMP22 * TMP23 *
      (-cI * (TMP1) + cI * (TMP2)) + 2.666666667 * (TMP33 * (-cI * (TMP1) + cI
      * (TMP2)))) + (+8. * (P2[0] * TMP23 * (-cI * (TMP2) + cI * (TMP1))) +
      P1[0] * 8. * TMP22 * (-cI * (TMP2) + cI * (TMP1)))) + 2.666666667 *
      (TMP22 * TMP23 * (-cI * (TMP1) + cI * (TMP2)))) + (P1[0] * 8. * P2[0] *
      (-cI * (TMP1) + cI * (TMP2)) + 2.666666667 * (TMP33 * (-cI * (TMP2) + cI
      * (TMP1)))));
  T3[3] = denom * (OM3 * (P3[0] * (P3[1] * (OM3 * 5.333333333 * TMP22 * TMP23 *
      (-cI * (TMP1) + cI * (TMP2)) + 2.666666667 * (TMP33 * (-cI * (TMP1) + cI
      * (TMP2)))) + (P1[1] * 4. * TMP22 * (-cI * (TMP2) + cI * (TMP1)) + 4. *
      (P2[1] * TMP23 * (-cI * (TMP2) + cI * (TMP1))))) + P3[1] * (P1[0] * 4. *
      TMP22 * (-cI * (TMP2) + cI * (TMP1)) + 4. * (P2[0] * TMP23 * (-cI *
      (TMP2) + cI * (TMP1))))) + (P1[0] * 4. * P2[1] * (-cI * (TMP1) + cI *
      (TMP2)) + 4. * (P1[1] * P2[0] * (-cI * (TMP1) + cI * (TMP2)))));
  T3[4] = denom * (OM3 * (P3[0] * (P3[2] * (OM3 * 5.333333333 * TMP22 * TMP23 *
      (-cI * (TMP1) + cI * (TMP2)) + 2.666666667 * (TMP33 * (-cI * (TMP1) + cI
      * (TMP2)))) + (P1[2] * 4. * TMP22 * (-cI * (TMP2) + cI * (TMP1)) + 4. *
      (P2[2] * TMP23 * (-cI * (TMP2) + cI * (TMP1))))) + P3[2] * (P1[0] * 4. *
      TMP22 * (-cI * (TMP2) + cI * (TMP1)) + 4. * (P2[0] * TMP23 * (-cI *
      (TMP2) + cI * (TMP1))))) + (P1[0] * 4. * P2[2] * (-cI * (TMP1) + cI *
      (TMP2)) + 4. * (P1[2] * P2[0] * (-cI * (TMP1) + cI * (TMP2)))));
  T3[5] = denom * (OM3 * (P3[0] * (P3[3] * (OM3 * 5.333333333 * TMP22 * TMP23 *
      (-cI * (TMP1) + cI * (TMP2)) + 2.666666667 * (TMP33 * (-cI * (TMP1) + cI
      * (TMP2)))) + (P1[3] * 4. * TMP22 * (-cI * (TMP2) + cI * (TMP1)) + 4. *
      (P2[3] * TMP23 * (-cI * (TMP2) + cI * (TMP1))))) + P3[3] * (P1[0] * 4. *
      TMP22 * (-cI * (TMP2) + cI * (TMP1)) + 4. * (P2[0] * TMP23 * (-cI *
      (TMP2) + cI * (TMP1))))) + (P1[0] * 4. * P2[3] * (-cI * (TMP1) + cI *
      (TMP2)) + 4. * (P1[3] * P2[0] * (-cI * (TMP1) + cI * (TMP2)))));
  T3[6] = denom * (OM3 * (P3[0] * (P3[1] * (OM3 * 5.333333333 * TMP22 * TMP23 *
      (-cI * (TMP1) + cI * (TMP2)) + 2.666666667 * (TMP33 * (-cI * (TMP1) + cI
      * (TMP2)))) + (P1[1] * 4. * TMP22 * (-cI * (TMP2) + cI * (TMP1)) + 4. *
      (P2[1] * TMP23 * (-cI * (TMP2) + cI * (TMP1))))) + P3[1] * (P1[0] * 4. *
      TMP22 * (-cI * (TMP2) + cI * (TMP1)) + 4. * (P2[0] * TMP23 * (-cI *
      (TMP2) + cI * (TMP1))))) + (P1[0] * 4. * P2[1] * (-cI * (TMP1) + cI *
      (TMP2)) + 4. * (P1[1] * P2[0] * (-cI * (TMP1) + cI * (TMP2)))));
  T3[7] = denom * (OM3 * (P3[1] * (P3[1] * (OM3 * 5.333333333 * TMP22 * TMP23 *
      (-cI * (TMP1) + cI * (TMP2)) + 2.666666667 * (TMP33 * (-cI * (TMP1) + cI
      * (TMP2)))) + (+8. * (P2[1] * TMP23 * (-cI * (TMP2) + cI * (TMP1))) +
      P1[1] * 8. * TMP22 * (-cI * (TMP2) + cI * (TMP1)))) + 2.666666667 *
      (TMP22 * TMP23 * (-cI * (TMP2) + cI * (TMP1)))) + (P1[1] * 8. * P2[1] *
      (-cI * (TMP1) + cI * (TMP2)) + 2.666666667 * (TMP33 * (-cI * (TMP1) + cI
      * (TMP2)))));
  T3[8] = denom * (OM3 * (P3[1] * (P3[2] * (OM3 * 5.333333333 * TMP22 * TMP23 *
      (-cI * (TMP1) + cI * (TMP2)) + 2.666666667 * (TMP33 * (-cI * (TMP1) + cI
      * (TMP2)))) + (P1[2] * 4. * TMP22 * (-cI * (TMP2) + cI * (TMP1)) + 4. *
      (P2[2] * TMP23 * (-cI * (TMP2) + cI * (TMP1))))) + P3[2] * (P1[1] * 4. *
      TMP22 * (-cI * (TMP2) + cI * (TMP1)) + 4. * (P2[1] * TMP23 * (-cI *
      (TMP2) + cI * (TMP1))))) + (P1[1] * 4. * P2[2] * (-cI * (TMP1) + cI *
      (TMP2)) + 4. * (P1[2] * P2[1] * (-cI * (TMP1) + cI * (TMP2)))));
  T3[9] = denom * (OM3 * (P3[1] * (P3[3] * (OM3 * 5.333333333 * TMP22 * TMP23 *
      (-cI * (TMP1) + cI * (TMP2)) + 2.666666667 * (TMP33 * (-cI * (TMP1) + cI
      * (TMP2)))) + (P1[3] * 4. * TMP22 * (-cI * (TMP2) + cI * (TMP1)) + 4. *
      (P2[3] * TMP23 * (-cI * (TMP2) + cI * (TMP1))))) + P3[3] * (P1[1] * 4. *
      TMP22 * (-cI * (TMP2) + cI * (TMP1)) + 4. * (P2[1] * TMP23 * (-cI *
      (TMP2) + cI * (TMP1))))) + (P1[1] * 4. * P2[3] * (-cI * (TMP1) + cI *
      (TMP2)) + 4. * (P1[3] * P2[1] * (-cI * (TMP1) + cI * (TMP2)))));
  T3[10] = denom * (OM3 * (P3[0] * (P3[2] * (OM3 * 5.333333333 * TMP22 * TMP23
      * (-cI * (TMP1) + cI * (TMP2)) + 2.666666667 * (TMP33 * (-cI * (TMP1) +
      cI * (TMP2)))) + (P1[2] * 4. * TMP22 * (-cI * (TMP2) + cI * (TMP1)) + 4.
      * (P2[2] * TMP23 * (-cI * (TMP2) + cI * (TMP1))))) + P3[2] * (P1[0] * 4.
      * TMP22 * (-cI * (TMP2) + cI * (TMP1)) + 4. * (P2[0] * TMP23 * (-cI *
      (TMP2) + cI * (TMP1))))) + (P1[0] * 4. * P2[2] * (-cI * (TMP1) + cI *
      (TMP2)) + 4. * (P1[2] * P2[0] * (-cI * (TMP1) + cI * (TMP2)))));
  T3[11] = denom * (OM3 * (P3[1] * (P3[2] * (OM3 * 5.333333333 * TMP22 * TMP23
      * (-cI * (TMP1) + cI * (TMP2)) + 2.666666667 * (TMP33 * (-cI * (TMP1) +
      cI * (TMP2)))) + (P1[2] * 4. * TMP22 * (-cI * (TMP2) + cI * (TMP1)) + 4.
      * (P2[2] * TMP23 * (-cI * (TMP2) + cI * (TMP1))))) + P3[2] * (P1[1] * 4.
      * TMP22 * (-cI * (TMP2) + cI * (TMP1)) + 4. * (P2[1] * TMP23 * (-cI *
      (TMP2) + cI * (TMP1))))) + (P1[1] * 4. * P2[2] * (-cI * (TMP1) + cI *
      (TMP2)) + 4. * (P1[2] * P2[1] * (-cI * (TMP1) + cI * (TMP2)))));
  T3[12] = denom * (OM3 * (P3[2] * (P3[2] * (OM3 * 5.333333333 * TMP22 * TMP23
      * (-cI * (TMP1) + cI * (TMP2)) + 2.666666667 * (TMP33 * (-cI * (TMP1) +
      cI * (TMP2)))) + (+8. * (P2[2] * TMP23 * (-cI * (TMP2) + cI * (TMP1))) +
      P1[2] * 8. * TMP22 * (-cI * (TMP2) + cI * (TMP1)))) + 2.666666667 *
      (TMP22 * TMP23 * (-cI * (TMP2) + cI * (TMP1)))) + (P1[2] * 8. * P2[2] *
      (-cI * (TMP1) + cI * (TMP2)) + 2.666666667 * (TMP33 * (-cI * (TMP1) + cI
      * (TMP2)))));
  T3[13] = denom * (OM3 * (P3[2] * (P3[3] * (OM3 * 5.333333333 * TMP22 * TMP23
      * (-cI * (TMP1) + cI * (TMP2)) + 2.666666667 * (TMP33 * (-cI * (TMP1) +
      cI * (TMP2)))) + (P1[3] * 4. * TMP22 * (-cI * (TMP2) + cI * (TMP1)) + 4.
      * (P2[3] * TMP23 * (-cI * (TMP2) + cI * (TMP1))))) + P3[3] * (P1[2] * 4.
      * TMP22 * (-cI * (TMP2) + cI * (TMP1)) + 4. * (P2[2] * TMP23 * (-cI *
      (TMP2) + cI * (TMP1))))) + (P1[2] * 4. * P2[3] * (-cI * (TMP1) + cI *
      (TMP2)) + 4. * (P1[3] * P2[2] * (-cI * (TMP1) + cI * (TMP2)))));
  T3[14] = denom * (OM3 * (P3[0] * (P3[3] * (OM3 * 5.333333333 * TMP22 * TMP23
      * (-cI * (TMP1) + cI * (TMP2)) + 2.666666667 * (TMP33 * (-cI * (TMP1) +
      cI * (TMP2)))) + (P1[3] * 4. * TMP22 * (-cI * (TMP2) + cI * (TMP1)) + 4.
      * (P2[3] * TMP23 * (-cI * (TMP2) + cI * (TMP1))))) + P3[3] * (P1[0] * 4.
      * TMP22 * (-cI * (TMP2) + cI * (TMP1)) + 4. * (P2[0] * TMP23 * (-cI *
      (TMP2) + cI * (TMP1))))) + (P1[0] * 4. * P2[3] * (-cI * (TMP1) + cI *
      (TMP2)) + 4. * (P1[3] * P2[0] * (-cI * (TMP1) + cI * (TMP2)))));
  T3[15] = denom * (OM3 * (P3[1] * (P3[3] * (OM3 * 5.333333333 * TMP22 * TMP23
      * (-cI * (TMP1) + cI * (TMP2)) + 2.666666667 * (TMP33 * (-cI * (TMP1) +
      cI * (TMP2)))) + (P1[3] * 4. * TMP22 * (-cI * (TMP2) + cI * (TMP1)) + 4.
      * (P2[3] * TMP23 * (-cI * (TMP2) + cI * (TMP1))))) + P3[3] * (P1[1] * 4.
      * TMP22 * (-cI * (TMP2) + cI * (TMP1)) + 4. * (P2[1] * TMP23 * (-cI *
      (TMP2) + cI * (TMP1))))) + (P1[1] * 4. * P2[3] * (-cI * (TMP1) + cI *
      (TMP2)) + 4. * (P1[3] * P2[1] * (-cI * (TMP1) + cI * (TMP2)))));
  T3[16] = denom * (OM3 * (P3[2] * (P3[3] * (OM3 * 5.333333333 * TMP22 * TMP23
      * (-cI * (TMP1) + cI * (TMP2)) + 2.666666667 * (TMP33 * (-cI * (TMP1) +
      cI * (TMP2)))) + (P1[3] * 4. * TMP22 * (-cI * (TMP2) + cI * (TMP1)) + 4.
      * (P2[3] * TMP23 * (-cI * (TMP2) + cI * (TMP1))))) + P3[3] * (P1[2] * 4.
      * TMP22 * (-cI * (TMP2) + cI * (TMP1)) + 4. * (P2[2] * TMP23 * (-cI *
      (TMP2) + cI * (TMP1))))) + (P1[2] * 4. * P2[3] * (-cI * (TMP1) + cI *
      (TMP2)) + 4. * (P1[3] * P2[2] * (-cI * (TMP1) + cI * (TMP2)))));
  T3[17] = denom * (OM3 * (P3[3] * (P3[3] * (OM3 * 5.333333333 * TMP22 * TMP23
      * (-cI * (TMP1) + cI * (TMP2)) + 2.666666667 * (TMP33 * (-cI * (TMP1) +
      cI * (TMP2)))) + (+8. * (P2[3] * TMP23 * (-cI * (TMP2) + cI * (TMP1))) +
      P1[3] * 8. * TMP22 * (-cI * (TMP2) + cI * (TMP1)))) + 2.666666667 *
      (TMP22 * TMP23 * (-cI * (TMP2) + cI * (TMP1)))) + (P1[3] * 8. * P2[3] *
      (-cI * (TMP1) + cI * (TMP2)) + 2.666666667 * (TMP33 * (-cI * (TMP1) + cI
      * (TMP2)))));
}

void VVT10_12_14_16_17_18_19_20_21_22_3(complex<double> V1[], complex<double>
    V2[], complex<double> COUP1, complex<double> COUP2, complex<double> COUP3,
    complex<double> COUP4, complex<double> COUP5, complex<double> COUP6,
    complex<double> COUP7, complex<double> COUP8, complex<double> COUP9,
    complex<double> COUP10, double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P3[4]; 
  double P2[4]; 
  double OM3; 
  double P1[4]; 
  complex<double> Ttmp[18]; 
  complex<double> denom; 
  int i; 
  VVT10_3(V1, V2, COUP1, M3, W3, T3); 
  VVT12_3(V1, V2, COUP2, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT14_3(V1, V2, COUP3, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT16_3(V1, V2, COUP4, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT17_3(V1, V2, COUP5, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT18_3(V1, V2, COUP6, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT19_3(V1, V2, COUP7, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT20_3(V1, V2, COUP8, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT21_3(V1, V2, COUP9, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT22_3(V1, V2, COUP10, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
}

void VVT21_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP17; 
  double P1[4]; 
  complex<double> TMP0; 
  double P2[4]; 
  complex<double> TMP9; 
  complex<double> TMP7; 
  complex<double> TMP20; 
  complex<double> TMP16; 
  complex<double> TMP21; 
  complex<double> TMP26; 
  complex<double> TMP25; 
  complex<double> TMP19; 
  complex<double> TMP3; 
  complex<double> TMP18; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  TMP25 = (P1[0] * - 1. * (P2[1] * T3[6] + P2[2] * T3[10] + P2[3] * T3[14] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[11] + P2[3] *
      T3[15] - P2[0] * T3[3]) + (P1[2] * (P2[1] * T3[8] + P2[2] * T3[12] +
      P2[3] * T3[16] - P2[0] * T3[4]) + P1[3] * (P2[1] * T3[9] + P2[2] * T3[13]
      + P2[3] * T3[17] - P2[0] * T3[5]))));
  TMP26 = (P1[0] * - 1. * (P2[1] * T3[3] + P2[2] * T3[4] + P2[3] * T3[5] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[8] + P2[3] * T3[9]
      - P2[0] * T3[6]) + (P1[2] * (P2[1] * T3[11] + P2[2] * T3[12] + P2[3] *
      T3[13] - P2[0] * T3[10]) + P1[3] * (P2[1] * T3[15] + P2[2] * T3[16] +
      P2[3] * T3[17] - P2[0] * T3[14]))));
  TMP20 = (P1[0] * - 1. * (V2[3] * T3[6] + V2[4] * T3[10] + V2[5] * T3[14] -
      V2[2] * T3[2]) + (P1[1] * (V2[3] * T3[7] + V2[4] * T3[11] + V2[5] *
      T3[15] - V2[2] * T3[3]) + (P1[2] * (V2[3] * T3[8] + V2[4] * T3[12] +
      V2[5] * T3[16] - V2[2] * T3[4]) + P1[3] * (V2[3] * T3[9] + V2[4] * T3[13]
      + V2[5] * T3[17] - V2[2] * T3[5]))));
  TMP21 = (P1[0] * - 1. * (V2[3] * T3[3] + V2[4] * T3[4] + V2[5] * T3[5] -
      V2[2] * T3[2]) + (P1[1] * (V2[3] * T3[7] + V2[4] * T3[8] + V2[5] * T3[9]
      - V2[2] * T3[6]) + (P1[2] * (V2[3] * T3[11] + V2[4] * T3[12] + V2[5] *
      T3[13] - V2[2] * T3[10]) + P1[3] * (V2[3] * T3[15] + V2[4] * T3[16] +
      V2[5] * T3[17] - V2[2] * T3[14]))));
  TMP0 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP19 = (P2[0] * - 1. * (V1[3] * T3[3] + V1[4] * T3[4] + V1[5] * T3[5] -
      V1[2] * T3[2]) + (P2[1] * (V1[3] * T3[7] + V1[4] * T3[8] + V1[5] * T3[9]
      - V1[2] * T3[6]) + (P2[2] * (V1[3] * T3[11] + V1[4] * T3[12] + V1[5] *
      T3[13] - V1[2] * T3[10]) + P2[3] * (V1[3] * T3[15] + V1[4] * T3[16] +
      V1[5] * T3[17] - V1[2] * T3[14]))));
  TMP18 = (P2[0] * - 1. * (V1[3] * T3[6] + V1[4] * T3[10] + V1[5] * T3[14] -
      V1[2] * T3[2]) + (P2[1] * (V1[3] * T3[7] + V1[4] * T3[11] + V1[5] *
      T3[15] - V1[2] * T3[3]) + (P2[2] * (V1[3] * T3[8] + V1[4] * T3[12] +
      V1[5] * T3[16] - V1[2] * T3[4]) + P2[3] * (V1[3] * T3[9] + V1[4] * T3[13]
      + V1[5] * T3[17] - V1[2] * T3[5]))));
  TMP9 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP17 = (V1[2] * - 1. * (V2[3] * T3[3] + V2[4] * T3[4] + V2[5] * T3[5] -
      V2[2] * T3[2]) + (V1[3] * (V2[3] * T3[7] + V2[4] * T3[8] + V2[5] * T3[9]
      - V2[2] * T3[6]) + (V1[4] * (V2[3] * T3[11] + V2[4] * T3[12] + V2[5] *
      T3[13] - V2[2] * T3[10]) + V1[5] * (V2[3] * T3[15] + V2[4] * T3[16] +
      V2[5] * T3[17] - V2[2] * T3[14]))));
  TMP16 = (V1[2] * - 1. * (V2[3] * T3[6] + V2[4] * T3[10] + V2[5] * T3[14] -
      V2[2] * T3[2]) + (V1[3] * (V2[3] * T3[7] + V2[4] * T3[11] + V2[5] *
      T3[15] - V2[2] * T3[3]) + (V1[4] * (V2[3] * T3[8] + V2[4] * T3[12] +
      V2[5] * T3[16] - V2[2] * T3[4]) + V1[5] * (V2[3] * T3[9] + V2[4] * T3[13]
      + V2[5] * T3[17] - V2[2] * T3[5]))));
  TMP7 = (V1[2] * P2[0] - V1[3] * P2[1] - V1[4] * P2[2] - V1[5] * P2[3]); 
  TMP3 = (V2[2] * P1[0] - V2[3] * P1[1] - V2[4] * P1[2] - V2[5] * P1[3]); 
  vertex = COUP * (TMP0 * - 1. * (+cI * (TMP25 + TMP26)) + (TMP3 * (+cI *
      (TMP18 + TMP19)) + (TMP7 * (+cI * (TMP20 + TMP21)) - TMP9 * (+cI * (TMP16
      + TMP17)))));
}


void VVT12_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> TMP22; 
  double P2[4]; 
  complex<double> TMP23; 
  double P3[4]; 
  complex<double> TMP49; 
  complex<double> denom; 
  double OM3; 
  complex<double> TMP50; 
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
  TMP50 = -1. * (P2[0] * (P3[1] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] *
      (V2[5] * V1[3] - V2[3] * V1[5]) + P3[3] * (V2[3] * V1[4] - V2[4] *
      V1[3]))) + (P2[1] * (P3[0] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] * (V2[4] * V1[2] - V2[2] *
      V1[4]))) + (P2[2] * (P3[0] * (V2[3] * V1[5] - V2[5] * V1[3]) + (P3[1] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] * (V2[2] * V1[3] - V2[3] *
      V1[2]))) + P2[3] * (P3[0] * (V2[4] * V1[3] - V2[3] * V1[4]) + (P3[1] *
      (V2[2] * V1[4] - V2[4] * V1[2]) + P3[2] * (V2[3] * V1[2] - V2[2] *
      V1[3]))))));
  TMP49 = -1. * (P1[0] * (P3[1] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] *
      (V2[5] * V1[3] - V2[3] * V1[5]) + P3[3] * (V2[3] * V1[4] - V2[4] *
      V1[3]))) + (P1[1] * (P3[0] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] * (V2[4] * V1[2] - V2[2] *
      V1[4]))) + (P1[2] * (P3[0] * (V2[3] * V1[5] - V2[5] * V1[3]) + (P3[1] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] * (V2[2] * V1[3] - V2[3] *
      V1[2]))) + P1[3] * (P3[0] * (V2[4] * V1[3] - V2[3] * V1[4]) + (P3[1] *
      (V2[2] * V1[4] - V2[4] * V1[2]) + P3[2] * (V2[3] * V1[2] - V2[2] *
      V1[3]))))));
  TMP22 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  TMP23 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * - 0.500000000 * cI * (OM3 * P3[0] * (TMP22 * (P3[1] * 4. *
      (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] * 4. * (V2[5] * V1[3] - V2[3] *
      V1[5]) + 4. * (P3[3] * (V2[3] * V1[4] - V2[4] * V1[3])))) + (TMP23 *
      (P3[1] * 4. * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] * 4. * (V2[3] *
      V1[5] - V2[5] * V1[3]) + 4. * (P3[3] * (V2[4] * V1[3] - V2[3] * V1[4]))))
      + 1.333333333 * (P3[0] * (TMP50 - TMP49)))) + (P1[0] * (P3[1] * 4. *
      (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] * 4. * (V2[5] * V1[3] - V2[3] *
      V1[5]) + 4. * (P3[3] * (V2[3] * V1[4] - V2[4] * V1[3])))) + (P2[0] *
      (P3[1] * 4. * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] * 4. * (V2[3] *
      V1[5] - V2[5] * V1[3]) + 4. * (P3[3] * (V2[4] * V1[3] - V2[3] * V1[4]))))
      + (-1.333333333 * (TMP50) + 1.333333333 * (TMP49)))));
  T3[6] = denom * - 0.500000000 * cI * (OM3 * (P3[0] * (TMP22 * (P3[0] * 2. *
      (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] * 2. * (V2[5] * V1[2] - V2[2] *
      V1[5]) + 2. * (P3[3] * (V2[2] * V1[4] - V2[4] * V1[2])))) + (TMP23 *
      (P3[0] * 2. * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] * 2. * (V2[2] *
      V1[5] - V2[5] * V1[2]) + 2. * (P3[3] * (V2[4] * V1[2] - V2[2] * V1[4]))))
      + 1.333333333 * (P3[1] * (TMP50 - TMP49)))) + P3[1] * (TMP22 * (P3[1] *
      2. * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] * 2. * (V2[5] * V1[3] -
      V2[3] * V1[5]) + 2. * (P3[3] * (V2[3] * V1[4] - V2[4] * V1[3])))) + TMP23
      * (P3[1] * 2. * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] * 2. * (V2[3] *
      V1[5] - V2[5] * V1[3]) + 2. * (P3[3] * (V2[4] * V1[3] - V2[3] *
      V1[4])))))) + (P3[2] * (V1[5] * (V2[2] * 2. * (P2[0] - P1[0]) + 2. *
      (V2[3] * (P2[1] - P1[1]))) + V2[5] * (V1[2] * 2. * (P1[0] - P2[0]) + 2. *
      (V1[3] * (P1[1] - P2[1])))) + (P3[3] * (V1[4] * (V2[2] * 2. * (P1[0] -
      P2[0]) + 2. * (V2[3] * (P1[1] - P2[1]))) + V2[4] * (V1[2] * 2. * (P2[0] -
      P1[0]) + 2. * (V1[3] * (P2[1] - P1[1])))) + (P3[0] * (V1[4] * 2. * V2[5]
      * (P2[0] - P1[0]) + 2. * (V1[5] * V2[4] * (P1[0] - P2[0]))) + P3[1] *
      (V1[4] * 2. * V2[5] * (P2[1] - P1[1]) + 2. * (V1[5] * V2[4] * (P1[1] -
      P2[1])))))));
  T3[10] = denom * - 0.500000000 * cI * (OM3 * (P3[0] * (TMP22 * (P3[0] * 2. *
      (V2[5] * V1[3] - V2[3] * V1[5]) + (P3[1] * 2. * (V2[2] * V1[5] - V2[5] *
      V1[2]) + 2. * (P3[3] * (V2[3] * V1[2] - V2[2] * V1[3])))) + (TMP23 *
      (P3[0] * 2. * (V2[3] * V1[5] - V2[5] * V1[3]) + (P3[1] * 2. * (V2[5] *
      V1[2] - V2[2] * V1[5]) + 2. * (P3[3] * (V2[2] * V1[3] - V2[3] * V1[2]))))
      + 1.333333333 * (P3[2] * (TMP50 - TMP49)))) + P3[2] * (TMP22 * (P3[1] *
      2. * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] * 2. * (V2[5] * V1[3] -
      V2[3] * V1[5]) + 2. * (P3[3] * (V2[3] * V1[4] - V2[4] * V1[3])))) + TMP23
      * (P3[1] * 2. * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] * 2. * (V2[3] *
      V1[5] - V2[5] * V1[3]) + 2. * (P3[3] * (V2[4] * V1[3] - V2[3] *
      V1[4])))))) + (P3[1] * (V1[5] * (V2[2] * 2. * (P1[0] - P2[0]) + 2. *
      (V2[4] * (P1[2] - P2[2]))) + V2[5] * (V1[2] * 2. * (P2[0] - P1[0]) + 2. *
      (V1[4] * (P2[2] - P1[2])))) + (P3[3] * (V1[3] * (V2[2] * 2. * (P2[0] -
      P1[0]) + 2. * (V2[4] * (P2[2] - P1[2]))) + V2[3] * (V1[2] * 2. * (P1[0] -
      P2[0]) + 2. * (V1[4] * (P1[2] - P2[2])))) + (P3[0] * (V1[3] * 2. * V2[5]
      * (P1[0] - P2[0]) + 2. * (V1[5] * V2[3] * (P2[0] - P1[0]))) + P3[2] *
      (V1[3] * 2. * V2[5] * (P1[2] - P2[2]) + 2. * (V1[5] * V2[3] * (P2[2] -
      P1[2])))))));
  T3[14] = denom * - 0.500000000 * cI * (OM3 * (P3[0] * (TMP22 * (P3[0] * 2. *
      (V2[3] * V1[4] - V2[4] * V1[3]) + (P3[1] * 2. * (V2[4] * V1[2] - V2[2] *
      V1[4]) + 2. * (P3[2] * (V2[2] * V1[3] - V2[3] * V1[2])))) + (TMP23 *
      (P3[0] * 2. * (V2[4] * V1[3] - V2[3] * V1[4]) + (P3[1] * 2. * (V2[2] *
      V1[4] - V2[4] * V1[2]) + 2. * (P3[2] * (V2[3] * V1[2] - V2[2] * V1[3]))))
      + 1.333333333 * (P3[3] * (TMP50 - TMP49)))) + P3[3] * (TMP22 * (P3[1] *
      2. * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] * 2. * (V2[5] * V1[3] -
      V2[3] * V1[5]) + 2. * (P3[3] * (V2[3] * V1[4] - V2[4] * V1[3])))) + TMP23
      * (P3[1] * 2. * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] * 2. * (V2[3] *
      V1[5] - V2[5] * V1[3]) + 2. * (P3[3] * (V2[4] * V1[3] - V2[3] *
      V1[4])))))) + (P3[1] * (V1[4] * (V2[2] * 2. * (P2[0] - P1[0]) + 2. *
      (V2[5] * (P2[3] - P1[3]))) + V2[4] * (V1[2] * 2. * (P1[0] - P2[0]) + 2. *
      (V1[5] * (P1[3] - P2[3])))) + (P3[2] * (V1[3] * (V2[2] * 2. * (P1[0] -
      P2[0]) + 2. * (V2[5] * (P1[3] - P2[3]))) + V2[3] * (V1[2] * 2. * (P2[0] -
      P1[0]) + 2. * (V1[5] * (P2[3] - P1[3])))) + (P3[0] * (V1[3] * 2. * V2[4]
      * (P2[0] - P1[0]) + 2. * (V1[4] * V2[3] * (P1[0] - P2[0]))) + P3[3] *
      (V1[3] * 2. * V2[4] * (P2[3] - P1[3]) + 2. * (V1[4] * V2[3] * (P1[3] -
      P2[3])))))));
  T3[3] = denom * 0.500000000 * cI * (OM3 * (P3[0] * (TMP22 * (P3[0] * 2. *
      (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] * 2. * (V2[2] * V1[5] - V2[5] *
      V1[2]) + 2. * (P3[3] * (V2[4] * V1[2] - V2[2] * V1[4])))) + (TMP23 *
      (P3[0] * 2. * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] * 2. * (V2[5] *
      V1[2] - V2[2] * V1[5]) + 2. * (P3[3] * (V2[2] * V1[4] - V2[4] * V1[2]))))
      + 1.333333333 * (P3[1] * (TMP49 - TMP50)))) + P3[1] * (TMP22 * (P3[1] *
      2. * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] * 2. * (V2[3] * V1[5] -
      V2[5] * V1[3]) + 2. * (P3[3] * (V2[4] * V1[3] - V2[3] * V1[4])))) + TMP23
      * (P3[1] * 2. * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] * 2. * (V2[5] *
      V1[3] - V2[3] * V1[5]) + 2. * (P3[3] * (V2[3] * V1[4] - V2[4] *
      V1[3])))))) + (P3[2] * (V1[5] * (V2[2] * 2. * (P1[0] - P2[0]) + 2. *
      (V2[3] * (P1[1] - P2[1]))) + V2[5] * (V1[2] * 2. * (P2[0] - P1[0]) + 2. *
      (V1[3] * (P2[1] - P1[1])))) + (P3[3] * (V1[4] * (V2[2] * 2. * (P2[0] -
      P1[0]) + 2. * (V2[3] * (P2[1] - P1[1]))) + V2[4] * (V1[2] * 2. * (P1[0] -
      P2[0]) + 2. * (V1[3] * (P1[1] - P2[1])))) + (P3[0] * (V1[4] * 2. * V2[5]
      * (P1[0] - P2[0]) + 2. * (V1[5] * V2[4] * (P2[0] - P1[0]))) + P3[1] *
      (V1[4] * 2. * V2[5] * (P1[1] - P2[1]) + 2. * (V1[5] * V2[4] * (P2[1] -
      P1[1])))))));
  T3[7] = denom * 0.500000000 * cI * (OM3 * P3[1] * (TMP22 * (P3[0] * 4. *
      (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] * 4. * (V2[2] * V1[5] - V2[5] *
      V1[2]) + 4. * (P3[3] * (V2[4] * V1[2] - V2[2] * V1[4])))) + (TMP23 *
      (P3[0] * 4. * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] * 4. * (V2[5] *
      V1[2] - V2[2] * V1[5]) + 4. * (P3[3] * (V2[2] * V1[4] - V2[4] * V1[2]))))
      + 1.333333333 * (P3[1] * (TMP49 - TMP50)))) + (P1[1] * (P3[0] * 4. *
      (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] * 4. * (V2[2] * V1[5] - V2[5] *
      V1[2]) + 4. * (P3[3] * (V2[4] * V1[2] - V2[2] * V1[4])))) + (P2[1] *
      (P3[0] * 4. * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] * 4. * (V2[5] *
      V1[2] - V2[2] * V1[5]) + 4. * (P3[3] * (V2[2] * V1[4] - V2[4] * V1[2]))))
      + (-1.333333333 * (TMP50) + 1.333333333 * (TMP49)))));
  T3[11] = denom * 0.500000000 * cI * (OM3 * (P3[1] * (TMP22 * (P3[0] * 2. *
      (V2[3] * V1[5] - V2[5] * V1[3]) + (P3[1] * 2. * (V2[5] * V1[2] - V2[2] *
      V1[5]) + 2. * (P3[3] * (V2[2] * V1[3] - V2[3] * V1[2])))) + (TMP23 *
      (P3[0] * 2. * (V2[5] * V1[3] - V2[3] * V1[5]) + (P3[1] * 2. * (V2[2] *
      V1[5] - V2[5] * V1[2]) + 2. * (P3[3] * (V2[3] * V1[2] - V2[2] * V1[3]))))
      + 1.333333333 * (P3[2] * (TMP49 - TMP50)))) + P3[2] * (TMP22 * (P3[0] *
      2. * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] * 2. * (V2[2] * V1[5] -
      V2[5] * V1[2]) + 2. * (P3[3] * (V2[4] * V1[2] - V2[2] * V1[4])))) + TMP23
      * (P3[0] * 2. * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] * 2. * (V2[5] *
      V1[2] - V2[2] * V1[5]) + 2. * (P3[3] * (V2[2] * V1[4] - V2[4] *
      V1[2])))))) + (P3[0] * (V1[5] * (V2[3] * 2. * (P1[1] - P2[1]) + 2. *
      (V2[4] * (P2[2] - P1[2]))) + V2[5] * (V1[3] * 2. * (P2[1] - P1[1]) + 2. *
      (V1[4] * (P1[2] - P2[2])))) + (P3[3] * (V1[2] * (V2[3] * 2. * (P2[1] -
      P1[1]) + 2. * (V2[4] * (P1[2] - P2[2]))) + V2[2] * (V1[3] * 2. * (P1[1] -
      P2[1]) + 2. * (V1[4] * (P2[2] - P1[2])))) + (P3[1] * (V1[2] * 2. * V2[5]
      * (P1[1] - P2[1]) + 2. * (V1[5] * V2[2] * (P2[1] - P1[1]))) + P3[2] *
      (V1[2] * 2. * V2[5] * (P2[2] - P1[2]) + 2. * (V1[5] * V2[2] * (P1[2] -
      P2[2])))))));
  T3[15] = denom * 0.500000000 * cI * (OM3 * (P3[1] * (TMP22 * (P3[0] * 2. *
      (V2[4] * V1[3] - V2[3] * V1[4]) + (P3[1] * 2. * (V2[2] * V1[4] - V2[4] *
      V1[2]) + 2. * (P3[2] * (V2[3] * V1[2] - V2[2] * V1[3])))) + (TMP23 *
      (P3[0] * 2. * (V2[3] * V1[4] - V2[4] * V1[3]) + (P3[1] * 2. * (V2[4] *
      V1[2] - V2[2] * V1[4]) + 2. * (P3[2] * (V2[2] * V1[3] - V2[3] * V1[2]))))
      + 1.333333333 * (P3[3] * (TMP49 - TMP50)))) + P3[3] * (TMP22 * (P3[0] *
      2. * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] * 2. * (V2[2] * V1[5] -
      V2[5] * V1[2]) + 2. * (P3[3] * (V2[4] * V1[2] - V2[2] * V1[4])))) + TMP23
      * (P3[0] * 2. * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] * 2. * (V2[5] *
      V1[2] - V2[2] * V1[5]) + 2. * (P3[3] * (V2[2] * V1[4] - V2[4] *
      V1[2])))))) + (P3[0] * (V1[4] * (V2[3] * 2. * (P2[1] - P1[1]) + 2. *
      (V2[5] * (P1[3] - P2[3]))) + V2[4] * (V1[3] * 2. * (P1[1] - P2[1]) + 2. *
      (V1[5] * (P2[3] - P1[3])))) + (P3[2] * (V1[2] * (V2[3] * 2. * (P1[1] -
      P2[1]) + 2. * (V2[5] * (P2[3] - P1[3]))) + V2[2] * (V1[3] * 2. * (P2[1] -
      P1[1]) + 2. * (V1[5] * (P1[3] - P2[3])))) + (P3[1] * (V1[2] * 2. * V2[4]
      * (P2[1] - P1[1]) + 2. * (V1[4] * V2[2] * (P1[1] - P2[1]))) + P3[3] *
      (V1[2] * 2. * V2[4] * (P1[3] - P2[3]) + 2. * (V1[4] * V2[2] * (P2[3] -
      P1[3])))))));
  T3[4] = denom * 0.500000000 * cI * (OM3 * (P3[0] * (TMP22 * (P3[0] * 2. *
      (V2[3] * V1[5] - V2[5] * V1[3]) + (P3[1] * 2. * (V2[5] * V1[2] - V2[2] *
      V1[5]) + 2. * (P3[3] * (V2[2] * V1[3] - V2[3] * V1[2])))) + (TMP23 *
      (P3[0] * 2. * (V2[5] * V1[3] - V2[3] * V1[5]) + (P3[1] * 2. * (V2[2] *
      V1[5] - V2[5] * V1[2]) + 2. * (P3[3] * (V2[3] * V1[2] - V2[2] * V1[3]))))
      + 1.333333333 * (P3[2] * (TMP49 - TMP50)))) + P3[2] * (TMP22 * (P3[1] *
      2. * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] * 2. * (V2[3] * V1[5] -
      V2[5] * V1[3]) + 2. * (P3[3] * (V2[4] * V1[3] - V2[3] * V1[4])))) + TMP23
      * (P3[1] * 2. * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] * 2. * (V2[5] *
      V1[3] - V2[3] * V1[5]) + 2. * (P3[3] * (V2[3] * V1[4] - V2[4] *
      V1[3])))))) + (P3[1] * (V1[5] * (V2[2] * 2. * (P2[0] - P1[0]) + 2. *
      (V2[4] * (P2[2] - P1[2]))) + V2[5] * (V1[2] * 2. * (P1[0] - P2[0]) + 2. *
      (V1[4] * (P1[2] - P2[2])))) + (P3[3] * (V1[3] * (V2[2] * 2. * (P1[0] -
      P2[0]) + 2. * (V2[4] * (P1[2] - P2[2]))) + V2[3] * (V1[2] * 2. * (P2[0] -
      P1[0]) + 2. * (V1[4] * (P2[2] - P1[2])))) + (P3[0] * (V1[3] * 2. * V2[5]
      * (P2[0] - P1[0]) + 2. * (V1[5] * V2[3] * (P1[0] - P2[0]))) + P3[2] *
      (V1[3] * 2. * V2[5] * (P2[2] - P1[2]) + 2. * (V1[5] * V2[3] * (P1[2] -
      P2[2])))))));
  T3[8] = denom * 0.500000000 * cI * (OM3 * (P3[1] * (TMP22 * (P3[0] * 2. *
      (V2[3] * V1[5] - V2[5] * V1[3]) + (P3[1] * 2. * (V2[5] * V1[2] - V2[2] *
      V1[5]) + 2. * (P3[3] * (V2[2] * V1[3] - V2[3] * V1[2])))) + (TMP23 *
      (P3[0] * 2. * (V2[5] * V1[3] - V2[3] * V1[5]) + (P3[1] * 2. * (V2[2] *
      V1[5] - V2[5] * V1[2]) + 2. * (P3[3] * (V2[3] * V1[2] - V2[2] * V1[3]))))
      + 1.333333333 * (P3[2] * (TMP49 - TMP50)))) + P3[2] * (TMP22 * (P3[0] *
      2. * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] * 2. * (V2[2] * V1[5] -
      V2[5] * V1[2]) + 2. * (P3[3] * (V2[4] * V1[2] - V2[2] * V1[4])))) + TMP23
      * (P3[0] * 2. * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] * 2. * (V2[5] *
      V1[2] - V2[2] * V1[5]) + 2. * (P3[3] * (V2[2] * V1[4] - V2[4] *
      V1[2])))))) + (P3[0] * (V1[5] * (V2[3] * 2. * (P1[1] - P2[1]) + 2. *
      (V2[4] * (P2[2] - P1[2]))) + V2[5] * (V1[3] * 2. * (P2[1] - P1[1]) + 2. *
      (V1[4] * (P1[2] - P2[2])))) + (P3[3] * (V1[2] * (V2[3] * 2. * (P2[1] -
      P1[1]) + 2. * (V2[4] * (P1[2] - P2[2]))) + V2[2] * (V1[3] * 2. * (P1[1] -
      P2[1]) + 2. * (V1[4] * (P2[2] - P1[2])))) + (P3[1] * (V1[2] * 2. * V2[5]
      * (P1[1] - P2[1]) + 2. * (V1[5] * V2[2] * (P2[1] - P1[1]))) + P3[2] *
      (V1[2] * 2. * V2[5] * (P2[2] - P1[2]) + 2. * (V1[5] * V2[2] * (P1[2] -
      P2[2])))))));
  T3[12] = denom * 0.500000000 * cI * (OM3 * P3[2] * (TMP22 * (P3[0] * 4. *
      (V2[3] * V1[5] - V2[5] * V1[3]) + (P3[1] * 4. * (V2[5] * V1[2] - V2[2] *
      V1[5]) + 4. * (P3[3] * (V2[2] * V1[3] - V2[3] * V1[2])))) + (TMP23 *
      (P3[0] * 4. * (V2[5] * V1[3] - V2[3] * V1[5]) + (P3[1] * 4. * (V2[2] *
      V1[5] - V2[5] * V1[2]) + 4. * (P3[3] * (V2[3] * V1[2] - V2[2] * V1[3]))))
      + 1.333333333 * (P3[2] * (TMP49 - TMP50)))) + (P1[2] * (P3[0] * 4. *
      (V2[3] * V1[5] - V2[5] * V1[3]) + (P3[1] * 4. * (V2[5] * V1[2] - V2[2] *
      V1[5]) + 4. * (P3[3] * (V2[2] * V1[3] - V2[3] * V1[2])))) + (P2[2] *
      (P3[0] * 4. * (V2[5] * V1[3] - V2[3] * V1[5]) + (P3[1] * 4. * (V2[2] *
      V1[5] - V2[5] * V1[2]) + 4. * (P3[3] * (V2[3] * V1[2] - V2[2] * V1[3]))))
      + (-1.333333333 * (TMP50) + 1.333333333 * (TMP49)))));
  T3[16] = denom * 0.500000000 * cI * (OM3 * (P3[2] * (TMP22 * (P3[0] * 2. *
      (V2[4] * V1[3] - V2[3] * V1[4]) + (P3[1] * 2. * (V2[2] * V1[4] - V2[4] *
      V1[2]) + 2. * (P3[2] * (V2[3] * V1[2] - V2[2] * V1[3])))) + (TMP23 *
      (P3[0] * 2. * (V2[3] * V1[4] - V2[4] * V1[3]) + (P3[1] * 2. * (V2[4] *
      V1[2] - V2[2] * V1[4]) + 2. * (P3[2] * (V2[2] * V1[3] - V2[3] * V1[2]))))
      + 1.333333333 * (P3[3] * (TMP49 - TMP50)))) + P3[3] * (TMP22 * (P3[0] *
      2. * (V2[3] * V1[5] - V2[5] * V1[3]) + (P3[1] * 2. * (V2[5] * V1[2] -
      V2[2] * V1[5]) + 2. * (P3[3] * (V2[2] * V1[3] - V2[3] * V1[2])))) + TMP23
      * (P3[0] * 2. * (V2[5] * V1[3] - V2[3] * V1[5]) + (P3[1] * 2. * (V2[2] *
      V1[5] - V2[5] * V1[2]) + 2. * (P3[3] * (V2[3] * V1[2] - V2[2] *
      V1[3])))))) + (P3[0] * (V1[3] * (V2[4] * 2. * (P1[2] - P2[2]) + 2. *
      (V2[5] * (P2[3] - P1[3]))) + V2[3] * (V1[4] * 2. * (P2[2] - P1[2]) + 2. *
      (V1[5] * (P1[3] - P2[3])))) + (P3[1] * (V1[2] * (V2[4] * 2. * (P2[2] -
      P1[2]) + 2. * (V2[5] * (P1[3] - P2[3]))) + V2[2] * (V1[4] * 2. * (P1[2] -
      P2[2]) + 2. * (V1[5] * (P2[3] - P1[3])))) + (P3[2] * (V1[2] * 2. * V2[3]
      * (P1[2] - P2[2]) + 2. * (V1[3] * V2[2] * (P2[2] - P1[2]))) + P3[3] *
      (V1[2] * 2. * V2[3] * (P2[3] - P1[3]) + 2. * (V1[3] * V2[2] * (P1[3] -
      P2[3])))))));
  T3[5] = denom * 0.500000000 * cI * (OM3 * (P3[0] * (TMP22 * (P3[0] * 2. *
      (V2[4] * V1[3] - V2[3] * V1[4]) + (P3[1] * 2. * (V2[2] * V1[4] - V2[4] *
      V1[2]) + 2. * (P3[2] * (V2[3] * V1[2] - V2[2] * V1[3])))) + (TMP23 *
      (P3[0] * 2. * (V2[3] * V1[4] - V2[4] * V1[3]) + (P3[1] * 2. * (V2[4] *
      V1[2] - V2[2] * V1[4]) + 2. * (P3[2] * (V2[2] * V1[3] - V2[3] * V1[2]))))
      + 1.333333333 * (P3[3] * (TMP49 - TMP50)))) + P3[3] * (TMP22 * (P3[1] *
      2. * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] * 2. * (V2[3] * V1[5] -
      V2[5] * V1[3]) + 2. * (P3[3] * (V2[4] * V1[3] - V2[3] * V1[4])))) + TMP23
      * (P3[1] * 2. * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] * 2. * (V2[5] *
      V1[3] - V2[3] * V1[5]) + 2. * (P3[3] * (V2[3] * V1[4] - V2[4] *
      V1[3])))))) + (P3[1] * (V1[4] * (V2[2] * 2. * (P1[0] - P2[0]) + 2. *
      (V2[5] * (P1[3] - P2[3]))) + V2[4] * (V1[2] * 2. * (P2[0] - P1[0]) + 2. *
      (V1[5] * (P2[3] - P1[3])))) + (P3[2] * (V1[3] * (V2[2] * 2. * (P2[0] -
      P1[0]) + 2. * (V2[5] * (P2[3] - P1[3]))) + V2[3] * (V1[2] * 2. * (P1[0] -
      P2[0]) + 2. * (V1[5] * (P1[3] - P2[3])))) + (P3[0] * (V1[3] * 2. * V2[4]
      * (P1[0] - P2[0]) + 2. * (V1[4] * V2[3] * (P2[0] - P1[0]))) + P3[3] *
      (V1[3] * 2. * V2[4] * (P1[3] - P2[3]) + 2. * (V1[4] * V2[3] * (P2[3] -
      P1[3])))))));
  T3[9] = denom * 0.500000000 * cI * (OM3 * (P3[1] * (TMP22 * (P3[0] * 2. *
      (V2[4] * V1[3] - V2[3] * V1[4]) + (P3[1] * 2. * (V2[2] * V1[4] - V2[4] *
      V1[2]) + 2. * (P3[2] * (V2[3] * V1[2] - V2[2] * V1[3])))) + (TMP23 *
      (P3[0] * 2. * (V2[3] * V1[4] - V2[4] * V1[3]) + (P3[1] * 2. * (V2[4] *
      V1[2] - V2[2] * V1[4]) + 2. * (P3[2] * (V2[2] * V1[3] - V2[3] * V1[2]))))
      + 1.333333333 * (P3[3] * (TMP49 - TMP50)))) + P3[3] * (TMP22 * (P3[0] *
      2. * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] * 2. * (V2[2] * V1[5] -
      V2[5] * V1[2]) + 2. * (P3[3] * (V2[4] * V1[2] - V2[2] * V1[4])))) + TMP23
      * (P3[0] * 2. * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] * 2. * (V2[5] *
      V1[2] - V2[2] * V1[5]) + 2. * (P3[3] * (V2[2] * V1[4] - V2[4] *
      V1[2])))))) + (P3[0] * (V1[4] * (V2[3] * 2. * (P2[1] - P1[1]) + 2. *
      (V2[5] * (P1[3] - P2[3]))) + V2[4] * (V1[3] * 2. * (P1[1] - P2[1]) + 2. *
      (V1[5] * (P2[3] - P1[3])))) + (P3[2] * (V1[2] * (V2[3] * 2. * (P1[1] -
      P2[1]) + 2. * (V2[5] * (P2[3] - P1[3]))) + V2[2] * (V1[3] * 2. * (P2[1] -
      P1[1]) + 2. * (V1[5] * (P1[3] - P2[3])))) + (P3[1] * (V1[2] * 2. * V2[4]
      * (P2[1] - P1[1]) + 2. * (V1[4] * V2[2] * (P1[1] - P2[1]))) + P3[3] *
      (V1[2] * 2. * V2[4] * (P1[3] - P2[3]) + 2. * (V1[4] * V2[2] * (P2[3] -
      P1[3])))))));
  T3[13] = denom * 0.500000000 * cI * (OM3 * (P3[2] * (TMP22 * (P3[0] * 2. *
      (V2[4] * V1[3] - V2[3] * V1[4]) + (P3[1] * 2. * (V2[2] * V1[4] - V2[4] *
      V1[2]) + 2. * (P3[2] * (V2[3] * V1[2] - V2[2] * V1[3])))) + (TMP23 *
      (P3[0] * 2. * (V2[3] * V1[4] - V2[4] * V1[3]) + (P3[1] * 2. * (V2[4] *
      V1[2] - V2[2] * V1[4]) + 2. * (P3[2] * (V2[2] * V1[3] - V2[3] * V1[2]))))
      + 1.333333333 * (P3[3] * (TMP49 - TMP50)))) + P3[3] * (TMP22 * (P3[0] *
      2. * (V2[3] * V1[5] - V2[5] * V1[3]) + (P3[1] * 2. * (V2[5] * V1[2] -
      V2[2] * V1[5]) + 2. * (P3[3] * (V2[2] * V1[3] - V2[3] * V1[2])))) + TMP23
      * (P3[0] * 2. * (V2[5] * V1[3] - V2[3] * V1[5]) + (P3[1] * 2. * (V2[2] *
      V1[5] - V2[5] * V1[2]) + 2. * (P3[3] * (V2[3] * V1[2] - V2[2] *
      V1[3])))))) + (P3[0] * (V1[3] * (V2[4] * 2. * (P1[2] - P2[2]) + 2. *
      (V2[5] * (P2[3] - P1[3]))) + V2[3] * (V1[4] * 2. * (P2[2] - P1[2]) + 2. *
      (V1[5] * (P1[3] - P2[3])))) + (P3[1] * (V1[2] * (V2[4] * 2. * (P2[2] -
      P1[2]) + 2. * (V2[5] * (P1[3] - P2[3]))) + V2[2] * (V1[4] * 2. * (P1[2] -
      P2[2]) + 2. * (V1[5] * (P2[3] - P1[3])))) + (P3[2] * (V1[2] * 2. * V2[3]
      * (P1[2] - P2[2]) + 2. * (V1[3] * V2[2] * (P2[2] - P1[2]))) + P3[3] *
      (V1[2] * 2. * V2[3] * (P2[3] - P1[3]) + 2. * (V1[3] * V2[2] * (P1[3] -
      P2[3])))))));
  T3[17] = denom * 0.500000000 * cI * (OM3 * P3[3] * (TMP22 * (P3[0] * 4. *
      (V2[4] * V1[3] - V2[3] * V1[4]) + (P3[1] * 4. * (V2[2] * V1[4] - V2[4] *
      V1[2]) + 4. * (P3[2] * (V2[3] * V1[2] - V2[2] * V1[3])))) + (TMP23 *
      (P3[0] * 4. * (V2[3] * V1[4] - V2[4] * V1[3]) + (P3[1] * 4. * (V2[4] *
      V1[2] - V2[2] * V1[4]) + 4. * (P3[2] * (V2[2] * V1[3] - V2[3] * V1[2]))))
      + 1.333333333 * (P3[3] * (TMP49 - TMP50)))) + (P1[3] * (P3[0] * 4. *
      (V2[4] * V1[3] - V2[3] * V1[4]) + (P3[1] * 4. * (V2[2] * V1[4] - V2[4] *
      V1[2]) + 4. * (P3[2] * (V2[3] * V1[2] - V2[2] * V1[3])))) + (P2[3] *
      (P3[0] * 4. * (V2[3] * V1[4] - V2[4] * V1[3]) + (P3[1] * 4. * (V2[4] *
      V1[2] - V2[2] * V1[4]) + 4. * (P3[2] * (V2[2] * V1[3] - V2[3] * V1[2]))))
      + (-1.333333333 * (TMP50) + 1.333333333 * (TMP49)))));
}


void VVT19_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP36; 
  double P3[4]; 
  complex<double> TMP5; 
  complex<double> TMP32; 
  complex<double> TMP9; 
  double P2[4]; 
  complex<double> TMP6; 
  complex<double> TMP33; 
  double OM3; 
  complex<double> TMP3; 
  complex<double> TMP34; 
  double P1[4]; 
  complex<double> TMP23; 
  complex<double> TMP7; 
  complex<double> denom; 
  complex<double> TMP35; 
  complex<double> TMP22; 
  complex<double> TMP0; 
  complex<double> TMP4; 
  complex<double> TMP8; 
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
  TMP22 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  TMP23 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP0 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP3 = (V2[2] * P1[0] - V2[3] * P1[1] - V2[4] * P1[2] - V2[5] * P1[3]); 
  TMP33 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP9 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP8 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP32 = (P1[0] * P1[0] - P1[1] * P1[1] - P1[2] * P1[2] - P1[3] * P1[3]); 
  TMP5 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  TMP4 = (V1[2] * P3[0] - V1[3] * P3[1] - V1[4] * P3[2] - V1[5] * P3[3]); 
  TMP7 = (V1[2] * P2[0] - V1[3] * P2[1] - V1[4] * P2[2] - V1[5] * P2[3]); 
  TMP6 = (V2[2] * P3[0] - V2[3] * P3[1] - V2[4] * P3[2] - V2[5] * P3[3]); 
  TMP36 = (V2[2] * P2[0] - V2[3] * P2[1] - V2[4] * P2[2] - V2[5] * P2[3]); 
  TMP35 = (V1[2] * P1[0] - V1[3] * P1[1] - V1[4] * P1[2] - V1[5] * P1[3]); 
  TMP34 = (P2[0] * P2[0] - P2[1] * P2[1] - P2[2] * P2[2] - P2[3] * P2[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * (OM3 * (P3[0] * (P3[0] * (OM3 * (TMP22 * (TMP22 * (TMP0 * -
      0.666666667 * (+cI * (TMP9 + TMP8)) + 0.666666667 * (TMP3 * (+cI * (TMP7
      + TMP4)))) + (TMP23 * 1.333333333 * (-cI * (TMP3 * TMP7) + cI * (TMP0 *
      TMP9)) + 0.666666667 * (TMP6 * (-cI * (TMP4 * TMP9) + cI * (TMP7 *
      TMP8))))) + TMP23 * (TMP23 * (TMP0 * - 0.666666667 * (+cI * (TMP9 +
      TMP5)) + 0.666666667 * (TMP7 * (+cI * (TMP3 + TMP6)))) + 0.666666667 *
      (TMP4 * (-cI * (TMP6 * TMP9) + cI * (TMP3 * TMP5))))) + (TMP0 * (TMP9 * -
      0.333333333 * (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (-0.333333333
      * cI * (TMP5 * TMP32 + TMP8 * TMP34))) + (TMP3 * (TMP7 * 0.333333333 *
      (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (+0.333333333 * cI * (TMP4
      * TMP34 + TMP5 * TMP35))) + (TMP36 * 0.333333333 * (-cI * (TMP4 * TMP9) +
      cI * (TMP7 * TMP8)) + 0.333333333 * (TMP6 * (-cI * (TMP9 * TMP35) + cI *
      (TMP7 * TMP32))))))) + (TMP9 * (TMP0 * (P1[0] * 2. * (-cI * (TMP22) + cI
      * (TMP23)) + 2. * (P2[0] * (-cI * (TMP23) + cI * (TMP22)))) + (TMP4 *
      (TMP6 * (+cI * (P1[0] + P2[0])) + cI * (V2[2] * TMP22)) + cI * (V1[2] *
      TMP6 * TMP23))) + (TMP3 * (TMP7 * (P1[0] * 2. * (-cI * (TMP23) + cI *
      (TMP22)) + 2. * (P2[0] * (-cI * (TMP22) + cI * (TMP23)))) + (TMP4 * - 1.
      * (+cI * (P1[0] * TMP5) + 2. * cI * (P2[0] * TMP22)) - cI * (V1[2] * TMP5
      * TMP23))) + (TMP8 * (P2[0] * (-cI * (TMP6 * TMP7) + 2. * cI * (TMP0 *
      TMP22)) - cI * (V2[2] * TMP7 * TMP22)) + 2. * (P1[0] * TMP23 * (-cI *
      (TMP6 * TMP7) + cI * (TMP0 * TMP5))))))) + (TMP22 * (TMP22 * (TMP0 * -
      0.333333333 * (+cI * (TMP9 + TMP8)) + 0.333333333 * (TMP3 * (+cI * (TMP7
      + TMP4)))) + (TMP23 * 0.666666667 * (-cI * (TMP3 * TMP7) + cI * (TMP0 *
      TMP9)) + 0.333333333 * (TMP6 * (-cI * (TMP4 * TMP9) + cI * (TMP7 *
      TMP8))))) + TMP23 * (TMP23 * (TMP0 * - 0.333333333 * (+cI * (TMP9 +
      TMP5)) + 0.333333333 * (TMP7 * (+cI * (TMP3 + TMP6)))) + 0.333333333 *
      (TMP4 * (-cI * (TMP6 * TMP9) + cI * (TMP3 * TMP5)))))) + (TMP0 * (TMP9 *
      (P1[0] * (-cI * (P1[0]) + 2. * cI * (P2[0])) + (-0.666666667 * cI *
      (TMP33) - cI * (P2[0] * P2[0]) + 0.333333333 * cI * (TMP32 + TMP34))) +
      (TMP5 * (-cI * (P1[0] * P1[0]) + 0.333333333 * cI * (TMP32)) + TMP8 *
      (-cI * (P2[0] * P2[0]) + 0.333333333 * cI * (TMP34)))) + (TMP3 * (TMP7 *
      (P1[0] * (-2. * cI * (P2[0]) + cI * (P1[0])) + (-0.333333333 * cI *
      (TMP32 + TMP34) + cI * (P2[0] * P2[0]) + 0.666666667 * cI * (TMP33))) +
      (TMP4 * (-0.333333333 * cI * (TMP34) + cI * (P2[0] * P2[0])) + TMP5 *
      (-0.333333333 * cI * (TMP35) + cI * (V1[2] * P1[0])))) + (TMP6 * (TMP7 *
      (-0.333333333 * cI * (TMP32) + cI * (P1[0] * P1[0])) + TMP9 * (-cI *
      (V1[2] * P1[0]) + 0.333333333 * cI * (TMP35))) + (TMP4 * TMP9 * (-cI *
      (V2[2] * P2[0]) + 0.333333333 * cI * (TMP36)) + TMP7 * TMP8 *
      (-0.333333333 * cI * (TMP36) + cI * (V2[2] * P2[0])))))));
  T3[6] = denom * (OM3 * (P3[0] * (P3[1] * (OM3 * (TMP22 * (TMP22 * (TMP0 * -
      0.666666667 * (+cI * (TMP9 + TMP8)) + 0.666666667 * (TMP3 * (+cI * (TMP7
      + TMP4)))) + (TMP23 * 1.333333333 * (-cI * (TMP3 * TMP7) + cI * (TMP0 *
      TMP9)) + 0.666666667 * (TMP6 * (-cI * (TMP4 * TMP9) + cI * (TMP7 *
      TMP8))))) + TMP23 * (TMP23 * (TMP0 * - 0.666666667 * (+cI * (TMP9 +
      TMP5)) + 0.666666667 * (TMP7 * (+cI * (TMP3 + TMP6)))) + 0.666666667 *
      (TMP4 * (-cI * (TMP6 * TMP9) + cI * (TMP3 * TMP5))))) + (TMP0 * (TMP9 * -
      0.333333333 * (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (-0.333333333
      * cI * (TMP5 * TMP32 + TMP8 * TMP34))) + (TMP3 * (TMP7 * 0.333333333 *
      (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (+0.333333333 * cI * (TMP4
      * TMP34 + TMP5 * TMP35))) + (TMP36 * 0.333333333 * (-cI * (TMP4 * TMP9) +
      cI * (TMP7 * TMP8)) + 0.333333333 * (TMP6 * (-cI * (TMP9 * TMP35) + cI *
      (TMP7 * TMP32))))))) + (TMP9 * (TMP0 * (P1[1] * (-cI * (TMP22) + cI *
      (TMP23)) + P2[1] * (-cI * (TMP23) + cI * (TMP22))) + (+0.500000000 * cI *
      (V1[3] * TMP6 * TMP23) + TMP4 * 0.500000000 * (+2. * (TMP6 * 0.500000000
      * (+cI * (P1[1] + P2[1]))) + cI * (V2[3] * TMP22)))) + (TMP3 * (TMP7 *
      (P1[1] * (-cI * (TMP23) + cI * (TMP22)) + P2[1] * (-cI * (TMP22) + cI *
      (TMP23))) + (-0.500000000 * cI * (V1[3] * TMP5 * TMP23) + TMP4 * -
      0.500000000 * (+cI * (P1[1] * TMP5) + 2. * cI * (P2[1] * TMP22)))) +
      (TMP8 * 0.500000000 * (+2. * (P2[1] * 0.500000000 * (-cI * (TMP6 * TMP7)
      + 2. * cI * (TMP0 * TMP22))) - cI * (V2[3] * TMP7 * TMP22)) + P1[1] *
      TMP23 * (-cI * (TMP6 * TMP7) + cI * (TMP0 * TMP5)))))) + P3[1] * (TMP9 *
      (TMP0 * (P1[0] * (-cI * (TMP22) + cI * (TMP23)) + P2[0] * (-cI * (TMP23)
      + cI * (TMP22))) + (+0.500000000 * cI * (V1[2] * TMP6 * TMP23) + TMP4 *
      0.500000000 * (+2. * (TMP6 * 0.500000000 * (+cI * (P1[0] + P2[0]))) + cI
      * (V2[2] * TMP22)))) + (TMP3 * (TMP7 * (P1[0] * (-cI * (TMP23) + cI *
      (TMP22)) + P2[0] * (-cI * (TMP22) + cI * (TMP23))) + (-0.500000000 * cI *
      (V1[2] * TMP5 * TMP23) + TMP4 * - 0.500000000 * (+cI * (P1[0] * TMP5) +
      2. * cI * (P2[0] * TMP22)))) + (TMP8 * 0.500000000 * (+2. * (P2[0] *
      0.500000000 * (-cI * (TMP6 * TMP7) + 2. * cI * (TMP0 * TMP22))) - cI *
      (V2[2] * TMP7 * TMP22)) + P1[0] * TMP23 * (-cI * (TMP6 * TMP7) + cI *
      (TMP0 * TMP5)))))) + (P1[0] * (P1[1] * (TMP0 * - 1. * (+cI * (TMP9 +
      TMP5)) + TMP7 * (+cI * (TMP3 + TMP6))) + (+0.500000000 * (V1[3] * (-cI *
      (TMP6 * TMP9) + cI * (TMP3 * TMP5))) + P2[1] * (-cI * (TMP3 * TMP7) + cI
      * (TMP0 * TMP9)))) + (P2[0] * (P2[1] * (TMP0 * - 1. * (+cI * (TMP9 +
      TMP8)) + TMP3 * (+cI * (TMP7 + TMP4))) + (+0.500000000 * (V2[3] * (-cI *
      (TMP4 * TMP9) + cI * (TMP7 * TMP8))) + P1[1] * (-cI * (TMP3 * TMP7) + cI
      * (TMP0 * TMP9)))) + (+0.500000000 * (P2[1] * V2[2] * (-cI * (TMP4 *
      TMP9) + cI * (TMP7 * TMP8))) + P1[1] * 0.500000000 * V1[2] * (-cI * (TMP6
      * TMP9) + cI * (TMP3 * TMP5))))));
  T3[10] = denom * (OM3 * (P3[0] * (P3[2] * (OM3 * (TMP22 * (TMP22 * (TMP0 * -
      0.666666667 * (+cI * (TMP9 + TMP8)) + 0.666666667 * (TMP3 * (+cI * (TMP7
      + TMP4)))) + (TMP23 * 1.333333333 * (-cI * (TMP3 * TMP7) + cI * (TMP0 *
      TMP9)) + 0.666666667 * (TMP6 * (-cI * (TMP4 * TMP9) + cI * (TMP7 *
      TMP8))))) + TMP23 * (TMP23 * (TMP0 * - 0.666666667 * (+cI * (TMP9 +
      TMP5)) + 0.666666667 * (TMP7 * (+cI * (TMP3 + TMP6)))) + 0.666666667 *
      (TMP4 * (-cI * (TMP6 * TMP9) + cI * (TMP3 * TMP5))))) + (TMP0 * (TMP9 * -
      0.333333333 * (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (-0.333333333
      * cI * (TMP5 * TMP32 + TMP8 * TMP34))) + (TMP3 * (TMP7 * 0.333333333 *
      (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (+0.333333333 * cI * (TMP4
      * TMP34 + TMP5 * TMP35))) + (TMP36 * 0.333333333 * (-cI * (TMP4 * TMP9) +
      cI * (TMP7 * TMP8)) + 0.333333333 * (TMP6 * (-cI * (TMP9 * TMP35) + cI *
      (TMP7 * TMP32))))))) + (TMP9 * (TMP0 * (P1[2] * (-cI * (TMP22) + cI *
      (TMP23)) + P2[2] * (-cI * (TMP23) + cI * (TMP22))) + (+0.500000000 * cI *
      (V1[4] * TMP6 * TMP23) + TMP4 * 0.500000000 * (+2. * (TMP6 * 0.500000000
      * (+cI * (P1[2] + P2[2]))) + cI * (V2[4] * TMP22)))) + (TMP3 * (TMP7 *
      (P1[2] * (-cI * (TMP23) + cI * (TMP22)) + P2[2] * (-cI * (TMP22) + cI *
      (TMP23))) + (-0.500000000 * cI * (V1[4] * TMP5 * TMP23) + TMP4 * -
      0.500000000 * (+cI * (P1[2] * TMP5) + 2. * cI * (P2[2] * TMP22)))) +
      (TMP8 * 0.500000000 * (+2. * (P2[2] * 0.500000000 * (-cI * (TMP6 * TMP7)
      + 2. * cI * (TMP0 * TMP22))) - cI * (V2[4] * TMP7 * TMP22)) + P1[2] *
      TMP23 * (-cI * (TMP6 * TMP7) + cI * (TMP0 * TMP5)))))) + P3[2] * (TMP9 *
      (TMP0 * (P1[0] * (-cI * (TMP22) + cI * (TMP23)) + P2[0] * (-cI * (TMP23)
      + cI * (TMP22))) + (+0.500000000 * cI * (V1[2] * TMP6 * TMP23) + TMP4 *
      0.500000000 * (+2. * (TMP6 * 0.500000000 * (+cI * (P1[0] + P2[0]))) + cI
      * (V2[2] * TMP22)))) + (TMP3 * (TMP7 * (P1[0] * (-cI * (TMP23) + cI *
      (TMP22)) + P2[0] * (-cI * (TMP22) + cI * (TMP23))) + (-0.500000000 * cI *
      (V1[2] * TMP5 * TMP23) + TMP4 * - 0.500000000 * (+cI * (P1[0] * TMP5) +
      2. * cI * (P2[0] * TMP22)))) + (TMP8 * 0.500000000 * (+2. * (P2[0] *
      0.500000000 * (-cI * (TMP6 * TMP7) + 2. * cI * (TMP0 * TMP22))) - cI *
      (V2[2] * TMP7 * TMP22)) + P1[0] * TMP23 * (-cI * (TMP6 * TMP7) + cI *
      (TMP0 * TMP5)))))) + (P1[0] * (P1[2] * (TMP0 * - 1. * (+cI * (TMP9 +
      TMP5)) + TMP7 * (+cI * (TMP3 + TMP6))) + (+0.500000000 * (V1[4] * (-cI *
      (TMP6 * TMP9) + cI * (TMP3 * TMP5))) + P2[2] * (-cI * (TMP3 * TMP7) + cI
      * (TMP0 * TMP9)))) + (P2[0] * (P2[2] * (TMP0 * - 1. * (+cI * (TMP9 +
      TMP8)) + TMP3 * (+cI * (TMP7 + TMP4))) + (+0.500000000 * (V2[4] * (-cI *
      (TMP4 * TMP9) + cI * (TMP7 * TMP8))) + P1[2] * (-cI * (TMP3 * TMP7) + cI
      * (TMP0 * TMP9)))) + (+0.500000000 * (P2[2] * V2[2] * (-cI * (TMP4 *
      TMP9) + cI * (TMP7 * TMP8))) + P1[2] * 0.500000000 * V1[2] * (-cI * (TMP6
      * TMP9) + cI * (TMP3 * TMP5))))));
  T3[14] = denom * (OM3 * (P3[0] * (P3[3] * (OM3 * (TMP22 * (TMP22 * (TMP0 * -
      0.666666667 * (+cI * (TMP9 + TMP8)) + 0.666666667 * (TMP3 * (+cI * (TMP7
      + TMP4)))) + (TMP23 * 1.333333333 * (-cI * (TMP3 * TMP7) + cI * (TMP0 *
      TMP9)) + 0.666666667 * (TMP6 * (-cI * (TMP4 * TMP9) + cI * (TMP7 *
      TMP8))))) + TMP23 * (TMP23 * (TMP0 * - 0.666666667 * (+cI * (TMP9 +
      TMP5)) + 0.666666667 * (TMP7 * (+cI * (TMP3 + TMP6)))) + 0.666666667 *
      (TMP4 * (-cI * (TMP6 * TMP9) + cI * (TMP3 * TMP5))))) + (TMP0 * (TMP9 * -
      0.333333333 * (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (-0.333333333
      * cI * (TMP5 * TMP32 + TMP8 * TMP34))) + (TMP3 * (TMP7 * 0.333333333 *
      (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (+0.333333333 * cI * (TMP4
      * TMP34 + TMP5 * TMP35))) + (TMP36 * 0.333333333 * (-cI * (TMP4 * TMP9) +
      cI * (TMP7 * TMP8)) + 0.333333333 * (TMP6 * (-cI * (TMP9 * TMP35) + cI *
      (TMP7 * TMP32))))))) + (TMP9 * (TMP0 * (P1[3] * (-cI * (TMP22) + cI *
      (TMP23)) + P2[3] * (-cI * (TMP23) + cI * (TMP22))) + (+0.500000000 * cI *
      (V1[5] * TMP6 * TMP23) + TMP4 * 0.500000000 * (+2. * (TMP6 * 0.500000000
      * (+cI * (P1[3] + P2[3]))) + cI * (V2[5] * TMP22)))) + (TMP3 * (TMP7 *
      (P1[3] * (-cI * (TMP23) + cI * (TMP22)) + P2[3] * (-cI * (TMP22) + cI *
      (TMP23))) + (-0.500000000 * cI * (V1[5] * TMP5 * TMP23) + TMP4 * -
      0.500000000 * (+cI * (P1[3] * TMP5) + 2. * cI * (P2[3] * TMP22)))) +
      (TMP8 * 0.500000000 * (+2. * (P2[3] * 0.500000000 * (-cI * (TMP6 * TMP7)
      + 2. * cI * (TMP0 * TMP22))) - cI * (V2[5] * TMP7 * TMP22)) + P1[3] *
      TMP23 * (-cI * (TMP6 * TMP7) + cI * (TMP0 * TMP5)))))) + P3[3] * (TMP9 *
      (TMP0 * (P1[0] * (-cI * (TMP22) + cI * (TMP23)) + P2[0] * (-cI * (TMP23)
      + cI * (TMP22))) + (+0.500000000 * cI * (V1[2] * TMP6 * TMP23) + TMP4 *
      0.500000000 * (+2. * (TMP6 * 0.500000000 * (+cI * (P1[0] + P2[0]))) + cI
      * (V2[2] * TMP22)))) + (TMP3 * (TMP7 * (P1[0] * (-cI * (TMP23) + cI *
      (TMP22)) + P2[0] * (-cI * (TMP22) + cI * (TMP23))) + (-0.500000000 * cI *
      (V1[2] * TMP5 * TMP23) + TMP4 * - 0.500000000 * (+cI * (P1[0] * TMP5) +
      2. * cI * (P2[0] * TMP22)))) + (TMP8 * 0.500000000 * (+2. * (P2[0] *
      0.500000000 * (-cI * (TMP6 * TMP7) + 2. * cI * (TMP0 * TMP22))) - cI *
      (V2[2] * TMP7 * TMP22)) + P1[0] * TMP23 * (-cI * (TMP6 * TMP7) + cI *
      (TMP0 * TMP5)))))) + (P1[0] * (P1[3] * (TMP0 * - 1. * (+cI * (TMP9 +
      TMP5)) + TMP7 * (+cI * (TMP3 + TMP6))) + (+0.500000000 * (V1[5] * (-cI *
      (TMP6 * TMP9) + cI * (TMP3 * TMP5))) + P2[3] * (-cI * (TMP3 * TMP7) + cI
      * (TMP0 * TMP9)))) + (P2[0] * (P2[3] * (TMP0 * - 1. * (+cI * (TMP9 +
      TMP8)) + TMP3 * (+cI * (TMP7 + TMP4))) + (+0.500000000 * (V2[5] * (-cI *
      (TMP4 * TMP9) + cI * (TMP7 * TMP8))) + P1[3] * (-cI * (TMP3 * TMP7) + cI
      * (TMP0 * TMP9)))) + (+0.500000000 * (P2[3] * V2[2] * (-cI * (TMP4 *
      TMP9) + cI * (TMP7 * TMP8))) + P1[3] * 0.500000000 * V1[2] * (-cI * (TMP6
      * TMP9) + cI * (TMP3 * TMP5))))));
  T3[3] = denom * (OM3 * (P3[0] * (P3[1] * (OM3 * (TMP22 * (TMP22 * (TMP0 * -
      0.666666667 * (+cI * (TMP9 + TMP8)) + 0.666666667 * (TMP3 * (+cI * (TMP7
      + TMP4)))) + (TMP23 * 1.333333333 * (-cI * (TMP3 * TMP7) + cI * (TMP0 *
      TMP9)) + 0.666666667 * (TMP6 * (-cI * (TMP4 * TMP9) + cI * (TMP7 *
      TMP8))))) + TMP23 * (TMP23 * (TMP0 * - 0.666666667 * (+cI * (TMP9 +
      TMP5)) + 0.666666667 * (TMP7 * (+cI * (TMP3 + TMP6)))) + 0.666666667 *
      (TMP4 * (-cI * (TMP6 * TMP9) + cI * (TMP3 * TMP5))))) + (TMP0 * (TMP9 * -
      0.333333333 * (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (-0.333333333
      * cI * (TMP5 * TMP32 + TMP8 * TMP34))) + (TMP3 * (TMP7 * 0.333333333 *
      (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (+0.333333333 * cI * (TMP4
      * TMP34 + TMP5 * TMP35))) + (TMP36 * 0.333333333 * (-cI * (TMP4 * TMP9) +
      cI * (TMP7 * TMP8)) + 0.333333333 * (TMP6 * (-cI * (TMP9 * TMP35) + cI *
      (TMP7 * TMP32))))))) + (TMP9 * (TMP0 * (P1[1] * (-cI * (TMP22) + cI *
      (TMP23)) + P2[1] * (-cI * (TMP23) + cI * (TMP22))) + (+0.500000000 * cI *
      (V1[3] * TMP6 * TMP23) + TMP4 * 0.500000000 * (+2. * (TMP6 * 0.500000000
      * (+cI * (P1[1] + P2[1]))) + cI * (V2[3] * TMP22)))) + (TMP3 * (TMP7 *
      (P1[1] * (-cI * (TMP23) + cI * (TMP22)) + P2[1] * (-cI * (TMP22) + cI *
      (TMP23))) + (-0.500000000 * cI * (V1[3] * TMP5 * TMP23) + TMP4 * -
      0.500000000 * (+cI * (P1[1] * TMP5) + 2. * cI * (P2[1] * TMP22)))) +
      (TMP8 * 0.500000000 * (+2. * (P2[1] * 0.500000000 * (-cI * (TMP6 * TMP7)
      + 2. * cI * (TMP0 * TMP22))) - cI * (V2[3] * TMP7 * TMP22)) + P1[1] *
      TMP23 * (-cI * (TMP6 * TMP7) + cI * (TMP0 * TMP5)))))) + P3[1] * (TMP9 *
      (TMP0 * (P1[0] * (-cI * (TMP22) + cI * (TMP23)) + P2[0] * (-cI * (TMP23)
      + cI * (TMP22))) + (+0.500000000 * cI * (V1[2] * TMP6 * TMP23) + TMP4 *
      0.500000000 * (+2. * (TMP6 * 0.500000000 * (+cI * (P1[0] + P2[0]))) + cI
      * (V2[2] * TMP22)))) + (TMP3 * (TMP7 * (P1[0] * (-cI * (TMP23) + cI *
      (TMP22)) + P2[0] * (-cI * (TMP22) + cI * (TMP23))) + (-0.500000000 * cI *
      (V1[2] * TMP5 * TMP23) + TMP4 * - 0.500000000 * (+cI * (P1[0] * TMP5) +
      2. * cI * (P2[0] * TMP22)))) + (TMP8 * 0.500000000 * (+2. * (P2[0] *
      0.500000000 * (-cI * (TMP6 * TMP7) + 2. * cI * (TMP0 * TMP22))) - cI *
      (V2[2] * TMP7 * TMP22)) + P1[0] * TMP23 * (-cI * (TMP6 * TMP7) + cI *
      (TMP0 * TMP5)))))) + (P1[0] * (P1[1] * (TMP0 * - 1. * (+cI * (TMP9 +
      TMP5)) + TMP7 * (+cI * (TMP3 + TMP6))) + (+0.500000000 * (V1[3] * (-cI *
      (TMP6 * TMP9) + cI * (TMP3 * TMP5))) + P2[1] * (-cI * (TMP3 * TMP7) + cI
      * (TMP0 * TMP9)))) + (P2[0] * (P2[1] * (TMP0 * - 1. * (+cI * (TMP9 +
      TMP8)) + TMP3 * (+cI * (TMP7 + TMP4))) + (+0.500000000 * (V2[3] * (-cI *
      (TMP4 * TMP9) + cI * (TMP7 * TMP8))) + P1[1] * (-cI * (TMP3 * TMP7) + cI
      * (TMP0 * TMP9)))) + (+0.500000000 * (P2[1] * V2[2] * (-cI * (TMP4 *
      TMP9) + cI * (TMP7 * TMP8))) + P1[1] * 0.500000000 * V1[2] * (-cI * (TMP6
      * TMP9) + cI * (TMP3 * TMP5))))));
  T3[7] = denom * (OM3 * (P3[1] * (P3[1] * (OM3 * (TMP22 * (TMP22 * (TMP0 * -
      0.666666667 * (+cI * (TMP9 + TMP8)) + 0.666666667 * (TMP3 * (+cI * (TMP7
      + TMP4)))) + (TMP23 * 1.333333333 * (-cI * (TMP3 * TMP7) + cI * (TMP0 *
      TMP9)) + 0.666666667 * (TMP6 * (-cI * (TMP4 * TMP9) + cI * (TMP7 *
      TMP8))))) + TMP23 * (TMP23 * (TMP0 * - 0.666666667 * (+cI * (TMP9 +
      TMP5)) + 0.666666667 * (TMP7 * (+cI * (TMP3 + TMP6)))) + 0.666666667 *
      (TMP4 * (-cI * (TMP6 * TMP9) + cI * (TMP3 * TMP5))))) + (TMP0 * (TMP9 * -
      0.333333333 * (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (-0.333333333
      * cI * (TMP5 * TMP32 + TMP8 * TMP34))) + (TMP3 * (TMP7 * 0.333333333 *
      (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (+0.333333333 * cI * (TMP4
      * TMP34 + TMP5 * TMP35))) + (TMP36 * 0.333333333 * (-cI * (TMP4 * TMP9) +
      cI * (TMP7 * TMP8)) + 0.333333333 * (TMP6 * (-cI * (TMP9 * TMP35) + cI *
      (TMP7 * TMP32))))))) + (TMP9 * (TMP0 * (P1[1] * 2. * (-cI * (TMP22) + cI
      * (TMP23)) + 2. * (P2[1] * (-cI * (TMP23) + cI * (TMP22)))) + (TMP4 *
      (TMP6 * (+cI * (P1[1] + P2[1])) + cI * (V2[3] * TMP22)) + cI * (V1[3] *
      TMP6 * TMP23))) + (TMP3 * (TMP7 * (P1[1] * 2. * (-cI * (TMP23) + cI *
      (TMP22)) + 2. * (P2[1] * (-cI * (TMP22) + cI * (TMP23)))) + (TMP4 * - 1.
      * (+cI * (P1[1] * TMP5) + 2. * cI * (P2[1] * TMP22)) - cI * (V1[3] * TMP5
      * TMP23))) + (TMP8 * (P2[1] * (-cI * (TMP6 * TMP7) + 2. * cI * (TMP0 *
      TMP22)) - cI * (V2[3] * TMP7 * TMP22)) + 2. * (P1[1] * TMP23 * (-cI *
      (TMP6 * TMP7) + cI * (TMP0 * TMP5))))))) + (TMP22 * (TMP22 * (TMP0 *
      0.333333333 * (+cI * (TMP9 + TMP8)) - 0.333333333 * (TMP3 * (+cI * (TMP7
      + TMP4)))) + (TMP23 * 0.666666667 * (-cI * (TMP0 * TMP9) + cI * (TMP3 *
      TMP7)) + 0.333333333 * (TMP6 * (-cI * (TMP7 * TMP8) + cI * (TMP4 *
      TMP9))))) + TMP23 * (TMP23 * (TMP0 * 0.333333333 * (+cI * (TMP9 + TMP5))
      - 0.333333333 * (TMP7 * (+cI * (TMP3 + TMP6)))) + 0.333333333 * (TMP4 *
      (-cI * (TMP3 * TMP5) + cI * (TMP6 * TMP9)))))) + (TMP0 * (TMP9 * (P1[1] *
      (-cI * (P1[1]) + 2. * cI * (P2[1])) + (-0.333333333 * cI * (TMP32 +
      TMP34) - cI * (P2[1] * P2[1]) + 0.666666667 * cI * (TMP33))) + (TMP5 * -
      1. * (+cI * (P1[1] * P1[1]) + 0.333333333 * cI * (TMP32)) - TMP8 * (+cI *
      (P2[1] * P2[1]) + 0.333333333 * cI * (TMP34)))) + (TMP3 * (TMP7 * (P1[1]
      * (-2. * cI * (P2[1]) + cI * (P1[1])) + (-0.666666667 * cI * (TMP33) + cI
      * (P2[1] * P2[1]) + 0.333333333 * cI * (TMP32 + TMP34))) + (TMP4 * (+cI *
      (P2[1] * P2[1]) + 0.333333333 * cI * (TMP34)) + TMP5 * (+cI * (V1[3] *
      P1[1]) + 0.333333333 * cI * (TMP35)))) + (TMP6 * (TMP7 * (+cI * (P1[1] *
      P1[1]) + 0.333333333 * cI * (TMP32)) - TMP9 * (+cI * (V1[3] * P1[1]) +
      0.333333333 * cI * (TMP35))) + (TMP4 * - TMP9 * (+cI * (V2[3] * P2[1]) +
      0.333333333 * cI * (TMP36)) + TMP7 * TMP8 * (+cI * (V2[3] * P2[1]) +
      0.333333333 * cI * (TMP36)))))));
  T3[11] = denom * (OM3 * (P3[1] * (P3[2] * (OM3 * (TMP22 * (TMP22 * (TMP0 * -
      0.666666667 * (+cI * (TMP9 + TMP8)) + 0.666666667 * (TMP3 * (+cI * (TMP7
      + TMP4)))) + (TMP23 * 1.333333333 * (-cI * (TMP3 * TMP7) + cI * (TMP0 *
      TMP9)) + 0.666666667 * (TMP6 * (-cI * (TMP4 * TMP9) + cI * (TMP7 *
      TMP8))))) + TMP23 * (TMP23 * (TMP0 * - 0.666666667 * (+cI * (TMP9 +
      TMP5)) + 0.666666667 * (TMP7 * (+cI * (TMP3 + TMP6)))) + 0.666666667 *
      (TMP4 * (-cI * (TMP6 * TMP9) + cI * (TMP3 * TMP5))))) + (TMP0 * (TMP9 * -
      0.333333333 * (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (-0.333333333
      * cI * (TMP5 * TMP32 + TMP8 * TMP34))) + (TMP3 * (TMP7 * 0.333333333 *
      (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (+0.333333333 * cI * (TMP4
      * TMP34 + TMP5 * TMP35))) + (TMP36 * 0.333333333 * (-cI * (TMP4 * TMP9) +
      cI * (TMP7 * TMP8)) + 0.333333333 * (TMP6 * (-cI * (TMP9 * TMP35) + cI *
      (TMP7 * TMP32))))))) + (TMP9 * (TMP0 * (P1[2] * (-cI * (TMP22) + cI *
      (TMP23)) + P2[2] * (-cI * (TMP23) + cI * (TMP22))) + (+0.500000000 * cI *
      (V1[4] * TMP6 * TMP23) + TMP4 * 0.500000000 * (+2. * (TMP6 * 0.500000000
      * (+cI * (P1[2] + P2[2]))) + cI * (V2[4] * TMP22)))) + (TMP3 * (TMP7 *
      (P1[2] * (-cI * (TMP23) + cI * (TMP22)) + P2[2] * (-cI * (TMP22) + cI *
      (TMP23))) + (-0.500000000 * cI * (V1[4] * TMP5 * TMP23) + TMP4 * -
      0.500000000 * (+cI * (P1[2] * TMP5) + 2. * cI * (P2[2] * TMP22)))) +
      (TMP8 * 0.500000000 * (+2. * (P2[2] * 0.500000000 * (-cI * (TMP6 * TMP7)
      + 2. * cI * (TMP0 * TMP22))) - cI * (V2[4] * TMP7 * TMP22)) + P1[2] *
      TMP23 * (-cI * (TMP6 * TMP7) + cI * (TMP0 * TMP5)))))) + P3[2] * (TMP9 *
      (TMP0 * (P1[1] * (-cI * (TMP22) + cI * (TMP23)) + P2[1] * (-cI * (TMP23)
      + cI * (TMP22))) + (+0.500000000 * cI * (V1[3] * TMP6 * TMP23) + TMP4 *
      0.500000000 * (+2. * (TMP6 * 0.500000000 * (+cI * (P1[1] + P2[1]))) + cI
      * (V2[3] * TMP22)))) + (TMP3 * (TMP7 * (P1[1] * (-cI * (TMP23) + cI *
      (TMP22)) + P2[1] * (-cI * (TMP22) + cI * (TMP23))) + (-0.500000000 * cI *
      (V1[3] * TMP5 * TMP23) + TMP4 * - 0.500000000 * (+cI * (P1[1] * TMP5) +
      2. * cI * (P2[1] * TMP22)))) + (TMP8 * 0.500000000 * (+2. * (P2[1] *
      0.500000000 * (-cI * (TMP6 * TMP7) + 2. * cI * (TMP0 * TMP22))) - cI *
      (V2[3] * TMP7 * TMP22)) + P1[1] * TMP23 * (-cI * (TMP6 * TMP7) + cI *
      (TMP0 * TMP5)))))) + (P1[1] * (P1[2] * (TMP0 * - 1. * (+cI * (TMP9 +
      TMP5)) + TMP7 * (+cI * (TMP3 + TMP6))) + (+0.500000000 * (V1[4] * (-cI *
      (TMP6 * TMP9) + cI * (TMP3 * TMP5))) + P2[2] * (-cI * (TMP3 * TMP7) + cI
      * (TMP0 * TMP9)))) + (P2[1] * (P2[2] * (TMP0 * - 1. * (+cI * (TMP9 +
      TMP8)) + TMP3 * (+cI * (TMP7 + TMP4))) + (+0.500000000 * (V2[4] * (-cI *
      (TMP4 * TMP9) + cI * (TMP7 * TMP8))) + P1[2] * (-cI * (TMP3 * TMP7) + cI
      * (TMP0 * TMP9)))) + (+0.500000000 * (P2[2] * V2[3] * (-cI * (TMP4 *
      TMP9) + cI * (TMP7 * TMP8))) + P1[2] * 0.500000000 * V1[3] * (-cI * (TMP6
      * TMP9) + cI * (TMP3 * TMP5))))));
  T3[15] = denom * (OM3 * (P3[1] * (P3[3] * (OM3 * (TMP22 * (TMP22 * (TMP0 * -
      0.666666667 * (+cI * (TMP9 + TMP8)) + 0.666666667 * (TMP3 * (+cI * (TMP7
      + TMP4)))) + (TMP23 * 1.333333333 * (-cI * (TMP3 * TMP7) + cI * (TMP0 *
      TMP9)) + 0.666666667 * (TMP6 * (-cI * (TMP4 * TMP9) + cI * (TMP7 *
      TMP8))))) + TMP23 * (TMP23 * (TMP0 * - 0.666666667 * (+cI * (TMP9 +
      TMP5)) + 0.666666667 * (TMP7 * (+cI * (TMP3 + TMP6)))) + 0.666666667 *
      (TMP4 * (-cI * (TMP6 * TMP9) + cI * (TMP3 * TMP5))))) + (TMP0 * (TMP9 * -
      0.333333333 * (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (-0.333333333
      * cI * (TMP5 * TMP32 + TMP8 * TMP34))) + (TMP3 * (TMP7 * 0.333333333 *
      (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (+0.333333333 * cI * (TMP4
      * TMP34 + TMP5 * TMP35))) + (TMP36 * 0.333333333 * (-cI * (TMP4 * TMP9) +
      cI * (TMP7 * TMP8)) + 0.333333333 * (TMP6 * (-cI * (TMP9 * TMP35) + cI *
      (TMP7 * TMP32))))))) + (TMP9 * (TMP0 * (P1[3] * (-cI * (TMP22) + cI *
      (TMP23)) + P2[3] * (-cI * (TMP23) + cI * (TMP22))) + (+0.500000000 * cI *
      (V1[5] * TMP6 * TMP23) + TMP4 * 0.500000000 * (+2. * (TMP6 * 0.500000000
      * (+cI * (P1[3] + P2[3]))) + cI * (V2[5] * TMP22)))) + (TMP3 * (TMP7 *
      (P1[3] * (-cI * (TMP23) + cI * (TMP22)) + P2[3] * (-cI * (TMP22) + cI *
      (TMP23))) + (-0.500000000 * cI * (V1[5] * TMP5 * TMP23) + TMP4 * -
      0.500000000 * (+cI * (P1[3] * TMP5) + 2. * cI * (P2[3] * TMP22)))) +
      (TMP8 * 0.500000000 * (+2. * (P2[3] * 0.500000000 * (-cI * (TMP6 * TMP7)
      + 2. * cI * (TMP0 * TMP22))) - cI * (V2[5] * TMP7 * TMP22)) + P1[3] *
      TMP23 * (-cI * (TMP6 * TMP7) + cI * (TMP0 * TMP5)))))) + P3[3] * (TMP9 *
      (TMP0 * (P1[1] * (-cI * (TMP22) + cI * (TMP23)) + P2[1] * (-cI * (TMP23)
      + cI * (TMP22))) + (+0.500000000 * cI * (V1[3] * TMP6 * TMP23) + TMP4 *
      0.500000000 * (+2. * (TMP6 * 0.500000000 * (+cI * (P1[1] + P2[1]))) + cI
      * (V2[3] * TMP22)))) + (TMP3 * (TMP7 * (P1[1] * (-cI * (TMP23) + cI *
      (TMP22)) + P2[1] * (-cI * (TMP22) + cI * (TMP23))) + (-0.500000000 * cI *
      (V1[3] * TMP5 * TMP23) + TMP4 * - 0.500000000 * (+cI * (P1[1] * TMP5) +
      2. * cI * (P2[1] * TMP22)))) + (TMP8 * 0.500000000 * (+2. * (P2[1] *
      0.500000000 * (-cI * (TMP6 * TMP7) + 2. * cI * (TMP0 * TMP22))) - cI *
      (V2[3] * TMP7 * TMP22)) + P1[1] * TMP23 * (-cI * (TMP6 * TMP7) + cI *
      (TMP0 * TMP5)))))) + (P1[1] * (P1[3] * (TMP0 * - 1. * (+cI * (TMP9 +
      TMP5)) + TMP7 * (+cI * (TMP3 + TMP6))) + (+0.500000000 * (V1[5] * (-cI *
      (TMP6 * TMP9) + cI * (TMP3 * TMP5))) + P2[3] * (-cI * (TMP3 * TMP7) + cI
      * (TMP0 * TMP9)))) + (P2[1] * (P2[3] * (TMP0 * - 1. * (+cI * (TMP9 +
      TMP8)) + TMP3 * (+cI * (TMP7 + TMP4))) + (+0.500000000 * (V2[5] * (-cI *
      (TMP4 * TMP9) + cI * (TMP7 * TMP8))) + P1[3] * (-cI * (TMP3 * TMP7) + cI
      * (TMP0 * TMP9)))) + (+0.500000000 * (P2[3] * V2[3] * (-cI * (TMP4 *
      TMP9) + cI * (TMP7 * TMP8))) + P1[3] * 0.500000000 * V1[3] * (-cI * (TMP6
      * TMP9) + cI * (TMP3 * TMP5))))));
  T3[4] = denom * (OM3 * (P3[0] * (P3[2] * (OM3 * (TMP22 * (TMP22 * (TMP0 * -
      0.666666667 * (+cI * (TMP9 + TMP8)) + 0.666666667 * (TMP3 * (+cI * (TMP7
      + TMP4)))) + (TMP23 * 1.333333333 * (-cI * (TMP3 * TMP7) + cI * (TMP0 *
      TMP9)) + 0.666666667 * (TMP6 * (-cI * (TMP4 * TMP9) + cI * (TMP7 *
      TMP8))))) + TMP23 * (TMP23 * (TMP0 * - 0.666666667 * (+cI * (TMP9 +
      TMP5)) + 0.666666667 * (TMP7 * (+cI * (TMP3 + TMP6)))) + 0.666666667 *
      (TMP4 * (-cI * (TMP6 * TMP9) + cI * (TMP3 * TMP5))))) + (TMP0 * (TMP9 * -
      0.333333333 * (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (-0.333333333
      * cI * (TMP5 * TMP32 + TMP8 * TMP34))) + (TMP3 * (TMP7 * 0.333333333 *
      (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (+0.333333333 * cI * (TMP4
      * TMP34 + TMP5 * TMP35))) + (TMP36 * 0.333333333 * (-cI * (TMP4 * TMP9) +
      cI * (TMP7 * TMP8)) + 0.333333333 * (TMP6 * (-cI * (TMP9 * TMP35) + cI *
      (TMP7 * TMP32))))))) + (TMP9 * (TMP0 * (P1[2] * (-cI * (TMP22) + cI *
      (TMP23)) + P2[2] * (-cI * (TMP23) + cI * (TMP22))) + (+0.500000000 * cI *
      (V1[4] * TMP6 * TMP23) + TMP4 * 0.500000000 * (+2. * (TMP6 * 0.500000000
      * (+cI * (P1[2] + P2[2]))) + cI * (V2[4] * TMP22)))) + (TMP3 * (TMP7 *
      (P1[2] * (-cI * (TMP23) + cI * (TMP22)) + P2[2] * (-cI * (TMP22) + cI *
      (TMP23))) + (-0.500000000 * cI * (V1[4] * TMP5 * TMP23) + TMP4 * -
      0.500000000 * (+cI * (P1[2] * TMP5) + 2. * cI * (P2[2] * TMP22)))) +
      (TMP8 * 0.500000000 * (+2. * (P2[2] * 0.500000000 * (-cI * (TMP6 * TMP7)
      + 2. * cI * (TMP0 * TMP22))) - cI * (V2[4] * TMP7 * TMP22)) + P1[2] *
      TMP23 * (-cI * (TMP6 * TMP7) + cI * (TMP0 * TMP5)))))) + P3[2] * (TMP9 *
      (TMP0 * (P1[0] * (-cI * (TMP22) + cI * (TMP23)) + P2[0] * (-cI * (TMP23)
      + cI * (TMP22))) + (+0.500000000 * cI * (V1[2] * TMP6 * TMP23) + TMP4 *
      0.500000000 * (+2. * (TMP6 * 0.500000000 * (+cI * (P1[0] + P2[0]))) + cI
      * (V2[2] * TMP22)))) + (TMP3 * (TMP7 * (P1[0] * (-cI * (TMP23) + cI *
      (TMP22)) + P2[0] * (-cI * (TMP22) + cI * (TMP23))) + (-0.500000000 * cI *
      (V1[2] * TMP5 * TMP23) + TMP4 * - 0.500000000 * (+cI * (P1[0] * TMP5) +
      2. * cI * (P2[0] * TMP22)))) + (TMP8 * 0.500000000 * (+2. * (P2[0] *
      0.500000000 * (-cI * (TMP6 * TMP7) + 2. * cI * (TMP0 * TMP22))) - cI *
      (V2[2] * TMP7 * TMP22)) + P1[0] * TMP23 * (-cI * (TMP6 * TMP7) + cI *
      (TMP0 * TMP5)))))) + (P1[0] * (P1[2] * (TMP0 * - 1. * (+cI * (TMP9 +
      TMP5)) + TMP7 * (+cI * (TMP3 + TMP6))) + (+0.500000000 * (V1[4] * (-cI *
      (TMP6 * TMP9) + cI * (TMP3 * TMP5))) + P2[2] * (-cI * (TMP3 * TMP7) + cI
      * (TMP0 * TMP9)))) + (P2[0] * (P2[2] * (TMP0 * - 1. * (+cI * (TMP9 +
      TMP8)) + TMP3 * (+cI * (TMP7 + TMP4))) + (+0.500000000 * (V2[4] * (-cI *
      (TMP4 * TMP9) + cI * (TMP7 * TMP8))) + P1[2] * (-cI * (TMP3 * TMP7) + cI
      * (TMP0 * TMP9)))) + (+0.500000000 * (P2[2] * V2[2] * (-cI * (TMP4 *
      TMP9) + cI * (TMP7 * TMP8))) + P1[2] * 0.500000000 * V1[2] * (-cI * (TMP6
      * TMP9) + cI * (TMP3 * TMP5))))));
  T3[8] = denom * (OM3 * (P3[1] * (P3[2] * (OM3 * (TMP22 * (TMP22 * (TMP0 * -
      0.666666667 * (+cI * (TMP9 + TMP8)) + 0.666666667 * (TMP3 * (+cI * (TMP7
      + TMP4)))) + (TMP23 * 1.333333333 * (-cI * (TMP3 * TMP7) + cI * (TMP0 *
      TMP9)) + 0.666666667 * (TMP6 * (-cI * (TMP4 * TMP9) + cI * (TMP7 *
      TMP8))))) + TMP23 * (TMP23 * (TMP0 * - 0.666666667 * (+cI * (TMP9 +
      TMP5)) + 0.666666667 * (TMP7 * (+cI * (TMP3 + TMP6)))) + 0.666666667 *
      (TMP4 * (-cI * (TMP6 * TMP9) + cI * (TMP3 * TMP5))))) + (TMP0 * (TMP9 * -
      0.333333333 * (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (-0.333333333
      * cI * (TMP5 * TMP32 + TMP8 * TMP34))) + (TMP3 * (TMP7 * 0.333333333 *
      (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (+0.333333333 * cI * (TMP4
      * TMP34 + TMP5 * TMP35))) + (TMP36 * 0.333333333 * (-cI * (TMP4 * TMP9) +
      cI * (TMP7 * TMP8)) + 0.333333333 * (TMP6 * (-cI * (TMP9 * TMP35) + cI *
      (TMP7 * TMP32))))))) + (TMP9 * (TMP0 * (P1[2] * (-cI * (TMP22) + cI *
      (TMP23)) + P2[2] * (-cI * (TMP23) + cI * (TMP22))) + (+0.500000000 * cI *
      (V1[4] * TMP6 * TMP23) + TMP4 * 0.500000000 * (+2. * (TMP6 * 0.500000000
      * (+cI * (P1[2] + P2[2]))) + cI * (V2[4] * TMP22)))) + (TMP3 * (TMP7 *
      (P1[2] * (-cI * (TMP23) + cI * (TMP22)) + P2[2] * (-cI * (TMP22) + cI *
      (TMP23))) + (-0.500000000 * cI * (V1[4] * TMP5 * TMP23) + TMP4 * -
      0.500000000 * (+cI * (P1[2] * TMP5) + 2. * cI * (P2[2] * TMP22)))) +
      (TMP8 * 0.500000000 * (+2. * (P2[2] * 0.500000000 * (-cI * (TMP6 * TMP7)
      + 2. * cI * (TMP0 * TMP22))) - cI * (V2[4] * TMP7 * TMP22)) + P1[2] *
      TMP23 * (-cI * (TMP6 * TMP7) + cI * (TMP0 * TMP5)))))) + P3[2] * (TMP9 *
      (TMP0 * (P1[1] * (-cI * (TMP22) + cI * (TMP23)) + P2[1] * (-cI * (TMP23)
      + cI * (TMP22))) + (+0.500000000 * cI * (V1[3] * TMP6 * TMP23) + TMP4 *
      0.500000000 * (+2. * (TMP6 * 0.500000000 * (+cI * (P1[1] + P2[1]))) + cI
      * (V2[3] * TMP22)))) + (TMP3 * (TMP7 * (P1[1] * (-cI * (TMP23) + cI *
      (TMP22)) + P2[1] * (-cI * (TMP22) + cI * (TMP23))) + (-0.500000000 * cI *
      (V1[3] * TMP5 * TMP23) + TMP4 * - 0.500000000 * (+cI * (P1[1] * TMP5) +
      2. * cI * (P2[1] * TMP22)))) + (TMP8 * 0.500000000 * (+2. * (P2[1] *
      0.500000000 * (-cI * (TMP6 * TMP7) + 2. * cI * (TMP0 * TMP22))) - cI *
      (V2[3] * TMP7 * TMP22)) + P1[1] * TMP23 * (-cI * (TMP6 * TMP7) + cI *
      (TMP0 * TMP5)))))) + (P1[1] * (P1[2] * (TMP0 * - 1. * (+cI * (TMP9 +
      TMP5)) + TMP7 * (+cI * (TMP3 + TMP6))) + (+0.500000000 * (V1[4] * (-cI *
      (TMP6 * TMP9) + cI * (TMP3 * TMP5))) + P2[2] * (-cI * (TMP3 * TMP7) + cI
      * (TMP0 * TMP9)))) + (P2[1] * (P2[2] * (TMP0 * - 1. * (+cI * (TMP9 +
      TMP8)) + TMP3 * (+cI * (TMP7 + TMP4))) + (+0.500000000 * (V2[4] * (-cI *
      (TMP4 * TMP9) + cI * (TMP7 * TMP8))) + P1[2] * (-cI * (TMP3 * TMP7) + cI
      * (TMP0 * TMP9)))) + (+0.500000000 * (P2[2] * V2[3] * (-cI * (TMP4 *
      TMP9) + cI * (TMP7 * TMP8))) + P1[2] * 0.500000000 * V1[3] * (-cI * (TMP6
      * TMP9) + cI * (TMP3 * TMP5))))));
  T3[12] = denom * (OM3 * (P3[2] * (P3[2] * (OM3 * (TMP22 * (TMP22 * (TMP0 * -
      0.666666667 * (+cI * (TMP9 + TMP8)) + 0.666666667 * (TMP3 * (+cI * (TMP7
      + TMP4)))) + (TMP23 * 1.333333333 * (-cI * (TMP3 * TMP7) + cI * (TMP0 *
      TMP9)) + 0.666666667 * (TMP6 * (-cI * (TMP4 * TMP9) + cI * (TMP7 *
      TMP8))))) + TMP23 * (TMP23 * (TMP0 * - 0.666666667 * (+cI * (TMP9 +
      TMP5)) + 0.666666667 * (TMP7 * (+cI * (TMP3 + TMP6)))) + 0.666666667 *
      (TMP4 * (-cI * (TMP6 * TMP9) + cI * (TMP3 * TMP5))))) + (TMP0 * (TMP9 * -
      0.333333333 * (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (-0.333333333
      * cI * (TMP5 * TMP32 + TMP8 * TMP34))) + (TMP3 * (TMP7 * 0.333333333 *
      (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (+0.333333333 * cI * (TMP4
      * TMP34 + TMP5 * TMP35))) + (TMP36 * 0.333333333 * (-cI * (TMP4 * TMP9) +
      cI * (TMP7 * TMP8)) + 0.333333333 * (TMP6 * (-cI * (TMP9 * TMP35) + cI *
      (TMP7 * TMP32))))))) + (TMP9 * (TMP0 * (P1[2] * 2. * (-cI * (TMP22) + cI
      * (TMP23)) + 2. * (P2[2] * (-cI * (TMP23) + cI * (TMP22)))) + (TMP4 *
      (TMP6 * (+cI * (P1[2] + P2[2])) + cI * (V2[4] * TMP22)) + cI * (V1[4] *
      TMP6 * TMP23))) + (TMP3 * (TMP7 * (P1[2] * 2. * (-cI * (TMP23) + cI *
      (TMP22)) + 2. * (P2[2] * (-cI * (TMP22) + cI * (TMP23)))) + (TMP4 * - 1.
      * (+cI * (P1[2] * TMP5) + 2. * cI * (P2[2] * TMP22)) - cI * (V1[4] * TMP5
      * TMP23))) + (TMP8 * (P2[2] * (-cI * (TMP6 * TMP7) + 2. * cI * (TMP0 *
      TMP22)) - cI * (V2[4] * TMP7 * TMP22)) + 2. * (P1[2] * TMP23 * (-cI *
      (TMP6 * TMP7) + cI * (TMP0 * TMP5))))))) + (TMP22 * (TMP22 * (TMP0 *
      0.333333333 * (+cI * (TMP9 + TMP8)) - 0.333333333 * (TMP3 * (+cI * (TMP7
      + TMP4)))) + (TMP23 * 0.666666667 * (-cI * (TMP0 * TMP9) + cI * (TMP3 *
      TMP7)) + 0.333333333 * (TMP6 * (-cI * (TMP7 * TMP8) + cI * (TMP4 *
      TMP9))))) + TMP23 * (TMP23 * (TMP0 * 0.333333333 * (+cI * (TMP9 + TMP5))
      - 0.333333333 * (TMP7 * (+cI * (TMP3 + TMP6)))) + 0.333333333 * (TMP4 *
      (-cI * (TMP3 * TMP5) + cI * (TMP6 * TMP9)))))) + (TMP0 * (TMP9 * (P1[2] *
      (-cI * (P1[2]) + 2. * cI * (P2[2])) + (-0.333333333 * cI * (TMP32 +
      TMP34) - cI * (P2[2] * P2[2]) + 0.666666667 * cI * (TMP33))) + (TMP5 * -
      1. * (+cI * (P1[2] * P1[2]) + 0.333333333 * cI * (TMP32)) - TMP8 * (+cI *
      (P2[2] * P2[2]) + 0.333333333 * cI * (TMP34)))) + (TMP3 * (TMP7 * (P1[2]
      * (-2. * cI * (P2[2]) + cI * (P1[2])) + (-0.666666667 * cI * (TMP33) + cI
      * (P2[2] * P2[2]) + 0.333333333 * cI * (TMP32 + TMP34))) + (TMP4 * (+cI *
      (P2[2] * P2[2]) + 0.333333333 * cI * (TMP34)) + TMP5 * (+cI * (V1[4] *
      P1[2]) + 0.333333333 * cI * (TMP35)))) + (TMP6 * (TMP7 * (+cI * (P1[2] *
      P1[2]) + 0.333333333 * cI * (TMP32)) - TMP9 * (+cI * (V1[4] * P1[2]) +
      0.333333333 * cI * (TMP35))) + (TMP4 * - TMP9 * (+cI * (V2[4] * P2[2]) +
      0.333333333 * cI * (TMP36)) + TMP7 * TMP8 * (+cI * (V2[4] * P2[2]) +
      0.333333333 * cI * (TMP36)))))));
  T3[16] = denom * (OM3 * (P3[2] * (P3[3] * (OM3 * (TMP22 * (TMP22 * (TMP0 * -
      0.666666667 * (+cI * (TMP9 + TMP8)) + 0.666666667 * (TMP3 * (+cI * (TMP7
      + TMP4)))) + (TMP23 * 1.333333333 * (-cI * (TMP3 * TMP7) + cI * (TMP0 *
      TMP9)) + 0.666666667 * (TMP6 * (-cI * (TMP4 * TMP9) + cI * (TMP7 *
      TMP8))))) + TMP23 * (TMP23 * (TMP0 * - 0.666666667 * (+cI * (TMP9 +
      TMP5)) + 0.666666667 * (TMP7 * (+cI * (TMP3 + TMP6)))) + 0.666666667 *
      (TMP4 * (-cI * (TMP6 * TMP9) + cI * (TMP3 * TMP5))))) + (TMP0 * (TMP9 * -
      0.333333333 * (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (-0.333333333
      * cI * (TMP5 * TMP32 + TMP8 * TMP34))) + (TMP3 * (TMP7 * 0.333333333 *
      (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (+0.333333333 * cI * (TMP4
      * TMP34 + TMP5 * TMP35))) + (TMP36 * 0.333333333 * (-cI * (TMP4 * TMP9) +
      cI * (TMP7 * TMP8)) + 0.333333333 * (TMP6 * (-cI * (TMP9 * TMP35) + cI *
      (TMP7 * TMP32))))))) + (TMP9 * (TMP0 * (P1[3] * (-cI * (TMP22) + cI *
      (TMP23)) + P2[3] * (-cI * (TMP23) + cI * (TMP22))) + (+0.500000000 * cI *
      (V1[5] * TMP6 * TMP23) + TMP4 * 0.500000000 * (+2. * (TMP6 * 0.500000000
      * (+cI * (P1[3] + P2[3]))) + cI * (V2[5] * TMP22)))) + (TMP3 * (TMP7 *
      (P1[3] * (-cI * (TMP23) + cI * (TMP22)) + P2[3] * (-cI * (TMP22) + cI *
      (TMP23))) + (-0.500000000 * cI * (V1[5] * TMP5 * TMP23) + TMP4 * -
      0.500000000 * (+cI * (P1[3] * TMP5) + 2. * cI * (P2[3] * TMP22)))) +
      (TMP8 * 0.500000000 * (+2. * (P2[3] * 0.500000000 * (-cI * (TMP6 * TMP7)
      + 2. * cI * (TMP0 * TMP22))) - cI * (V2[5] * TMP7 * TMP22)) + P1[3] *
      TMP23 * (-cI * (TMP6 * TMP7) + cI * (TMP0 * TMP5)))))) + P3[3] * (TMP9 *
      (TMP0 * (P1[2] * (-cI * (TMP22) + cI * (TMP23)) + P2[2] * (-cI * (TMP23)
      + cI * (TMP22))) + (+0.500000000 * cI * (V1[4] * TMP6 * TMP23) + TMP4 *
      0.500000000 * (+2. * (TMP6 * 0.500000000 * (+cI * (P1[2] + P2[2]))) + cI
      * (V2[4] * TMP22)))) + (TMP3 * (TMP7 * (P1[2] * (-cI * (TMP23) + cI *
      (TMP22)) + P2[2] * (-cI * (TMP22) + cI * (TMP23))) + (-0.500000000 * cI *
      (V1[4] * TMP5 * TMP23) + TMP4 * - 0.500000000 * (+cI * (P1[2] * TMP5) +
      2. * cI * (P2[2] * TMP22)))) + (TMP8 * 0.500000000 * (+2. * (P2[2] *
      0.500000000 * (-cI * (TMP6 * TMP7) + 2. * cI * (TMP0 * TMP22))) - cI *
      (V2[4] * TMP7 * TMP22)) + P1[2] * TMP23 * (-cI * (TMP6 * TMP7) + cI *
      (TMP0 * TMP5)))))) + (P1[2] * (P1[3] * (TMP0 * - 1. * (+cI * (TMP9 +
      TMP5)) + TMP7 * (+cI * (TMP3 + TMP6))) + (+0.500000000 * (V1[5] * (-cI *
      (TMP6 * TMP9) + cI * (TMP3 * TMP5))) + P2[3] * (-cI * (TMP3 * TMP7) + cI
      * (TMP0 * TMP9)))) + (P2[2] * (P2[3] * (TMP0 * - 1. * (+cI * (TMP9 +
      TMP8)) + TMP3 * (+cI * (TMP7 + TMP4))) + (+0.500000000 * (V2[5] * (-cI *
      (TMP4 * TMP9) + cI * (TMP7 * TMP8))) + P1[3] * (-cI * (TMP3 * TMP7) + cI
      * (TMP0 * TMP9)))) + (+0.500000000 * (P2[3] * V2[4] * (-cI * (TMP4 *
      TMP9) + cI * (TMP7 * TMP8))) + P1[3] * 0.500000000 * V1[4] * (-cI * (TMP6
      * TMP9) + cI * (TMP3 * TMP5))))));
  T3[5] = denom * (OM3 * (P3[0] * (P3[3] * (OM3 * (TMP22 * (TMP22 * (TMP0 * -
      0.666666667 * (+cI * (TMP9 + TMP8)) + 0.666666667 * (TMP3 * (+cI * (TMP7
      + TMP4)))) + (TMP23 * 1.333333333 * (-cI * (TMP3 * TMP7) + cI * (TMP0 *
      TMP9)) + 0.666666667 * (TMP6 * (-cI * (TMP4 * TMP9) + cI * (TMP7 *
      TMP8))))) + TMP23 * (TMP23 * (TMP0 * - 0.666666667 * (+cI * (TMP9 +
      TMP5)) + 0.666666667 * (TMP7 * (+cI * (TMP3 + TMP6)))) + 0.666666667 *
      (TMP4 * (-cI * (TMP6 * TMP9) + cI * (TMP3 * TMP5))))) + (TMP0 * (TMP9 * -
      0.333333333 * (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (-0.333333333
      * cI * (TMP5 * TMP32 + TMP8 * TMP34))) + (TMP3 * (TMP7 * 0.333333333 *
      (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (+0.333333333 * cI * (TMP4
      * TMP34 + TMP5 * TMP35))) + (TMP36 * 0.333333333 * (-cI * (TMP4 * TMP9) +
      cI * (TMP7 * TMP8)) + 0.333333333 * (TMP6 * (-cI * (TMP9 * TMP35) + cI *
      (TMP7 * TMP32))))))) + (TMP9 * (TMP0 * (P1[3] * (-cI * (TMP22) + cI *
      (TMP23)) + P2[3] * (-cI * (TMP23) + cI * (TMP22))) + (+0.500000000 * cI *
      (V1[5] * TMP6 * TMP23) + TMP4 * 0.500000000 * (+2. * (TMP6 * 0.500000000
      * (+cI * (P1[3] + P2[3]))) + cI * (V2[5] * TMP22)))) + (TMP3 * (TMP7 *
      (P1[3] * (-cI * (TMP23) + cI * (TMP22)) + P2[3] * (-cI * (TMP22) + cI *
      (TMP23))) + (-0.500000000 * cI * (V1[5] * TMP5 * TMP23) + TMP4 * -
      0.500000000 * (+cI * (P1[3] * TMP5) + 2. * cI * (P2[3] * TMP22)))) +
      (TMP8 * 0.500000000 * (+2. * (P2[3] * 0.500000000 * (-cI * (TMP6 * TMP7)
      + 2. * cI * (TMP0 * TMP22))) - cI * (V2[5] * TMP7 * TMP22)) + P1[3] *
      TMP23 * (-cI * (TMP6 * TMP7) + cI * (TMP0 * TMP5)))))) + P3[3] * (TMP9 *
      (TMP0 * (P1[0] * (-cI * (TMP22) + cI * (TMP23)) + P2[0] * (-cI * (TMP23)
      + cI * (TMP22))) + (+0.500000000 * cI * (V1[2] * TMP6 * TMP23) + TMP4 *
      0.500000000 * (+2. * (TMP6 * 0.500000000 * (+cI * (P1[0] + P2[0]))) + cI
      * (V2[2] * TMP22)))) + (TMP3 * (TMP7 * (P1[0] * (-cI * (TMP23) + cI *
      (TMP22)) + P2[0] * (-cI * (TMP22) + cI * (TMP23))) + (-0.500000000 * cI *
      (V1[2] * TMP5 * TMP23) + TMP4 * - 0.500000000 * (+cI * (P1[0] * TMP5) +
      2. * cI * (P2[0] * TMP22)))) + (TMP8 * 0.500000000 * (+2. * (P2[0] *
      0.500000000 * (-cI * (TMP6 * TMP7) + 2. * cI * (TMP0 * TMP22))) - cI *
      (V2[2] * TMP7 * TMP22)) + P1[0] * TMP23 * (-cI * (TMP6 * TMP7) + cI *
      (TMP0 * TMP5)))))) + (P1[0] * (P1[3] * (TMP0 * - 1. * (+cI * (TMP9 +
      TMP5)) + TMP7 * (+cI * (TMP3 + TMP6))) + (+0.500000000 * (V1[5] * (-cI *
      (TMP6 * TMP9) + cI * (TMP3 * TMP5))) + P2[3] * (-cI * (TMP3 * TMP7) + cI
      * (TMP0 * TMP9)))) + (P2[0] * (P2[3] * (TMP0 * - 1. * (+cI * (TMP9 +
      TMP8)) + TMP3 * (+cI * (TMP7 + TMP4))) + (+0.500000000 * (V2[5] * (-cI *
      (TMP4 * TMP9) + cI * (TMP7 * TMP8))) + P1[3] * (-cI * (TMP3 * TMP7) + cI
      * (TMP0 * TMP9)))) + (+0.500000000 * (P2[3] * V2[2] * (-cI * (TMP4 *
      TMP9) + cI * (TMP7 * TMP8))) + P1[3] * 0.500000000 * V1[2] * (-cI * (TMP6
      * TMP9) + cI * (TMP3 * TMP5))))));
  T3[9] = denom * (OM3 * (P3[1] * (P3[3] * (OM3 * (TMP22 * (TMP22 * (TMP0 * -
      0.666666667 * (+cI * (TMP9 + TMP8)) + 0.666666667 * (TMP3 * (+cI * (TMP7
      + TMP4)))) + (TMP23 * 1.333333333 * (-cI * (TMP3 * TMP7) + cI * (TMP0 *
      TMP9)) + 0.666666667 * (TMP6 * (-cI * (TMP4 * TMP9) + cI * (TMP7 *
      TMP8))))) + TMP23 * (TMP23 * (TMP0 * - 0.666666667 * (+cI * (TMP9 +
      TMP5)) + 0.666666667 * (TMP7 * (+cI * (TMP3 + TMP6)))) + 0.666666667 *
      (TMP4 * (-cI * (TMP6 * TMP9) + cI * (TMP3 * TMP5))))) + (TMP0 * (TMP9 * -
      0.333333333 * (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (-0.333333333
      * cI * (TMP5 * TMP32 + TMP8 * TMP34))) + (TMP3 * (TMP7 * 0.333333333 *
      (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (+0.333333333 * cI * (TMP4
      * TMP34 + TMP5 * TMP35))) + (TMP36 * 0.333333333 * (-cI * (TMP4 * TMP9) +
      cI * (TMP7 * TMP8)) + 0.333333333 * (TMP6 * (-cI * (TMP9 * TMP35) + cI *
      (TMP7 * TMP32))))))) + (TMP9 * (TMP0 * (P1[3] * (-cI * (TMP22) + cI *
      (TMP23)) + P2[3] * (-cI * (TMP23) + cI * (TMP22))) + (+0.500000000 * cI *
      (V1[5] * TMP6 * TMP23) + TMP4 * 0.500000000 * (+2. * (TMP6 * 0.500000000
      * (+cI * (P1[3] + P2[3]))) + cI * (V2[5] * TMP22)))) + (TMP3 * (TMP7 *
      (P1[3] * (-cI * (TMP23) + cI * (TMP22)) + P2[3] * (-cI * (TMP22) + cI *
      (TMP23))) + (-0.500000000 * cI * (V1[5] * TMP5 * TMP23) + TMP4 * -
      0.500000000 * (+cI * (P1[3] * TMP5) + 2. * cI * (P2[3] * TMP22)))) +
      (TMP8 * 0.500000000 * (+2. * (P2[3] * 0.500000000 * (-cI * (TMP6 * TMP7)
      + 2. * cI * (TMP0 * TMP22))) - cI * (V2[5] * TMP7 * TMP22)) + P1[3] *
      TMP23 * (-cI * (TMP6 * TMP7) + cI * (TMP0 * TMP5)))))) + P3[3] * (TMP9 *
      (TMP0 * (P1[1] * (-cI * (TMP22) + cI * (TMP23)) + P2[1] * (-cI * (TMP23)
      + cI * (TMP22))) + (+0.500000000 * cI * (V1[3] * TMP6 * TMP23) + TMP4 *
      0.500000000 * (+2. * (TMP6 * 0.500000000 * (+cI * (P1[1] + P2[1]))) + cI
      * (V2[3] * TMP22)))) + (TMP3 * (TMP7 * (P1[1] * (-cI * (TMP23) + cI *
      (TMP22)) + P2[1] * (-cI * (TMP22) + cI * (TMP23))) + (-0.500000000 * cI *
      (V1[3] * TMP5 * TMP23) + TMP4 * - 0.500000000 * (+cI * (P1[1] * TMP5) +
      2. * cI * (P2[1] * TMP22)))) + (TMP8 * 0.500000000 * (+2. * (P2[1] *
      0.500000000 * (-cI * (TMP6 * TMP7) + 2. * cI * (TMP0 * TMP22))) - cI *
      (V2[3] * TMP7 * TMP22)) + P1[1] * TMP23 * (-cI * (TMP6 * TMP7) + cI *
      (TMP0 * TMP5)))))) + (P1[1] * (P1[3] * (TMP0 * - 1. * (+cI * (TMP9 +
      TMP5)) + TMP7 * (+cI * (TMP3 + TMP6))) + (+0.500000000 * (V1[5] * (-cI *
      (TMP6 * TMP9) + cI * (TMP3 * TMP5))) + P2[3] * (-cI * (TMP3 * TMP7) + cI
      * (TMP0 * TMP9)))) + (P2[1] * (P2[3] * (TMP0 * - 1. * (+cI * (TMP9 +
      TMP8)) + TMP3 * (+cI * (TMP7 + TMP4))) + (+0.500000000 * (V2[5] * (-cI *
      (TMP4 * TMP9) + cI * (TMP7 * TMP8))) + P1[3] * (-cI * (TMP3 * TMP7) + cI
      * (TMP0 * TMP9)))) + (+0.500000000 * (P2[3] * V2[3] * (-cI * (TMP4 *
      TMP9) + cI * (TMP7 * TMP8))) + P1[3] * 0.500000000 * V1[3] * (-cI * (TMP6
      * TMP9) + cI * (TMP3 * TMP5))))));
  T3[13] = denom * (OM3 * (P3[2] * (P3[3] * (OM3 * (TMP22 * (TMP22 * (TMP0 * -
      0.666666667 * (+cI * (TMP9 + TMP8)) + 0.666666667 * (TMP3 * (+cI * (TMP7
      + TMP4)))) + (TMP23 * 1.333333333 * (-cI * (TMP3 * TMP7) + cI * (TMP0 *
      TMP9)) + 0.666666667 * (TMP6 * (-cI * (TMP4 * TMP9) + cI * (TMP7 *
      TMP8))))) + TMP23 * (TMP23 * (TMP0 * - 0.666666667 * (+cI * (TMP9 +
      TMP5)) + 0.666666667 * (TMP7 * (+cI * (TMP3 + TMP6)))) + 0.666666667 *
      (TMP4 * (-cI * (TMP6 * TMP9) + cI * (TMP3 * TMP5))))) + (TMP0 * (TMP9 * -
      0.333333333 * (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (-0.333333333
      * cI * (TMP5 * TMP32 + TMP8 * TMP34))) + (TMP3 * (TMP7 * 0.333333333 *
      (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (+0.333333333 * cI * (TMP4
      * TMP34 + TMP5 * TMP35))) + (TMP36 * 0.333333333 * (-cI * (TMP4 * TMP9) +
      cI * (TMP7 * TMP8)) + 0.333333333 * (TMP6 * (-cI * (TMP9 * TMP35) + cI *
      (TMP7 * TMP32))))))) + (TMP9 * (TMP0 * (P1[3] * (-cI * (TMP22) + cI *
      (TMP23)) + P2[3] * (-cI * (TMP23) + cI * (TMP22))) + (+0.500000000 * cI *
      (V1[5] * TMP6 * TMP23) + TMP4 * 0.500000000 * (+2. * (TMP6 * 0.500000000
      * (+cI * (P1[3] + P2[3]))) + cI * (V2[5] * TMP22)))) + (TMP3 * (TMP7 *
      (P1[3] * (-cI * (TMP23) + cI * (TMP22)) + P2[3] * (-cI * (TMP22) + cI *
      (TMP23))) + (-0.500000000 * cI * (V1[5] * TMP5 * TMP23) + TMP4 * -
      0.500000000 * (+cI * (P1[3] * TMP5) + 2. * cI * (P2[3] * TMP22)))) +
      (TMP8 * 0.500000000 * (+2. * (P2[3] * 0.500000000 * (-cI * (TMP6 * TMP7)
      + 2. * cI * (TMP0 * TMP22))) - cI * (V2[5] * TMP7 * TMP22)) + P1[3] *
      TMP23 * (-cI * (TMP6 * TMP7) + cI * (TMP0 * TMP5)))))) + P3[3] * (TMP9 *
      (TMP0 * (P1[2] * (-cI * (TMP22) + cI * (TMP23)) + P2[2] * (-cI * (TMP23)
      + cI * (TMP22))) + (+0.500000000 * cI * (V1[4] * TMP6 * TMP23) + TMP4 *
      0.500000000 * (+2. * (TMP6 * 0.500000000 * (+cI * (P1[2] + P2[2]))) + cI
      * (V2[4] * TMP22)))) + (TMP3 * (TMP7 * (P1[2] * (-cI * (TMP23) + cI *
      (TMP22)) + P2[2] * (-cI * (TMP22) + cI * (TMP23))) + (-0.500000000 * cI *
      (V1[4] * TMP5 * TMP23) + TMP4 * - 0.500000000 * (+cI * (P1[2] * TMP5) +
      2. * cI * (P2[2] * TMP22)))) + (TMP8 * 0.500000000 * (+2. * (P2[2] *
      0.500000000 * (-cI * (TMP6 * TMP7) + 2. * cI * (TMP0 * TMP22))) - cI *
      (V2[4] * TMP7 * TMP22)) + P1[2] * TMP23 * (-cI * (TMP6 * TMP7) + cI *
      (TMP0 * TMP5)))))) + (P1[2] * (P1[3] * (TMP0 * - 1. * (+cI * (TMP9 +
      TMP5)) + TMP7 * (+cI * (TMP3 + TMP6))) + (+0.500000000 * (V1[5] * (-cI *
      (TMP6 * TMP9) + cI * (TMP3 * TMP5))) + P2[3] * (-cI * (TMP3 * TMP7) + cI
      * (TMP0 * TMP9)))) + (P2[2] * (P2[3] * (TMP0 * - 1. * (+cI * (TMP9 +
      TMP8)) + TMP3 * (+cI * (TMP7 + TMP4))) + (+0.500000000 * (V2[5] * (-cI *
      (TMP4 * TMP9) + cI * (TMP7 * TMP8))) + P1[3] * (-cI * (TMP3 * TMP7) + cI
      * (TMP0 * TMP9)))) + (+0.500000000 * (P2[3] * V2[4] * (-cI * (TMP4 *
      TMP9) + cI * (TMP7 * TMP8))) + P1[3] * 0.500000000 * V1[4] * (-cI * (TMP6
      * TMP9) + cI * (TMP3 * TMP5))))));
  T3[17] = denom * (OM3 * (P3[3] * (P3[3] * (OM3 * (TMP22 * (TMP22 * (TMP0 * -
      0.666666667 * (+cI * (TMP9 + TMP8)) + 0.666666667 * (TMP3 * (+cI * (TMP7
      + TMP4)))) + (TMP23 * 1.333333333 * (-cI * (TMP3 * TMP7) + cI * (TMP0 *
      TMP9)) + 0.666666667 * (TMP6 * (-cI * (TMP4 * TMP9) + cI * (TMP7 *
      TMP8))))) + TMP23 * (TMP23 * (TMP0 * - 0.666666667 * (+cI * (TMP9 +
      TMP5)) + 0.666666667 * (TMP7 * (+cI * (TMP3 + TMP6)))) + 0.666666667 *
      (TMP4 * (-cI * (TMP6 * TMP9) + cI * (TMP3 * TMP5))))) + (TMP0 * (TMP9 * -
      0.333333333 * (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (-0.333333333
      * cI * (TMP5 * TMP32 + TMP8 * TMP34))) + (TMP3 * (TMP7 * 0.333333333 *
      (-2. * cI * (TMP33) + cI * (TMP32 + TMP34)) + (+0.333333333 * cI * (TMP4
      * TMP34 + TMP5 * TMP35))) + (TMP36 * 0.333333333 * (-cI * (TMP4 * TMP9) +
      cI * (TMP7 * TMP8)) + 0.333333333 * (TMP6 * (-cI * (TMP9 * TMP35) + cI *
      (TMP7 * TMP32))))))) + (TMP9 * (TMP0 * (P1[3] * 2. * (-cI * (TMP22) + cI
      * (TMP23)) + 2. * (P2[3] * (-cI * (TMP23) + cI * (TMP22)))) + (TMP4 *
      (TMP6 * (+cI * (P1[3] + P2[3])) + cI * (V2[5] * TMP22)) + cI * (V1[5] *
      TMP6 * TMP23))) + (TMP3 * (TMP7 * (P1[3] * 2. * (-cI * (TMP23) + cI *
      (TMP22)) + 2. * (P2[3] * (-cI * (TMP22) + cI * (TMP23)))) + (TMP4 * - 1.
      * (+cI * (P1[3] * TMP5) + 2. * cI * (P2[3] * TMP22)) - cI * (V1[5] * TMP5
      * TMP23))) + (TMP8 * (P2[3] * (-cI * (TMP6 * TMP7) + 2. * cI * (TMP0 *
      TMP22)) - cI * (V2[5] * TMP7 * TMP22)) + 2. * (P1[3] * TMP23 * (-cI *
      (TMP6 * TMP7) + cI * (TMP0 * TMP5))))))) + (TMP22 * (TMP22 * (TMP0 *
      0.333333333 * (+cI * (TMP9 + TMP8)) - 0.333333333 * (TMP3 * (+cI * (TMP7
      + TMP4)))) + (TMP23 * 0.666666667 * (-cI * (TMP0 * TMP9) + cI * (TMP3 *
      TMP7)) + 0.333333333 * (TMP6 * (-cI * (TMP7 * TMP8) + cI * (TMP4 *
      TMP9))))) + TMP23 * (TMP23 * (TMP0 * 0.333333333 * (+cI * (TMP9 + TMP5))
      - 0.333333333 * (TMP7 * (+cI * (TMP3 + TMP6)))) + 0.333333333 * (TMP4 *
      (-cI * (TMP3 * TMP5) + cI * (TMP6 * TMP9)))))) + (TMP0 * (TMP9 * (P1[3] *
      (-cI * (P1[3]) + 2. * cI * (P2[3])) + (-0.333333333 * cI * (TMP32 +
      TMP34) - cI * (P2[3] * P2[3]) + 0.666666667 * cI * (TMP33))) + (TMP5 * -
      1. * (+cI * (P1[3] * P1[3]) + 0.333333333 * cI * (TMP32)) - TMP8 * (+cI *
      (P2[3] * P2[3]) + 0.333333333 * cI * (TMP34)))) + (TMP3 * (TMP7 * (P1[3]
      * (-2. * cI * (P2[3]) + cI * (P1[3])) + (-0.666666667 * cI * (TMP33) + cI
      * (P2[3] * P2[3]) + 0.333333333 * cI * (TMP32 + TMP34))) + (TMP4 * (+cI *
      (P2[3] * P2[3]) + 0.333333333 * cI * (TMP34)) + TMP5 * (+cI * (V1[5] *
      P1[3]) + 0.333333333 * cI * (TMP35)))) + (TMP6 * (TMP7 * (+cI * (P1[3] *
      P1[3]) + 0.333333333 * cI * (TMP32)) - TMP9 * (+cI * (V1[5] * P1[3]) +
      0.333333333 * cI * (TMP35))) + (TMP4 * - TMP9 * (+cI * (V2[5] * P2[3]) +
      0.333333333 * cI * (TMP36)) + TMP7 * TMP8 * (+cI * (V2[5] * P2[3]) +
      0.333333333 * cI * (TMP36)))))));
}


void VVT17_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> TMP0; 
  double P2[4]; 
  complex<double> TMP7; 
  complex<double> TMP26; 
  complex<double> TMP25; 
  complex<double> TMP9; 
  complex<double> TMP3; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  TMP25 = (P1[0] * - 1. * (P2[1] * T3[6] + P2[2] * T3[10] + P2[3] * T3[14] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[11] + P2[3] *
      T3[15] - P2[0] * T3[3]) + (P1[2] * (P2[1] * T3[8] + P2[2] * T3[12] +
      P2[3] * T3[16] - P2[0] * T3[4]) + P1[3] * (P2[1] * T3[9] + P2[2] * T3[13]
      + P2[3] * T3[17] - P2[0] * T3[5]))));
  TMP26 = (P1[0] * - 1. * (P2[1] * T3[3] + P2[2] * T3[4] + P2[3] * T3[5] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[8] + P2[3] * T3[9]
      - P2[0] * T3[6]) + (P1[2] * (P2[1] * T3[11] + P2[2] * T3[12] + P2[3] *
      T3[13] - P2[0] * T3[10]) + P1[3] * (P2[1] * T3[15] + P2[2] * T3[16] +
      P2[3] * T3[17] - P2[0] * T3[14]))));
  TMP9 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP7 = (V1[2] * P2[0] - V1[3] * P2[1] - V1[4] * P2[2] - V1[5] * P2[3]); 
  TMP0 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP3 = (V2[2] * P1[0] - V2[3] * P1[1] - V2[4] * P1[2] - V2[5] * P1[3]); 
  vertex = COUP * (TMP0 * TMP9 * (+cI * (TMP25 + TMP26)) - TMP3 * TMP7 * (+cI *
      (TMP25 + TMP26)));
}


void VVT22_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP5; 
  double P1[4]; 
  complex<double> TMP22; 
  complex<double> TMP0; 
  double P2[4]; 
  complex<double> TMP23; 
  complex<double> TMP7; 
  double P3[4]; 
  complex<double> TMP6; 
  complex<double> denom; 
  double OM3; 
  complex<double> TMP4; 
  complex<double> TMP33; 
  complex<double> TMP3; 
  complex<double> TMP8; 
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
  TMP22 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  TMP23 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP33 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP8 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP5 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  TMP4 = (V1[2] * P3[0] - V1[3] * P3[1] - V1[4] * P3[2] - V1[5] * P3[3]); 
  TMP7 = (V1[2] * P2[0] - V1[3] * P2[1] - V1[4] * P2[2] - V1[5] * P2[3]); 
  TMP6 = (V2[2] * P3[0] - V2[3] * P3[1] - V2[4] * P3[2] - V2[5] * P3[3]); 
  TMP0 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP3 = (V2[2] * P1[0] - V2[3] * P1[1] - V2[4] * P1[2] - V2[5] * P1[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * (TMP4 * (OM3 * (TMP6 * (P3[0] * (P3[0] * (OM3 * (TMP22 *
      1.333333333 * (-cI * (TMP8) + cI * (TMP23)) + 1.333333333 * (TMP5 * (-cI
      * (TMP23) + cI * (TMP8)))) + 0.666666667 * cI * (TMP33)) + (P1[0] * 2. *
      (-cI * (TMP22) + cI * (TMP5)) + 2. * (P2[0] * (-cI * (TMP23) + cI *
      (TMP8))))) + (TMP22 * 0.666666667 * (-cI * (TMP8) + cI * (TMP23)) +
      0.666666667 * (TMP5 * (-cI * (TMP23) + cI * (TMP8))))) + P3[0] * TMP5 *
      (V2[2] * 2. * (-cI * (TMP8) + cI * (TMP23)) - 0.666666667 * cI * (P3[0] *
      TMP3))) + (TMP5 * 2. * (-cI * (V2[2] * P1[0]) + 0.333333333 * cI *
      (TMP3)) + 2. * (TMP6 * (-0.333333333 * cI * (TMP33) + cI * (P1[0] *
      P2[0]))))) + TMP8 * (TMP6 * (OM3 * P3[0] * (V1[2] * 2. * (-cI * (TMP5) +
      cI * (TMP22)) - 0.666666667 * cI * (P3[0] * TMP7)) + (-2. * cI * (V1[2] *
      P2[0]) + 0.666666667 * cI * (TMP7))) + TMP5 * (TMP0 * 0.666666667 * (+cI
      * (P3[0] * P3[0] * OM3) + - 1.000000000 * cI) + 2. * cI * (V2[2] *
      V1[2]))));
  T3[3] = denom * (OM3 * (TMP4 * (TMP6 * (P3[0] * (P3[1] * (OM3 * (TMP22 *
      1.333333333 * (-cI * (TMP8) + cI * (TMP23)) + 1.333333333 * (TMP5 * (-cI
      * (TMP23) + cI * (TMP8)))) + 0.666666667 * cI * (TMP33)) + (P1[1] * (-cI
      * (TMP22) + cI * (TMP5)) + P2[1] * (-cI * (TMP23) + cI * (TMP8)))) +
      P3[1] * (P1[0] * (-cI * (TMP22) + cI * (TMP5)) + P2[0] * (-cI * (TMP23) +
      cI * (TMP8)))) + TMP5 * (P3[0] * (V2[3] * (-cI * (TMP8) + cI * (TMP23)) -
      0.666666667 * cI * (P3[1] * TMP3)) + P3[1] * V2[2] * (-cI * (TMP8) + cI *
      (TMP23)))) + TMP8 * (TMP6 * (P3[0] * (V1[3] * (-cI * (TMP5) + cI *
      (TMP22)) - 0.666666667 * cI * (P3[1] * TMP7)) + P3[1] * V1[2] * (-cI *
      (TMP5) + cI * (TMP22))) + 0.666666667 * cI * (TMP0 * P3[0] * P3[1] *
      TMP5))) + (TMP4 * (P1[0] * (-cI * (V2[3] * TMP5) + cI * (P2[1] * TMP6)) +
      P1[1] * (-cI * (V2[2] * TMP5) + cI * (P2[0] * TMP6))) + TMP8 * (TMP5 *
      (+cI * (V2[2] * V1[3] + V2[3] * V1[2])) - TMP6 * (+cI * (V1[2] * P2[1] +
      V1[3] * P2[0])))));
  T3[4] = denom * (OM3 * (TMP4 * (TMP6 * (P3[0] * (P3[2] * (OM3 * (TMP22 *
      1.333333333 * (-cI * (TMP8) + cI * (TMP23)) + 1.333333333 * (TMP5 * (-cI
      * (TMP23) + cI * (TMP8)))) + 0.666666667 * cI * (TMP33)) + (P1[2] * (-cI
      * (TMP22) + cI * (TMP5)) + P2[2] * (-cI * (TMP23) + cI * (TMP8)))) +
      P3[2] * (P1[0] * (-cI * (TMP22) + cI * (TMP5)) + P2[0] * (-cI * (TMP23) +
      cI * (TMP8)))) + TMP5 * (P3[0] * (V2[4] * (-cI * (TMP8) + cI * (TMP23)) -
      0.666666667 * cI * (P3[2] * TMP3)) + P3[2] * V2[2] * (-cI * (TMP8) + cI *
      (TMP23)))) + TMP8 * (TMP6 * (P3[0] * (V1[4] * (-cI * (TMP5) + cI *
      (TMP22)) - 0.666666667 * cI * (P3[2] * TMP7)) + P3[2] * V1[2] * (-cI *
      (TMP5) + cI * (TMP22))) + 0.666666667 * cI * (TMP0 * P3[0] * P3[2] *
      TMP5))) + (TMP4 * (P1[0] * (-cI * (V2[4] * TMP5) + cI * (P2[2] * TMP6)) +
      P1[2] * (-cI * (V2[2] * TMP5) + cI * (P2[0] * TMP6))) + TMP8 * (TMP5 *
      (+cI * (V2[2] * V1[4] + V2[4] * V1[2])) - TMP6 * (+cI * (V1[2] * P2[2] +
      V1[4] * P2[0])))));
  T3[5] = denom * (OM3 * (TMP4 * (TMP6 * (P3[0] * (P3[3] * (OM3 * (TMP22 *
      1.333333333 * (-cI * (TMP8) + cI * (TMP23)) + 1.333333333 * (TMP5 * (-cI
      * (TMP23) + cI * (TMP8)))) + 0.666666667 * cI * (TMP33)) + (P1[3] * (-cI
      * (TMP22) + cI * (TMP5)) + P2[3] * (-cI * (TMP23) + cI * (TMP8)))) +
      P3[3] * (P1[0] * (-cI * (TMP22) + cI * (TMP5)) + P2[0] * (-cI * (TMP23) +
      cI * (TMP8)))) + TMP5 * (P3[0] * (V2[5] * (-cI * (TMP8) + cI * (TMP23)) -
      0.666666667 * cI * (P3[3] * TMP3)) + P3[3] * V2[2] * (-cI * (TMP8) + cI *
      (TMP23)))) + TMP8 * (TMP6 * (P3[0] * (V1[5] * (-cI * (TMP5) + cI *
      (TMP22)) - 0.666666667 * cI * (P3[3] * TMP7)) + P3[3] * V1[2] * (-cI *
      (TMP5) + cI * (TMP22))) + 0.666666667 * cI * (TMP0 * P3[0] * P3[3] *
      TMP5))) + (TMP4 * (P1[0] * (-cI * (V2[5] * TMP5) + cI * (P2[3] * TMP6)) +
      P1[3] * (-cI * (V2[2] * TMP5) + cI * (P2[0] * TMP6))) + TMP8 * (TMP5 *
      (+cI * (V2[2] * V1[5] + V2[5] * V1[2])) - TMP6 * (+cI * (V1[2] * P2[3] +
      V1[5] * P2[0])))));
  T3[6] = denom * (OM3 * (TMP4 * (TMP6 * (P3[0] * (P3[1] * (OM3 * (TMP22 *
      1.333333333 * (-cI * (TMP8) + cI * (TMP23)) + 1.333333333 * (TMP5 * (-cI
      * (TMP23) + cI * (TMP8)))) + 0.666666667 * cI * (TMP33)) + (P1[1] * (-cI
      * (TMP22) + cI * (TMP5)) + P2[1] * (-cI * (TMP23) + cI * (TMP8)))) +
      P3[1] * (P1[0] * (-cI * (TMP22) + cI * (TMP5)) + P2[0] * (-cI * (TMP23) +
      cI * (TMP8)))) + TMP5 * (P3[0] * (V2[3] * (-cI * (TMP8) + cI * (TMP23)) -
      0.666666667 * cI * (P3[1] * TMP3)) + P3[1] * V2[2] * (-cI * (TMP8) + cI *
      (TMP23)))) + TMP8 * (TMP6 * (P3[0] * (V1[3] * (-cI * (TMP5) + cI *
      (TMP22)) - 0.666666667 * cI * (P3[1] * TMP7)) + P3[1] * V1[2] * (-cI *
      (TMP5) + cI * (TMP22))) + 0.666666667 * cI * (TMP0 * P3[0] * P3[1] *
      TMP5))) + (TMP4 * (P1[0] * (-cI * (V2[3] * TMP5) + cI * (P2[1] * TMP6)) +
      P1[1] * (-cI * (V2[2] * TMP5) + cI * (P2[0] * TMP6))) + TMP8 * (TMP5 *
      (+cI * (V2[3] * V1[2] + V2[2] * V1[3])) - TMP6 * (+cI * (V1[3] * P2[0] +
      V1[2] * P2[1])))));
  T3[7] = denom * (TMP4 * (OM3 * (TMP6 * (P3[1] * (P3[1] * (OM3 * (TMP22 *
      1.333333333 * (-cI * (TMP8) + cI * (TMP23)) + 1.333333333 * (TMP5 * (-cI
      * (TMP23) + cI * (TMP8)))) + 0.666666667 * cI * (TMP33)) + (P1[1] * 2. *
      (-cI * (TMP22) + cI * (TMP5)) + 2. * (P2[1] * (-cI * (TMP23) + cI *
      (TMP8))))) + (TMP22 * 0.666666667 * (-cI * (TMP23) + cI * (TMP8)) +
      0.666666667 * (TMP5 * (-cI * (TMP8) + cI * (TMP23))))) + P3[1] * TMP5 *
      (V2[3] * 2. * (-cI * (TMP8) + cI * (TMP23)) - 0.666666667 * cI * (P3[1] *
      TMP3))) + (TMP5 * - 2. * (+cI * (V2[3] * P1[1]) + 0.333333333 * cI *
      (TMP3)) + 2. * (TMP6 * (+cI * (P1[1] * P2[1]) + 0.333333333 * cI *
      (TMP33))))) + TMP8 * (TMP6 * (OM3 * P3[1] * (V1[3] * 2. * (-cI * (TMP5) +
      cI * (TMP22)) - 0.666666667 * cI * (P3[1] * TMP7)) + (-0.666666667 * cI *
      (TMP7) - 2. * cI * (V1[3] * P2[1]))) + TMP5 * (TMP0 * 0.666666667 * (+cI
      * (P3[1] * P3[1] * OM3) + 1.000000000 * cI) + 2. * cI * (V2[3] *
      V1[3]))));
  T3[8] = denom * (OM3 * (TMP4 * (TMP6 * (P3[1] * (P3[2] * (OM3 * (TMP22 *
      1.333333333 * (-cI * (TMP8) + cI * (TMP23)) + 1.333333333 * (TMP5 * (-cI
      * (TMP23) + cI * (TMP8)))) + 0.666666667 * cI * (TMP33)) + (P1[2] * (-cI
      * (TMP22) + cI * (TMP5)) + P2[2] * (-cI * (TMP23) + cI * (TMP8)))) +
      P3[2] * (P1[1] * (-cI * (TMP22) + cI * (TMP5)) + P2[1] * (-cI * (TMP23) +
      cI * (TMP8)))) + TMP5 * (P3[1] * (V2[4] * (-cI * (TMP8) + cI * (TMP23)) -
      0.666666667 * cI * (P3[2] * TMP3)) + P3[2] * V2[3] * (-cI * (TMP8) + cI *
      (TMP23)))) + TMP8 * (TMP6 * (P3[1] * (V1[4] * (-cI * (TMP5) + cI *
      (TMP22)) - 0.666666667 * cI * (P3[2] * TMP7)) + P3[2] * V1[3] * (-cI *
      (TMP5) + cI * (TMP22))) + 0.666666667 * cI * (TMP0 * P3[1] * P3[2] *
      TMP5))) + (TMP4 * (P1[1] * (-cI * (V2[4] * TMP5) + cI * (P2[2] * TMP6)) +
      P1[2] * (-cI * (V2[3] * TMP5) + cI * (P2[1] * TMP6))) + TMP8 * (TMP5 *
      (+cI * (V2[3] * V1[4] + V2[4] * V1[3])) - TMP6 * (+cI * (V1[3] * P2[2] +
      V1[4] * P2[1])))));
  T3[9] = denom * (OM3 * (TMP4 * (TMP6 * (P3[1] * (P3[3] * (OM3 * (TMP22 *
      1.333333333 * (-cI * (TMP8) + cI * (TMP23)) + 1.333333333 * (TMP5 * (-cI
      * (TMP23) + cI * (TMP8)))) + 0.666666667 * cI * (TMP33)) + (P1[3] * (-cI
      * (TMP22) + cI * (TMP5)) + P2[3] * (-cI * (TMP23) + cI * (TMP8)))) +
      P3[3] * (P1[1] * (-cI * (TMP22) + cI * (TMP5)) + P2[1] * (-cI * (TMP23) +
      cI * (TMP8)))) + TMP5 * (P3[1] * (V2[5] * (-cI * (TMP8) + cI * (TMP23)) -
      0.666666667 * cI * (P3[3] * TMP3)) + P3[3] * V2[3] * (-cI * (TMP8) + cI *
      (TMP23)))) + TMP8 * (TMP6 * (P3[1] * (V1[5] * (-cI * (TMP5) + cI *
      (TMP22)) - 0.666666667 * cI * (P3[3] * TMP7)) + P3[3] * V1[3] * (-cI *
      (TMP5) + cI * (TMP22))) + 0.666666667 * cI * (TMP0 * P3[1] * P3[3] *
      TMP5))) + (TMP4 * (P1[1] * (-cI * (V2[5] * TMP5) + cI * (P2[3] * TMP6)) +
      P1[3] * (-cI * (V2[3] * TMP5) + cI * (P2[1] * TMP6))) + TMP8 * (TMP5 *
      (+cI * (V2[3] * V1[5] + V2[5] * V1[3])) - TMP6 * (+cI * (V1[3] * P2[3] +
      V1[5] * P2[1])))));
  T3[10] = denom * (OM3 * (TMP4 * (TMP6 * (P3[0] * (P3[2] * (OM3 * (TMP22 *
      1.333333333 * (-cI * (TMP8) + cI * (TMP23)) + 1.333333333 * (TMP5 * (-cI
      * (TMP23) + cI * (TMP8)))) + 0.666666667 * cI * (TMP33)) + (P1[2] * (-cI
      * (TMP22) + cI * (TMP5)) + P2[2] * (-cI * (TMP23) + cI * (TMP8)))) +
      P3[2] * (P1[0] * (-cI * (TMP22) + cI * (TMP5)) + P2[0] * (-cI * (TMP23) +
      cI * (TMP8)))) + TMP5 * (P3[0] * (V2[4] * (-cI * (TMP8) + cI * (TMP23)) -
      0.666666667 * cI * (P3[2] * TMP3)) + P3[2] * V2[2] * (-cI * (TMP8) + cI *
      (TMP23)))) + TMP8 * (TMP6 * (P3[0] * (V1[4] * (-cI * (TMP5) + cI *
      (TMP22)) - 0.666666667 * cI * (P3[2] * TMP7)) + P3[2] * V1[2] * (-cI *
      (TMP5) + cI * (TMP22))) + 0.666666667 * cI * (TMP0 * P3[0] * P3[2] *
      TMP5))) + (TMP4 * (P1[0] * (-cI * (V2[4] * TMP5) + cI * (P2[2] * TMP6)) +
      P1[2] * (-cI * (V2[2] * TMP5) + cI * (P2[0] * TMP6))) + TMP8 * (TMP5 *
      (+cI * (V2[4] * V1[2] + V2[2] * V1[4])) - TMP6 * (+cI * (V1[4] * P2[0] +
      V1[2] * P2[2])))));
  T3[11] = denom * (OM3 * (TMP4 * (TMP6 * (P3[1] * (P3[2] * (OM3 * (TMP22 *
      1.333333333 * (-cI * (TMP8) + cI * (TMP23)) + 1.333333333 * (TMP5 * (-cI
      * (TMP23) + cI * (TMP8)))) + 0.666666667 * cI * (TMP33)) + (P1[2] * (-cI
      * (TMP22) + cI * (TMP5)) + P2[2] * (-cI * (TMP23) + cI * (TMP8)))) +
      P3[2] * (P1[1] * (-cI * (TMP22) + cI * (TMP5)) + P2[1] * (-cI * (TMP23) +
      cI * (TMP8)))) + TMP5 * (P3[1] * (V2[4] * (-cI * (TMP8) + cI * (TMP23)) -
      0.666666667 * cI * (P3[2] * TMP3)) + P3[2] * V2[3] * (-cI * (TMP8) + cI *
      (TMP23)))) + TMP8 * (TMP6 * (P3[1] * (V1[4] * (-cI * (TMP5) + cI *
      (TMP22)) - 0.666666667 * cI * (P3[2] * TMP7)) + P3[2] * V1[3] * (-cI *
      (TMP5) + cI * (TMP22))) + 0.666666667 * cI * (TMP0 * P3[1] * P3[2] *
      TMP5))) + (TMP4 * (P1[1] * (-cI * (V2[4] * TMP5) + cI * (P2[2] * TMP6)) +
      P1[2] * (-cI * (V2[3] * TMP5) + cI * (P2[1] * TMP6))) + TMP8 * (TMP5 *
      (+cI * (V2[4] * V1[3] + V2[3] * V1[4])) - TMP6 * (+cI * (V1[4] * P2[1] +
      V1[3] * P2[2])))));
  T3[12] = denom * (TMP4 * (OM3 * (TMP6 * (P3[2] * (P3[2] * (OM3 * (TMP22 *
      1.333333333 * (-cI * (TMP8) + cI * (TMP23)) + 1.333333333 * (TMP5 * (-cI
      * (TMP23) + cI * (TMP8)))) + 0.666666667 * cI * (TMP33)) + (P1[2] * 2. *
      (-cI * (TMP22) + cI * (TMP5)) + 2. * (P2[2] * (-cI * (TMP23) + cI *
      (TMP8))))) + (TMP22 * 0.666666667 * (-cI * (TMP23) + cI * (TMP8)) +
      0.666666667 * (TMP5 * (-cI * (TMP8) + cI * (TMP23))))) + P3[2] * TMP5 *
      (V2[4] * 2. * (-cI * (TMP8) + cI * (TMP23)) - 0.666666667 * cI * (P3[2] *
      TMP3))) + (TMP5 * - 2. * (+cI * (V2[4] * P1[2]) + 0.333333333 * cI *
      (TMP3)) + 2. * (TMP6 * (+cI * (P1[2] * P2[2]) + 0.333333333 * cI *
      (TMP33))))) + TMP8 * (TMP6 * (OM3 * P3[2] * (V1[4] * 2. * (-cI * (TMP5) +
      cI * (TMP22)) - 0.666666667 * cI * (P3[2] * TMP7)) + (-0.666666667 * cI *
      (TMP7) - 2. * cI * (V1[4] * P2[2]))) + TMP5 * (TMP0 * 0.666666667 * (+cI
      * (P3[2] * P3[2] * OM3) + 1.000000000 * cI) + 2. * cI * (V2[4] *
      V1[4]))));
  T3[13] = denom * (OM3 * (TMP4 * (TMP6 * (P3[2] * (P3[3] * (OM3 * (TMP22 *
      1.333333333 * (-cI * (TMP8) + cI * (TMP23)) + 1.333333333 * (TMP5 * (-cI
      * (TMP23) + cI * (TMP8)))) + 0.666666667 * cI * (TMP33)) + (P1[3] * (-cI
      * (TMP22) + cI * (TMP5)) + P2[3] * (-cI * (TMP23) + cI * (TMP8)))) +
      P3[3] * (P1[2] * (-cI * (TMP22) + cI * (TMP5)) + P2[2] * (-cI * (TMP23) +
      cI * (TMP8)))) + TMP5 * (P3[2] * (V2[5] * (-cI * (TMP8) + cI * (TMP23)) -
      0.666666667 * cI * (P3[3] * TMP3)) + P3[3] * V2[4] * (-cI * (TMP8) + cI *
      (TMP23)))) + TMP8 * (TMP6 * (P3[2] * (V1[5] * (-cI * (TMP5) + cI *
      (TMP22)) - 0.666666667 * cI * (P3[3] * TMP7)) + P3[3] * V1[4] * (-cI *
      (TMP5) + cI * (TMP22))) + 0.666666667 * cI * (TMP0 * P3[2] * P3[3] *
      TMP5))) + (TMP4 * (P1[2] * (-cI * (V2[5] * TMP5) + cI * (P2[3] * TMP6)) +
      P1[3] * (-cI * (V2[4] * TMP5) + cI * (P2[2] * TMP6))) + TMP8 * (TMP5 *
      (+cI * (V2[4] * V1[5] + V2[5] * V1[4])) - TMP6 * (+cI * (V1[4] * P2[3] +
      V1[5] * P2[2])))));
  T3[14] = denom * (OM3 * (TMP4 * (TMP6 * (P3[0] * (P3[3] * (OM3 * (TMP22 *
      1.333333333 * (-cI * (TMP8) + cI * (TMP23)) + 1.333333333 * (TMP5 * (-cI
      * (TMP23) + cI * (TMP8)))) + 0.666666667 * cI * (TMP33)) + (P1[3] * (-cI
      * (TMP22) + cI * (TMP5)) + P2[3] * (-cI * (TMP23) + cI * (TMP8)))) +
      P3[3] * (P1[0] * (-cI * (TMP22) + cI * (TMP5)) + P2[0] * (-cI * (TMP23) +
      cI * (TMP8)))) + TMP5 * (P3[0] * (V2[5] * (-cI * (TMP8) + cI * (TMP23)) -
      0.666666667 * cI * (P3[3] * TMP3)) + P3[3] * V2[2] * (-cI * (TMP8) + cI *
      (TMP23)))) + TMP8 * (TMP6 * (P3[0] * (V1[5] * (-cI * (TMP5) + cI *
      (TMP22)) - 0.666666667 * cI * (P3[3] * TMP7)) + P3[3] * V1[2] * (-cI *
      (TMP5) + cI * (TMP22))) + 0.666666667 * cI * (TMP0 * P3[0] * P3[3] *
      TMP5))) + (TMP4 * (P1[0] * (-cI * (V2[5] * TMP5) + cI * (P2[3] * TMP6)) +
      P1[3] * (-cI * (V2[2] * TMP5) + cI * (P2[0] * TMP6))) + TMP8 * (TMP5 *
      (+cI * (V2[5] * V1[2] + V2[2] * V1[5])) - TMP6 * (+cI * (V1[5] * P2[0] +
      V1[2] * P2[3])))));
  T3[15] = denom * (OM3 * (TMP4 * (TMP6 * (P3[1] * (P3[3] * (OM3 * (TMP22 *
      1.333333333 * (-cI * (TMP8) + cI * (TMP23)) + 1.333333333 * (TMP5 * (-cI
      * (TMP23) + cI * (TMP8)))) + 0.666666667 * cI * (TMP33)) + (P1[3] * (-cI
      * (TMP22) + cI * (TMP5)) + P2[3] * (-cI * (TMP23) + cI * (TMP8)))) +
      P3[3] * (P1[1] * (-cI * (TMP22) + cI * (TMP5)) + P2[1] * (-cI * (TMP23) +
      cI * (TMP8)))) + TMP5 * (P3[1] * (V2[5] * (-cI * (TMP8) + cI * (TMP23)) -
      0.666666667 * cI * (P3[3] * TMP3)) + P3[3] * V2[3] * (-cI * (TMP8) + cI *
      (TMP23)))) + TMP8 * (TMP6 * (P3[1] * (V1[5] * (-cI * (TMP5) + cI *
      (TMP22)) - 0.666666667 * cI * (P3[3] * TMP7)) + P3[3] * V1[3] * (-cI *
      (TMP5) + cI * (TMP22))) + 0.666666667 * cI * (TMP0 * P3[1] * P3[3] *
      TMP5))) + (TMP4 * (P1[1] * (-cI * (V2[5] * TMP5) + cI * (P2[3] * TMP6)) +
      P1[3] * (-cI * (V2[3] * TMP5) + cI * (P2[1] * TMP6))) + TMP8 * (TMP5 *
      (+cI * (V2[5] * V1[3] + V2[3] * V1[5])) - TMP6 * (+cI * (V1[5] * P2[1] +
      V1[3] * P2[3])))));
  T3[16] = denom * (OM3 * (TMP4 * (TMP6 * (P3[2] * (P3[3] * (OM3 * (TMP22 *
      1.333333333 * (-cI * (TMP8) + cI * (TMP23)) + 1.333333333 * (TMP5 * (-cI
      * (TMP23) + cI * (TMP8)))) + 0.666666667 * cI * (TMP33)) + (P1[3] * (-cI
      * (TMP22) + cI * (TMP5)) + P2[3] * (-cI * (TMP23) + cI * (TMP8)))) +
      P3[3] * (P1[2] * (-cI * (TMP22) + cI * (TMP5)) + P2[2] * (-cI * (TMP23) +
      cI * (TMP8)))) + TMP5 * (P3[2] * (V2[5] * (-cI * (TMP8) + cI * (TMP23)) -
      0.666666667 * cI * (P3[3] * TMP3)) + P3[3] * V2[4] * (-cI * (TMP8) + cI *
      (TMP23)))) + TMP8 * (TMP6 * (P3[2] * (V1[5] * (-cI * (TMP5) + cI *
      (TMP22)) - 0.666666667 * cI * (P3[3] * TMP7)) + P3[3] * V1[4] * (-cI *
      (TMP5) + cI * (TMP22))) + 0.666666667 * cI * (TMP0 * P3[2] * P3[3] *
      TMP5))) + (TMP4 * (P1[2] * (-cI * (V2[5] * TMP5) + cI * (P2[3] * TMP6)) +
      P1[3] * (-cI * (V2[4] * TMP5) + cI * (P2[2] * TMP6))) + TMP8 * (TMP5 *
      (+cI * (V2[5] * V1[4] + V2[4] * V1[5])) - TMP6 * (+cI * (V1[5] * P2[2] +
      V1[4] * P2[3])))));
  T3[17] = denom * (TMP4 * (OM3 * (TMP6 * (P3[3] * (P3[3] * (OM3 * (TMP22 *
      1.333333333 * (-cI * (TMP8) + cI * (TMP23)) + 1.333333333 * (TMP5 * (-cI
      * (TMP23) + cI * (TMP8)))) + 0.666666667 * cI * (TMP33)) + (P1[3] * 2. *
      (-cI * (TMP22) + cI * (TMP5)) + 2. * (P2[3] * (-cI * (TMP23) + cI *
      (TMP8))))) + (TMP22 * 0.666666667 * (-cI * (TMP23) + cI * (TMP8)) +
      0.666666667 * (TMP5 * (-cI * (TMP8) + cI * (TMP23))))) + P3[3] * TMP5 *
      (V2[5] * 2. * (-cI * (TMP8) + cI * (TMP23)) - 0.666666667 * cI * (P3[3] *
      TMP3))) + (TMP5 * - 2. * (+cI * (V2[5] * P1[3]) + 0.333333333 * cI *
      (TMP3)) + 2. * (TMP6 * (+cI * (P1[3] * P2[3]) + 0.333333333 * cI *
      (TMP33))))) + TMP8 * (TMP6 * (OM3 * P3[3] * (V1[5] * 2. * (-cI * (TMP5) +
      cI * (TMP22)) - 0.666666667 * cI * (P3[3] * TMP7)) + (-0.666666667 * cI *
      (TMP7) - 2. * cI * (V1[5] * P2[3]))) + TMP5 * (TMP0 * 0.666666667 * (+cI
      * (P3[3] * P3[3] * OM3) + 1.000000000 * cI) + 2. * cI * (V2[5] *
      V1[5]))));
}


void FFV44_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  complex<double> TMP11; 
  complex<double> TMP10; 
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
  TMP11 = (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) +
      F1[5] * (F2[2] * (+cI * (P3[2]) - P3[1]) + F2[3] * (P3[0] + P3[3])));
  TMP10 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
      F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * - cI * (OM3 * - P3[0] * (TMP10 + 2. * (TMP11)) + (F2[4] *
      F1[2] + F2[5] * F1[3] + 2. * (F2[2] * F1[4] + F2[3] * F1[5])));
  V3[3] = denom * - cI * (OM3 * - P3[1] * (TMP10 + 2. * (TMP11)) + (+2. *
      (F2[3] * F1[4] + F2[2] * F1[5]) - F2[5] * F1[2] - F2[4] * F1[3]));
  V3[4] = denom * - cI * (OM3 * - P3[2] * (TMP10 + 2. * (TMP11)) + (-cI *
      (F2[5] * F1[2]) + cI * (F2[4] * F1[3]) - 2. * cI * (F2[2] * F1[5]) + 2. *
      cI * (F2[3] * F1[4])));
  V3[5] = denom * - cI * (OM3 * - P3[3] * (TMP10 + 2. * (TMP11)) + (F2[5] *
      F1[3] + 2. * (F2[2] * F1[4]) - 2. * (F2[3] * F1[5]) - F2[4] * F1[2]));
}


void VVT16_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> TMP0; 
  double P2[4]; 
  complex<double> TMP26; 
  complex<double> TMP25; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  TMP25 = (P1[0] * - 1. * (P2[1] * T3[6] + P2[2] * T3[10] + P2[3] * T3[14] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[11] + P2[3] *
      T3[15] - P2[0] * T3[3]) + (P1[2] * (P2[1] * T3[8] + P2[2] * T3[12] +
      P2[3] * T3[16] - P2[0] * T3[4]) + P1[3] * (P2[1] * T3[9] + P2[2] * T3[13]
      + P2[3] * T3[17] - P2[0] * T3[5]))));
  TMP26 = (P1[0] * - 1. * (P2[1] * T3[3] + P2[2] * T3[4] + P2[3] * T3[5] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[8] + P2[3] * T3[9]
      - P2[0] * T3[6]) + (P1[2] * (P2[1] * T3[11] + P2[2] * T3[12] + P2[3] *
      T3[13] - P2[0] * T3[10]) + P1[3] * (P2[1] * T3[15] + P2[2] * T3[16] +
      P2[3] * T3[17] - P2[0] * T3[14]))));
  TMP0 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  vertex = COUP * - TMP0 * (+cI * (TMP25 + TMP26)); 
}


void VVT17_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP33; 
  double P1[4]; 
  complex<double> TMP22; 
  complex<double> TMP0; 
  double P2[4]; 
  complex<double> TMP23; 
  complex<double> TMP7; 
  double P3[4]; 
  complex<double> denom; 
  double OM3; 
  complex<double> TMP9; 
  complex<double> TMP3; 
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
  TMP22 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  TMP23 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP9 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP33 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP7 = (V1[2] * P2[0] - V1[3] * P2[1] - V1[4] * P2[2] - V1[5] * P2[3]); 
  TMP0 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP3 = (V2[2] * P1[0] - V2[3] * P1[1] - V2[4] * P1[2] - V2[5] * P1[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * (OM3 * (P3[0] * (P3[0] * (OM3 * 1.333333333 * TMP22 * TMP23 *
      (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7)) + 0.666666667 * (TMP33 * (-cI
      * (TMP0 * TMP9) + cI * (TMP3 * TMP7)))) + (P1[0] * 2. * TMP22 * (-cI *
      (TMP3 * TMP7) + cI * (TMP0 * TMP9)) + 2. * (P2[0] * TMP23 * (-cI * (TMP3
      * TMP7) + cI * (TMP0 * TMP9))))) + 0.666666667 * (TMP22 * TMP23 * (-cI *
      (TMP0 * TMP9) + cI * (TMP3 * TMP7)))) + (TMP0 * 2. * TMP9 * (-cI * (P1[0]
      * P2[0]) + 0.333333333 * cI * (TMP33)) + 2. * (TMP3 * TMP7 *
      (-0.333333333 * cI * (TMP33) + cI * (P1[0] * P2[0])))));
  T3[3] = denom * (OM3 * (P3[0] * (P3[1] * (OM3 * 1.333333333 * TMP22 * TMP23 *
      (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7)) + 0.666666667 * (TMP33 * (-cI
      * (TMP0 * TMP9) + cI * (TMP3 * TMP7)))) + (P1[1] * TMP22 * (-cI * (TMP3 *
      TMP7) + cI * (TMP0 * TMP9)) + P2[1] * TMP23 * (-cI * (TMP3 * TMP7) + cI *
      (TMP0 * TMP9)))) + P3[1] * (P1[0] * TMP22 * (-cI * (TMP3 * TMP7) + cI *
      (TMP0 * TMP9)) + P2[0] * TMP23 * (-cI * (TMP3 * TMP7) + cI * (TMP0 *
      TMP9)))) + (P1[0] * P2[1] * (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7)) +
      P1[1] * P2[0] * (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7))));
  T3[4] = denom * (OM3 * (P3[0] * (P3[2] * (OM3 * 1.333333333 * TMP22 * TMP23 *
      (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7)) + 0.666666667 * (TMP33 * (-cI
      * (TMP0 * TMP9) + cI * (TMP3 * TMP7)))) + (P1[2] * TMP22 * (-cI * (TMP3 *
      TMP7) + cI * (TMP0 * TMP9)) + P2[2] * TMP23 * (-cI * (TMP3 * TMP7) + cI *
      (TMP0 * TMP9)))) + P3[2] * (P1[0] * TMP22 * (-cI * (TMP3 * TMP7) + cI *
      (TMP0 * TMP9)) + P2[0] * TMP23 * (-cI * (TMP3 * TMP7) + cI * (TMP0 *
      TMP9)))) + (P1[0] * P2[2] * (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7)) +
      P1[2] * P2[0] * (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7))));
  T3[5] = denom * (OM3 * (P3[0] * (P3[3] * (OM3 * 1.333333333 * TMP22 * TMP23 *
      (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7)) + 0.666666667 * (TMP33 * (-cI
      * (TMP0 * TMP9) + cI * (TMP3 * TMP7)))) + (P1[3] * TMP22 * (-cI * (TMP3 *
      TMP7) + cI * (TMP0 * TMP9)) + P2[3] * TMP23 * (-cI * (TMP3 * TMP7) + cI *
      (TMP0 * TMP9)))) + P3[3] * (P1[0] * TMP22 * (-cI * (TMP3 * TMP7) + cI *
      (TMP0 * TMP9)) + P2[0] * TMP23 * (-cI * (TMP3 * TMP7) + cI * (TMP0 *
      TMP9)))) + (P1[0] * P2[3] * (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7)) +
      P1[3] * P2[0] * (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7))));
  T3[6] = denom * (OM3 * (P3[0] * (P3[1] * (OM3 * 1.333333333 * TMP22 * TMP23 *
      (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7)) + 0.666666667 * (TMP33 * (-cI
      * (TMP0 * TMP9) + cI * (TMP3 * TMP7)))) + (P1[1] * TMP22 * (-cI * (TMP3 *
      TMP7) + cI * (TMP0 * TMP9)) + P2[1] * TMP23 * (-cI * (TMP3 * TMP7) + cI *
      (TMP0 * TMP9)))) + P3[1] * (P1[0] * TMP22 * (-cI * (TMP3 * TMP7) + cI *
      (TMP0 * TMP9)) + P2[0] * TMP23 * (-cI * (TMP3 * TMP7) + cI * (TMP0 *
      TMP9)))) + (P1[0] * P2[1] * (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7)) +
      P1[1] * P2[0] * (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7))));
  T3[7] = denom * (OM3 * (P3[1] * (P3[1] * (OM3 * 1.333333333 * TMP22 * TMP23 *
      (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7)) + 0.666666667 * (TMP33 * (-cI
      * (TMP0 * TMP9) + cI * (TMP3 * TMP7)))) + (P1[1] * 2. * TMP22 * (-cI *
      (TMP3 * TMP7) + cI * (TMP0 * TMP9)) + 2. * (P2[1] * TMP23 * (-cI * (TMP3
      * TMP7) + cI * (TMP0 * TMP9))))) + 0.666666667 * (TMP22 * TMP23 * (-cI *
      (TMP3 * TMP7) + cI * (TMP0 * TMP9)))) + (TMP0 * - 2. * TMP9 * (+cI *
      (P1[1] * P2[1]) + 0.333333333 * cI * (TMP33)) + 2. * (TMP3 * TMP7 * (+cI
      * (P1[1] * P2[1]) + 0.333333333 * cI * (TMP33)))));
  T3[8] = denom * (OM3 * (P3[1] * (P3[2] * (OM3 * 1.333333333 * TMP22 * TMP23 *
      (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7)) + 0.666666667 * (TMP33 * (-cI
      * (TMP0 * TMP9) + cI * (TMP3 * TMP7)))) + (P1[2] * TMP22 * (-cI * (TMP3 *
      TMP7) + cI * (TMP0 * TMP9)) + P2[2] * TMP23 * (-cI * (TMP3 * TMP7) + cI *
      (TMP0 * TMP9)))) + P3[2] * (P1[1] * TMP22 * (-cI * (TMP3 * TMP7) + cI *
      (TMP0 * TMP9)) + P2[1] * TMP23 * (-cI * (TMP3 * TMP7) + cI * (TMP0 *
      TMP9)))) + (P1[1] * P2[2] * (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7)) +
      P1[2] * P2[1] * (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7))));
  T3[9] = denom * (OM3 * (P3[1] * (P3[3] * (OM3 * 1.333333333 * TMP22 * TMP23 *
      (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7)) + 0.666666667 * (TMP33 * (-cI
      * (TMP0 * TMP9) + cI * (TMP3 * TMP7)))) + (P1[3] * TMP22 * (-cI * (TMP3 *
      TMP7) + cI * (TMP0 * TMP9)) + P2[3] * TMP23 * (-cI * (TMP3 * TMP7) + cI *
      (TMP0 * TMP9)))) + P3[3] * (P1[1] * TMP22 * (-cI * (TMP3 * TMP7) + cI *
      (TMP0 * TMP9)) + P2[1] * TMP23 * (-cI * (TMP3 * TMP7) + cI * (TMP0 *
      TMP9)))) + (P1[1] * P2[3] * (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7)) +
      P1[3] * P2[1] * (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7))));
  T3[10] = denom * (OM3 * (P3[0] * (P3[2] * (OM3 * 1.333333333 * TMP22 * TMP23
      * (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7)) + 0.666666667 * (TMP33 *
      (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7)))) + (P1[2] * TMP22 * (-cI *
      (TMP3 * TMP7) + cI * (TMP0 * TMP9)) + P2[2] * TMP23 * (-cI * (TMP3 *
      TMP7) + cI * (TMP0 * TMP9)))) + P3[2] * (P1[0] * TMP22 * (-cI * (TMP3 *
      TMP7) + cI * (TMP0 * TMP9)) + P2[0] * TMP23 * (-cI * (TMP3 * TMP7) + cI *
      (TMP0 * TMP9)))) + (P1[0] * P2[2] * (-cI * (TMP0 * TMP9) + cI * (TMP3 *
      TMP7)) + P1[2] * P2[0] * (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7))));
  T3[11] = denom * (OM3 * (P3[1] * (P3[2] * (OM3 * 1.333333333 * TMP22 * TMP23
      * (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7)) + 0.666666667 * (TMP33 *
      (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7)))) + (P1[2] * TMP22 * (-cI *
      (TMP3 * TMP7) + cI * (TMP0 * TMP9)) + P2[2] * TMP23 * (-cI * (TMP3 *
      TMP7) + cI * (TMP0 * TMP9)))) + P3[2] * (P1[1] * TMP22 * (-cI * (TMP3 *
      TMP7) + cI * (TMP0 * TMP9)) + P2[1] * TMP23 * (-cI * (TMP3 * TMP7) + cI *
      (TMP0 * TMP9)))) + (P1[1] * P2[2] * (-cI * (TMP0 * TMP9) + cI * (TMP3 *
      TMP7)) + P1[2] * P2[1] * (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7))));
  T3[12] = denom * (OM3 * (P3[2] * (P3[2] * (OM3 * 1.333333333 * TMP22 * TMP23
      * (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7)) + 0.666666667 * (TMP33 *
      (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7)))) + (P1[2] * 2. * TMP22 * (-cI
      * (TMP3 * TMP7) + cI * (TMP0 * TMP9)) + 2. * (P2[2] * TMP23 * (-cI *
      (TMP3 * TMP7) + cI * (TMP0 * TMP9))))) + 0.666666667 * (TMP22 * TMP23 *
      (-cI * (TMP3 * TMP7) + cI * (TMP0 * TMP9)))) + (TMP0 * - 2. * TMP9 * (+cI
      * (P1[2] * P2[2]) + 0.333333333 * cI * (TMP33)) + 2. * (TMP3 * TMP7 *
      (+cI * (P1[2] * P2[2]) + 0.333333333 * cI * (TMP33)))));
  T3[13] = denom * (OM3 * (P3[2] * (P3[3] * (OM3 * 1.333333333 * TMP22 * TMP23
      * (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7)) + 0.666666667 * (TMP33 *
      (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7)))) + (P1[3] * TMP22 * (-cI *
      (TMP3 * TMP7) + cI * (TMP0 * TMP9)) + P2[3] * TMP23 * (-cI * (TMP3 *
      TMP7) + cI * (TMP0 * TMP9)))) + P3[3] * (P1[2] * TMP22 * (-cI * (TMP3 *
      TMP7) + cI * (TMP0 * TMP9)) + P2[2] * TMP23 * (-cI * (TMP3 * TMP7) + cI *
      (TMP0 * TMP9)))) + (P1[2] * P2[3] * (-cI * (TMP0 * TMP9) + cI * (TMP3 *
      TMP7)) + P1[3] * P2[2] * (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7))));
  T3[14] = denom * (OM3 * (P3[0] * (P3[3] * (OM3 * 1.333333333 * TMP22 * TMP23
      * (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7)) + 0.666666667 * (TMP33 *
      (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7)))) + (P1[3] * TMP22 * (-cI *
      (TMP3 * TMP7) + cI * (TMP0 * TMP9)) + P2[3] * TMP23 * (-cI * (TMP3 *
      TMP7) + cI * (TMP0 * TMP9)))) + P3[3] * (P1[0] * TMP22 * (-cI * (TMP3 *
      TMP7) + cI * (TMP0 * TMP9)) + P2[0] * TMP23 * (-cI * (TMP3 * TMP7) + cI *
      (TMP0 * TMP9)))) + (P1[0] * P2[3] * (-cI * (TMP0 * TMP9) + cI * (TMP3 *
      TMP7)) + P1[3] * P2[0] * (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7))));
  T3[15] = denom * (OM3 * (P3[1] * (P3[3] * (OM3 * 1.333333333 * TMP22 * TMP23
      * (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7)) + 0.666666667 * (TMP33 *
      (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7)))) + (P1[3] * TMP22 * (-cI *
      (TMP3 * TMP7) + cI * (TMP0 * TMP9)) + P2[3] * TMP23 * (-cI * (TMP3 *
      TMP7) + cI * (TMP0 * TMP9)))) + P3[3] * (P1[1] * TMP22 * (-cI * (TMP3 *
      TMP7) + cI * (TMP0 * TMP9)) + P2[1] * TMP23 * (-cI * (TMP3 * TMP7) + cI *
      (TMP0 * TMP9)))) + (P1[1] * P2[3] * (-cI * (TMP0 * TMP9) + cI * (TMP3 *
      TMP7)) + P1[3] * P2[1] * (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7))));
  T3[16] = denom * (OM3 * (P3[2] * (P3[3] * (OM3 * 1.333333333 * TMP22 * TMP23
      * (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7)) + 0.666666667 * (TMP33 *
      (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7)))) + (P1[3] * TMP22 * (-cI *
      (TMP3 * TMP7) + cI * (TMP0 * TMP9)) + P2[3] * TMP23 * (-cI * (TMP3 *
      TMP7) + cI * (TMP0 * TMP9)))) + P3[3] * (P1[2] * TMP22 * (-cI * (TMP3 *
      TMP7) + cI * (TMP0 * TMP9)) + P2[2] * TMP23 * (-cI * (TMP3 * TMP7) + cI *
      (TMP0 * TMP9)))) + (P1[2] * P2[3] * (-cI * (TMP0 * TMP9) + cI * (TMP3 *
      TMP7)) + P1[3] * P2[2] * (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7))));
  T3[17] = denom * (OM3 * (P3[3] * (P3[3] * (OM3 * 1.333333333 * TMP22 * TMP23
      * (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7)) + 0.666666667 * (TMP33 *
      (-cI * (TMP0 * TMP9) + cI * (TMP3 * TMP7)))) + (P1[3] * 2. * TMP22 * (-cI
      * (TMP3 * TMP7) + cI * (TMP0 * TMP9)) + 2. * (P2[3] * TMP23 * (-cI *
      (TMP3 * TMP7) + cI * (TMP0 * TMP9))))) + 0.666666667 * (TMP22 * TMP23 *
      (-cI * (TMP3 * TMP7) + cI * (TMP0 * TMP9)))) + (TMP0 * - 2. * TMP9 * (+cI
      * (P1[3] * P2[3]) + 0.333333333 * cI * (TMP33)) + 2. * (TMP3 * TMP7 *
      (+cI * (P1[3] * P2[3]) + 0.333333333 * cI * (TMP33)))));
}


void VVT13_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> TMP22; 
  double P2[4]; 
  complex<double> TMP23; 
  double P3[4]; 
  complex<double> TMP6; 
  complex<double> denom; 
  double OM3; 
  complex<double> TMP4; 
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
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +V1[0] + V2[0]; 
  T3[1] = +V1[1] + V2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP51 = -1. * (P1[0] * (P2[1] * (V2[4] * P3[3] - V2[5] * P3[2]) + (P2[2] *
      (V2[5] * P3[1] - V2[3] * P3[3]) + P2[3] * (V2[3] * P3[2] - V2[4] *
      P3[1]))) + (P1[1] * (P2[0] * (V2[5] * P3[2] - V2[4] * P3[3]) + (P2[2] *
      (V2[2] * P3[3] - V2[5] * P3[0]) + P2[3] * (V2[4] * P3[0] - V2[2] *
      P3[2]))) + (P1[2] * (P2[0] * (V2[3] * P3[3] - V2[5] * P3[1]) + (P2[1] *
      (V2[5] * P3[0] - V2[2] * P3[3]) + P2[3] * (V2[2] * P3[1] - V2[3] *
      P3[0]))) + P1[3] * (P2[0] * (V2[4] * P3[1] - V2[3] * P3[2]) + (P2[1] *
      (V2[2] * P3[2] - V2[4] * P3[0]) + P2[2] * (V2[3] * P3[0] - V2[2] *
      P3[1]))))));
  TMP52 = -1. * (P1[0] * (P2[1] * (V1[5] * P3[2] - V1[4] * P3[3]) + (P2[2] *
      (V1[3] * P3[3] - V1[5] * P3[1]) + P2[3] * (V1[4] * P3[1] - V1[3] *
      P3[2]))) + (P1[1] * (P2[0] * (V1[4] * P3[3] - V1[5] * P3[2]) + (P2[2] *
      (V1[5] * P3[0] - V1[2] * P3[3]) + P2[3] * (V1[2] * P3[2] - V1[4] *
      P3[0]))) + (P1[2] * (P2[0] * (V1[5] * P3[1] - V1[3] * P3[3]) + (P2[1] *
      (V1[2] * P3[3] - V1[5] * P3[0]) + P2[3] * (V1[3] * P3[0] - V1[2] *
      P3[1]))) + P1[3] * (P2[0] * (V1[3] * P3[2] - V1[4] * P3[1]) + (P2[1] *
      (V1[4] * P3[0] - V1[2] * P3[2]) + P2[2] * (V1[2] * P3[1] - V1[3] *
      P3[0]))))));
  TMP22 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  TMP23 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP4 = (V1[2] * P3[0] - V1[3] * P3[1] - V1[4] * P3[2] - V1[5] * P3[3]); 
  TMP6 = (V2[2] * P3[0] - V2[3] * P3[1] - V2[4] * P3[2] - V2[5] * P3[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * - 0.500000000 * cI * (TMP4 * (OM3 * P3[0] * (TMP23 * (P2[1] *
      4. * (V2[5] * P3[2] - V2[4] * P3[3]) + (P2[2] * 4. * (V2[3] * P3[3] -
      V2[5] * P3[1]) + 4. * (P2[3] * (V2[4] * P3[1] - V2[3] * P3[2])))) -
      1.333333333 * (P3[0] * TMP51)) + (P1[0] * (P2[1] * 4. * (V2[4] * P3[3] -
      V2[5] * P3[2]) + (P2[2] * 4. * (V2[5] * P3[1] - V2[3] * P3[3]) + 4. *
      (P2[3] * (V2[3] * P3[2] - V2[4] * P3[1])))) + 1.333333333 * (TMP51))) +
      TMP6 * (OM3 * P3[0] * (TMP22 * (P1[1] * 4. * (V1[5] * P3[2] - V1[4] *
      P3[3]) + (P1[2] * 4. * (V1[3] * P3[3] - V1[5] * P3[1]) + 4. * (P1[3] *
      (V1[4] * P3[1] - V1[3] * P3[2])))) - 1.333333333 * (P3[0] * TMP52)) +
      (P2[0] * (P1[1] * 4. * (V1[4] * P3[3] - V1[5] * P3[2]) + (P1[2] * 4. *
      (V1[5] * P3[1] - V1[3] * P3[3]) + 4. * (P1[3] * (V1[3] * P3[2] - V1[4] *
      P3[1])))) + 1.333333333 * (TMP52))));
  T3[6] = denom * - 0.500000000 * cI * (OM3 * (P3[0] * (TMP4 * (TMP23 * (P2[0]
      * 2. * (V2[5] * P3[2] - V2[4] * P3[3]) + (P2[2] * 2. * (V2[2] * P3[3] -
      V2[5] * P3[0]) + 2. * (P2[3] * (V2[4] * P3[0] - V2[2] * P3[2])))) -
      1.333333333 * (P3[1] * TMP51)) + TMP6 * (TMP22 * (P1[0] * 2. * (V1[5] *
      P3[2] - V1[4] * P3[3]) + (P1[2] * 2. * (V1[2] * P3[3] - V1[5] * P3[0]) +
      2. * (P1[3] * (V1[4] * P3[0] - V1[2] * P3[2])))) - 1.333333333 * (P3[1] *
      TMP52))) + P3[1] * (TMP22 * TMP6 * (P1[1] * 2. * (V1[5] * P3[2] - V1[4] *
      P3[3]) + (P1[2] * 2. * (V1[3] * P3[3] - V1[5] * P3[1]) + 2. * (P1[3] *
      (V1[4] * P3[1] - V1[3] * P3[2])))) + TMP23 * TMP4 * (P2[1] * 2. * (V2[5]
      * P3[2] - V2[4] * P3[3]) + (P2[2] * 2. * (V2[3] * P3[3] - V2[5] * P3[1])
      + 2. * (P2[3] * (V2[4] * P3[1] - V2[3] * P3[2])))))) + (TMP4 * (P1[0] *
      (P2[0] * 2. * (V2[4] * P3[3] - V2[5] * P3[2]) + (P2[2] * 2. * (V2[5] *
      P3[0] - V2[2] * P3[3]) + 2. * (P2[3] * (V2[2] * P3[2] - V2[4] * P3[0]))))
      + P1[1] * (P2[1] * 2. * (V2[4] * P3[3] - V2[5] * P3[2]) + (P2[2] * 2. *
      (V2[5] * P3[1] - V2[3] * P3[3]) + 2. * (P2[3] * (V2[3] * P3[2] - V2[4] *
      P3[1]))))) + TMP6 * (P2[0] * (P1[0] * 2. * (V1[4] * P3[3] - V1[5] *
      P3[2]) + (P1[2] * 2. * (V1[5] * P3[0] - V1[2] * P3[3]) + 2. * (P1[3] *
      (V1[2] * P3[2] - V1[4] * P3[0])))) + P2[1] * (P1[1] * 2. * (V1[4] * P3[3]
      - V1[5] * P3[2]) + (P1[2] * 2. * (V1[5] * P3[1] - V1[3] * P3[3]) + 2. *
      (P1[3] * (V1[3] * P3[2] - V1[4] * P3[1])))))));
  T3[10] = denom * - 0.500000000 * cI * (OM3 * (P3[0] * (TMP4 * (TMP23 * (P2[0]
      * 2. * (V2[3] * P3[3] - V2[5] * P3[1]) + (P2[1] * 2. * (V2[5] * P3[0] -
      V2[2] * P3[3]) + 2. * (P2[3] * (V2[2] * P3[1] - V2[3] * P3[0])))) -
      1.333333333 * (P3[2] * TMP51)) + TMP6 * (TMP22 * (P1[0] * 2. * (V1[3] *
      P3[3] - V1[5] * P3[1]) + (P1[1] * 2. * (V1[5] * P3[0] - V1[2] * P3[3]) +
      2. * (P1[3] * (V1[2] * P3[1] - V1[3] * P3[0])))) - 1.333333333 * (P3[2] *
      TMP52))) + P3[2] * (TMP22 * TMP6 * (P1[1] * 2. * (V1[5] * P3[2] - V1[4] *
      P3[3]) + (P1[2] * 2. * (V1[3] * P3[3] - V1[5] * P3[1]) + 2. * (P1[3] *
      (V1[4] * P3[1] - V1[3] * P3[2])))) + TMP23 * TMP4 * (P2[1] * 2. * (V2[5]
      * P3[2] - V2[4] * P3[3]) + (P2[2] * 2. * (V2[3] * P3[3] - V2[5] * P3[1])
      + 2. * (P2[3] * (V2[4] * P3[1] - V2[3] * P3[2])))))) + (TMP4 * (P1[0] *
      (P2[0] * 2. * (V2[5] * P3[1] - V2[3] * P3[3]) + (P2[1] * 2. * (V2[2] *
      P3[3] - V2[5] * P3[0]) + 2. * (P2[3] * (V2[3] * P3[0] - V2[2] * P3[1]))))
      + P1[2] * (P2[1] * 2. * (V2[4] * P3[3] - V2[5] * P3[2]) + (P2[2] * 2. *
      (V2[5] * P3[1] - V2[3] * P3[3]) + 2. * (P2[3] * (V2[3] * P3[2] - V2[4] *
      P3[1]))))) + TMP6 * (P2[0] * (P1[0] * 2. * (V1[5] * P3[1] - V1[3] *
      P3[3]) + (P1[1] * 2. * (V1[2] * P3[3] - V1[5] * P3[0]) + 2. * (P1[3] *
      (V1[3] * P3[0] - V1[2] * P3[1])))) + P2[2] * (P1[1] * 2. * (V1[4] * P3[3]
      - V1[5] * P3[2]) + (P1[2] * 2. * (V1[5] * P3[1] - V1[3] * P3[3]) + 2. *
      (P1[3] * (V1[3] * P3[2] - V1[4] * P3[1])))))));
  T3[14] = denom * - 0.500000000 * cI * (OM3 * (P3[0] * (TMP4 * (TMP23 * (P2[0]
      * 2. * (V2[4] * P3[1] - V2[3] * P3[2]) + (P2[1] * 2. * (V2[2] * P3[2] -
      V2[4] * P3[0]) + 2. * (P2[2] * (V2[3] * P3[0] - V2[2] * P3[1])))) -
      1.333333333 * (P3[3] * TMP51)) + TMP6 * (TMP22 * (P1[0] * 2. * (V1[4] *
      P3[1] - V1[3] * P3[2]) + (P1[1] * 2. * (V1[2] * P3[2] - V1[4] * P3[0]) +
      2. * (P1[2] * (V1[3] * P3[0] - V1[2] * P3[1])))) - 1.333333333 * (P3[3] *
      TMP52))) + P3[3] * (TMP22 * TMP6 * (P1[1] * 2. * (V1[5] * P3[2] - V1[4] *
      P3[3]) + (P1[2] * 2. * (V1[3] * P3[3] - V1[5] * P3[1]) + 2. * (P1[3] *
      (V1[4] * P3[1] - V1[3] * P3[2])))) + TMP23 * TMP4 * (P2[1] * 2. * (V2[5]
      * P3[2] - V2[4] * P3[3]) + (P2[2] * 2. * (V2[3] * P3[3] - V2[5] * P3[1])
      + 2. * (P2[3] * (V2[4] * P3[1] - V2[3] * P3[2])))))) + (TMP4 * (P1[0] *
      (P2[0] * 2. * (V2[3] * P3[2] - V2[4] * P3[1]) + (P2[1] * 2. * (V2[4] *
      P3[0] - V2[2] * P3[2]) + 2. * (P2[2] * (V2[2] * P3[1] - V2[3] * P3[0]))))
      + P1[3] * (P2[1] * 2. * (V2[4] * P3[3] - V2[5] * P3[2]) + (P2[2] * 2. *
      (V2[5] * P3[1] - V2[3] * P3[3]) + 2. * (P2[3] * (V2[3] * P3[2] - V2[4] *
      P3[1]))))) + TMP6 * (P2[0] * (P1[0] * 2. * (V1[3] * P3[2] - V1[4] *
      P3[1]) + (P1[1] * 2. * (V1[4] * P3[0] - V1[2] * P3[2]) + 2. * (P1[2] *
      (V1[2] * P3[1] - V1[3] * P3[0])))) + P2[3] * (P1[1] * 2. * (V1[4] * P3[3]
      - V1[5] * P3[2]) + (P1[2] * 2. * (V1[5] * P3[1] - V1[3] * P3[3]) + 2. *
      (P1[3] * (V1[3] * P3[2] - V1[4] * P3[1])))))));
  T3[3] = denom * 0.500000000 * cI * (OM3 * (P3[0] * (TMP4 * (TMP23 * (P2[0] *
      2. * (V2[4] * P3[3] - V2[5] * P3[2]) + (P2[2] * 2. * (V2[5] * P3[0] -
      V2[2] * P3[3]) + 2. * (P2[3] * (V2[2] * P3[2] - V2[4] * P3[0])))) +
      1.333333333 * (P3[1] * TMP51)) + TMP6 * (TMP22 * (P1[0] * 2. * (V1[4] *
      P3[3] - V1[5] * P3[2]) + (P1[2] * 2. * (V1[5] * P3[0] - V1[2] * P3[3]) +
      2. * (P1[3] * (V1[2] * P3[2] - V1[4] * P3[0])))) + 1.333333333 * (P3[1] *
      TMP52))) + P3[1] * (TMP22 * TMP6 * (P1[1] * 2. * (V1[4] * P3[3] - V1[5] *
      P3[2]) + (P1[2] * 2. * (V1[5] * P3[1] - V1[3] * P3[3]) + 2. * (P1[3] *
      (V1[3] * P3[2] - V1[4] * P3[1])))) + TMP23 * TMP4 * (P2[1] * 2. * (V2[4]
      * P3[3] - V2[5] * P3[2]) + (P2[2] * 2. * (V2[5] * P3[1] - V2[3] * P3[3])
      + 2. * (P2[3] * (V2[3] * P3[2] - V2[4] * P3[1])))))) + (TMP4 * (P1[0] *
      (P2[0] * 2. * (V2[5] * P3[2] - V2[4] * P3[3]) + (P2[2] * 2. * (V2[2] *
      P3[3] - V2[5] * P3[0]) + 2. * (P2[3] * (V2[4] * P3[0] - V2[2] * P3[2]))))
      + P1[1] * (P2[1] * 2. * (V2[5] * P3[2] - V2[4] * P3[3]) + (P2[2] * 2. *
      (V2[3] * P3[3] - V2[5] * P3[1]) + 2. * (P2[3] * (V2[4] * P3[1] - V2[3] *
      P3[2]))))) + TMP6 * (P2[0] * (P1[0] * 2. * (V1[5] * P3[2] - V1[4] *
      P3[3]) + (P1[2] * 2. * (V1[2] * P3[3] - V1[5] * P3[0]) + 2. * (P1[3] *
      (V1[4] * P3[0] - V1[2] * P3[2])))) + P2[1] * (P1[1] * 2. * (V1[5] * P3[2]
      - V1[4] * P3[3]) + (P1[2] * 2. * (V1[3] * P3[3] - V1[5] * P3[1]) + 2. *
      (P1[3] * (V1[4] * P3[1] - V1[3] * P3[2])))))));
  T3[7] = denom * 0.500000000 * cI * (TMP4 * (OM3 * P3[1] * (TMP23 * (P2[0] *
      4. * (V2[4] * P3[3] - V2[5] * P3[2]) + (P2[2] * 4. * (V2[5] * P3[0] -
      V2[2] * P3[3]) + 4. * (P2[3] * (V2[2] * P3[2] - V2[4] * P3[0])))) +
      1.333333333 * (P3[1] * TMP51)) + (P1[1] * (P2[0] * 4. * (V2[5] * P3[2] -
      V2[4] * P3[3]) + (P2[2] * 4. * (V2[2] * P3[3] - V2[5] * P3[0]) + 4. *
      (P2[3] * (V2[4] * P3[0] - V2[2] * P3[2])))) + 1.333333333 * (TMP51))) +
      TMP6 * (OM3 * P3[1] * (TMP22 * (P1[0] * 4. * (V1[4] * P3[3] - V1[5] *
      P3[2]) + (P1[2] * 4. * (V1[5] * P3[0] - V1[2] * P3[3]) + 4. * (P1[3] *
      (V1[2] * P3[2] - V1[4] * P3[0])))) + 1.333333333 * (P3[1] * TMP52)) +
      (P2[1] * (P1[0] * 4. * (V1[5] * P3[2] - V1[4] * P3[3]) + (P1[2] * 4. *
      (V1[2] * P3[3] - V1[5] * P3[0]) + 4. * (P1[3] * (V1[4] * P3[0] - V1[2] *
      P3[2])))) + 1.333333333 * (TMP52))));
  T3[11] = denom * 0.500000000 * cI * (OM3 * (P3[1] * (TMP4 * (TMP23 * (P2[0] *
      2. * (V2[5] * P3[1] - V2[3] * P3[3]) + (P2[1] * 2. * (V2[2] * P3[3] -
      V2[5] * P3[0]) + 2. * (P2[3] * (V2[3] * P3[0] - V2[2] * P3[1])))) +
      1.333333333 * (P3[2] * TMP51)) + TMP6 * (TMP22 * (P1[0] * 2. * (V1[5] *
      P3[1] - V1[3] * P3[3]) + (P1[1] * 2. * (V1[2] * P3[3] - V1[5] * P3[0]) +
      2. * (P1[3] * (V1[3] * P3[0] - V1[2] * P3[1])))) + 1.333333333 * (P3[2] *
      TMP52))) + P3[2] * (TMP22 * TMP6 * (P1[0] * 2. * (V1[4] * P3[3] - V1[5] *
      P3[2]) + (P1[2] * 2. * (V1[5] * P3[0] - V1[2] * P3[3]) + 2. * (P1[3] *
      (V1[2] * P3[2] - V1[4] * P3[0])))) + TMP23 * TMP4 * (P2[0] * 2. * (V2[4]
      * P3[3] - V2[5] * P3[2]) + (P2[2] * 2. * (V2[5] * P3[0] - V2[2] * P3[3])
      + 2. * (P2[3] * (V2[2] * P3[2] - V2[4] * P3[0])))))) + (TMP4 * (P1[1] *
      (P2[0] * 2. * (V2[3] * P3[3] - V2[5] * P3[1]) + (P2[1] * 2. * (V2[5] *
      P3[0] - V2[2] * P3[3]) + 2. * (P2[3] * (V2[2] * P3[1] - V2[3] * P3[0]))))
      + P1[2] * (P2[0] * 2. * (V2[5] * P3[2] - V2[4] * P3[3]) + (P2[2] * 2. *
      (V2[2] * P3[3] - V2[5] * P3[0]) + 2. * (P2[3] * (V2[4] * P3[0] - V2[2] *
      P3[2]))))) + TMP6 * (P2[1] * (P1[0] * 2. * (V1[3] * P3[3] - V1[5] *
      P3[1]) + (P1[1] * 2. * (V1[5] * P3[0] - V1[2] * P3[3]) + 2. * (P1[3] *
      (V1[2] * P3[1] - V1[3] * P3[0])))) + P2[2] * (P1[0] * 2. * (V1[5] * P3[2]
      - V1[4] * P3[3]) + (P1[2] * 2. * (V1[2] * P3[3] - V1[5] * P3[0]) + 2. *
      (P1[3] * (V1[4] * P3[0] - V1[2] * P3[2])))))));
  T3[15] = denom * 0.500000000 * cI * (OM3 * (P3[1] * (TMP4 * (TMP23 * (P2[0] *
      2. * (V2[3] * P3[2] - V2[4] * P3[1]) + (P2[1] * 2. * (V2[4] * P3[0] -
      V2[2] * P3[2]) + 2. * (P2[2] * (V2[2] * P3[1] - V2[3] * P3[0])))) +
      1.333333333 * (P3[3] * TMP51)) + TMP6 * (TMP22 * (P1[0] * 2. * (V1[3] *
      P3[2] - V1[4] * P3[1]) + (P1[1] * 2. * (V1[4] * P3[0] - V1[2] * P3[2]) +
      2. * (P1[2] * (V1[2] * P3[1] - V1[3] * P3[0])))) + 1.333333333 * (P3[3] *
      TMP52))) + P3[3] * (TMP22 * TMP6 * (P1[0] * 2. * (V1[4] * P3[3] - V1[5] *
      P3[2]) + (P1[2] * 2. * (V1[5] * P3[0] - V1[2] * P3[3]) + 2. * (P1[3] *
      (V1[2] * P3[2] - V1[4] * P3[0])))) + TMP23 * TMP4 * (P2[0] * 2. * (V2[4]
      * P3[3] - V2[5] * P3[2]) + (P2[2] * 2. * (V2[5] * P3[0] - V2[2] * P3[3])
      + 2. * (P2[3] * (V2[2] * P3[2] - V2[4] * P3[0])))))) + (TMP4 * (P1[1] *
      (P2[0] * 2. * (V2[4] * P3[1] - V2[3] * P3[2]) + (P2[1] * 2. * (V2[2] *
      P3[2] - V2[4] * P3[0]) + 2. * (P2[2] * (V2[3] * P3[0] - V2[2] * P3[1]))))
      + P1[3] * (P2[0] * 2. * (V2[5] * P3[2] - V2[4] * P3[3]) + (P2[2] * 2. *
      (V2[2] * P3[3] - V2[5] * P3[0]) + 2. * (P2[3] * (V2[4] * P3[0] - V2[2] *
      P3[2]))))) + TMP6 * (P2[1] * (P1[0] * 2. * (V1[4] * P3[1] - V1[3] *
      P3[2]) + (P1[1] * 2. * (V1[2] * P3[2] - V1[4] * P3[0]) + 2. * (P1[2] *
      (V1[3] * P3[0] - V1[2] * P3[1])))) + P2[3] * (P1[0] * 2. * (V1[5] * P3[2]
      - V1[4] * P3[3]) + (P1[2] * 2. * (V1[2] * P3[3] - V1[5] * P3[0]) + 2. *
      (P1[3] * (V1[4] * P3[0] - V1[2] * P3[2])))))));
  T3[4] = denom * 0.500000000 * cI * (OM3 * (P3[0] * (TMP4 * (TMP23 * (P2[0] *
      2. * (V2[5] * P3[1] - V2[3] * P3[3]) + (P2[1] * 2. * (V2[2] * P3[3] -
      V2[5] * P3[0]) + 2. * (P2[3] * (V2[3] * P3[0] - V2[2] * P3[1])))) +
      1.333333333 * (P3[2] * TMP51)) + TMP6 * (TMP22 * (P1[0] * 2. * (V1[5] *
      P3[1] - V1[3] * P3[3]) + (P1[1] * 2. * (V1[2] * P3[3] - V1[5] * P3[0]) +
      2. * (P1[3] * (V1[3] * P3[0] - V1[2] * P3[1])))) + 1.333333333 * (P3[2] *
      TMP52))) + P3[2] * (TMP22 * TMP6 * (P1[1] * 2. * (V1[4] * P3[3] - V1[5] *
      P3[2]) + (P1[2] * 2. * (V1[5] * P3[1] - V1[3] * P3[3]) + 2. * (P1[3] *
      (V1[3] * P3[2] - V1[4] * P3[1])))) + TMP23 * TMP4 * (P2[1] * 2. * (V2[4]
      * P3[3] - V2[5] * P3[2]) + (P2[2] * 2. * (V2[5] * P3[1] - V2[3] * P3[3])
      + 2. * (P2[3] * (V2[3] * P3[2] - V2[4] * P3[1])))))) + (TMP4 * (P1[0] *
      (P2[0] * 2. * (V2[3] * P3[3] - V2[5] * P3[1]) + (P2[1] * 2. * (V2[5] *
      P3[0] - V2[2] * P3[3]) + 2. * (P2[3] * (V2[2] * P3[1] - V2[3] * P3[0]))))
      + P1[2] * (P2[1] * 2. * (V2[5] * P3[2] - V2[4] * P3[3]) + (P2[2] * 2. *
      (V2[3] * P3[3] - V2[5] * P3[1]) + 2. * (P2[3] * (V2[4] * P3[1] - V2[3] *
      P3[2]))))) + TMP6 * (P2[0] * (P1[0] * 2. * (V1[3] * P3[3] - V1[5] *
      P3[1]) + (P1[1] * 2. * (V1[5] * P3[0] - V1[2] * P3[3]) + 2. * (P1[3] *
      (V1[2] * P3[1] - V1[3] * P3[0])))) + P2[2] * (P1[1] * 2. * (V1[5] * P3[2]
      - V1[4] * P3[3]) + (P1[2] * 2. * (V1[3] * P3[3] - V1[5] * P3[1]) + 2. *
      (P1[3] * (V1[4] * P3[1] - V1[3] * P3[2])))))));
  T3[8] = denom * 0.500000000 * cI * (OM3 * (P3[1] * (TMP4 * (TMP23 * (P2[0] *
      2. * (V2[5] * P3[1] - V2[3] * P3[3]) + (P2[1] * 2. * (V2[2] * P3[3] -
      V2[5] * P3[0]) + 2. * (P2[3] * (V2[3] * P3[0] - V2[2] * P3[1])))) +
      1.333333333 * (P3[2] * TMP51)) + TMP6 * (TMP22 * (P1[0] * 2. * (V1[5] *
      P3[1] - V1[3] * P3[3]) + (P1[1] * 2. * (V1[2] * P3[3] - V1[5] * P3[0]) +
      2. * (P1[3] * (V1[3] * P3[0] - V1[2] * P3[1])))) + 1.333333333 * (P3[2] *
      TMP52))) + P3[2] * (TMP22 * TMP6 * (P1[0] * 2. * (V1[4] * P3[3] - V1[5] *
      P3[2]) + (P1[2] * 2. * (V1[5] * P3[0] - V1[2] * P3[3]) + 2. * (P1[3] *
      (V1[2] * P3[2] - V1[4] * P3[0])))) + TMP23 * TMP4 * (P2[0] * 2. * (V2[4]
      * P3[3] - V2[5] * P3[2]) + (P2[2] * 2. * (V2[5] * P3[0] - V2[2] * P3[3])
      + 2. * (P2[3] * (V2[2] * P3[2] - V2[4] * P3[0])))))) + (TMP4 * (P1[1] *
      (P2[0] * 2. * (V2[3] * P3[3] - V2[5] * P3[1]) + (P2[1] * 2. * (V2[5] *
      P3[0] - V2[2] * P3[3]) + 2. * (P2[3] * (V2[2] * P3[1] - V2[3] * P3[0]))))
      + P1[2] * (P2[0] * 2. * (V2[5] * P3[2] - V2[4] * P3[3]) + (P2[2] * 2. *
      (V2[2] * P3[3] - V2[5] * P3[0]) + 2. * (P2[3] * (V2[4] * P3[0] - V2[2] *
      P3[2]))))) + TMP6 * (P2[1] * (P1[0] * 2. * (V1[3] * P3[3] - V1[5] *
      P3[1]) + (P1[1] * 2. * (V1[5] * P3[0] - V1[2] * P3[3]) + 2. * (P1[3] *
      (V1[2] * P3[1] - V1[3] * P3[0])))) + P2[2] * (P1[0] * 2. * (V1[5] * P3[2]
      - V1[4] * P3[3]) + (P1[2] * 2. * (V1[2] * P3[3] - V1[5] * P3[0]) + 2. *
      (P1[3] * (V1[4] * P3[0] - V1[2] * P3[2])))))));
  T3[12] = denom * 0.500000000 * cI * (TMP4 * (OM3 * P3[2] * (TMP23 * (P2[0] *
      4. * (V2[5] * P3[1] - V2[3] * P3[3]) + (P2[1] * 4. * (V2[2] * P3[3] -
      V2[5] * P3[0]) + 4. * (P2[3] * (V2[3] * P3[0] - V2[2] * P3[1])))) +
      1.333333333 * (P3[2] * TMP51)) + (P1[2] * (P2[0] * 4. * (V2[3] * P3[3] -
      V2[5] * P3[1]) + (P2[1] * 4. * (V2[5] * P3[0] - V2[2] * P3[3]) + 4. *
      (P2[3] * (V2[2] * P3[1] - V2[3] * P3[0])))) + 1.333333333 * (TMP51))) +
      TMP6 * (OM3 * P3[2] * (TMP22 * (P1[0] * 4. * (V1[5] * P3[1] - V1[3] *
      P3[3]) + (P1[1] * 4. * (V1[2] * P3[3] - V1[5] * P3[0]) + 4. * (P1[3] *
      (V1[3] * P3[0] - V1[2] * P3[1])))) + 1.333333333 * (P3[2] * TMP52)) +
      (P2[2] * (P1[0] * 4. * (V1[3] * P3[3] - V1[5] * P3[1]) + (P1[1] * 4. *
      (V1[5] * P3[0] - V1[2] * P3[3]) + 4. * (P1[3] * (V1[2] * P3[1] - V1[3] *
      P3[0])))) + 1.333333333 * (TMP52))));
  T3[16] = denom * 0.500000000 * cI * (OM3 * (P3[2] * (TMP4 * (TMP23 * (P2[0] *
      2. * (V2[3] * P3[2] - V2[4] * P3[1]) + (P2[1] * 2. * (V2[4] * P3[0] -
      V2[2] * P3[2]) + 2. * (P2[2] * (V2[2] * P3[1] - V2[3] * P3[0])))) +
      1.333333333 * (P3[3] * TMP51)) + TMP6 * (TMP22 * (P1[0] * 2. * (V1[3] *
      P3[2] - V1[4] * P3[1]) + (P1[1] * 2. * (V1[4] * P3[0] - V1[2] * P3[2]) +
      2. * (P1[2] * (V1[2] * P3[1] - V1[3] * P3[0])))) + 1.333333333 * (P3[3] *
      TMP52))) + P3[3] * (TMP22 * TMP6 * (P1[0] * 2. * (V1[5] * P3[1] - V1[3] *
      P3[3]) + (P1[1] * 2. * (V1[2] * P3[3] - V1[5] * P3[0]) + 2. * (P1[3] *
      (V1[3] * P3[0] - V1[2] * P3[1])))) + TMP23 * TMP4 * (P2[0] * 2. * (V2[5]
      * P3[1] - V2[3] * P3[3]) + (P2[1] * 2. * (V2[2] * P3[3] - V2[5] * P3[0])
      + 2. * (P2[3] * (V2[3] * P3[0] - V2[2] * P3[1])))))) + (TMP4 * (P1[2] *
      (P2[0] * 2. * (V2[4] * P3[1] - V2[3] * P3[2]) + (P2[1] * 2. * (V2[2] *
      P3[2] - V2[4] * P3[0]) + 2. * (P2[2] * (V2[3] * P3[0] - V2[2] * P3[1]))))
      + P1[3] * (P2[0] * 2. * (V2[3] * P3[3] - V2[5] * P3[1]) + (P2[1] * 2. *
      (V2[5] * P3[0] - V2[2] * P3[3]) + 2. * (P2[3] * (V2[2] * P3[1] - V2[3] *
      P3[0]))))) + TMP6 * (P2[2] * (P1[0] * 2. * (V1[4] * P3[1] - V1[3] *
      P3[2]) + (P1[1] * 2. * (V1[2] * P3[2] - V1[4] * P3[0]) + 2. * (P1[2] *
      (V1[3] * P3[0] - V1[2] * P3[1])))) + P2[3] * (P1[0] * 2. * (V1[3] * P3[3]
      - V1[5] * P3[1]) + (P1[1] * 2. * (V1[5] * P3[0] - V1[2] * P3[3]) + 2. *
      (P1[3] * (V1[2] * P3[1] - V1[3] * P3[0])))))));
  T3[5] = denom * 0.500000000 * cI * (OM3 * (P3[0] * (TMP4 * (TMP23 * (P2[0] *
      2. * (V2[3] * P3[2] - V2[4] * P3[1]) + (P2[1] * 2. * (V2[4] * P3[0] -
      V2[2] * P3[2]) + 2. * (P2[2] * (V2[2] * P3[1] - V2[3] * P3[0])))) +
      1.333333333 * (P3[3] * TMP51)) + TMP6 * (TMP22 * (P1[0] * 2. * (V1[3] *
      P3[2] - V1[4] * P3[1]) + (P1[1] * 2. * (V1[4] * P3[0] - V1[2] * P3[2]) +
      2. * (P1[2] * (V1[2] * P3[1] - V1[3] * P3[0])))) + 1.333333333 * (P3[3] *
      TMP52))) + P3[3] * (TMP22 * TMP6 * (P1[1] * 2. * (V1[4] * P3[3] - V1[5] *
      P3[2]) + (P1[2] * 2. * (V1[5] * P3[1] - V1[3] * P3[3]) + 2. * (P1[3] *
      (V1[3] * P3[2] - V1[4] * P3[1])))) + TMP23 * TMP4 * (P2[1] * 2. * (V2[4]
      * P3[3] - V2[5] * P3[2]) + (P2[2] * 2. * (V2[5] * P3[1] - V2[3] * P3[3])
      + 2. * (P2[3] * (V2[3] * P3[2] - V2[4] * P3[1])))))) + (TMP4 * (P1[0] *
      (P2[0] * 2. * (V2[4] * P3[1] - V2[3] * P3[2]) + (P2[1] * 2. * (V2[2] *
      P3[2] - V2[4] * P3[0]) + 2. * (P2[2] * (V2[3] * P3[0] - V2[2] * P3[1]))))
      + P1[3] * (P2[1] * 2. * (V2[5] * P3[2] - V2[4] * P3[3]) + (P2[2] * 2. *
      (V2[3] * P3[3] - V2[5] * P3[1]) + 2. * (P2[3] * (V2[4] * P3[1] - V2[3] *
      P3[2]))))) + TMP6 * (P2[0] * (P1[0] * 2. * (V1[4] * P3[1] - V1[3] *
      P3[2]) + (P1[1] * 2. * (V1[2] * P3[2] - V1[4] * P3[0]) + 2. * (P1[2] *
      (V1[3] * P3[0] - V1[2] * P3[1])))) + P2[3] * (P1[1] * 2. * (V1[5] * P3[2]
      - V1[4] * P3[3]) + (P1[2] * 2. * (V1[3] * P3[3] - V1[5] * P3[1]) + 2. *
      (P1[3] * (V1[4] * P3[1] - V1[3] * P3[2])))))));
  T3[9] = denom * 0.500000000 * cI * (OM3 * (P3[1] * (TMP4 * (TMP23 * (P2[0] *
      2. * (V2[3] * P3[2] - V2[4] * P3[1]) + (P2[1] * 2. * (V2[4] * P3[0] -
      V2[2] * P3[2]) + 2. * (P2[2] * (V2[2] * P3[1] - V2[3] * P3[0])))) +
      1.333333333 * (P3[3] * TMP51)) + TMP6 * (TMP22 * (P1[0] * 2. * (V1[3] *
      P3[2] - V1[4] * P3[1]) + (P1[1] * 2. * (V1[4] * P3[0] - V1[2] * P3[2]) +
      2. * (P1[2] * (V1[2] * P3[1] - V1[3] * P3[0])))) + 1.333333333 * (P3[3] *
      TMP52))) + P3[3] * (TMP22 * TMP6 * (P1[0] * 2. * (V1[4] * P3[3] - V1[5] *
      P3[2]) + (P1[2] * 2. * (V1[5] * P3[0] - V1[2] * P3[3]) + 2. * (P1[3] *
      (V1[2] * P3[2] - V1[4] * P3[0])))) + TMP23 * TMP4 * (P2[0] * 2. * (V2[4]
      * P3[3] - V2[5] * P3[2]) + (P2[2] * 2. * (V2[5] * P3[0] - V2[2] * P3[3])
      + 2. * (P2[3] * (V2[2] * P3[2] - V2[4] * P3[0])))))) + (TMP4 * (P1[1] *
      (P2[0] * 2. * (V2[4] * P3[1] - V2[3] * P3[2]) + (P2[1] * 2. * (V2[2] *
      P3[2] - V2[4] * P3[0]) + 2. * (P2[2] * (V2[3] * P3[0] - V2[2] * P3[1]))))
      + P1[3] * (P2[0] * 2. * (V2[5] * P3[2] - V2[4] * P3[3]) + (P2[2] * 2. *
      (V2[2] * P3[3] - V2[5] * P3[0]) + 2. * (P2[3] * (V2[4] * P3[0] - V2[2] *
      P3[2]))))) + TMP6 * (P2[1] * (P1[0] * 2. * (V1[4] * P3[1] - V1[3] *
      P3[2]) + (P1[1] * 2. * (V1[2] * P3[2] - V1[4] * P3[0]) + 2. * (P1[2] *
      (V1[3] * P3[0] - V1[2] * P3[1])))) + P2[3] * (P1[0] * 2. * (V1[5] * P3[2]
      - V1[4] * P3[3]) + (P1[2] * 2. * (V1[2] * P3[3] - V1[5] * P3[0]) + 2. *
      (P1[3] * (V1[4] * P3[0] - V1[2] * P3[2])))))));
  T3[13] = denom * 0.500000000 * cI * (OM3 * (P3[2] * (TMP4 * (TMP23 * (P2[0] *
      2. * (V2[3] * P3[2] - V2[4] * P3[1]) + (P2[1] * 2. * (V2[4] * P3[0] -
      V2[2] * P3[2]) + 2. * (P2[2] * (V2[2] * P3[1] - V2[3] * P3[0])))) +
      1.333333333 * (P3[3] * TMP51)) + TMP6 * (TMP22 * (P1[0] * 2. * (V1[3] *
      P3[2] - V1[4] * P3[1]) + (P1[1] * 2. * (V1[4] * P3[0] - V1[2] * P3[2]) +
      2. * (P1[2] * (V1[2] * P3[1] - V1[3] * P3[0])))) + 1.333333333 * (P3[3] *
      TMP52))) + P3[3] * (TMP22 * TMP6 * (P1[0] * 2. * (V1[5] * P3[1] - V1[3] *
      P3[3]) + (P1[1] * 2. * (V1[2] * P3[3] - V1[5] * P3[0]) + 2. * (P1[3] *
      (V1[3] * P3[0] - V1[2] * P3[1])))) + TMP23 * TMP4 * (P2[0] * 2. * (V2[5]
      * P3[1] - V2[3] * P3[3]) + (P2[1] * 2. * (V2[2] * P3[3] - V2[5] * P3[0])
      + 2. * (P2[3] * (V2[3] * P3[0] - V2[2] * P3[1])))))) + (TMP4 * (P1[2] *
      (P2[0] * 2. * (V2[4] * P3[1] - V2[3] * P3[2]) + (P2[1] * 2. * (V2[2] *
      P3[2] - V2[4] * P3[0]) + 2. * (P2[2] * (V2[3] * P3[0] - V2[2] * P3[1]))))
      + P1[3] * (P2[0] * 2. * (V2[3] * P3[3] - V2[5] * P3[1]) + (P2[1] * 2. *
      (V2[5] * P3[0] - V2[2] * P3[3]) + 2. * (P2[3] * (V2[2] * P3[1] - V2[3] *
      P3[0]))))) + TMP6 * (P2[2] * (P1[0] * 2. * (V1[4] * P3[1] - V1[3] *
      P3[2]) + (P1[1] * 2. * (V1[2] * P3[2] - V1[4] * P3[0]) + 2. * (P1[2] *
      (V1[3] * P3[0] - V1[2] * P3[1])))) + P2[3] * (P1[0] * 2. * (V1[3] * P3[3]
      - V1[5] * P3[1]) + (P1[1] * 2. * (V1[5] * P3[0] - V1[2] * P3[3]) + 2. *
      (P1[3] * (V1[2] * P3[1] - V1[3] * P3[0])))))));
  T3[17] = denom * 0.500000000 * cI * (TMP4 * (OM3 * P3[3] * (TMP23 * (P2[0] *
      4. * (V2[3] * P3[2] - V2[4] * P3[1]) + (P2[1] * 4. * (V2[4] * P3[0] -
      V2[2] * P3[2]) + 4. * (P2[2] * (V2[2] * P3[1] - V2[3] * P3[0])))) +
      1.333333333 * (P3[3] * TMP51)) + (P1[3] * (P2[0] * 4. * (V2[4] * P3[1] -
      V2[3] * P3[2]) + (P2[1] * 4. * (V2[2] * P3[2] - V2[4] * P3[0]) + 4. *
      (P2[2] * (V2[3] * P3[0] - V2[2] * P3[1])))) + 1.333333333 * (TMP51))) +
      TMP6 * (OM3 * P3[3] * (TMP22 * (P1[0] * 4. * (V1[3] * P3[2] - V1[4] *
      P3[1]) + (P1[1] * 4. * (V1[4] * P3[0] - V1[2] * P3[2]) + 4. * (P1[2] *
      (V1[2] * P3[1] - V1[3] * P3[0])))) + 1.333333333 * (P3[3] * TMP52)) +
      (P2[3] * (P1[0] * 4. * (V1[4] * P3[1] - V1[3] * P3[2]) + (P1[1] * 4. *
      (V1[2] * P3[2] - V1[4] * P3[0]) + 4. * (P1[2] * (V1[3] * P3[0] - V1[2] *
      P3[1])))) + 1.333333333 * (TMP52))));
}


void VVT20_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP17; 
  complex<double> TMP16; 
  TMP17 = (V1[2] * - 1. * (V2[3] * T3[3] + V2[4] * T3[4] + V2[5] * T3[5] -
      V2[2] * T3[2]) + (V1[3] * (V2[3] * T3[7] + V2[4] * T3[8] + V2[5] * T3[9]
      - V2[2] * T3[6]) + (V1[4] * (V2[3] * T3[11] + V2[4] * T3[12] + V2[5] *
      T3[13] - V2[2] * T3[10]) + V1[5] * (V2[3] * T3[15] + V2[4] * T3[16] +
      V2[5] * T3[17] - V2[2] * T3[14]))));
  TMP16 = (V1[2] * - 1. * (V2[3] * T3[6] + V2[4] * T3[10] + V2[5] * T3[14] -
      V2[2] * T3[2]) + (V1[3] * (V2[3] * T3[7] + V2[4] * T3[11] + V2[5] *
      T3[15] - V2[2] * T3[3]) + (V1[4] * (V2[3] * T3[8] + V2[4] * T3[12] +
      V2[5] * T3[16] - V2[2] * T3[4]) + V1[5] * (V2[3] * T3[9] + V2[4] * T3[13]
      + V2[5] * T3[17] - V2[2] * T3[5]))));
  vertex = COUP * (-cI * (TMP16 + TMP17)); 
}


void VVT15_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP39; 
  complex<double> TMP37; 
  double P1[4]; 
  complex<double> TMP44; 
  double P2[4]; 
  double P3[4]; 
  complex<double> TMP6; 
  complex<double> TMP40; 
  complex<double> TMP4; 
  complex<double> TMP41; 
  complex<double> TMP42; 
  complex<double> TMP43; 
  complex<double> TMP38; 
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
  TMP42 = (P1[0] * (P2[0] * (P3[1] * (V1[5] * T3[10] - V1[4] * T3[14]) + (P3[2]
      * (V1[3] * T3[14] - V1[5] * T3[6]) + P3[3] * (V1[4] * T3[6] - V1[3] *
      T3[10]))) + (P2[1] * (P3[0] * (V1[4] * T3[14] - V1[5] * T3[10]) + (P3[2]
      * (V1[5] * T3[2] - V1[2] * T3[14]) + P3[3] * (V1[2] * T3[10] - V1[4] *
      T3[2]))) + (P2[2] * (P3[0] * (V1[5] * T3[6] - V1[3] * T3[14]) + (P3[1] *
      (V1[2] * T3[14] - V1[5] * T3[2]) + P3[3] * (V1[3] * T3[2] - V1[2] *
      T3[6]))) + P2[3] * (P3[0] * (V1[3] * T3[10] - V1[4] * T3[6]) + (P3[1] *
      (V1[4] * T3[2] - V1[2] * T3[10]) + P3[2] * (V1[2] * T3[6] - V1[3] *
      T3[2])))))) + (P1[1] * (P2[0] * (P3[1] * (V1[4] * T3[15] - V1[5] *
      T3[11]) + (P3[2] * (V1[5] * T3[7] - V1[3] * T3[15]) + P3[3] * (V1[3] *
      T3[11] - V1[4] * T3[7]))) + (P2[1] * (P3[0] * (V1[5] * T3[11] - V1[4] *
      T3[15]) + (P3[2] * (V1[2] * T3[15] - V1[5] * T3[3]) + P3[3] * (V1[4] *
      T3[3] - V1[2] * T3[11]))) + (P2[2] * (P3[0] * (V1[3] * T3[15] - V1[5] *
      T3[7]) + (P3[1] * (V1[5] * T3[3] - V1[2] * T3[15]) + P3[3] * (V1[2] *
      T3[7] - V1[3] * T3[3]))) + P2[3] * (P3[0] * (V1[4] * T3[7] - V1[3] *
      T3[11]) + (P3[1] * (V1[2] * T3[11] - V1[4] * T3[3]) + P3[2] * (V1[3] *
      T3[3] - V1[2] * T3[7])))))) + (P1[2] * (P2[0] * (P3[1] * (V1[4] * T3[16]
      - V1[5] * T3[12]) + (P3[2] * (V1[5] * T3[8] - V1[3] * T3[16]) + P3[3] *
      (V1[3] * T3[12] - V1[4] * T3[8]))) + (P2[1] * (P3[0] * (V1[5] * T3[12] -
      V1[4] * T3[16]) + (P3[2] * (V1[2] * T3[16] - V1[5] * T3[4]) + P3[3] *
      (V1[4] * T3[4] - V1[2] * T3[12]))) + (P2[2] * (P3[0] * (V1[3] * T3[16] -
      V1[5] * T3[8]) + (P3[1] * (V1[5] * T3[4] - V1[2] * T3[16]) + P3[3] *
      (V1[2] * T3[8] - V1[3] * T3[4]))) + P2[3] * (P3[0] * (V1[4] * T3[8] -
      V1[3] * T3[12]) + (P3[1] * (V1[2] * T3[12] - V1[4] * T3[4]) + P3[2] *
      (V1[3] * T3[4] - V1[2] * T3[8])))))) + P1[3] * (P2[0] * (P3[1] * (V1[4] *
      T3[17] - V1[5] * T3[13]) + (P3[2] * (V1[5] * T3[9] - V1[3] * T3[17]) +
      P3[3] * (V1[3] * T3[13] - V1[4] * T3[9]))) + (P2[1] * (P3[0] * (V1[5] *
      T3[13] - V1[4] * T3[17]) + (P3[2] * (V1[2] * T3[17] - V1[5] * T3[5]) +
      P3[3] * (V1[4] * T3[5] - V1[2] * T3[13]))) + (P2[2] * (P3[0] * (V1[3] *
      T3[17] - V1[5] * T3[9]) + (P3[1] * (V1[5] * T3[5] - V1[2] * T3[17]) +
      P3[3] * (V1[2] * T3[9] - V1[3] * T3[5]))) + P2[3] * (P3[0] * (V1[4] *
      T3[9] - V1[3] * T3[13]) + (P3[1] * (V1[2] * T3[13] - V1[4] * T3[5]) +
      P3[2] * (V1[3] * T3[5] - V1[2] * T3[9])))))))));
  TMP43 = (P2[0] * (P1[0] * (P3[1] * (V2[5] * T3[10] - V2[4] * T3[14]) + (P3[2]
      * (V2[3] * T3[14] - V2[5] * T3[6]) + P3[3] * (V2[4] * T3[6] - V2[3] *
      T3[10]))) + (P1[1] * (P3[0] * (V2[4] * T3[14] - V2[5] * T3[10]) + (P3[2]
      * (V2[5] * T3[2] - V2[2] * T3[14]) + P3[3] * (V2[2] * T3[10] - V2[4] *
      T3[2]))) + (P1[2] * (P3[0] * (V2[5] * T3[6] - V2[3] * T3[14]) + (P3[1] *
      (V2[2] * T3[14] - V2[5] * T3[2]) + P3[3] * (V2[3] * T3[2] - V2[2] *
      T3[6]))) + P1[3] * (P3[0] * (V2[3] * T3[10] - V2[4] * T3[6]) + (P3[1] *
      (V2[4] * T3[2] - V2[2] * T3[10]) + P3[2] * (V2[2] * T3[6] - V2[3] *
      T3[2])))))) + (P2[1] * (P1[0] * (P3[1] * (V2[4] * T3[15] - V2[5] *
      T3[11]) + (P3[2] * (V2[5] * T3[7] - V2[3] * T3[15]) + P3[3] * (V2[3] *
      T3[11] - V2[4] * T3[7]))) + (P1[1] * (P3[0] * (V2[5] * T3[11] - V2[4] *
      T3[15]) + (P3[2] * (V2[2] * T3[15] - V2[5] * T3[3]) + P3[3] * (V2[4] *
      T3[3] - V2[2] * T3[11]))) + (P1[2] * (P3[0] * (V2[3] * T3[15] - V2[5] *
      T3[7]) + (P3[1] * (V2[5] * T3[3] - V2[2] * T3[15]) + P3[3] * (V2[2] *
      T3[7] - V2[3] * T3[3]))) + P1[3] * (P3[0] * (V2[4] * T3[7] - V2[3] *
      T3[11]) + (P3[1] * (V2[2] * T3[11] - V2[4] * T3[3]) + P3[2] * (V2[3] *
      T3[3] - V2[2] * T3[7])))))) + (P2[2] * (P1[0] * (P3[1] * (V2[4] * T3[16]
      - V2[5] * T3[12]) + (P3[2] * (V2[5] * T3[8] - V2[3] * T3[16]) + P3[3] *
      (V2[3] * T3[12] - V2[4] * T3[8]))) + (P1[1] * (P3[0] * (V2[5] * T3[12] -
      V2[4] * T3[16]) + (P3[2] * (V2[2] * T3[16] - V2[5] * T3[4]) + P3[3] *
      (V2[4] * T3[4] - V2[2] * T3[12]))) + (P1[2] * (P3[0] * (V2[3] * T3[16] -
      V2[5] * T3[8]) + (P3[1] * (V2[5] * T3[4] - V2[2] * T3[16]) + P3[3] *
      (V2[2] * T3[8] - V2[3] * T3[4]))) + P1[3] * (P3[0] * (V2[4] * T3[8] -
      V2[3] * T3[12]) + (P3[1] * (V2[2] * T3[12] - V2[4] * T3[4]) + P3[2] *
      (V2[3] * T3[4] - V2[2] * T3[8])))))) + P2[3] * (P1[0] * (P3[1] * (V2[4] *
      T3[17] - V2[5] * T3[13]) + (P3[2] * (V2[5] * T3[9] - V2[3] * T3[17]) +
      P3[3] * (V2[3] * T3[13] - V2[4] * T3[9]))) + (P1[1] * (P3[0] * (V2[5] *
      T3[13] - V2[4] * T3[17]) + (P3[2] * (V2[2] * T3[17] - V2[5] * T3[5]) +
      P3[3] * (V2[4] * T3[5] - V2[2] * T3[13]))) + (P1[2] * (P3[0] * (V2[3] *
      T3[17] - V2[5] * T3[9]) + (P3[1] * (V2[5] * T3[5] - V2[2] * T3[17]) +
      P3[3] * (V2[2] * T3[9] - V2[3] * T3[5]))) + P1[3] * (P3[0] * (V2[4] *
      T3[9] - V2[3] * T3[13]) + (P3[1] * (V2[2] * T3[13] - V2[4] * T3[5]) +
      P3[2] * (V2[3] * T3[5] - V2[2] * T3[9])))))))));
  TMP40 = (P2[0] * (P1[0] * (P3[1] * (V1[5] * T3[4] - V1[4] * T3[5]) + (P3[2] *
      (V1[3] * T3[5] - V1[5] * T3[3]) + P3[3] * (V1[4] * T3[3] - V1[3] *
      T3[4]))) + (P1[1] * (P3[0] * (V1[4] * T3[5] - V1[5] * T3[4]) + (P3[2] *
      (V1[5] * T3[2] - V1[2] * T3[5]) + P3[3] * (V1[2] * T3[4] - V1[4] *
      T3[2]))) + (P1[2] * (P3[0] * (V1[5] * T3[3] - V1[3] * T3[5]) + (P3[1] *
      (V1[2] * T3[5] - V1[5] * T3[2]) + P3[3] * (V1[3] * T3[2] - V1[2] *
      T3[3]))) + P1[3] * (P3[0] * (V1[3] * T3[4] - V1[4] * T3[3]) + (P3[1] *
      (V1[4] * T3[2] - V1[2] * T3[4]) + P3[2] * (V1[2] * T3[3] - V1[3] *
      T3[2])))))) + (P2[1] * (P1[0] * (P3[1] * (V1[4] * T3[9] - V1[5] * T3[8])
      + (P3[2] * (V1[5] * T3[7] - V1[3] * T3[9]) + P3[3] * (V1[3] * T3[8] -
      V1[4] * T3[7]))) + (P1[1] * (P3[0] * (V1[5] * T3[8] - V1[4] * T3[9]) +
      (P3[2] * (V1[2] * T3[9] - V1[5] * T3[6]) + P3[3] * (V1[4] * T3[6] - V1[2]
      * T3[8]))) + (P1[2] * (P3[0] * (V1[3] * T3[9] - V1[5] * T3[7]) + (P3[1] *
      (V1[5] * T3[6] - V1[2] * T3[9]) + P3[3] * (V1[2] * T3[7] - V1[3] *
      T3[6]))) + P1[3] * (P3[0] * (V1[4] * T3[7] - V1[3] * T3[8]) + (P3[1] *
      (V1[2] * T3[8] - V1[4] * T3[6]) + P3[2] * (V1[3] * T3[6] - V1[2] *
      T3[7])))))) + (P2[2] * (P1[0] * (P3[1] * (V1[4] * T3[13] - V1[5] *
      T3[12]) + (P3[2] * (V1[5] * T3[11] - V1[3] * T3[13]) + P3[3] * (V1[3] *
      T3[12] - V1[4] * T3[11]))) + (P1[1] * (P3[0] * (V1[5] * T3[12] - V1[4] *
      T3[13]) + (P3[2] * (V1[2] * T3[13] - V1[5] * T3[10]) + P3[3] * (V1[4] *
      T3[10] - V1[2] * T3[12]))) + (P1[2] * (P3[0] * (V1[3] * T3[13] - V1[5] *
      T3[11]) + (P3[1] * (V1[5] * T3[10] - V1[2] * T3[13]) + P3[3] * (V1[2] *
      T3[11] - V1[3] * T3[10]))) + P1[3] * (P3[0] * (V1[4] * T3[11] - V1[3] *
      T3[12]) + (P3[1] * (V1[2] * T3[12] - V1[4] * T3[10]) + P3[2] * (V1[3] *
      T3[10] - V1[2] * T3[11])))))) + P2[3] * (P1[0] * (P3[1] * (V1[4] * T3[17]
      - V1[5] * T3[16]) + (P3[2] * (V1[5] * T3[15] - V1[3] * T3[17]) + P3[3] *
      (V1[3] * T3[16] - V1[4] * T3[15]))) + (P1[1] * (P3[0] * (V1[5] * T3[16] -
      V1[4] * T3[17]) + (P3[2] * (V1[2] * T3[17] - V1[5] * T3[14]) + P3[3] *
      (V1[4] * T3[14] - V1[2] * T3[16]))) + (P1[2] * (P3[0] * (V1[3] * T3[17] -
      V1[5] * T3[15]) + (P3[1] * (V1[5] * T3[14] - V1[2] * T3[17]) + P3[3] *
      (V1[2] * T3[15] - V1[3] * T3[14]))) + P1[3] * (P3[0] * (V1[4] * T3[15] -
      V1[3] * T3[16]) + (P3[1] * (V1[2] * T3[16] - V1[4] * T3[14]) + P3[2] *
      (V1[3] * T3[14] - V1[2] * T3[15])))))))));
  TMP41 = (P1[0] * (P2[0] * (P3[1] * (V2[5] * T3[10] - V2[4] * T3[14]) + (P3[2]
      * (V2[3] * T3[14] - V2[5] * T3[6]) + P3[3] * (V2[4] * T3[6] - V2[3] *
      T3[10]))) + (P2[1] * (P3[0] * (V2[4] * T3[14] - V2[5] * T3[10]) + (P3[2]
      * (V2[5] * T3[2] - V2[2] * T3[14]) + P3[3] * (V2[2] * T3[10] - V2[4] *
      T3[2]))) + (P2[2] * (P3[0] * (V2[5] * T3[6] - V2[3] * T3[14]) + (P3[1] *
      (V2[2] * T3[14] - V2[5] * T3[2]) + P3[3] * (V2[3] * T3[2] - V2[2] *
      T3[6]))) + P2[3] * (P3[0] * (V2[3] * T3[10] - V2[4] * T3[6]) + (P3[1] *
      (V2[4] * T3[2] - V2[2] * T3[10]) + P3[2] * (V2[2] * T3[6] - V2[3] *
      T3[2])))))) + (P1[1] * (P2[0] * (P3[1] * (V2[4] * T3[15] - V2[5] *
      T3[11]) + (P3[2] * (V2[5] * T3[7] - V2[3] * T3[15]) + P3[3] * (V2[3] *
      T3[11] - V2[4] * T3[7]))) + (P2[1] * (P3[0] * (V2[5] * T3[11] - V2[4] *
      T3[15]) + (P3[2] * (V2[2] * T3[15] - V2[5] * T3[3]) + P3[3] * (V2[4] *
      T3[3] - V2[2] * T3[11]))) + (P2[2] * (P3[0] * (V2[3] * T3[15] - V2[5] *
      T3[7]) + (P3[1] * (V2[5] * T3[3] - V2[2] * T3[15]) + P3[3] * (V2[2] *
      T3[7] - V2[3] * T3[3]))) + P2[3] * (P3[0] * (V2[4] * T3[7] - V2[3] *
      T3[11]) + (P3[1] * (V2[2] * T3[11] - V2[4] * T3[3]) + P3[2] * (V2[3] *
      T3[3] - V2[2] * T3[7])))))) + (P1[2] * (P2[0] * (P3[1] * (V2[4] * T3[16]
      - V2[5] * T3[12]) + (P3[2] * (V2[5] * T3[8] - V2[3] * T3[16]) + P3[3] *
      (V2[3] * T3[12] - V2[4] * T3[8]))) + (P2[1] * (P3[0] * (V2[5] * T3[12] -
      V2[4] * T3[16]) + (P3[2] * (V2[2] * T3[16] - V2[5] * T3[4]) + P3[3] *
      (V2[4] * T3[4] - V2[2] * T3[12]))) + (P2[2] * (P3[0] * (V2[3] * T3[16] -
      V2[5] * T3[8]) + (P3[1] * (V2[5] * T3[4] - V2[2] * T3[16]) + P3[3] *
      (V2[2] * T3[8] - V2[3] * T3[4]))) + P2[3] * (P3[0] * (V2[4] * T3[8] -
      V2[3] * T3[12]) + (P3[1] * (V2[2] * T3[12] - V2[4] * T3[4]) + P3[2] *
      (V2[3] * T3[4] - V2[2] * T3[8])))))) + P1[3] * (P2[0] * (P3[1] * (V2[4] *
      T3[17] - V2[5] * T3[13]) + (P3[2] * (V2[5] * T3[9] - V2[3] * T3[17]) +
      P3[3] * (V2[3] * T3[13] - V2[4] * T3[9]))) + (P2[1] * (P3[0] * (V2[5] *
      T3[13] - V2[4] * T3[17]) + (P3[2] * (V2[2] * T3[17] - V2[5] * T3[5]) +
      P3[3] * (V2[4] * T3[5] - V2[2] * T3[13]))) + (P2[2] * (P3[0] * (V2[3] *
      T3[17] - V2[5] * T3[9]) + (P3[1] * (V2[5] * T3[5] - V2[2] * T3[17]) +
      P3[3] * (V2[2] * T3[9] - V2[3] * T3[5]))) + P2[3] * (P3[0] * (V2[4] *
      T3[9] - V2[3] * T3[13]) + (P3[1] * (V2[2] * T3[13] - V2[4] * T3[5]) +
      P3[2] * (V2[3] * T3[5] - V2[2] * T3[9])))))))));
  TMP44 = (P2[0] * (P1[0] * (P3[1] * (V1[5] * T3[10] - V1[4] * T3[14]) + (P3[2]
      * (V1[3] * T3[14] - V1[5] * T3[6]) + P3[3] * (V1[4] * T3[6] - V1[3] *
      T3[10]))) + (P1[1] * (P3[0] * (V1[4] * T3[14] - V1[5] * T3[10]) + (P3[2]
      * (V1[5] * T3[2] - V1[2] * T3[14]) + P3[3] * (V1[2] * T3[10] - V1[4] *
      T3[2]))) + (P1[2] * (P3[0] * (V1[5] * T3[6] - V1[3] * T3[14]) + (P3[1] *
      (V1[2] * T3[14] - V1[5] * T3[2]) + P3[3] * (V1[3] * T3[2] - V1[2] *
      T3[6]))) + P1[3] * (P3[0] * (V1[3] * T3[10] - V1[4] * T3[6]) + (P3[1] *
      (V1[4] * T3[2] - V1[2] * T3[10]) + P3[2] * (V1[2] * T3[6] - V1[3] *
      T3[2])))))) + (P2[1] * (P1[0] * (P3[1] * (V1[4] * T3[15] - V1[5] *
      T3[11]) + (P3[2] * (V1[5] * T3[7] - V1[3] * T3[15]) + P3[3] * (V1[3] *
      T3[11] - V1[4] * T3[7]))) + (P1[1] * (P3[0] * (V1[5] * T3[11] - V1[4] *
      T3[15]) + (P3[2] * (V1[2] * T3[15] - V1[5] * T3[3]) + P3[3] * (V1[4] *
      T3[3] - V1[2] * T3[11]))) + (P1[2] * (P3[0] * (V1[3] * T3[15] - V1[5] *
      T3[7]) + (P3[1] * (V1[5] * T3[3] - V1[2] * T3[15]) + P3[3] * (V1[2] *
      T3[7] - V1[3] * T3[3]))) + P1[3] * (P3[0] * (V1[4] * T3[7] - V1[3] *
      T3[11]) + (P3[1] * (V1[2] * T3[11] - V1[4] * T3[3]) + P3[2] * (V1[3] *
      T3[3] - V1[2] * T3[7])))))) + (P2[2] * (P1[0] * (P3[1] * (V1[4] * T3[16]
      - V1[5] * T3[12]) + (P3[2] * (V1[5] * T3[8] - V1[3] * T3[16]) + P3[3] *
      (V1[3] * T3[12] - V1[4] * T3[8]))) + (P1[1] * (P3[0] * (V1[5] * T3[12] -
      V1[4] * T3[16]) + (P3[2] * (V1[2] * T3[16] - V1[5] * T3[4]) + P3[3] *
      (V1[4] * T3[4] - V1[2] * T3[12]))) + (P1[2] * (P3[0] * (V1[3] * T3[16] -
      V1[5] * T3[8]) + (P3[1] * (V1[5] * T3[4] - V1[2] * T3[16]) + P3[3] *
      (V1[2] * T3[8] - V1[3] * T3[4]))) + P1[3] * (P3[0] * (V1[4] * T3[8] -
      V1[3] * T3[12]) + (P3[1] * (V1[2] * T3[12] - V1[4] * T3[4]) + P3[2] *
      (V1[3] * T3[4] - V1[2] * T3[8])))))) + P2[3] * (P1[0] * (P3[1] * (V1[4] *
      T3[17] - V1[5] * T3[13]) + (P3[2] * (V1[5] * T3[9] - V1[3] * T3[17]) +
      P3[3] * (V1[3] * T3[13] - V1[4] * T3[9]))) + (P1[1] * (P3[0] * (V1[5] *
      T3[13] - V1[4] * T3[17]) + (P3[2] * (V1[2] * T3[17] - V1[5] * T3[5]) +
      P3[3] * (V1[4] * T3[5] - V1[2] * T3[13]))) + (P1[2] * (P3[0] * (V1[3] *
      T3[17] - V1[5] * T3[9]) + (P3[1] * (V1[5] * T3[5] - V1[2] * T3[17]) +
      P3[3] * (V1[2] * T3[9] - V1[3] * T3[5]))) + P1[3] * (P3[0] * (V1[4] *
      T3[9] - V1[3] * T3[13]) + (P3[1] * (V1[2] * T3[13] - V1[4] * T3[5]) +
      P3[2] * (V1[3] * T3[5] - V1[2] * T3[9])))))))));
  TMP39 = (P2[0] * (P1[0] * (P3[1] * (V2[5] * T3[4] - V2[4] * T3[5]) + (P3[2] *
      (V2[3] * T3[5] - V2[5] * T3[3]) + P3[3] * (V2[4] * T3[3] - V2[3] *
      T3[4]))) + (P1[1] * (P3[0] * (V2[4] * T3[5] - V2[5] * T3[4]) + (P3[2] *
      (V2[5] * T3[2] - V2[2] * T3[5]) + P3[3] * (V2[2] * T3[4] - V2[4] *
      T3[2]))) + (P1[2] * (P3[0] * (V2[5] * T3[3] - V2[3] * T3[5]) + (P3[1] *
      (V2[2] * T3[5] - V2[5] * T3[2]) + P3[3] * (V2[3] * T3[2] - V2[2] *
      T3[3]))) + P1[3] * (P3[0] * (V2[3] * T3[4] - V2[4] * T3[3]) + (P3[1] *
      (V2[4] * T3[2] - V2[2] * T3[4]) + P3[2] * (V2[2] * T3[3] - V2[3] *
      T3[2])))))) + (P2[1] * (P1[0] * (P3[1] * (V2[4] * T3[9] - V2[5] * T3[8])
      + (P3[2] * (V2[5] * T3[7] - V2[3] * T3[9]) + P3[3] * (V2[3] * T3[8] -
      V2[4] * T3[7]))) + (P1[1] * (P3[0] * (V2[5] * T3[8] - V2[4] * T3[9]) +
      (P3[2] * (V2[2] * T3[9] - V2[5] * T3[6]) + P3[3] * (V2[4] * T3[6] - V2[2]
      * T3[8]))) + (P1[2] * (P3[0] * (V2[3] * T3[9] - V2[5] * T3[7]) + (P3[1] *
      (V2[5] * T3[6] - V2[2] * T3[9]) + P3[3] * (V2[2] * T3[7] - V2[3] *
      T3[6]))) + P1[3] * (P3[0] * (V2[4] * T3[7] - V2[3] * T3[8]) + (P3[1] *
      (V2[2] * T3[8] - V2[4] * T3[6]) + P3[2] * (V2[3] * T3[6] - V2[2] *
      T3[7])))))) + (P2[2] * (P1[0] * (P3[1] * (V2[4] * T3[13] - V2[5] *
      T3[12]) + (P3[2] * (V2[5] * T3[11] - V2[3] * T3[13]) + P3[3] * (V2[3] *
      T3[12] - V2[4] * T3[11]))) + (P1[1] * (P3[0] * (V2[5] * T3[12] - V2[4] *
      T3[13]) + (P3[2] * (V2[2] * T3[13] - V2[5] * T3[10]) + P3[3] * (V2[4] *
      T3[10] - V2[2] * T3[12]))) + (P1[2] * (P3[0] * (V2[3] * T3[13] - V2[5] *
      T3[11]) + (P3[1] * (V2[5] * T3[10] - V2[2] * T3[13]) + P3[3] * (V2[2] *
      T3[11] - V2[3] * T3[10]))) + P1[3] * (P3[0] * (V2[4] * T3[11] - V2[3] *
      T3[12]) + (P3[1] * (V2[2] * T3[12] - V2[4] * T3[10]) + P3[2] * (V2[3] *
      T3[10] - V2[2] * T3[11])))))) + P2[3] * (P1[0] * (P3[1] * (V2[4] * T3[17]
      - V2[5] * T3[16]) + (P3[2] * (V2[5] * T3[15] - V2[3] * T3[17]) + P3[3] *
      (V2[3] * T3[16] - V2[4] * T3[15]))) + (P1[1] * (P3[0] * (V2[5] * T3[16] -
      V2[4] * T3[17]) + (P3[2] * (V2[2] * T3[17] - V2[5] * T3[14]) + P3[3] *
      (V2[4] * T3[14] - V2[2] * T3[16]))) + (P1[2] * (P3[0] * (V2[3] * T3[17] -
      V2[5] * T3[15]) + (P3[1] * (V2[5] * T3[14] - V2[2] * T3[17]) + P3[3] *
      (V2[2] * T3[15] - V2[3] * T3[14]))) + P1[3] * (P3[0] * (V2[4] * T3[15] -
      V2[3] * T3[16]) + (P3[1] * (V2[2] * T3[16] - V2[4] * T3[14]) + P3[2] *
      (V2[3] * T3[14] - V2[2] * T3[15])))))))));
  TMP38 = (P1[0] * (P2[0] * (P3[1] * (V1[5] * T3[4] - V1[4] * T3[5]) + (P3[2] *
      (V1[3] * T3[5] - V1[5] * T3[3]) + P3[3] * (V1[4] * T3[3] - V1[3] *
      T3[4]))) + (P2[1] * (P3[0] * (V1[4] * T3[5] - V1[5] * T3[4]) + (P3[2] *
      (V1[5] * T3[2] - V1[2] * T3[5]) + P3[3] * (V1[2] * T3[4] - V1[4] *
      T3[2]))) + (P2[2] * (P3[0] * (V1[5] * T3[3] - V1[3] * T3[5]) + (P3[1] *
      (V1[2] * T3[5] - V1[5] * T3[2]) + P3[3] * (V1[3] * T3[2] - V1[2] *
      T3[3]))) + P2[3] * (P3[0] * (V1[3] * T3[4] - V1[4] * T3[3]) + (P3[1] *
      (V1[4] * T3[2] - V1[2] * T3[4]) + P3[2] * (V1[2] * T3[3] - V1[3] *
      T3[2])))))) + (P1[1] * (P2[0] * (P3[1] * (V1[4] * T3[9] - V1[5] * T3[8])
      + (P3[2] * (V1[5] * T3[7] - V1[3] * T3[9]) + P3[3] * (V1[3] * T3[8] -
      V1[4] * T3[7]))) + (P2[1] * (P3[0] * (V1[5] * T3[8] - V1[4] * T3[9]) +
      (P3[2] * (V1[2] * T3[9] - V1[5] * T3[6]) + P3[3] * (V1[4] * T3[6] - V1[2]
      * T3[8]))) + (P2[2] * (P3[0] * (V1[3] * T3[9] - V1[5] * T3[7]) + (P3[1] *
      (V1[5] * T3[6] - V1[2] * T3[9]) + P3[3] * (V1[2] * T3[7] - V1[3] *
      T3[6]))) + P2[3] * (P3[0] * (V1[4] * T3[7] - V1[3] * T3[8]) + (P3[1] *
      (V1[2] * T3[8] - V1[4] * T3[6]) + P3[2] * (V1[3] * T3[6] - V1[2] *
      T3[7])))))) + (P1[2] * (P2[0] * (P3[1] * (V1[4] * T3[13] - V1[5] *
      T3[12]) + (P3[2] * (V1[5] * T3[11] - V1[3] * T3[13]) + P3[3] * (V1[3] *
      T3[12] - V1[4] * T3[11]))) + (P2[1] * (P3[0] * (V1[5] * T3[12] - V1[4] *
      T3[13]) + (P3[2] * (V1[2] * T3[13] - V1[5] * T3[10]) + P3[3] * (V1[4] *
      T3[10] - V1[2] * T3[12]))) + (P2[2] * (P3[0] * (V1[3] * T3[13] - V1[5] *
      T3[11]) + (P3[1] * (V1[5] * T3[10] - V1[2] * T3[13]) + P3[3] * (V1[2] *
      T3[11] - V1[3] * T3[10]))) + P2[3] * (P3[0] * (V1[4] * T3[11] - V1[3] *
      T3[12]) + (P3[1] * (V1[2] * T3[12] - V1[4] * T3[10]) + P3[2] * (V1[3] *
      T3[10] - V1[2] * T3[11])))))) + P1[3] * (P2[0] * (P3[1] * (V1[4] * T3[17]
      - V1[5] * T3[16]) + (P3[2] * (V1[5] * T3[15] - V1[3] * T3[17]) + P3[3] *
      (V1[3] * T3[16] - V1[4] * T3[15]))) + (P2[1] * (P3[0] * (V1[5] * T3[16] -
      V1[4] * T3[17]) + (P3[2] * (V1[2] * T3[17] - V1[5] * T3[14]) + P3[3] *
      (V1[4] * T3[14] - V1[2] * T3[16]))) + (P2[2] * (P3[0] * (V1[3] * T3[17] -
      V1[5] * T3[15]) + (P3[1] * (V1[5] * T3[14] - V1[2] * T3[17]) + P3[3] *
      (V1[2] * T3[15] - V1[3] * T3[14]))) + P2[3] * (P3[0] * (V1[4] * T3[15] -
      V1[3] * T3[16]) + (P3[1] * (V1[2] * T3[16] - V1[4] * T3[14]) + P3[2] *
      (V1[3] * T3[14] - V1[2] * T3[15])))))))));
  TMP4 = (V1[2] * P3[0] - V1[3] * P3[1] - V1[4] * P3[2] - V1[5] * P3[3]); 
  TMP6 = (V2[2] * P3[0] - V2[3] * P3[1] - V2[4] * P3[2] - V2[5] * P3[3]); 
  TMP37 = (P1[0] * (P2[0] * (P3[1] * (V2[5] * T3[4] - V2[4] * T3[5]) + (P3[2] *
      (V2[3] * T3[5] - V2[5] * T3[3]) + P3[3] * (V2[4] * T3[3] - V2[3] *
      T3[4]))) + (P2[1] * (P3[0] * (V2[4] * T3[5] - V2[5] * T3[4]) + (P3[2] *
      (V2[5] * T3[2] - V2[2] * T3[5]) + P3[3] * (V2[2] * T3[4] - V2[4] *
      T3[2]))) + (P2[2] * (P3[0] * (V2[5] * T3[3] - V2[3] * T3[5]) + (P3[1] *
      (V2[2] * T3[5] - V2[5] * T3[2]) + P3[3] * (V2[3] * T3[2] - V2[2] *
      T3[3]))) + P2[3] * (P3[0] * (V2[3] * T3[4] - V2[4] * T3[3]) + (P3[1] *
      (V2[4] * T3[2] - V2[2] * T3[4]) + P3[2] * (V2[2] * T3[3] - V2[3] *
      T3[2])))))) + (P1[1] * (P2[0] * (P3[1] * (V2[4] * T3[9] - V2[5] * T3[8])
      + (P3[2] * (V2[5] * T3[7] - V2[3] * T3[9]) + P3[3] * (V2[3] * T3[8] -
      V2[4] * T3[7]))) + (P2[1] * (P3[0] * (V2[5] * T3[8] - V2[4] * T3[9]) +
      (P3[2] * (V2[2] * T3[9] - V2[5] * T3[6]) + P3[3] * (V2[4] * T3[6] - V2[2]
      * T3[8]))) + (P2[2] * (P3[0] * (V2[3] * T3[9] - V2[5] * T3[7]) + (P3[1] *
      (V2[5] * T3[6] - V2[2] * T3[9]) + P3[3] * (V2[2] * T3[7] - V2[3] *
      T3[6]))) + P2[3] * (P3[0] * (V2[4] * T3[7] - V2[3] * T3[8]) + (P3[1] *
      (V2[2] * T3[8] - V2[4] * T3[6]) + P3[2] * (V2[3] * T3[6] - V2[2] *
      T3[7])))))) + (P1[2] * (P2[0] * (P3[1] * (V2[4] * T3[13] - V2[5] *
      T3[12]) + (P3[2] * (V2[5] * T3[11] - V2[3] * T3[13]) + P3[3] * (V2[3] *
      T3[12] - V2[4] * T3[11]))) + (P2[1] * (P3[0] * (V2[5] * T3[12] - V2[4] *
      T3[13]) + (P3[2] * (V2[2] * T3[13] - V2[5] * T3[10]) + P3[3] * (V2[4] *
      T3[10] - V2[2] * T3[12]))) + (P2[2] * (P3[0] * (V2[3] * T3[13] - V2[5] *
      T3[11]) + (P3[1] * (V2[5] * T3[10] - V2[2] * T3[13]) + P3[3] * (V2[2] *
      T3[11] - V2[3] * T3[10]))) + P2[3] * (P3[0] * (V2[4] * T3[11] - V2[3] *
      T3[12]) + (P3[1] * (V2[2] * T3[12] - V2[4] * T3[10]) + P3[2] * (V2[3] *
      T3[10] - V2[2] * T3[11])))))) + P1[3] * (P2[0] * (P3[1] * (V2[4] * T3[17]
      - V2[5] * T3[16]) + (P3[2] * (V2[5] * T3[15] - V2[3] * T3[17]) + P3[3] *
      (V2[3] * T3[16] - V2[4] * T3[15]))) + (P2[1] * (P3[0] * (V2[5] * T3[16] -
      V2[4] * T3[17]) + (P3[2] * (V2[2] * T3[17] - V2[5] * T3[14]) + P3[3] *
      (V2[4] * T3[14] - V2[2] * T3[16]))) + (P2[2] * (P3[0] * (V2[3] * T3[17] -
      V2[5] * T3[15]) + (P3[1] * (V2[5] * T3[14] - V2[2] * T3[17]) + P3[3] *
      (V2[2] * T3[15] - V2[3] * T3[14]))) + P2[3] * (P3[0] * (V2[4] * T3[15] -
      V2[3] * T3[16]) + (P3[1] * (V2[2] * T3[16] - V2[4] * T3[14]) + P3[2] *
      (V2[3] * T3[14] - V2[2] * T3[15])))))))));
  vertex = COUP * (TMP4 * - 1. * (+cI * (TMP37 + TMP39 + TMP41 + TMP43)) - TMP6
      * (+cI * (TMP38 + TMP40 + TMP42 + TMP44)));
}


void VVT14_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP53; 
  double P1[4]; 
  complex<double> TMP22; 
  double P2[4]; 
  complex<double> TMP23; 
  double P3[4]; 
  complex<double> TMP6; 
  complex<double> TMP54; 
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
  TMP53 = -1. * (P1[0] * (P2[1] * (V1[4] * P3[3] - V1[5] * P3[2]) + (P2[2] *
      (V1[5] * P3[1] - V1[3] * P3[3]) + P2[3] * (V1[3] * P3[2] - V1[4] *
      P3[1]))) + (P1[1] * (P2[0] * (V1[5] * P3[2] - V1[4] * P3[3]) + (P2[2] *
      (V1[2] * P3[3] - V1[5] * P3[0]) + P2[3] * (V1[4] * P3[0] - V1[2] *
      P3[2]))) + (P1[2] * (P2[0] * (V1[3] * P3[3] - V1[5] * P3[1]) + (P2[1] *
      (V1[5] * P3[0] - V1[2] * P3[3]) + P2[3] * (V1[2] * P3[1] - V1[3] *
      P3[0]))) + P1[3] * (P2[0] * (V1[4] * P3[1] - V1[3] * P3[2]) + (P2[1] *
      (V1[2] * P3[2] - V1[4] * P3[0]) + P2[2] * (V1[3] * P3[0] - V1[2] *
      P3[1]))))));
  TMP54 = -1. * (P1[0] * (P2[1] * (V2[5] * P3[2] - V2[4] * P3[3]) + (P2[2] *
      (V2[3] * P3[3] - V2[5] * P3[1]) + P2[3] * (V2[4] * P3[1] - V2[3] *
      P3[2]))) + (P1[1] * (P2[0] * (V2[4] * P3[3] - V2[5] * P3[2]) + (P2[2] *
      (V2[5] * P3[0] - V2[2] * P3[3]) + P2[3] * (V2[2] * P3[2] - V2[4] *
      P3[0]))) + (P1[2] * (P2[0] * (V2[5] * P3[1] - V2[3] * P3[3]) + (P2[1] *
      (V2[2] * P3[3] - V2[5] * P3[0]) + P2[3] * (V2[3] * P3[0] - V2[2] *
      P3[1]))) + P1[3] * (P2[0] * (V2[3] * P3[2] - V2[4] * P3[1]) + (P2[1] *
      (V2[4] * P3[0] - V2[2] * P3[2]) + P2[2] * (V2[2] * P3[1] - V2[3] *
      P3[0]))))));
  TMP22 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  TMP23 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP4 = (V1[2] * P3[0] - V1[3] * P3[1] - V1[4] * P3[2] - V1[5] * P3[3]); 
  TMP6 = (V2[2] * P3[0] - V2[3] * P3[1] - V2[4] * P3[2] - V2[5] * P3[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * - 0.500000000 * cI * (TMP4 * (OM3 * P3[0] * (TMP22 * (P1[1] *
      4. * (V2[5] * P3[2] - V2[4] * P3[3]) + (P1[2] * 4. * (V2[3] * P3[3] -
      V2[5] * P3[1]) + 4. * (P1[3] * (V2[4] * P3[1] - V2[3] * P3[2])))) -
      1.333333333 * (P3[0] * TMP54)) + (P2[0] * (P1[1] * 4. * (V2[4] * P3[3] -
      V2[5] * P3[2]) + (P1[2] * 4. * (V2[5] * P3[1] - V2[3] * P3[3]) + 4. *
      (P1[3] * (V2[3] * P3[2] - V2[4] * P3[1])))) + 1.333333333 * (TMP54))) +
      TMP6 * (OM3 * P3[0] * (TMP23 * (P2[1] * 4. * (V1[5] * P3[2] - V1[4] *
      P3[3]) + (P2[2] * 4. * (V1[3] * P3[3] - V1[5] * P3[1]) + 4. * (P2[3] *
      (V1[4] * P3[1] - V1[3] * P3[2])))) - 1.333333333 * (P3[0] * TMP53)) +
      (P1[0] * (P2[1] * 4. * (V1[4] * P3[3] - V1[5] * P3[2]) + (P2[2] * 4. *
      (V1[5] * P3[1] - V1[3] * P3[3]) + 4. * (P2[3] * (V1[3] * P3[2] - V1[4] *
      P3[1])))) + 1.333333333 * (TMP53))));
  T3[6] = denom * - 0.500000000 * cI * (OM3 * (P3[0] * (TMP4 * (TMP22 * (P1[0]
      * 2. * (V2[5] * P3[2] - V2[4] * P3[3]) + (P1[2] * 2. * (V2[2] * P3[3] -
      V2[5] * P3[0]) + 2. * (P1[3] * (V2[4] * P3[0] - V2[2] * P3[2])))) -
      1.333333333 * (P3[1] * TMP54)) + TMP6 * (TMP23 * (P2[0] * 2. * (V1[5] *
      P3[2] - V1[4] * P3[3]) + (P2[2] * 2. * (V1[2] * P3[3] - V1[5] * P3[0]) +
      2. * (P2[3] * (V1[4] * P3[0] - V1[2] * P3[2])))) - 1.333333333 * (P3[1] *
      TMP53))) + P3[1] * (TMP22 * TMP4 * (P1[1] * 2. * (V2[5] * P3[2] - V2[4] *
      P3[3]) + (P1[2] * 2. * (V2[3] * P3[3] - V2[5] * P3[1]) + 2. * (P1[3] *
      (V2[4] * P3[1] - V2[3] * P3[2])))) + TMP23 * TMP6 * (P2[1] * 2. * (V1[5]
      * P3[2] - V1[4] * P3[3]) + (P2[2] * 2. * (V1[3] * P3[3] - V1[5] * P3[1])
      + 2. * (P2[3] * (V1[4] * P3[1] - V1[3] * P3[2])))))) + (TMP4 * (P2[0] *
      (P1[0] * 2. * (V2[4] * P3[3] - V2[5] * P3[2]) + (P1[2] * 2. * (V2[5] *
      P3[0] - V2[2] * P3[3]) + 2. * (P1[3] * (V2[2] * P3[2] - V2[4] * P3[0]))))
      + P2[1] * (P1[1] * 2. * (V2[4] * P3[3] - V2[5] * P3[2]) + (P1[2] * 2. *
      (V2[5] * P3[1] - V2[3] * P3[3]) + 2. * (P1[3] * (V2[3] * P3[2] - V2[4] *
      P3[1]))))) + TMP6 * (P1[0] * (P2[0] * 2. * (V1[4] * P3[3] - V1[5] *
      P3[2]) + (P2[2] * 2. * (V1[5] * P3[0] - V1[2] * P3[3]) + 2. * (P2[3] *
      (V1[2] * P3[2] - V1[4] * P3[0])))) + P1[1] * (P2[1] * 2. * (V1[4] * P3[3]
      - V1[5] * P3[2]) + (P2[2] * 2. * (V1[5] * P3[1] - V1[3] * P3[3]) + 2. *
      (P2[3] * (V1[3] * P3[2] - V1[4] * P3[1])))))));
  T3[10] = denom * - 0.500000000 * cI * (OM3 * (P3[0] * (TMP4 * (TMP22 * (P1[0]
      * 2. * (V2[3] * P3[3] - V2[5] * P3[1]) + (P1[1] * 2. * (V2[5] * P3[0] -
      V2[2] * P3[3]) + 2. * (P1[3] * (V2[2] * P3[1] - V2[3] * P3[0])))) -
      1.333333333 * (P3[2] * TMP54)) + TMP6 * (TMP23 * (P2[0] * 2. * (V1[3] *
      P3[3] - V1[5] * P3[1]) + (P2[1] * 2. * (V1[5] * P3[0] - V1[2] * P3[3]) +
      2. * (P2[3] * (V1[2] * P3[1] - V1[3] * P3[0])))) - 1.333333333 * (P3[2] *
      TMP53))) + P3[2] * (TMP22 * TMP4 * (P1[1] * 2. * (V2[5] * P3[2] - V2[4] *
      P3[3]) + (P1[2] * 2. * (V2[3] * P3[3] - V2[5] * P3[1]) + 2. * (P1[3] *
      (V2[4] * P3[1] - V2[3] * P3[2])))) + TMP23 * TMP6 * (P2[1] * 2. * (V1[5]
      * P3[2] - V1[4] * P3[3]) + (P2[2] * 2. * (V1[3] * P3[3] - V1[5] * P3[1])
      + 2. * (P2[3] * (V1[4] * P3[1] - V1[3] * P3[2])))))) + (TMP4 * (P2[0] *
      (P1[0] * 2. * (V2[5] * P3[1] - V2[3] * P3[3]) + (P1[1] * 2. * (V2[2] *
      P3[3] - V2[5] * P3[0]) + 2. * (P1[3] * (V2[3] * P3[0] - V2[2] * P3[1]))))
      + P2[2] * (P1[1] * 2. * (V2[4] * P3[3] - V2[5] * P3[2]) + (P1[2] * 2. *
      (V2[5] * P3[1] - V2[3] * P3[3]) + 2. * (P1[3] * (V2[3] * P3[2] - V2[4] *
      P3[1]))))) + TMP6 * (P1[0] * (P2[0] * 2. * (V1[5] * P3[1] - V1[3] *
      P3[3]) + (P2[1] * 2. * (V1[2] * P3[3] - V1[5] * P3[0]) + 2. * (P2[3] *
      (V1[3] * P3[0] - V1[2] * P3[1])))) + P1[2] * (P2[1] * 2. * (V1[4] * P3[3]
      - V1[5] * P3[2]) + (P2[2] * 2. * (V1[5] * P3[1] - V1[3] * P3[3]) + 2. *
      (P2[3] * (V1[3] * P3[2] - V1[4] * P3[1])))))));
  T3[14] = denom * - 0.500000000 * cI * (OM3 * (P3[0] * (TMP4 * (TMP22 * (P1[0]
      * 2. * (V2[4] * P3[1] - V2[3] * P3[2]) + (P1[1] * 2. * (V2[2] * P3[2] -
      V2[4] * P3[0]) + 2. * (P1[2] * (V2[3] * P3[0] - V2[2] * P3[1])))) -
      1.333333333 * (P3[3] * TMP54)) + TMP6 * (TMP23 * (P2[0] * 2. * (V1[4] *
      P3[1] - V1[3] * P3[2]) + (P2[1] * 2. * (V1[2] * P3[2] - V1[4] * P3[0]) +
      2. * (P2[2] * (V1[3] * P3[0] - V1[2] * P3[1])))) - 1.333333333 * (P3[3] *
      TMP53))) + P3[3] * (TMP22 * TMP4 * (P1[1] * 2. * (V2[5] * P3[2] - V2[4] *
      P3[3]) + (P1[2] * 2. * (V2[3] * P3[3] - V2[5] * P3[1]) + 2. * (P1[3] *
      (V2[4] * P3[1] - V2[3] * P3[2])))) + TMP23 * TMP6 * (P2[1] * 2. * (V1[5]
      * P3[2] - V1[4] * P3[3]) + (P2[2] * 2. * (V1[3] * P3[3] - V1[5] * P3[1])
      + 2. * (P2[3] * (V1[4] * P3[1] - V1[3] * P3[2])))))) + (TMP4 * (P2[0] *
      (P1[0] * 2. * (V2[3] * P3[2] - V2[4] * P3[1]) + (P1[1] * 2. * (V2[4] *
      P3[0] - V2[2] * P3[2]) + 2. * (P1[2] * (V2[2] * P3[1] - V2[3] * P3[0]))))
      + P2[3] * (P1[1] * 2. * (V2[4] * P3[3] - V2[5] * P3[2]) + (P1[2] * 2. *
      (V2[5] * P3[1] - V2[3] * P3[3]) + 2. * (P1[3] * (V2[3] * P3[2] - V2[4] *
      P3[1]))))) + TMP6 * (P1[0] * (P2[0] * 2. * (V1[3] * P3[2] - V1[4] *
      P3[1]) + (P2[1] * 2. * (V1[4] * P3[0] - V1[2] * P3[2]) + 2. * (P2[2] *
      (V1[2] * P3[1] - V1[3] * P3[0])))) + P1[3] * (P2[1] * 2. * (V1[4] * P3[3]
      - V1[5] * P3[2]) + (P2[2] * 2. * (V1[5] * P3[1] - V1[3] * P3[3]) + 2. *
      (P2[3] * (V1[3] * P3[2] - V1[4] * P3[1])))))));
  T3[3] = denom * 0.500000000 * cI * (OM3 * (P3[0] * (TMP4 * (TMP22 * (P1[0] *
      2. * (V2[4] * P3[3] - V2[5] * P3[2]) + (P1[2] * 2. * (V2[5] * P3[0] -
      V2[2] * P3[3]) + 2. * (P1[3] * (V2[2] * P3[2] - V2[4] * P3[0])))) +
      1.333333333 * (P3[1] * TMP54)) + TMP6 * (TMP23 * (P2[0] * 2. * (V1[4] *
      P3[3] - V1[5] * P3[2]) + (P2[2] * 2. * (V1[5] * P3[0] - V1[2] * P3[3]) +
      2. * (P2[3] * (V1[2] * P3[2] - V1[4] * P3[0])))) + 1.333333333 * (P3[1] *
      TMP53))) + P3[1] * (TMP22 * TMP4 * (P1[1] * 2. * (V2[4] * P3[3] - V2[5] *
      P3[2]) + (P1[2] * 2. * (V2[5] * P3[1] - V2[3] * P3[3]) + 2. * (P1[3] *
      (V2[3] * P3[2] - V2[4] * P3[1])))) + TMP23 * TMP6 * (P2[1] * 2. * (V1[4]
      * P3[3] - V1[5] * P3[2]) + (P2[2] * 2. * (V1[5] * P3[1] - V1[3] * P3[3])
      + 2. * (P2[3] * (V1[3] * P3[2] - V1[4] * P3[1])))))) + (TMP4 * (P2[0] *
      (P1[0] * 2. * (V2[5] * P3[2] - V2[4] * P3[3]) + (P1[2] * 2. * (V2[2] *
      P3[3] - V2[5] * P3[0]) + 2. * (P1[3] * (V2[4] * P3[0] - V2[2] * P3[2]))))
      + P2[1] * (P1[1] * 2. * (V2[5] * P3[2] - V2[4] * P3[3]) + (P1[2] * 2. *
      (V2[3] * P3[3] - V2[5] * P3[1]) + 2. * (P1[3] * (V2[4] * P3[1] - V2[3] *
      P3[2]))))) + TMP6 * (P1[0] * (P2[0] * 2. * (V1[5] * P3[2] - V1[4] *
      P3[3]) + (P2[2] * 2. * (V1[2] * P3[3] - V1[5] * P3[0]) + 2. * (P2[3] *
      (V1[4] * P3[0] - V1[2] * P3[2])))) + P1[1] * (P2[1] * 2. * (V1[5] * P3[2]
      - V1[4] * P3[3]) + (P2[2] * 2. * (V1[3] * P3[3] - V1[5] * P3[1]) + 2. *
      (P2[3] * (V1[4] * P3[1] - V1[3] * P3[2])))))));
  T3[7] = denom * 0.500000000 * cI * (TMP4 * (OM3 * P3[1] * (TMP22 * (P1[0] *
      4. * (V2[4] * P3[3] - V2[5] * P3[2]) + (P1[2] * 4. * (V2[5] * P3[0] -
      V2[2] * P3[3]) + 4. * (P1[3] * (V2[2] * P3[2] - V2[4] * P3[0])))) +
      1.333333333 * (P3[1] * TMP54)) + (P2[1] * (P1[0] * 4. * (V2[5] * P3[2] -
      V2[4] * P3[3]) + (P1[2] * 4. * (V2[2] * P3[3] - V2[5] * P3[0]) + 4. *
      (P1[3] * (V2[4] * P3[0] - V2[2] * P3[2])))) + 1.333333333 * (TMP54))) +
      TMP6 * (OM3 * P3[1] * (TMP23 * (P2[0] * 4. * (V1[4] * P3[3] - V1[5] *
      P3[2]) + (P2[2] * 4. * (V1[5] * P3[0] - V1[2] * P3[3]) + 4. * (P2[3] *
      (V1[2] * P3[2] - V1[4] * P3[0])))) + 1.333333333 * (P3[1] * TMP53)) +
      (P1[1] * (P2[0] * 4. * (V1[5] * P3[2] - V1[4] * P3[3]) + (P2[2] * 4. *
      (V1[2] * P3[3] - V1[5] * P3[0]) + 4. * (P2[3] * (V1[4] * P3[0] - V1[2] *
      P3[2])))) + 1.333333333 * (TMP53))));
  T3[11] = denom * 0.500000000 * cI * (OM3 * (P3[1] * (TMP4 * (TMP22 * (P1[0] *
      2. * (V2[5] * P3[1] - V2[3] * P3[3]) + (P1[1] * 2. * (V2[2] * P3[3] -
      V2[5] * P3[0]) + 2. * (P1[3] * (V2[3] * P3[0] - V2[2] * P3[1])))) +
      1.333333333 * (P3[2] * TMP54)) + TMP6 * (TMP23 * (P2[0] * 2. * (V1[5] *
      P3[1] - V1[3] * P3[3]) + (P2[1] * 2. * (V1[2] * P3[3] - V1[5] * P3[0]) +
      2. * (P2[3] * (V1[3] * P3[0] - V1[2] * P3[1])))) + 1.333333333 * (P3[2] *
      TMP53))) + P3[2] * (TMP22 * TMP4 * (P1[0] * 2. * (V2[4] * P3[3] - V2[5] *
      P3[2]) + (P1[2] * 2. * (V2[5] * P3[0] - V2[2] * P3[3]) + 2. * (P1[3] *
      (V2[2] * P3[2] - V2[4] * P3[0])))) + TMP23 * TMP6 * (P2[0] * 2. * (V1[4]
      * P3[3] - V1[5] * P3[2]) + (P2[2] * 2. * (V1[5] * P3[0] - V1[2] * P3[3])
      + 2. * (P2[3] * (V1[2] * P3[2] - V1[4] * P3[0])))))) + (TMP4 * (P2[1] *
      (P1[0] * 2. * (V2[3] * P3[3] - V2[5] * P3[1]) + (P1[1] * 2. * (V2[5] *
      P3[0] - V2[2] * P3[3]) + 2. * (P1[3] * (V2[2] * P3[1] - V2[3] * P3[0]))))
      + P2[2] * (P1[0] * 2. * (V2[5] * P3[2] - V2[4] * P3[3]) + (P1[2] * 2. *
      (V2[2] * P3[3] - V2[5] * P3[0]) + 2. * (P1[3] * (V2[4] * P3[0] - V2[2] *
      P3[2]))))) + TMP6 * (P1[1] * (P2[0] * 2. * (V1[3] * P3[3] - V1[5] *
      P3[1]) + (P2[1] * 2. * (V1[5] * P3[0] - V1[2] * P3[3]) + 2. * (P2[3] *
      (V1[2] * P3[1] - V1[3] * P3[0])))) + P1[2] * (P2[0] * 2. * (V1[5] * P3[2]
      - V1[4] * P3[3]) + (P2[2] * 2. * (V1[2] * P3[3] - V1[5] * P3[0]) + 2. *
      (P2[3] * (V1[4] * P3[0] - V1[2] * P3[2])))))));
  T3[15] = denom * 0.500000000 * cI * (OM3 * (P3[1] * (TMP4 * (TMP22 * (P1[0] *
      2. * (V2[3] * P3[2] - V2[4] * P3[1]) + (P1[1] * 2. * (V2[4] * P3[0] -
      V2[2] * P3[2]) + 2. * (P1[2] * (V2[2] * P3[1] - V2[3] * P3[0])))) +
      1.333333333 * (P3[3] * TMP54)) + TMP6 * (TMP23 * (P2[0] * 2. * (V1[3] *
      P3[2] - V1[4] * P3[1]) + (P2[1] * 2. * (V1[4] * P3[0] - V1[2] * P3[2]) +
      2. * (P2[2] * (V1[2] * P3[1] - V1[3] * P3[0])))) + 1.333333333 * (P3[3] *
      TMP53))) + P3[3] * (TMP22 * TMP4 * (P1[0] * 2. * (V2[4] * P3[3] - V2[5] *
      P3[2]) + (P1[2] * 2. * (V2[5] * P3[0] - V2[2] * P3[3]) + 2. * (P1[3] *
      (V2[2] * P3[2] - V2[4] * P3[0])))) + TMP23 * TMP6 * (P2[0] * 2. * (V1[4]
      * P3[3] - V1[5] * P3[2]) + (P2[2] * 2. * (V1[5] * P3[0] - V1[2] * P3[3])
      + 2. * (P2[3] * (V1[2] * P3[2] - V1[4] * P3[0])))))) + (TMP4 * (P2[1] *
      (P1[0] * 2. * (V2[4] * P3[1] - V2[3] * P3[2]) + (P1[1] * 2. * (V2[2] *
      P3[2] - V2[4] * P3[0]) + 2. * (P1[2] * (V2[3] * P3[0] - V2[2] * P3[1]))))
      + P2[3] * (P1[0] * 2. * (V2[5] * P3[2] - V2[4] * P3[3]) + (P1[2] * 2. *
      (V2[2] * P3[3] - V2[5] * P3[0]) + 2. * (P1[3] * (V2[4] * P3[0] - V2[2] *
      P3[2]))))) + TMP6 * (P1[1] * (P2[0] * 2. * (V1[4] * P3[1] - V1[3] *
      P3[2]) + (P2[1] * 2. * (V1[2] * P3[2] - V1[4] * P3[0]) + 2. * (P2[2] *
      (V1[3] * P3[0] - V1[2] * P3[1])))) + P1[3] * (P2[0] * 2. * (V1[5] * P3[2]
      - V1[4] * P3[3]) + (P2[2] * 2. * (V1[2] * P3[3] - V1[5] * P3[0]) + 2. *
      (P2[3] * (V1[4] * P3[0] - V1[2] * P3[2])))))));
  T3[4] = denom * 0.500000000 * cI * (OM3 * (P3[0] * (TMP4 * (TMP22 * (P1[0] *
      2. * (V2[5] * P3[1] - V2[3] * P3[3]) + (P1[1] * 2. * (V2[2] * P3[3] -
      V2[5] * P3[0]) + 2. * (P1[3] * (V2[3] * P3[0] - V2[2] * P3[1])))) +
      1.333333333 * (P3[2] * TMP54)) + TMP6 * (TMP23 * (P2[0] * 2. * (V1[5] *
      P3[1] - V1[3] * P3[3]) + (P2[1] * 2. * (V1[2] * P3[3] - V1[5] * P3[0]) +
      2. * (P2[3] * (V1[3] * P3[0] - V1[2] * P3[1])))) + 1.333333333 * (P3[2] *
      TMP53))) + P3[2] * (TMP22 * TMP4 * (P1[1] * 2. * (V2[4] * P3[3] - V2[5] *
      P3[2]) + (P1[2] * 2. * (V2[5] * P3[1] - V2[3] * P3[3]) + 2. * (P1[3] *
      (V2[3] * P3[2] - V2[4] * P3[1])))) + TMP23 * TMP6 * (P2[1] * 2. * (V1[4]
      * P3[3] - V1[5] * P3[2]) + (P2[2] * 2. * (V1[5] * P3[1] - V1[3] * P3[3])
      + 2. * (P2[3] * (V1[3] * P3[2] - V1[4] * P3[1])))))) + (TMP4 * (P2[0] *
      (P1[0] * 2. * (V2[3] * P3[3] - V2[5] * P3[1]) + (P1[1] * 2. * (V2[5] *
      P3[0] - V2[2] * P3[3]) + 2. * (P1[3] * (V2[2] * P3[1] - V2[3] * P3[0]))))
      + P2[2] * (P1[1] * 2. * (V2[5] * P3[2] - V2[4] * P3[3]) + (P1[2] * 2. *
      (V2[3] * P3[3] - V2[5] * P3[1]) + 2. * (P1[3] * (V2[4] * P3[1] - V2[3] *
      P3[2]))))) + TMP6 * (P1[0] * (P2[0] * 2. * (V1[3] * P3[3] - V1[5] *
      P3[1]) + (P2[1] * 2. * (V1[5] * P3[0] - V1[2] * P3[3]) + 2. * (P2[3] *
      (V1[2] * P3[1] - V1[3] * P3[0])))) + P1[2] * (P2[1] * 2. * (V1[5] * P3[2]
      - V1[4] * P3[3]) + (P2[2] * 2. * (V1[3] * P3[3] - V1[5] * P3[1]) + 2. *
      (P2[3] * (V1[4] * P3[1] - V1[3] * P3[2])))))));
  T3[8] = denom * 0.500000000 * cI * (OM3 * (P3[1] * (TMP4 * (TMP22 * (P1[0] *
      2. * (V2[5] * P3[1] - V2[3] * P3[3]) + (P1[1] * 2. * (V2[2] * P3[3] -
      V2[5] * P3[0]) + 2. * (P1[3] * (V2[3] * P3[0] - V2[2] * P3[1])))) +
      1.333333333 * (P3[2] * TMP54)) + TMP6 * (TMP23 * (P2[0] * 2. * (V1[5] *
      P3[1] - V1[3] * P3[3]) + (P2[1] * 2. * (V1[2] * P3[3] - V1[5] * P3[0]) +
      2. * (P2[3] * (V1[3] * P3[0] - V1[2] * P3[1])))) + 1.333333333 * (P3[2] *
      TMP53))) + P3[2] * (TMP22 * TMP4 * (P1[0] * 2. * (V2[4] * P3[3] - V2[5] *
      P3[2]) + (P1[2] * 2. * (V2[5] * P3[0] - V2[2] * P3[3]) + 2. * (P1[3] *
      (V2[2] * P3[2] - V2[4] * P3[0])))) + TMP23 * TMP6 * (P2[0] * 2. * (V1[4]
      * P3[3] - V1[5] * P3[2]) + (P2[2] * 2. * (V1[5] * P3[0] - V1[2] * P3[3])
      + 2. * (P2[3] * (V1[2] * P3[2] - V1[4] * P3[0])))))) + (TMP4 * (P2[1] *
      (P1[0] * 2. * (V2[3] * P3[3] - V2[5] * P3[1]) + (P1[1] * 2. * (V2[5] *
      P3[0] - V2[2] * P3[3]) + 2. * (P1[3] * (V2[2] * P3[1] - V2[3] * P3[0]))))
      + P2[2] * (P1[0] * 2. * (V2[5] * P3[2] - V2[4] * P3[3]) + (P1[2] * 2. *
      (V2[2] * P3[3] - V2[5] * P3[0]) + 2. * (P1[3] * (V2[4] * P3[0] - V2[2] *
      P3[2]))))) + TMP6 * (P1[1] * (P2[0] * 2. * (V1[3] * P3[3] - V1[5] *
      P3[1]) + (P2[1] * 2. * (V1[5] * P3[0] - V1[2] * P3[3]) + 2. * (P2[3] *
      (V1[2] * P3[1] - V1[3] * P3[0])))) + P1[2] * (P2[0] * 2. * (V1[5] * P3[2]
      - V1[4] * P3[3]) + (P2[2] * 2. * (V1[2] * P3[3] - V1[5] * P3[0]) + 2. *
      (P2[3] * (V1[4] * P3[0] - V1[2] * P3[2])))))));
  T3[12] = denom * 0.500000000 * cI * (TMP4 * (OM3 * P3[2] * (TMP22 * (P1[0] *
      4. * (V2[5] * P3[1] - V2[3] * P3[3]) + (P1[1] * 4. * (V2[2] * P3[3] -
      V2[5] * P3[0]) + 4. * (P1[3] * (V2[3] * P3[0] - V2[2] * P3[1])))) +
      1.333333333 * (P3[2] * TMP54)) + (P2[2] * (P1[0] * 4. * (V2[3] * P3[3] -
      V2[5] * P3[1]) + (P1[1] * 4. * (V2[5] * P3[0] - V2[2] * P3[3]) + 4. *
      (P1[3] * (V2[2] * P3[1] - V2[3] * P3[0])))) + 1.333333333 * (TMP54))) +
      TMP6 * (OM3 * P3[2] * (TMP23 * (P2[0] * 4. * (V1[5] * P3[1] - V1[3] *
      P3[3]) + (P2[1] * 4. * (V1[2] * P3[3] - V1[5] * P3[0]) + 4. * (P2[3] *
      (V1[3] * P3[0] - V1[2] * P3[1])))) + 1.333333333 * (P3[2] * TMP53)) +
      (P1[2] * (P2[0] * 4. * (V1[3] * P3[3] - V1[5] * P3[1]) + (P2[1] * 4. *
      (V1[5] * P3[0] - V1[2] * P3[3]) + 4. * (P2[3] * (V1[2] * P3[1] - V1[3] *
      P3[0])))) + 1.333333333 * (TMP53))));
  T3[16] = denom * 0.500000000 * cI * (OM3 * (P3[2] * (TMP4 * (TMP22 * (P1[0] *
      2. * (V2[3] * P3[2] - V2[4] * P3[1]) + (P1[1] * 2. * (V2[4] * P3[0] -
      V2[2] * P3[2]) + 2. * (P1[2] * (V2[2] * P3[1] - V2[3] * P3[0])))) +
      1.333333333 * (P3[3] * TMP54)) + TMP6 * (TMP23 * (P2[0] * 2. * (V1[3] *
      P3[2] - V1[4] * P3[1]) + (P2[1] * 2. * (V1[4] * P3[0] - V1[2] * P3[2]) +
      2. * (P2[2] * (V1[2] * P3[1] - V1[3] * P3[0])))) + 1.333333333 * (P3[3] *
      TMP53))) + P3[3] * (TMP22 * TMP4 * (P1[0] * 2. * (V2[5] * P3[1] - V2[3] *
      P3[3]) + (P1[1] * 2. * (V2[2] * P3[3] - V2[5] * P3[0]) + 2. * (P1[3] *
      (V2[3] * P3[0] - V2[2] * P3[1])))) + TMP23 * TMP6 * (P2[0] * 2. * (V1[5]
      * P3[1] - V1[3] * P3[3]) + (P2[1] * 2. * (V1[2] * P3[3] - V1[5] * P3[0])
      + 2. * (P2[3] * (V1[3] * P3[0] - V1[2] * P3[1])))))) + (TMP4 * (P2[2] *
      (P1[0] * 2. * (V2[4] * P3[1] - V2[3] * P3[2]) + (P1[1] * 2. * (V2[2] *
      P3[2] - V2[4] * P3[0]) + 2. * (P1[2] * (V2[3] * P3[0] - V2[2] * P3[1]))))
      + P2[3] * (P1[0] * 2. * (V2[3] * P3[3] - V2[5] * P3[1]) + (P1[1] * 2. *
      (V2[5] * P3[0] - V2[2] * P3[3]) + 2. * (P1[3] * (V2[2] * P3[1] - V2[3] *
      P3[0]))))) + TMP6 * (P1[2] * (P2[0] * 2. * (V1[4] * P3[1] - V1[3] *
      P3[2]) + (P2[1] * 2. * (V1[2] * P3[2] - V1[4] * P3[0]) + 2. * (P2[2] *
      (V1[3] * P3[0] - V1[2] * P3[1])))) + P1[3] * (P2[0] * 2. * (V1[3] * P3[3]
      - V1[5] * P3[1]) + (P2[1] * 2. * (V1[5] * P3[0] - V1[2] * P3[3]) + 2. *
      (P2[3] * (V1[2] * P3[1] - V1[3] * P3[0])))))));
  T3[5] = denom * 0.500000000 * cI * (OM3 * (P3[0] * (TMP4 * (TMP22 * (P1[0] *
      2. * (V2[3] * P3[2] - V2[4] * P3[1]) + (P1[1] * 2. * (V2[4] * P3[0] -
      V2[2] * P3[2]) + 2. * (P1[2] * (V2[2] * P3[1] - V2[3] * P3[0])))) +
      1.333333333 * (P3[3] * TMP54)) + TMP6 * (TMP23 * (P2[0] * 2. * (V1[3] *
      P3[2] - V1[4] * P3[1]) + (P2[1] * 2. * (V1[4] * P3[0] - V1[2] * P3[2]) +
      2. * (P2[2] * (V1[2] * P3[1] - V1[3] * P3[0])))) + 1.333333333 * (P3[3] *
      TMP53))) + P3[3] * (TMP22 * TMP4 * (P1[1] * 2. * (V2[4] * P3[3] - V2[5] *
      P3[2]) + (P1[2] * 2. * (V2[5] * P3[1] - V2[3] * P3[3]) + 2. * (P1[3] *
      (V2[3] * P3[2] - V2[4] * P3[1])))) + TMP23 * TMP6 * (P2[1] * 2. * (V1[4]
      * P3[3] - V1[5] * P3[2]) + (P2[2] * 2. * (V1[5] * P3[1] - V1[3] * P3[3])
      + 2. * (P2[3] * (V1[3] * P3[2] - V1[4] * P3[1])))))) + (TMP4 * (P2[0] *
      (P1[0] * 2. * (V2[4] * P3[1] - V2[3] * P3[2]) + (P1[1] * 2. * (V2[2] *
      P3[2] - V2[4] * P3[0]) + 2. * (P1[2] * (V2[3] * P3[0] - V2[2] * P3[1]))))
      + P2[3] * (P1[1] * 2. * (V2[5] * P3[2] - V2[4] * P3[3]) + (P1[2] * 2. *
      (V2[3] * P3[3] - V2[5] * P3[1]) + 2. * (P1[3] * (V2[4] * P3[1] - V2[3] *
      P3[2]))))) + TMP6 * (P1[0] * (P2[0] * 2. * (V1[4] * P3[1] - V1[3] *
      P3[2]) + (P2[1] * 2. * (V1[2] * P3[2] - V1[4] * P3[0]) + 2. * (P2[2] *
      (V1[3] * P3[0] - V1[2] * P3[1])))) + P1[3] * (P2[1] * 2. * (V1[5] * P3[2]
      - V1[4] * P3[3]) + (P2[2] * 2. * (V1[3] * P3[3] - V1[5] * P3[1]) + 2. *
      (P2[3] * (V1[4] * P3[1] - V1[3] * P3[2])))))));
  T3[9] = denom * 0.500000000 * cI * (OM3 * (P3[1] * (TMP4 * (TMP22 * (P1[0] *
      2. * (V2[3] * P3[2] - V2[4] * P3[1]) + (P1[1] * 2. * (V2[4] * P3[0] -
      V2[2] * P3[2]) + 2. * (P1[2] * (V2[2] * P3[1] - V2[3] * P3[0])))) +
      1.333333333 * (P3[3] * TMP54)) + TMP6 * (TMP23 * (P2[0] * 2. * (V1[3] *
      P3[2] - V1[4] * P3[1]) + (P2[1] * 2. * (V1[4] * P3[0] - V1[2] * P3[2]) +
      2. * (P2[2] * (V1[2] * P3[1] - V1[3] * P3[0])))) + 1.333333333 * (P3[3] *
      TMP53))) + P3[3] * (TMP22 * TMP4 * (P1[0] * 2. * (V2[4] * P3[3] - V2[5] *
      P3[2]) + (P1[2] * 2. * (V2[5] * P3[0] - V2[2] * P3[3]) + 2. * (P1[3] *
      (V2[2] * P3[2] - V2[4] * P3[0])))) + TMP23 * TMP6 * (P2[0] * 2. * (V1[4]
      * P3[3] - V1[5] * P3[2]) + (P2[2] * 2. * (V1[5] * P3[0] - V1[2] * P3[3])
      + 2. * (P2[3] * (V1[2] * P3[2] - V1[4] * P3[0])))))) + (TMP4 * (P2[1] *
      (P1[0] * 2. * (V2[4] * P3[1] - V2[3] * P3[2]) + (P1[1] * 2. * (V2[2] *
      P3[2] - V2[4] * P3[0]) + 2. * (P1[2] * (V2[3] * P3[0] - V2[2] * P3[1]))))
      + P2[3] * (P1[0] * 2. * (V2[5] * P3[2] - V2[4] * P3[3]) + (P1[2] * 2. *
      (V2[2] * P3[3] - V2[5] * P3[0]) + 2. * (P1[3] * (V2[4] * P3[0] - V2[2] *
      P3[2]))))) + TMP6 * (P1[1] * (P2[0] * 2. * (V1[4] * P3[1] - V1[3] *
      P3[2]) + (P2[1] * 2. * (V1[2] * P3[2] - V1[4] * P3[0]) + 2. * (P2[2] *
      (V1[3] * P3[0] - V1[2] * P3[1])))) + P1[3] * (P2[0] * 2. * (V1[5] * P3[2]
      - V1[4] * P3[3]) + (P2[2] * 2. * (V1[2] * P3[3] - V1[5] * P3[0]) + 2. *
      (P2[3] * (V1[4] * P3[0] - V1[2] * P3[2])))))));
  T3[13] = denom * 0.500000000 * cI * (OM3 * (P3[2] * (TMP4 * (TMP22 * (P1[0] *
      2. * (V2[3] * P3[2] - V2[4] * P3[1]) + (P1[1] * 2. * (V2[4] * P3[0] -
      V2[2] * P3[2]) + 2. * (P1[2] * (V2[2] * P3[1] - V2[3] * P3[0])))) +
      1.333333333 * (P3[3] * TMP54)) + TMP6 * (TMP23 * (P2[0] * 2. * (V1[3] *
      P3[2] - V1[4] * P3[1]) + (P2[1] * 2. * (V1[4] * P3[0] - V1[2] * P3[2]) +
      2. * (P2[2] * (V1[2] * P3[1] - V1[3] * P3[0])))) + 1.333333333 * (P3[3] *
      TMP53))) + P3[3] * (TMP22 * TMP4 * (P1[0] * 2. * (V2[5] * P3[1] - V2[3] *
      P3[3]) + (P1[1] * 2. * (V2[2] * P3[3] - V2[5] * P3[0]) + 2. * (P1[3] *
      (V2[3] * P3[0] - V2[2] * P3[1])))) + TMP23 * TMP6 * (P2[0] * 2. * (V1[5]
      * P3[1] - V1[3] * P3[3]) + (P2[1] * 2. * (V1[2] * P3[3] - V1[5] * P3[0])
      + 2. * (P2[3] * (V1[3] * P3[0] - V1[2] * P3[1])))))) + (TMP4 * (P2[2] *
      (P1[0] * 2. * (V2[4] * P3[1] - V2[3] * P3[2]) + (P1[1] * 2. * (V2[2] *
      P3[2] - V2[4] * P3[0]) + 2. * (P1[2] * (V2[3] * P3[0] - V2[2] * P3[1]))))
      + P2[3] * (P1[0] * 2. * (V2[3] * P3[3] - V2[5] * P3[1]) + (P1[1] * 2. *
      (V2[5] * P3[0] - V2[2] * P3[3]) + 2. * (P1[3] * (V2[2] * P3[1] - V2[3] *
      P3[0]))))) + TMP6 * (P1[2] * (P2[0] * 2. * (V1[4] * P3[1] - V1[3] *
      P3[2]) + (P2[1] * 2. * (V1[2] * P3[2] - V1[4] * P3[0]) + 2. * (P2[2] *
      (V1[3] * P3[0] - V1[2] * P3[1])))) + P1[3] * (P2[0] * 2. * (V1[3] * P3[3]
      - V1[5] * P3[1]) + (P2[1] * 2. * (V1[5] * P3[0] - V1[2] * P3[3]) + 2. *
      (P2[3] * (V1[2] * P3[1] - V1[3] * P3[0])))))));
  T3[17] = denom * 0.500000000 * cI * (TMP4 * (OM3 * P3[3] * (TMP22 * (P1[0] *
      4. * (V2[3] * P3[2] - V2[4] * P3[1]) + (P1[1] * 4. * (V2[4] * P3[0] -
      V2[2] * P3[2]) + 4. * (P1[2] * (V2[2] * P3[1] - V2[3] * P3[0])))) +
      1.333333333 * (P3[3] * TMP54)) + (P2[3] * (P1[0] * 4. * (V2[4] * P3[1] -
      V2[3] * P3[2]) + (P1[1] * 4. * (V2[2] * P3[2] - V2[4] * P3[0]) + 4. *
      (P1[2] * (V2[3] * P3[0] - V2[2] * P3[1])))) + 1.333333333 * (TMP54))) +
      TMP6 * (OM3 * P3[3] * (TMP23 * (P2[0] * 4. * (V1[3] * P3[2] - V1[4] *
      P3[1]) + (P2[1] * 4. * (V1[4] * P3[0] - V1[2] * P3[2]) + 4. * (P2[2] *
      (V1[2] * P3[1] - V1[3] * P3[0])))) + 1.333333333 * (P3[3] * TMP53)) +
      (P1[3] * (P2[0] * 4. * (V1[4] * P3[1] - V1[3] * P3[2]) + (P2[1] * 4. *
      (V1[2] * P3[2] - V1[4] * P3[0]) + 4. * (P2[2] * (V1[3] * P3[0] - V1[2] *
      P3[1])))) + 1.333333333 * (TMP53))));
}


}  // end namespace $(namespace)s_HEF_UF

