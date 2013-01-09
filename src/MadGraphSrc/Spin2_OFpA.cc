//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph 5 v. 1.5.5, 2012-11-18
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#include "Spin2_OFpA.h"
#include "HelAmps_HZZ_Unitary_spin2pA.h"
#include "read_slha.h"

using namespace MG5_HZZ_Unitary_spin2pA; 

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: g g > xg > e- e+ mu- mu+ a / h zp WEIGHTED=10

//--------------------------------------------------------------------------
// Initialize process.

void Spin2_OFpA::initProc(string param_card_name) 
{
  // Instantiate the model class and set parameters that stay fixed during run
  pars = Parameters_HEF_UFO::getInstance(); 
  SLHAReader slha(param_card_name); 
  pars->setIndependentParameters(slha); 
  pars->setIndependentCouplings(); 
  // Set external particle masses for this matrix element
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->Me); 
  mME.push_back(pars->Me); 
  mME.push_back(pars->MM); 
  mME.push_back(pars->MM); 
  mME.push_back(pars->ZERO); 
  jamp2[0] = new double[1]; 
}

//--------------------------------------------------------------------------
// Update process.

void Spin2_OFpA::updateProc(SLHAReader slha) 
{
  pars->setIndependentParameters(slha); 
  pars->setIndependentCouplings(); 
  
  // Set external particle masses for this matrix element
  mME[0]=(pars->ZERO);
  mME[1]=(pars->ZERO);
  mME[2]=(pars->Me);
  mME[3]=(pars->Me);
  mME[4]=(pars->MM);
  mME[5]=(pars->MM);
  mME[6]=(pars->ZERO);
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void Spin2_OFpA::sigmaKin() 
{
  // Set the parameters which change event by event
  pars->setDependentParameters(); 
  pars->setDependentCouplings();

  // Reset color flows
  for(int i = 0; i < 1; i++ )
    jamp2[0][i] = 0.; 

  // Local variables and constants
  const int ncomb = 128; 
  static bool goodhel[ncomb] = {ncomb * false}; 
  static int ntry = 0, sum_hel = 0, ngood = 0; 
  static int igood[ncomb]; 
  static int jhel; 
//   std::complex<double> * * wfs; 
  double t[nprocesses]; 
  // Helicities for the process
  static const int helicities[ncomb][nexternal] = {{-1, -1, -1, -1, -1, -1,
      -1}, {-1, -1, -1, -1, -1, -1, 1}, {-1, -1, -1, -1, -1, 1, -1}, {-1, -1,
      -1, -1, -1, 1, 1}, {-1, -1, -1, -1, 1, -1, -1}, {-1, -1, -1, -1, 1, -1,
      1}, {-1, -1, -1, -1, 1, 1, -1}, {-1, -1, -1, -1, 1, 1, 1}, {-1, -1, -1,
      1, -1, -1, -1}, {-1, -1, -1, 1, -1, -1, 1}, {-1, -1, -1, 1, -1, 1, -1},
      {-1, -1, -1, 1, -1, 1, 1}, {-1, -1, -1, 1, 1, -1, -1}, {-1, -1, -1, 1, 1,
      -1, 1}, {-1, -1, -1, 1, 1, 1, -1}, {-1, -1, -1, 1, 1, 1, 1}, {-1, -1, 1,
      -1, -1, -1, -1}, {-1, -1, 1, -1, -1, -1, 1}, {-1, -1, 1, -1, -1, 1, -1},
      {-1, -1, 1, -1, -1, 1, 1}, {-1, -1, 1, -1, 1, -1, -1}, {-1, -1, 1, -1, 1,
      -1, 1}, {-1, -1, 1, -1, 1, 1, -1}, {-1, -1, 1, -1, 1, 1, 1}, {-1, -1, 1,
      1, -1, -1, -1}, {-1, -1, 1, 1, -1, -1, 1}, {-1, -1, 1, 1, -1, 1, -1},
      {-1, -1, 1, 1, -1, 1, 1}, {-1, -1, 1, 1, 1, -1, -1}, {-1, -1, 1, 1, 1,
      -1, 1}, {-1, -1, 1, 1, 1, 1, -1}, {-1, -1, 1, 1, 1, 1, 1}, {-1, 1, -1,
      -1, -1, -1, -1}, {-1, 1, -1, -1, -1, -1, 1}, {-1, 1, -1, -1, -1, 1, -1},
      {-1, 1, -1, -1, -1, 1, 1}, {-1, 1, -1, -1, 1, -1, -1}, {-1, 1, -1, -1, 1,
      -1, 1}, {-1, 1, -1, -1, 1, 1, -1}, {-1, 1, -1, -1, 1, 1, 1}, {-1, 1, -1,
      1, -1, -1, -1}, {-1, 1, -1, 1, -1, -1, 1}, {-1, 1, -1, 1, -1, 1, -1},
      {-1, 1, -1, 1, -1, 1, 1}, {-1, 1, -1, 1, 1, -1, -1}, {-1, 1, -1, 1, 1,
      -1, 1}, {-1, 1, -1, 1, 1, 1, -1}, {-1, 1, -1, 1, 1, 1, 1}, {-1, 1, 1, -1,
      -1, -1, -1}, {-1, 1, 1, -1, -1, -1, 1}, {-1, 1, 1, -1, -1, 1, -1}, {-1,
      1, 1, -1, -1, 1, 1}, {-1, 1, 1, -1, 1, -1, -1}, {-1, 1, 1, -1, 1, -1, 1},
      {-1, 1, 1, -1, 1, 1, -1}, {-1, 1, 1, -1, 1, 1, 1}, {-1, 1, 1, 1, -1, -1,
      -1}, {-1, 1, 1, 1, -1, -1, 1}, {-1, 1, 1, 1, -1, 1, -1}, {-1, 1, 1, 1,
      -1, 1, 1}, {-1, 1, 1, 1, 1, -1, -1}, {-1, 1, 1, 1, 1, -1, 1}, {-1, 1, 1,
      1, 1, 1, -1}, {-1, 1, 1, 1, 1, 1, 1}, {1, -1, -1, -1, -1, -1, -1}, {1,
      -1, -1, -1, -1, -1, 1}, {1, -1, -1, -1, -1, 1, -1}, {1, -1, -1, -1, -1,
      1, 1}, {1, -1, -1, -1, 1, -1, -1}, {1, -1, -1, -1, 1, -1, 1}, {1, -1, -1,
      -1, 1, 1, -1}, {1, -1, -1, -1, 1, 1, 1}, {1, -1, -1, 1, -1, -1, -1}, {1,
      -1, -1, 1, -1, -1, 1}, {1, -1, -1, 1, -1, 1, -1}, {1, -1, -1, 1, -1, 1,
      1}, {1, -1, -1, 1, 1, -1, -1}, {1, -1, -1, 1, 1, -1, 1}, {1, -1, -1, 1,
      1, 1, -1}, {1, -1, -1, 1, 1, 1, 1}, {1, -1, 1, -1, -1, -1, -1}, {1, -1,
      1, -1, -1, -1, 1}, {1, -1, 1, -1, -1, 1, -1}, {1, -1, 1, -1, -1, 1, 1},
      {1, -1, 1, -1, 1, -1, -1}, {1, -1, 1, -1, 1, -1, 1}, {1, -1, 1, -1, 1, 1,
      -1}, {1, -1, 1, -1, 1, 1, 1}, {1, -1, 1, 1, -1, -1, -1}, {1, -1, 1, 1,
      -1, -1, 1}, {1, -1, 1, 1, -1, 1, -1}, {1, -1, 1, 1, -1, 1, 1}, {1, -1, 1,
      1, 1, -1, -1}, {1, -1, 1, 1, 1, -1, 1}, {1, -1, 1, 1, 1, 1, -1}, {1, -1,
      1, 1, 1, 1, 1}, {1, 1, -1, -1, -1, -1, -1}, {1, 1, -1, -1, -1, -1, 1},
      {1, 1, -1, -1, -1, 1, -1}, {1, 1, -1, -1, -1, 1, 1}, {1, 1, -1, -1, 1,
      -1, -1}, {1, 1, -1, -1, 1, -1, 1}, {1, 1, -1, -1, 1, 1, -1}, {1, 1, -1,
      -1, 1, 1, 1}, {1, 1, -1, 1, -1, -1, -1}, {1, 1, -1, 1, -1, -1, 1}, {1, 1,
      -1, 1, -1, 1, -1}, {1, 1, -1, 1, -1, 1, 1}, {1, 1, -1, 1, 1, -1, -1}, {1,
      1, -1, 1, 1, -1, 1}, {1, 1, -1, 1, 1, 1, -1}, {1, 1, -1, 1, 1, 1, 1}, {1,
      1, 1, -1, -1, -1, -1}, {1, 1, 1, -1, -1, -1, 1}, {1, 1, 1, -1, -1, 1,
      -1}, {1, 1, 1, -1, -1, 1, 1}, {1, 1, 1, -1, 1, -1, -1}, {1, 1, 1, -1, 1,
      -1, 1}, {1, 1, 1, -1, 1, 1, -1}, {1, 1, 1, -1, 1, 1, 1}, {1, 1, 1, 1, -1,
      -1, -1}, {1, 1, 1, 1, -1, -1, 1}, {1, 1, 1, 1, -1, 1, -1}, {1, 1, 1, 1,
      -1, 1, 1}, {1, 1, 1, 1, 1, -1, -1}, {1, 1, 1, 1, 1, -1, 1}, {1, 1, 1, 1,
      1, 1, -1}, {1, 1, 1, 1, 1, 1, 1}};
  // Denominators: spins, colors and identical particles
  const int denominators[nprocesses] = {256}; 

  ntry = ntry + 1; 

  // Reset the matrix elements
  for(int i = 0; i < nprocesses; i++ )
  {
    matrix_element[i] = 0.; 
  }
  // Define permutation
  int perm[nexternal]; 
  for(int i = 0; i < nexternal; i++ )
  {
    perm[i] = i; 
  }

  if (sum_hel == 0 || ntry < 10)
  {
    // Calculate the matrix element for all helicities
    for(int ihel = 0; ihel < ncomb; ihel++ )
    {
      if (goodhel[ihel] || ntry < 2)
      {
        calculate_wavefunctions(perm, helicities[ihel]); 
        t[0] = matrix_gg_xg_emepmummupa_no_hzp(); 

        double tsum = 0; 
        for(int iproc = 0; iproc < nprocesses; iproc++ )
        {
          matrix_element[iproc] += t[iproc]; 
          tsum += t[iproc]; 
        }
        // Store which helicities give non-zero result
        if (tsum != 0. && !goodhel[ihel])
        {
          goodhel[ihel] = true; 
          ngood++; 
          igood[ngood] = ihel; 
        }
      }
    }
    jhel = 0; 
    sum_hel = min(sum_hel, ngood); 
  }
  else
  {
    // Only use the "good" helicities
    for(int j = 0; j < sum_hel; j++ )
    {
      jhel++; 
      if (jhel >= ngood)
        jhel = 0; 
      double hwgt = double(ngood)/double(sum_hel); 
      int ihel = igood[jhel]; 
      calculate_wavefunctions(perm, helicities[ihel]); 
      t[0] = matrix_gg_xg_emepmummupa_no_hzp(); 

      for(int iproc = 0; iproc < nprocesses; iproc++ )
      {
        matrix_element[iproc] += t[iproc] * hwgt; 
      }
    }
  }

  for (int i = 0; i < nprocesses; i++ )
    matrix_element[i] /= denominators[i]; 



}

//--------------------------------------------------------------------------
// Evaluate |M|^2, including incoming flavour dependence.

double Spin2_OFpA::sigmaHat() 
{
  // Select between the different processes
  if(id1 == 21 && id2 == 21)
  {
    // Add matrix elements for processes with beams (21, 21)
    return matrix_element[0]; 
  }
  else
  {
    // Return 0 if not correct initial state assignment
    return 0.; 
  }
}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

void Spin2_OFpA::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate wavefunctions for all processes
//   int i, j; 

  // Calculate all wavefunctions
  vxxxxx(p[perm[0]], mME[0], hel[0], -1, w[0]); 
  vxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]); 
  oxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]); 
  ixxxxx(p[perm[3]], mME[3], hel[3], -1, w[3]); 
  oxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]); 
  ixxxxx(p[perm[5]], mME[5], hel[5], -1, w[5]); 
  vxxxxx(p[perm[6]], mME[6], hel[6], +1, w[6]); 
  VVT1_10_11_12_13_3_5_7_8_9_3(w[0], w[1], pars->Unitary_GC_51, pars->Unitary_GC_37,
      pars->Unitary_GC_45, pars->Unitary_GC_29, pars->Unitary_GC_33, pars->Unitary_GC_56, pars->Unitary_GC_27,
      pars->Unitary_GC_49, pars->Unitary_GC_41, pars->Unitary_GC_47, pars->MXG, pars->WXG, w[7]);
  FFV2_4_3(w[3], w[2], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MZ, pars->WZ, w[8]); 
  FFV1_1(w[4], w[6], pars->Unitary_GC_7, pars->MM, pars->ZERO, w[9]); 
  VVT10_11_12_13_2_3_6_7_8_9_1(w[8], w[7], pars->Unitary_GC_40, pars->Unitary_GC_46,
      pars->Unitary_GC_32, pars->Unitary_GC_36, pars->Unitary_GC_55, pars->Unitary_GC_57, pars->Unitary_GC_28,
      pars->Unitary_GC_50, pars->Unitary_GC_44, pars->Unitary_GC_48, pars->MZ, pars->WZ, w[10]);
  VVT4_3(w[0], w[1], pars->Unitary_GC_28, pars->MXG, pars->WXG, w[11]); 
  VVT10_11_12_13_2_3_6_7_8_9_1(w[8], w[11], pars->Unitary_GC_40, pars->Unitary_GC_46,
      pars->Unitary_GC_32, pars->Unitary_GC_36, pars->Unitary_GC_55, pars->Unitary_GC_57, pars->Unitary_GC_28,
      pars->Unitary_GC_50, pars->Unitary_GC_44, pars->Unitary_GC_48, pars->MZ, pars->WZ, w[12]);
  FFV1_2(w[5], w[6], pars->Unitary_GC_7, pars->MM, pars->ZERO, w[13]); 
  FFV1_1(w[2], w[6], pars->Unitary_GC_7, pars->Me, pars->ZERO, w[14]); 
  FFV2_4_3(w[5], w[4], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MZ, pars->WZ, w[15]); 
  VVT10_11_12_13_2_3_6_7_8_9_1(w[15], w[7], pars->Unitary_GC_40, pars->Unitary_GC_46,
      pars->Unitary_GC_32, pars->Unitary_GC_36, pars->Unitary_GC_55, pars->Unitary_GC_57, pars->Unitary_GC_28,
      pars->Unitary_GC_50, pars->Unitary_GC_44, pars->Unitary_GC_48, pars->MZ, pars->WZ, w[16]);
  VVT10_11_12_13_2_3_6_7_8_9_1(w[15], w[11], pars->Unitary_GC_40, pars->Unitary_GC_46,
      pars->Unitary_GC_32, pars->Unitary_GC_36, pars->Unitary_GC_55, pars->Unitary_GC_57, pars->Unitary_GC_28,
      pars->Unitary_GC_50, pars->Unitary_GC_44, pars->Unitary_GC_48, pars->MZ, pars->WZ, w[17]);
  FFV1_2(w[3], w[6], pars->Unitary_GC_7, pars->Me, pars->ZERO, w[18]); 

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV2_4_0(w[5], w[9], w[10], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[0]); 
  FFV2_4_0(w[5], w[9], w[12], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[1]); 
  FFV2_4_0(w[13], w[4], w[10], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[2]); 
  FFV2_4_0(w[13], w[4], w[12], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[3]); 
  FFV2_4_0(w[3], w[14], w[16], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[4]); 
  FFV2_4_0(w[3], w[14], w[17], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[5]); 
  FFV2_4_0(w[18], w[2], w[16], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[6]); 
  FFV2_4_0(w[18], w[2], w[17], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[7]); 

}
double Spin2_OFpA::matrix_gg_xg_emepmummupa_no_hzp() 
{
  int i, j; 
  // Local variables
//   const int ngraphs = 8; 
  const int ncolor = 1; 
  std::complex<double> ztemp; 
  std::complex<double> jamp[ncolor]; 
  // The color matrix;
  static const double denom[ncolor] = {1}; 
  static const double cf[ncolor][ncolor] = {{2}}; 

  // Calculate color flows
  jamp[0] = +2. * (+amp[0] + amp[1] + amp[2] + amp[3] + amp[4] + amp[5] +
      amp[6] + amp[7]);

  // Sum and square the color flows to get the matrix element
  double matrix = 0; 
  for(i = 0; i < ncolor; i++ )
  {
    ztemp = 0.; 
    for(j = 0; j < ncolor; j++ )
      ztemp = ztemp + cf[i][j] * jamp[j]; 
    matrix = matrix + real(ztemp * conj(jamp[i]))/denom[i]; 
  }

  // Store the leading color flows for choice of color
  for(i = 0; i < ncolor; i++ )
    jamp2[0][i] += real(jamp[i] * conj(jamp[i])); 

  return matrix; 
}



