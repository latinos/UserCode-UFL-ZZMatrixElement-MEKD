/*************************************************************************
*  Authors:   MEKD fans
*  More info: http://mekd.ihepa.ufl.edu
*  Contact:   mekd@phys.ufl.edu
*************************************************************************/
#ifndef MEKD_MEKD_h
#define MEKD_MEKD_h

// C++ includes
#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>

#if (defined MEKD_STANDALONE && defined MEKD_with_ROOT) || !defined MEKD_STANDALONE
// ROOT includes
#include "TROOT.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include "TTree.h"
#include "TTreeIndex.h"
#include "TLorentzVector.h"
#endif

// MadGraph-based ME calculator
#include "MEKD_MG.h"


//////////////////////////////////////////////////////////////////////////
///  MEKD interface class.
///
///  Provides neessary interface to the MadGraph-based ME calculator
///  and computes MEs and KDs for the process specified by the user.
///
//////////////////////////////////////////////////////////////////////////

using namespace std;


class MEKD {
public:
    /// Constructor.
    ///
    /// \param collisionEnergy the sqrt(s) value in TeV (double, DEFAULT = 8).
    /// \param PDFName the name of the parton density functions to be used (string, DEFAULT = "CTEQ6L", NONE = "").
    ///
    MEKD(double collisionEnergy = 8, std::string PDFName = "CTEQ6L") {m_collisionEnergy = collisionEnergy; m_PDFName = PDFName;}

    /// Compute KDs and MEs for process A and process B out of the 4-momenta of 4 leptons (lepton ordering does not matter).
    ///
    /// Supported process names: "ZZ", "SMHiggs", "Higgs0M" (pseudo-scalar), "Graviton2PM" (minimal couplings graviton).
    ///
    /// \param[in]  processA, processB                   names of the processes X = A, B for which the KDs and MEs are computed (string, REQUIRED).
    /// \param[in]  lept1P, lept2P, lept3P, lept4P       the input arrays with 4-momentum (E,px,py,pz) values of leptons N=1..4 (double*, REQUIRED).
    /// \param[in]  lept1Id, lept2Id, lept3Id, lept4Id   the input IDs (PDG) of leptons N=1..4 (int, REQUIRED).
    /// \param[out] kd           the computed KD value for discrimination of processes A and B (double).
    /// \param[out] me2processA  the computed |ME|^2 for process A (double).
    /// \param[out] me2processB  the computed |ME|^2 for process B (double).
    /// \return     The error code of the computation: 0 = NO_ERR, 1 = ERR_SQRTS, 2 = ERR_PDFS, 3 = ERR_PROCESS
    ///
    int     computeKD(string processA, string processB,     // names of the processes
                      double lept1P[], int lept1Id,         // 4-momentum (E,px,py,pz) and id (PDG) of lepton 1
                      double lept2P[], int lept2Id,         // 4-momentum (E,px,py,pz) and id (PDG) of lepton 2
                      double lept3P[], int lept3Id,         // 4-momentum (E,px,py,pz) and id (PDG) of lepton 3
                      double lept4P[], int lept4Id,         // 4-momentum (E,px,py,pz) and id (PDG) of lepton 4
                      double& kd,                           // retun KD
                      double& me2processA,                  // retun |ME|^2 for process A
                      double& me2processB);                 // retun |ME|^2 for process B
private:
    // properties
    double      m_collisionEnergy;      // c.m. collision energy sqrt(s) in TeV
    std::string m_PDFName;              // name of the parton density functions to be used. Supported: CTEQ6l;
    std::string m_processA;             // name of the process A (background, signal hypotheses, etc.) Supported: Custom, SMHiggs, CPoddScalar, CPevenScalar, Spin2particle, ZZ
    std::string m_processB;             // name of the process B (background, signal hypotheses, etc.) Supported: Custom, SMHiggs, CPoddScalar, CPevenScalar, Spin2particle, ZZ
    bool        m_usePDF;               // flag to use PDFs (true) or not (false)
    bool        m_runBackgroundME;      // flat to run the ME for ZZ process (true) or not (false)
    enum        ERRCodes   {NO_ERR, ERR_SQRT, ERR_PDFS, ERR_PROCESS, NUM_ERRORS};
    // methods
    int     setProcessNames(string processA, string processB); // sanity check for input process names, translation to the the names supported by MEKD_MG.h
    int     processParameters(); // sanity check for internal paramters
    
#if (defined MEKD_STANDALONE && defined MEKD_with_ROOT) || !defined MEKD_STANDALONE
//------------------------------------------------------------------------
// ROOT-compatabile members
//------------------------------------------------------------------------
public:
    /// Compute KDs and MEs for process A and process B out of the 4-momenta of 4 leptons (lepton ordering does not matter).
    /// The overloaded method that supports input parameters of ROOT types TString and TLorentzVector.
    ///
    /// Supported process names: "ZZ", "SMHiggs", "Higgs0M" (pseudo-scalar), "Graviton2PM" (minimal couplings graviton).
    ///
    /// \param[in]  processA, processB                   names of the processes X = A, B for which the KDs and MEs are computed (TString, REQUIRED).
    /// \param[in]  lept1P, lept2P, lept3P, lept4P       the input arrays with 4-momentum (E,px,py,pz) values of leptons N=1..4 (TLorentzVector, REQUIRED).
    /// \param[in]  lept1Id, lept2Id, lept3Id, lept4Id   the input IDs (PDG) of leptons N=1..4 (int, REQUIRED).
    /// \param[out] kd           the computed KD value for discrimination of processes A and B (double).
    /// \param[out] me2processA  the computed |ME|^2 for process A (double).
    /// \param[out] me2processB  the computed |ME|^2 for process B (double).
    /// \return     The error code of the computation: 0 = NO_ERR, 1 = ERR_SQRTS, 2 = ERR_PDFS, 3 = ERR_PROCESS    
    ///
    int     computeKD(TString processA, TString processB,   // names of the processes
                      TLorentzVector lept1P, int lept1Id,   // 4-momentum and id (PDG) of lepton 1
                      TLorentzVector lept2P, int lept2Id,   // 4-momentum and id (PDG) of lepton 2
                      TLorentzVector lept3P, int lept3Id,   // 4-momentum and id (PDG) of lepton 3
                      TLorentzVector lept4P, int lept4Id,   // 4-momentum and id (PDG) of lepton 4
                      double& kd,                           // retun KD
                      double& me2processA,                  // retun |ME|^2 for process A
                      double& me2processB);                 // retun |ME|^2 for process B
#endif
};

//////////////////////////////////////////////////////////////////////////

#endif

