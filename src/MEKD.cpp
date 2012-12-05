/*************************************************************************
*  Authors:   MEKD fans
*  More info: http://mekd.ihepa.ufl.edu
*  Contact:   mekd@phys.ufl.edu
*************************************************************************/
#ifndef MEKD_MEKD_cpp
#define MEKD_MEKD_cpp

/// MEKD header
#include "../interface/MEKD.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////
//  MEKD class member implementation
//
//  Provides neessary interface to the MadGraph-based ME calculator
//  and computes MEs and KDs for the process specified by the user.
//
//////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------------
// MEKD::processParameters - sanity check for internal parameters
//------------------------------------------------------------------------
int MEKD::processParameters() {
    // check if the PDF name is supported and set PDF flag
    if (m_PDFName!="CTEQ6L" && m_PDFName!="" && m_PDFName!="no PDFs") return ERR_PDFS;
    m_usePDF = (m_PDFName=="CTEQ6L");
    // check if sqrt(s) is 7 or 8 TeV
    if (m_collisionEnergy!=7 && m_collisionEnergy!=8) return ERR_SQRT;
    return 0;
}

//------------------------------------------------------------------------
// MEKD::setProcessNames - sanity check and setting of process names
//------------------------------------------------------------------------
int MEKD::setProcessNames(string processA, string processB) {
    // processes A and B should be different
    if( processA == processB) return ERR_PROCESS;
    // check if processA is supported, translation of namings
    if     ( processA=="ZZ")                                 {m_processA="ZZ";}
    else if( processA=="Higgs"   || processA=="SMHiggs")     {m_processA="SMHiggs";}
    else if( processA=="CP-odd"  || processA=="Higgs0M")     {m_processA="CPoddScalar";}
    else if( processA=="CP-even" || processA=="Higgs0P")     {m_processA="CPevenScalar";}
    else if( processA=="Spin2PM" || processA=="Graviton2PM") {m_processA="Spin2particle";}
    else if( processA=="Custom")                             {m_processA="Custom";}
    else return ERR_PROCESS;
    // check if processB is supported, translation of namings
    if     ( processB=="ZZ")                                 {m_processB="ZZ";}
    else if( processB=="Higgs"   || processB=="SMHiggs")     {m_processB="SMHiggs";}
    else if( processB=="CP-odd"  || processB=="Higgs0M")     {m_processB="CPoddScalar";}
    else if( processB=="CP-even" || processB=="Higgs0P")     {m_processB="CPevenScalar";}
    else if( processB=="Spin2PM" || processB=="Graviton2PM") {m_processB="Spin2particle";}
    else if( processB=="Custom")                             {m_processB="Custom";}
    else return ERR_PROCESS;
    return 0;
}

//------------------------------------------------------------------------
// MEKD::computeKD - compute KD and MEs for input prcesses A and B
//------------------------------------------------------------------------
int MEKD::computeKD(string processA, string processB,
                    double lept1P[], int lept1Id, double lept2P[], int lept2Id,
                    double lept3P[], int lept3Id, double lept4P[], int lept4Id,
                    double& kd, double& me2processA, double& me2processB){
    // sanity check for input process names and internal parameters
    if (unsigned int err = setProcessNames(processA, processB)) return err;
    if (unsigned int err = processParameters()) return err;
    
    // Madgraph-based ME calculator
    MEKD_MG meKD_MG;

    // parameterise meKD_MG
	meKD_MG.Sqrt_s = m_collisionEnergy*1000.; // translate TeV in GeV
    meKD_MG.Test_Model = m_processA;
    meKD_MG.Use_PDF_w_pT0 = m_usePDF;

    // set input lepton kinematics
    meKD_MG.pl1 = lept1P;
    meKD_MG.id1 = lept1Id;
    meKD_MG.pl2 = lept2P;
    meKD_MG.id2 = lept2Id;
    meKD_MG.pl3 = lept3P;
    meKD_MG.id3 = lept3Id;
    meKD_MG.pl4 = lept4P;
    meKD_MG.id4 = lept4Id;

    // compute ME for process A only (e.g. signal 1)
    meKD_MG.Run_MEKD_MG(m_processA);
    // get ME for process A
    me2processA = meKD_MG.Signal_ME;
    // compute ME for process B only (e.g. signal 2 or background)
    meKD_MG.Run_MEKD_MG(m_processB);
    // get ME for process B
    me2processB = meKD_MG.Signal_ME;
    // build Kinematic Discriminant KD as a ratio of logs of MEs
    kd = log(me2processA / me2processB);
    return 0;
}

#if (defined(MEKD_STANDALONE) && defined(MEKD_with_ROOT)) || !(defined(MEKD_STANDALONE))
//------------------------------------------------------------------------
// MEKD::computeKD - compute KD and MEs for input prcesses A and B
// (ROOT-compatabile overloads)
//------------------------------------------------------------------------
int MEKD::computeKD(TString processA, TString processB,
                     TLorentzVector lept1P, int lept1Id, TLorentzVector lept2P, int lept2Id,
                     TLorentzVector lept3P, int lept3Id, TLorentzVector lept4P, int lept4Id,
                     double& kd, double& me2processA, double& me2processB){
    // prepare 4-momenta in the required format
    double lept1P4[4] = {lept1P.E(), lept1P.Px(), lept1P.Py(), lept1P.Pz()};
    double lept2P4[4] = {lept2P.E(), lept2P.Px(), lept2P.Py(), lept2P.Pz()};
    double lept3P4[4] = {lept3P.E(), lept3P.Px(), lept3P.Py(), lept3P.Pz()};
    double lept4P4[4] = {lept4P.E(), lept4P.Px(), lept4P.Py(), lept4P.Pz()};
    
    return (computeKD(processA.Data(), processB.Data(), lept1P4, lept1Id, lept2P4, lept2Id, lept3P4, lept3Id, lept4P4, lept4Id, kd, me2processA, me2processB));
}

//////////////////////////////////////////////////////////////////////////
#endif

#endif

