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
///  MEKD class member implementation
///
///  Provides neessary interface to the MadGraph-based ME calculator
///  and computes MEs and KDs for the process specified by the user.
///
//////////////////////////////////////////////////////////////////////////


///------------------------------------------------------------------------
/// MEKD::MEKD - a default constructor
///------------------------------------------------------------------------
MEKD::MEKD(double collisionEnergy, string PDFName)
{
	m_collisionEnergy = collisionEnergy;
	m_PDFName = PDFName;
	
	ME_ZZ = 0;
	ME_SMHiggs = 0;
	ME_CPoddScalar = 0;
	ME_Spin1 = 0;
	ME_Spin2 = 0;
}


///------------------------------------------------------------------------
/// MEKD::processParameters - sanity check for internal parameters
///------------------------------------------------------------------------
int MEKD::processParameters() {
    /// check if the PDF name is supported and set PDF flag
    if (m_PDFName!="CTEQ6L" && m_PDFName!="" && m_PDFName!="no PDFs") return ERR_PDFS;
    m_usePDF = (m_PDFName=="CTEQ6L");
    /// check if sqrt(s) is 7 or 8 TeV
	if (m_collisionEnergy!=7 && m_collisionEnergy!=8) cerr << "WARNING! You have set energy to be " << m_collisionEnergy << " TeV\n";
	
    return 0;
}


///------------------------------------------------------------------------
/// MEKD::setProcessNames - sanity check and setting of process names
///------------------------------------------------------------------------
int MEKD::setProcessNames(string processA, string processB) {
    /// processes A and B should be different
    if( processA == processB ) return ERR_PROCESS;
    /// check if processA is supported, translation of namings
    if     ( processA=="ZZ")                                 {m_processA="ZZ";}
    else if( processA=="Higgs"   || processA=="SMHiggs")     {m_processA="SMHiggs";}
    else if( processA=="CP-odd"  || processA=="Higgs0M")     {m_processA="CPoddScalar";}
    else if( processA=="CP-even" || processA=="Higgs0P")     {m_processA="CPevenScalar";}
    else if( processA=="Spin2PM" || processA=="Graviton2PM") {m_processA="Spin2particle";}
    else if( processA=="Custom")                             {m_processA="Custom";}
    else return ERR_PROCESS;
    /// check if processB is supported, translation of namings
    if     ( processB=="ZZ")                                 {m_processB="ZZ";}
    else if( processB=="Higgs"   || processB=="SMHiggs")     {m_processB="SMHiggs";}
    else if( processB=="CP-odd"  || processB=="Higgs0M")     {m_processB="CPoddScalar";}
    else if( processB=="CP-even" || processB=="Higgs0P")     {m_processB="CPevenScalar";}
    else if( processB=="Spin2PM" || processB=="Graviton2PM") {m_processB="Spin2particle";}
    else if( processB=="Custom")                             {m_processB="Custom";}
    else return ERR_PROCESS;
	
    return 0;
}


///------------------------------------------------------------------------
/// MEKD::computeKD - compute KD from precalculated MEs for input prcesses A and B
///------------------------------------------------------------------------
int MEKD::computeKD( string processA, string processB,
					double& kd,
					double& me2processA,
					double& me2processB )
{
	/// Sanity check for input process names
	if( (buffer_int=setProcessNames(processA, processB)) != 0 ) return buffer_int;
	
	/// Looking for the precalculated MEs
	if( ME_ZZ == 0 ) { cerr << "ERROR! The requested process has not been precalculated.\n" ; return ERR_PROCESS; }
	else if( m_processA=="ZZ" ) me2processA = ME_ZZ;
	else if( m_processA=="SMHiggs" ) me2processA = ME_SMHiggs;
	else if( m_processA=="CPoddScalar" ) me2processA = ME_CPoddScalar;
// 	else if( m_processA=="Spin1particle" ) me2processA = ME_Spin1;
	else if( m_processA=="Spin2particle" && ME_Spin2!=0 ) me2processA = ME_Spin2;
	else { cerr << "ERROR! The requested process has not been precalculated.\n" ; return ERR_PROCESS; }
	
	/// Looking for the precalculated MEs
	if( ME_ZZ == 0 ) { cerr << "ERROR! The requested process has not been precalculated.\n" ; return ERR_PROCESS; }
	else if( m_processB=="ZZ" ) me2processB = ME_ZZ;
	else if( m_processB=="SMHiggs" ) me2processB = ME_SMHiggs;
	else if( m_processB=="CPoddScalar" ) me2processB = ME_CPoddScalar;
// 	else if( m_processB=="Spin1particle" ) me2processB = ME_Spin1;
	else if( m_processB=="Spin2particle" ) me2processB = ME_Spin2;
	else { cerr << "ERROR! The requested process has not been precalculated.\n" ; return ERR_PROCESS; }
	
	
	/// Build Kinematic Discriminant (KD) as a ratio of logs of MEs
	kd = log(me2processA / me2processB);
	
	return 0;
}


///------------------------------------------------------------------------
/// MEKD::computeKD - compute KD and MEs for input prcesses A and B
///------------------------------------------------------------------------
int MEKD::computeKD(string processA, string processB,
                    double lept1P[], int lept1Id, double lept2P[], int lept2Id,
                    double lept3P[], int lept3Id, double lept4P[], int lept4Id,
                    double& kd, double& me2processA, double& me2processB)
{
    /// Sanity check for input process names and internal parameters
    if( (buffer_int=setProcessNames(processA, processB)) != 0 ) return buffer_int;
    if( (buffer_int=processParameters()) != 0 ) return buffer_int;

    /// Parameterize MEKD_MG_Calc
	MEKD_MG_Calc.Sqrt_s = m_collisionEnergy*1000;	// translate TeV in GeV
	MEKD_MG_Calc.Use_PDF_w_pT0 = m_usePDF;
	MEKD_MG_Calc.Test_Model = m_processA;

    /// Set input-lepton kinematics
    MEKD_MG_Calc.pl1 = lept1P;
    MEKD_MG_Calc.id1 = lept1Id;
    MEKD_MG_Calc.pl2 = lept2P;
    MEKD_MG_Calc.id2 = lept2Id;
    MEKD_MG_Calc.pl3 = lept3P;
    MEKD_MG_Calc.id3 = lept3Id;
    MEKD_MG_Calc.pl4 = lept4P;
    MEKD_MG_Calc.id4 = lept4Id;

    /// Compute ME for process A only (e.g. signal 1)
    buffer_int = MEKD_MG_Calc.Run_MEKD_MG( m_processA );
    /// Get ME for process A
    me2processA = MEKD_MG_Calc.Signal_ME;
	
    /// Compute ME for process B only (e.g. signal 2 or background)
    buffer_int = MEKD_MG_Calc.Run_MEKD_MG( m_processB );
    /// Get ME for process B
    me2processB = MEKD_MG_Calc.Signal_ME;
    /// Build Kinematic Discriminant (KD) as a ratio of logs of MEs
    kd = log(me2processA / me2processB);
	
	return buffer_int;
}


///------------------------------------------------------------------------
/// MEKD::computeMEs - compute MEs for a multiple reuse
///------------------------------------------------------------------------
int MEKD::computeMEs( double lept1P[], int lept1Id, double lept2P[], int lept2Id,
					double lept3P[], int lept3Id, double lept4P[], int lept4Id )
{
	/// Sanity check for internal parameters
	if( (buffer_int=processParameters()) != 0 ) return buffer_int;
	
	/// Parameterize MEKD_MG_Calc
	MEKD_MG_Calc.Sqrt_s = m_collisionEnergy*1000;	// translate TeV in GeV
	MEKD_MG_Calc.Use_PDF_w_pT0 = m_usePDF;
	if( MEKD_MG_Calc.Test_Models.size()==0 )	// Fills in intersting models to compute only once
	{
		MEKD_MG_Calc.Test_Models.push_back( "ZZ" );
		MEKD_MG_Calc.Test_Models.push_back( "SMHiggs" );
		MEKD_MG_Calc.Test_Models.push_back( "CPoddScalar" );
// 		MEKD_MG_Calc.Test_Models.push_back( "Spin1particle" );
		MEKD_MG_Calc.Test_Models.push_back( "Spin2particle" );
	}
	else if( MEKD_MG_Calc.Test_Models.size()!=4 ) return ERR_PROCESS;
	
	/// Set input-lepton kinematics
	MEKD_MG_Calc.pl1 = lept1P;
	MEKD_MG_Calc.id1 = lept1Id;
	MEKD_MG_Calc.pl2 = lept2P;
	MEKD_MG_Calc.id2 = lept2Id;
	MEKD_MG_Calc.pl3 = lept3P;
	MEKD_MG_Calc.id3 = lept3Id;
	MEKD_MG_Calc.pl4 = lept4P;
	MEKD_MG_Calc.id4 = lept4Id;
	
	/// Compute MEs
	buffer_int = MEKD_MG_Calc.Run_MEKD_MG();
	
	/// ME value readouts
	ME_ZZ = MEKD_MG_Calc.Signal_MEs[0];
	ME_SMHiggs = MEKD_MG_Calc.Signal_MEs[1];
	ME_CPoddScalar = MEKD_MG_Calc.Signal_MEs[2];
// 	ME_Spin1 = MEKD_MG_Calc.Signal_MEs[3];
	ME_Spin2 = MEKD_MG_Calc.Signal_MEs[3];
	
	return buffer_int;
}


#if (defined(MEKD_STANDALONE) && defined(MEKD_with_ROOT)) || !(defined(MEKD_STANDALONE))
///------------------------------------------------------------------------
/// (ROOT-compatabile overloads)
///------------------------------------------------------------------------


///------------------------------------------------------------------------
/// MEKD::computeKD - compute KD and MEs for input prcesses A and B
///------------------------------------------------------------------------
int MEKD::computeKD( TString processA, TString processB,
                     TLorentzVector lept1P, int lept1Id, TLorentzVector lept2P, int lept2Id,
                     TLorentzVector lept3P, int lept3Id, TLorentzVector lept4P, int lept4Id,
                     double& kd, double& me2processA, double& me2processB )
{
	/// prepare 4-momenta in the required format
	lept1P_i[0] = lept1P.E();
	lept1P_i[1] = lept1P.Px();
	lept1P_i[2] = lept1P.Py();
	lept1P_i[3] = lept1P.Pz();
	
	lept2P_i[0] = lept2P.E();
	lept2P_i[1] = lept2P.Px();
	lept2P_i[2] = lept2P.Py();
	lept2P_i[3] = lept2P.Pz();
	
	lept3P_i[0] = lept3P.E();
	lept3P_i[1] = lept3P.Px();
	lept3P_i[2] = lept3P.Py();
	lept3P_i[3] = lept3P.Pz();
	
	lept4P_i[0] = lept4P.E();
	lept4P_i[1] = lept4P.Px();
	lept4P_i[2] = lept4P.Py();
	lept4P_i[3] = lept4P.Pz();
	
	return computeKD( processA.Data(), processB.Data(), lept1P_i, lept1Id, lept2P_i, lept2Id, lept3P_i, lept3Id, lept4P_i, lept4Id, kd, me2processA, me2processB );
}


///------------------------------------------------------------------------
/// MEKD::computeMEs - compute MEs for a multiple reuse
///------------------------------------------------------------------------
int MEKD::computeMEs( TLorentzVector lept1P, int lept1Id, TLorentzVector lept2P, int lept2Id,
					TLorentzVector lept3P, int lept3Id, TLorentzVector lept4P, int lept4Id )
{
	/// prepare 4-momenta in the required format
	lept1P_i[0] = lept1P.E();
	lept1P_i[1] = lept1P.Px();
	lept1P_i[2] = lept1P.Py();
	lept1P_i[3] = lept1P.Pz();
	
	lept2P_i[0] = lept2P.E();
	lept2P_i[1] = lept2P.Px();
	lept2P_i[2] = lept2P.Py();
	lept2P_i[3] = lept2P.Pz();
	
	lept3P_i[0] = lept3P.E();
	lept3P_i[1] = lept3P.Px();
	lept3P_i[2] = lept3P.Py();
	lept3P_i[3] = lept3P.Pz();
	
	lept4P_i[0] = lept4P.E();
	lept4P_i[1] = lept4P.Px();
	lept4P_i[2] = lept4P.Py();
	lept4P_i[3] = lept4P.Pz();
	
	return computeMEs( lept1P_i, lept1Id, lept2P_i, lept2Id, lept3P_i, lept3Id, lept4P_i, lept4Id );
}

//////////////////////////////////////////////////////////////////////////
#endif	// end of (defined(MEKD_STANDALONE) && defined(MEKD_with_ROOT)) || !(defined(MEKD_STANDALONE))


#endif	// end of MEKD_MEKD_cpp
