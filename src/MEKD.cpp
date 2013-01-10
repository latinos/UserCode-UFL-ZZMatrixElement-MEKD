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
///  Provides necessary interface to the MadGraph-based ME calculator
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
	
	MEKD_MG_Calc.Sqrt_s = m_collisionEnergy*1000;	// translate TeV to GeV
	
	ME_ZZ = 0;
	ME_SMHiggs = 0;
	ME_CPoddScalar = 0;
	ME_Spin1 = 0;
	ME_Spin2 = 0;
	
	four_particle_IDs_i.resize( 4, 0 );
	four_particle_Ps_i.resize( 4, NULL );
}


///------------------------------------------------------------------------
/// MEKD::processParameters - sanity check for internal parameters
///------------------------------------------------------------------------
/// Might be merged to the constructor if these flags are not used as is(was)
///------------------------------------------------------------------------
int MEKD::processParameters()
{
    /// Check if the PDF name is supported and set PDF flag
    if (m_PDFName!="CTEQ6L" && m_PDFName!="" && m_PDFName!="no PDFs") return ERR_PDFS;
    m_usePDF = (m_PDFName=="CTEQ6L");
	
	MEKD_MG_Calc.Use_PDF_w_pT0 = m_usePDF;
	
    /// Check if sqrt(s) is 7 or 8 TeV
	if (m_collisionEnergy!=7 && m_collisionEnergy!=8) cerr << "WARNING! You have set energy to be " << m_collisionEnergy << " TeV\n";
	
    return 0;
}


///------------------------------------------------------------------------
/// MEKD::setProcessName - sanity check and setting of process names
///------------------------------------------------------------------------
int MEKD::setProcessName(string process)
{
    /// Check if process is supported, translation of namings
    if     ( process=="ZZ" )                                 {m_process="ZZ"; }
    else if( process=="Higgs"   || process=="SMHiggs" )     {m_process="SMHiggs"; }
    else if( process=="CP-odd"  || process=="Higgs0M" )     {m_process="CPoddScalar"; }
    else if( process=="CP-even" || process=="Higgs0P" )     {m_process="CPevenScalar"; }
    else if( process=="Spin2PM" || process=="Graviton2PM" ) {m_process="Spin2particle"; }
    else if( process=="Custom" )                             {m_process="Custom"; }
    else return ERR_PROCESS;
	
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
int MEKD::computeKD( string processA, string processB,
                     double lept1P[], int lept1Id, double lept2P[], int lept2Id,
                     double lept3P[], int lept3Id, double lept4P[], int lept4Id,
                     double& kd, double& me2processA, double& me2processB )
{
	/// Load internal containers
	four_particle_Ps_i[0] = lept1P;
	four_particle_Ps_i[1] = lept2P;
	four_particle_Ps_i[2] = lept3P;
	four_particle_Ps_i[3] = lept4P;
	
	four_particle_IDs_i[0] = lept1Id;
	four_particle_IDs_i[1] = lept2Id;
	four_particle_IDs_i[2] = lept3Id;
	four_particle_IDs_i[3] = lept4Id;
	
	return computeKD( processA, processB, four_particle_Ps_i, four_particle_IDs_i, kd, me2processA, me2processB );
}


///------------------------------------------------------------------------
/// MEKD::computeKD - compute KD and MEs for the input processes A and B
///------------------------------------------------------------------------
int MEKD::computeKD( string processA, string processB,
                     vector<double*> input_Ps, vector<int> input_IDs,
                     double& kd, double& me2processA, double& me2processB )
{
	/// Checks input for compatibility
	if( input_Ps.size() != input_IDs.size() ) return ERR_INPUT;
	if( input_Ps.size() != 4 && input_Ps.size() != 5 ) return ERR_INPUT;
	
    /// Sanity check for input process names and internal parameters
    if( (buffer_int=setProcessNames(processA, processB)) != 0 ) return buffer_int;
    if( (buffer_int=processParameters()) != 0 ) return buffer_int;

	/// Set input-particle kinematics
	MEKD_MG_Calc.p1 = input_Ps[0];
	MEKD_MG_Calc.id1 = input_IDs[0];
	MEKD_MG_Calc.p2 = input_Ps[1];
	MEKD_MG_Calc.id2 = input_IDs[1];
	MEKD_MG_Calc.p3 = input_Ps[2];
	MEKD_MG_Calc.id3 = input_IDs[2];
	MEKD_MG_Calc.p4 = input_Ps[3];
	MEKD_MG_Calc.id4 = input_IDs[3];
	if( input_IDs.size() == 5 )
	{
		MEKD_MG_Calc.p5 = input_Ps[4];
		MEKD_MG_Calc.id5 = input_IDs[4];
	}

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
/// MEKD::computeME - compute ME for the input process
///------------------------------------------------------------------------
int MEKD::computeME( string processName,
					vector<double*> input_Ps, vector<int> input_IDs,
					double& me2process )
{
	/// Checks input for compatibility
	if( input_Ps.size() != input_IDs.size() ) return ERR_INPUT;
	if( input_Ps.size() != 4 && input_Ps.size() != 5 ) return ERR_INPUT;
	
    /// Sanity check for input process names and internal parameters
    if( (buffer_int=setProcessName(processName)) != 0 ) return buffer_int;
    if( (buffer_int=processParameters()) != 0 ) return buffer_int;

	/// Set input-particle kinematics
	MEKD_MG_Calc.p1 = input_Ps[0];
	MEKD_MG_Calc.id1 = input_IDs[0];
	MEKD_MG_Calc.p2 = input_Ps[1];
	MEKD_MG_Calc.id2 = input_IDs[1];
	MEKD_MG_Calc.p3 = input_Ps[2];
	MEKD_MG_Calc.id3 = input_IDs[2];
	MEKD_MG_Calc.p4 = input_Ps[3];
	MEKD_MG_Calc.id4 = input_IDs[3];
	if( input_IDs.size() == 5 )
	{
		MEKD_MG_Calc.p5 = input_Ps[4];
		MEKD_MG_Calc.id5 = input_IDs[4];
	}
	
	/// Compute ME for the process (e.g. signal 1)
	buffer_int = MEKD_MG_Calc.Run_MEKD_MG( m_process );
	/// Get ME for the process
	me2process = MEKD_MG_Calc.Signal_ME;
	
	return buffer_int;
}


///------------------------------------------------------------------------
/// MEKD::computeMEs - compute MEs for a multiple reuse
///------------------------------------------------------------------------
int MEKD::computeMEs( double lept1P[], int lept1Id, double lept2P[], int lept2Id,
					double lept3P[], int lept3Id, double lept4P[], int lept4Id )
{
	/// Load internal containers
	four_particle_Ps_i[0] = lept1P;
	four_particle_Ps_i[1] = lept2P;
	four_particle_Ps_i[2] = lept3P;
	four_particle_Ps_i[3] = lept4P;
	
	four_particle_IDs_i[0] = lept1Id;
	four_particle_IDs_i[1] = lept2Id;
	four_particle_IDs_i[2] = lept3Id;
	four_particle_IDs_i[3] = lept4Id;
	
	return computeMEs( four_particle_Ps_i, four_particle_IDs_i );
}


///------------------------------------------------------------------------
/// MEKD::computeMEs - compute MEs for a multiple reuse
///------------------------------------------------------------------------
int MEKD::computeMEs( vector<double*> input_Ps, vector<int> input_IDs )
{
	/// Checks input for compatibility
	if( input_Ps.size() != input_IDs.size() ) return ERR_INPUT;
	if( input_Ps.size() != 4 && input_Ps.size() != 5 ) return ERR_INPUT;
	
	/// Sanity check for internal parameters
	if( (buffer_int=processParameters()) != 0 ) return buffer_int;
	
	/// Parameterize MEKD_MG_Calc
	if( MEKD_MG_Calc.Test_Models.size()==0 )	// Fills in intersting models to compute only once
	{
		MEKD_MG_Calc.Test_Models.push_back( "ZZ" );
		MEKD_MG_Calc.Test_Models.push_back( "SMHiggs" );
		MEKD_MG_Calc.Test_Models.push_back( "CPoddScalar" );
// 		MEKD_MG_Calc.Test_Models.push_back( "Spin1particle" );
		MEKD_MG_Calc.Test_Models.push_back( "Spin2particle" );
	}
	else if( MEKD_MG_Calc.Test_Models.size()!=4 ) return ERR_PROCESS;
	
	/// Set input-particle kinematics
	MEKD_MG_Calc.p1 = input_Ps[0];
	MEKD_MG_Calc.id1 = input_IDs[0];
	MEKD_MG_Calc.p2 = input_Ps[1];
	MEKD_MG_Calc.id2 = input_IDs[1];
	MEKD_MG_Calc.p3 = input_Ps[2];
	MEKD_MG_Calc.id3 = input_IDs[2];
	MEKD_MG_Calc.p4 = input_Ps[3];
	MEKD_MG_Calc.id4 = input_IDs[3];
	if( input_IDs.size() == 5 )
	{
		MEKD_MG_Calc.p5 = input_Ps[4];
		MEKD_MG_Calc.id5 = input_IDs[4];
	}
	
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
/// MEKD::computeKD - compute KD and MEs for the input processes A and B
///------------------------------------------------------------------------
int MEKD::computeKD( TString processA, TString processB,
                     TLorentzVector lept1P, int lept1Id, TLorentzVector lept2P, int lept2Id,
                     TLorentzVector lept3P, int lept3Id, TLorentzVector lept4P, int lept4Id,
                     double& kd, double& me2processA, double& me2processB )
{
	/// Prepare 4-momenta in the required format
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
	
	/// Load internal containers
	four_particle_Ps_i[0] = lept1P_i;
	four_particle_Ps_i[1] = lept2P_i;
	four_particle_Ps_i[2] = lept3P_i;
	four_particle_Ps_i[3] = lept4P_i;
	
	four_particle_IDs_i[0] = lept1Id;
	four_particle_IDs_i[1] = lept2Id;
	four_particle_IDs_i[2] = lept3Id;
	four_particle_IDs_i[3] = lept4Id;
	
	return computeKD( (string) processA.Data(), (string) processB.Data(), four_particle_Ps_i, four_particle_IDs_i, kd, me2processA, me2processB );
}


///------------------------------------------------------------------------
/// MEKD::computeKD - compute KD and MEs for the input processes A and B
///------------------------------------------------------------------------
int MEKD::computeKD( TString processA, TString processB,
					vector<TLorentzVector> input_Ps, vector<int> input_IDs,
					double& kd, double& me2processA, double& me2processB )
{
	/// Resize internal vector<double*> if needed
	if( input_Ps_i.size() != input_Ps.size() ) input_Ps_i.resize( input_Ps.size(), new double[4] );
	
	/// Put vector<TLorentzVector> into internal containers
	for( buffer_uint=0; buffer_uint < input_Ps_i.size(); buffer_uint++ )
	{
		input_Ps_i[buffer_uint][0] = input_Ps[buffer_uint].E();
		input_Ps_i[buffer_uint][1] = input_Ps[buffer_uint].Px();
		input_Ps_i[buffer_uint][2] = input_Ps[buffer_uint].Py();
		input_Ps_i[buffer_uint][3] = input_Ps[buffer_uint].Pz();
	}
	
	return computeKD( (string) processA.Data(), (string) processB.Data(), input_Ps_i, input_IDs, kd, me2processA, me2processB );
}


///------------------------------------------------------------------------
/// MEKD::computeME - compute ME for the input process
///------------------------------------------------------------------------
int MEKD::computeME( TString processName,
					vector<TLorentzVector> input_Ps, vector<int> input_IDs,
					double& me2process )
{
	/// Resize internal vector<double*> if needed
	if( input_Ps_i.size() != input_Ps.size() ) input_Ps_i.resize( input_Ps.size(), new double[4] );
	
	/// Put vector<TLorentzVector> into internal containers
	for( buffer_uint=0; buffer_uint < input_Ps_i.size(); buffer_uint++ )
	{
		input_Ps_i[buffer_uint][0] = input_Ps[buffer_uint].E();
		input_Ps_i[buffer_uint][1] = input_Ps[buffer_uint].Px();
		input_Ps_i[buffer_uint][2] = input_Ps[buffer_uint].Py();
		input_Ps_i[buffer_uint][3] = input_Ps[buffer_uint].Pz();
	}
	
	return computeME( (string) processName.Data(), input_Ps_i, input_IDs, me2process );
}


///------------------------------------------------------------------------
/// MEKD::computeMEs - compute MEs for a multiple reuse
///------------------------------------------------------------------------
int MEKD::computeMEs( TLorentzVector lept1P, int lept1Id, TLorentzVector lept2P, int lept2Id,
					TLorentzVector lept3P, int lept3Id, TLorentzVector lept4P, int lept4Id )
{
	/// Prepare 4-momenta in the required format
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
	
	/// Load internal containers
	four_particle_Ps_i[0] = lept1P_i;
	four_particle_Ps_i[1] = lept2P_i;
	four_particle_Ps_i[2] = lept3P_i;
	four_particle_Ps_i[3] = lept4P_i;
	
	four_particle_IDs_i[0] = lept1Id;
	four_particle_IDs_i[1] = lept2Id;
	four_particle_IDs_i[2] = lept3Id;
	four_particle_IDs_i[3] = lept4Id;
	
	return computeMEs( four_particle_Ps_i, four_particle_IDs_i );
}


///------------------------------------------------------------------------
/// MEKD::computeMEs - compute MEs for a multiple reuse
///------------------------------------------------------------------------
int MEKD::computeMEs( vector<TLorentzVector> input_Ps, vector<int> input_IDs )
{
	/// Resize internal vector<double*> if needed
	if( input_Ps_i.size() != input_Ps.size() ) input_Ps_i.resize( input_Ps.size(), new double[4] );
	
	/// Put vector<TLorentzVector> into internal containers
	for( buffer_uint=0; buffer_uint < input_Ps_i.size(); buffer_uint++ )
	{
		input_Ps_i[buffer_uint][0] = input_Ps[buffer_uint].E();
		input_Ps_i[buffer_uint][1] = input_Ps[buffer_uint].Px();
		input_Ps_i[buffer_uint][2] = input_Ps[buffer_uint].Py();
		input_Ps_i[buffer_uint][3] = input_Ps[buffer_uint].Pz();
	}
	
	return computeMEs( input_Ps_i, input_IDs );
}

//////////////////////////////////////////////////////////////////////////
#endif	// end of (defined(MEKD_STANDALONE) && defined(MEKD_with_ROOT)) || !(defined(MEKD_STANDALONE))


#endif	// end of MEKD_MEKD_cpp
