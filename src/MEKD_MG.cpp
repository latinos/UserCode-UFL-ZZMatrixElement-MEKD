#ifndef MEKD_MG_cpp
#define MEKD_MG_cpp

/// C++ libraries
#include <iostream>
// #include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>	// for sorting
#include <cmath>


/// CMSSW includes
#ifndef MEKD_STANDALONE
#include "FWCore/ParameterSet/interface/FileInPath.h"
#endif


/// MEs
#include "MadGraphSrc/BKG_DN_OF.h"
#include "MadGraphSrc/BKG_DN_SF.h"
#include "MadGraphSrc/BKG_UP_OF.h"
#include "MadGraphSrc/BKG_UP_SF.h"

#include "MadGraphSrc/Spin0_OF.h"
#include "MadGraphSrc/Spin0_SF.h"

#include "MadGraphSrc/Spin1_OF.h"
#include "MadGraphSrc/Spin1_SF.h"

#include "MadGraphSrc/Spin2_OF.h"
#include "MadGraphSrc/Spin2_SF.h"



extern "C"
{
#include "Extra_code/MEKD_CalcHEP_PDF.h"
#include "PDFTables/pdt.h"
}

#include "Extra_code/MEKD_CalcHEP_Extra_functions.h"
#include "Extra_code/MEKD_MG_Boosts.h"
#include "higgs_properties/hggeffective.h"
#include "MadGraphSrc/read_slha.h"

#include "../interface/MEKD_MG.h"

using namespace std;


/// Part of pdfreader
extern "C" pdtStr pdtSg, pdtSd, pdtSu, pdtSs, pdtSc,
	pdtSad, pdtSau, pdtSas, pdtSac;

// #define PDTFILE "PDFTables/cteq6l.pdt" // CalCHEP reads a table for CTEQ6L. You can change PDF set as you want.


BKG_DN_OF ME_Background_DownType_OF;
BKG_DN_SF ME_Background_DownType_SF;
BKG_UP_OF ME_Background_UpType_OF;
BKG_UP_SF ME_Background_UpType_SF;

Spin0_OF ME_Signal_Spin0_OF;
Spin0_SF ME_Signal_Spin0_SF;

Spin1_OF ME_Signal_Spin1_OF;
Spin1_SF ME_Signal_Spin1_SF;

Spin2_OF ME_Signal_Spin2_OF;
Spin2_SF ME_Signal_Spin2_SF;




MEKD_MG::MEKD_MG()
{
	Set_Default_MEKD_MG_Parameters();
	
	/// Cross-cheking MEs for consistency
	if( ME_Background_DownType_SF.nprocesses!=2 ) { cout << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Background_DownType_OF.nprocesses!=2 ) { cout << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Background_UpType_SF.nprocesses!=2 ) { cout << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Background_UpType_OF.nprocesses!=2 ) { cout << "Problem in ME class detected. Exiting.\n"; exit(1); }
	
	if( ME_Signal_Spin0_OF.nprocesses!=1 ) { cout << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_Spin0_SF.nprocesses!=1 ) { cout << "Problem in ME class detected. Exiting.\n"; exit(1); }
	
	if( ME_Signal_Spin1_OF.nprocesses!=1 ) { cout << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_Spin1_SF.nprocesses!=1 ) { cout << "Problem in ME class detected. Exiting.\n"; exit(1); }
	
	if( ME_Signal_Spin2_OF.nprocesses!=1 ) { cout << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_Spin2_SF.nprocesses!=1 ) { cout << "Problem in ME class detected. Exiting.\n"; exit(1); }
	
	p_set.push_back( new double[4] );
	p_set.push_back( new double[4] );
	p_set.push_back( new double[4] );
	p_set.push_back( new double[4] );
	p_set.push_back( new double[4] );
	p_set.push_back( new double[4] );
	
	pl1 = new double[4];
	pl2 = new double[4];
	pl3 = new double[4];
	pl4 = new double[4];
	
	id1=0; id2=0; id3=0; id4=0;
	
	id_set.push_back( 0 );
	id_set.push_back( 0 );
	id_set.push_back( 0 );
	id_set.push_back( 0 );
	
	Parameters_Are_Loaded = false;
}



MEKD_MG::~MEKD_MG()
{
	if( Parameters_Are_Loaded ) Unload_pdfreader();
	
	p_set.clear();
	id_set.clear();
}



int MEKD_MG::Load_Parameters()
{
	Set_Of_Model_Parameters.read_slha_file( Parameter_file );
	
	/// Initializing parameters
	ME_Background_UpType_SF.initProc( Parameter_file );
	ME_Background_UpType_OF.initProc( Parameter_file );
	ME_Background_DownType_SF.initProc( Parameter_file );
	ME_Background_DownType_OF.initProc( Parameter_file );
	
	ME_Signal_Spin0_SF.initProc( Parameter_file );
	ME_Signal_Spin0_OF.initProc( Parameter_file );
	ME_Signal_Spin1_SF.initProc( Parameter_file );
	ME_Signal_Spin1_OF.initProc( Parameter_file );
	ME_Signal_Spin2_SF.initProc( Parameter_file );
	ME_Signal_Spin2_OF.initProc( Parameter_file );
	
	m_d_temp = Set_Of_Model_Parameters.get_block_entry( "mass", 1, 0 );
	m_u_temp = Set_Of_Model_Parameters.get_block_entry( "mass", 2, 0 );
	m_s_temp = Set_Of_Model_Parameters.get_block_entry( "mass", 3, 0 );
	m_c_temp = Set_Of_Model_Parameters.get_block_entry( "mass", 4, 0 );
	m_e_temp = Set_Of_Model_Parameters.get_block_entry( "mass", 11, 0 );
	m_mu_temp = Set_Of_Model_Parameters.get_block_entry( "mass", 13, 0 );
	m_Z_temp = Set_Of_Model_Parameters.get_block_entry( "mass", 23, 9.11876e+01 );
	
	v_expectation = 1/sqrt( sqrt(2)*Set_Of_Model_Parameters.get_block_entry( "sminputs", 2, 1.166370e-05 ) );
	hZZ_coupling = 2*m_Z_temp*m_Z_temp/v_expectation;
	
	
	Load_pdfreader( const_cast<char*>(PDF_file.c_str()) );
	
	Parameters_Are_Loaded = true;
	return 0;
}



int MEKD_MG::Reload_Parameters()
{
	if( !Parameters_Are_Loaded ) return 1;
	
	Set_Of_Model_Parameters.read_slha_file( static_cast<string>(Parameter_file) );
	
	m_d_temp = Set_Of_Model_Parameters.get_block_entry( "mass", 1, 0 );
	m_u_temp = Set_Of_Model_Parameters.get_block_entry( "mass", 2, 0 );
	m_s_temp = Set_Of_Model_Parameters.get_block_entry( "mass", 3, 0 );
	m_c_temp = Set_Of_Model_Parameters.get_block_entry( "mass", 4, 0 );
	m_e_temp = Set_Of_Model_Parameters.get_block_entry( "mass", 11, 0 );
	m_mu_temp = Set_Of_Model_Parameters.get_block_entry( "mass", 13, 0 );
	m_Z_temp = Set_Of_Model_Parameters.get_block_entry( "mass", 23, 9.11876e+01 );
	
	v_expectation = 1/sqrt( sqrt(2)*Set_Of_Model_Parameters.get_block_entry( "sminputs", 2, 1.166370e-05 ) );
	hZZ_coupling = 2*m_Z_temp*m_Z_temp/v_expectation;
	
	
	Unload_pdfreader();
	Load_pdfreader( const_cast<char*>(PDF_file.c_str()) );
	
	return 0;
}



void MEKD_MG::Set_Default_MEKD_MG_Parameters()
{
	Boost_To_CM = true;	// for a boosted data
	Debug_Mode = false;	// Enable debugging mode
	Force_g3_running = false;
	Overwrite_e_and_mu_masses = false;	/// switch for manual m_e, m_mu masses
	Use_mh_eq_m4l = true;	// Set mh to m4l for every event
	Use_Higgs_width = true;	//	if false, width is fixed to =1
	Use_PDF_w_pT0 = true;	// Use PDFs in the pT=0 frame. If true, Boost_To_CM is ignored
	Vary_signal_couplings = true;	// Allow couplings to change with mass
	
	ContributionCoeff_d = 0;	//42	/// the value has no effect if PDF is used but the variable is always used
	ContributionCoeff_u = 1;	//217
	ContributionCoeff_s = 0;	//5
	ContributionCoeff_c = 0;	//3
// 	GG=0;	// Assign QCD coupling, force g3 running if needed
	Sqrt_s = 8000;	//Max energy, collision energy
	
	Electron_mass = 0;	//0.0005109989, for enabled overwriting
	Higgs_mass = 125;	// Works only if Use_mh_eq_m4l=false
	Higgs_width = 5.753088e-03;	// Practically not used, for future implementations
	Muon_mass = 0;	//0.10565837, for enabled overwriting
	Proton_mass = 0.93827205;	// Always used if needed
	
	Final_state = "2e2m";	// Final state, for the moment: 4e, 4mu, 2e2mu
	Test_Model = "1";	// 0 or Custon; 1 or SMHiggs; 2 or CPoddScalar; 3 or CPevenScalar

#ifndef MEKD_STANDALONE
	string inputParameterFile = "ZZMatrixElement/MEKD/src/Cards/param_card.dat";	// Location where a parameter card is stored
	string inputPDFFile = "ZZMatrixElement/MEKD/src/PDFTables/cteq6l.pdt";	// PDF/PDT table file
	edm::FileInPath parameterFileWithFullPath(inputParameterFile);
	edm::FileInPath pdfFileWithFullPath(inputPDFFile);
	Parameter_file = parameterFileWithFullPath.fullPath();
	PDF_file = pdfFileWithFullPath.fullPath();
#else
	Parameter_file = "src/Cards/param_card.dat";	// Location where a parameter card is stored
	PDF_file = "src/PDFTables/cteq6l.pdt";	// PDF/PDT table file
#endif
}



int MEKD_MG::Run_MEKD_MG()
{
	if( !Parameters_Are_Loaded ) Load_Parameters();
	if( Arrange_Internal_pls() == 1 ) { cerr << "Particle id error. Exiting.\n"; exit(1); }
	
	double CollisionE;
	
	PDFx1 = 0;
	PDFx2 = 0;
	Background_ME = 0;
	Signal_ME = 0;
	
	if( Overwrite_e_and_mu_masses )
	{
		Set_Of_Model_Parameters.set_block_entry( "mass", 11, Electron_mass );
		Set_Of_Model_Parameters.set_block_entry( "mass", 13, Muon_mass );
		m_e_temp = Electron_mass;
		m_mu_temp = Muon_mass;
	}
	
	if( Final_state=="4e" )
		{ ml1=Set_Of_Model_Parameters.get_block_entry( "mass", 11, Electron_mass ); ml2=ml1; ml3=ml1; ml4=ml1; }
	if( Final_state=="4m" || Final_state=="4mu" )
		{ ml1=Set_Of_Model_Parameters.get_block_entry( "mass", 13, Electron_mass ); ml2=ml1; ml3=ml1; ml4=ml1; }
	if( Final_state=="2e2m" || Final_state=="2e2mu" )
		{ ml1=Set_Of_Model_Parameters.get_block_entry( "mass", 11, Electron_mass ); ml2=ml1; ml3=Set_Of_Model_Parameters.get_block_entry( "mass", 13, Electron_mass ); ml4=ml3; }
	
	
	/// No boosting setup for initial partons
	for( int i=0; i<4; i++ )
	{
		p_set[0][i] = 0;
		p_set[1][i] = 0;
		p_set[2][i] = pl1_internal[i];
		p_set[3][i] = pl2_internal[i];
		p_set[4][i] = pl3_internal[i];
		p_set[5][i] = pl4_internal[i];
	}
	
	PDFx1 = ( (pl1_internal[0]+pl2_internal[0]+pl3_internal[0]+pl4_internal[0])+(pl1_internal[3]+pl2_internal[3]+pl3_internal[3]+pl4_internal[3]) )/Sqrt_s;
	PDFx2 = ( (pl1_internal[0]+pl2_internal[0]+pl3_internal[0]+pl4_internal[0])-(pl1_internal[3]+pl2_internal[3]+pl3_internal[3]+pl4_internal[3]) )/Sqrt_s;
	
	// 0 mass approximation
	p_set[0][0] = 0.5*PDFx1*Sqrt_s;
	p_set[0][1] = 0;
	p_set[0][2] = 0;
	p_set[0][3] = 0.5*PDFx1*Sqrt_s;	// to be updated

	p_set[1][0] = 0.5*PDFx2*Sqrt_s;
	p_set[1][1] = 0;
	p_set[1][2] = 0;
	p_set[1][3] = -0.5*PDFx2*Sqrt_s;	// to be updated
	
	
	/// Calculate values needed for the PDF in the pT=0 frame
	if( Use_PDF_w_pT0 )
	{
		p_set[0][0] = Sqrt_s/2;
		p_set[1][0] = Sqrt_s/2;
// 		p_set[0][3] = sqrt( p_set[0][0]*p_set[0][0] - Proton_mass*Proton_mass );
// 		p_set[1][3] = -p_set[0][3];
// 		/// Using forward-beam approximation
// 		Boost_4p_and_2p_2_pT0( ml1, p_set[2], ml2, p_set[3], ml3, p_set[4], ml4, p_set[5], Proton_mass, p_set[0], Proton_mass, p_set[1] );
		Boost_4p_2_pT0( ml1, p_set[2], ml2, p_set[3], ml3, p_set[4], ml4, p_set[5] );
		
		PDFx1 = ( (p_set[2][0]+p_set[3][0]+p_set[4][0]+p_set[5][0]) + (p_set[2][3]+p_set[3][3]+p_set[4][3]+p_set[5][3]) )/( 2*p_set[0][0] );
		PDFx2 = ( (p_set[2][0]+p_set[3][0]+p_set[4][0]+p_set[5][0]) - (p_set[2][3]+p_set[3][3]+p_set[4][3]+p_set[5][3]) )/( 2*p_set[0][0] );
		/// Setting up partons
		p_set[0][0]*= PDFx1;
		p_set[0][1] = 0;
		p_set[0][2] = 0;
		p_set[0][3] = 0;	// to be updated
		p_set[1][0]*= PDFx2;
		p_set[1][1] = 0;
		p_set[1][2] = 0;
		p_set[1][3] = 0;	// to be updated
		if( Debug_Mode ) { printf( "Coefficients for PDF ( x1, x2 ): ( %.10E, %.10E )\n", PDFx1, PDFx2 ); }
	}
	
	
	/// If flag is true, boost to CM frame if PDF is NOT included.
	if( Boost_To_CM && !Use_PDF_w_pT0 ) 
	{
		Boost2CM( ml1, p_set[2], ml2, p_set[3], ml3, p_set[4], ml4, p_set[5]);
		CollisionE = p_set[2][0]+p_set[3][0]+p_set[4][0]+p_set[5][0];
		p_set[0][0] = CollisionE/2;
		p_set[1][0] = CollisionE/2;
	}
	
	
	Mass_4l = sqrt( (pl1[0]+pl2[0]+pl3[0]+pl4[0])*(pl1[0]+pl2[0]+pl3[0]+pl4[0])
		- (pl1[1]+pl2[1]+pl3[1]+pl4[1])*(pl1[1]+pl2[1]+pl3[1]+pl4[1])
		- (pl1[2]+pl2[2]+pl3[2]+pl4[2])*(pl1[2]+pl2[2]+pl3[2]+pl4[2])
		- (pl1[3]+pl2[3]+pl3[3]+pl4[3])*(pl1[3]+pl2[3]+pl3[3]+pl4[3]) );
	
	
	if( !Use_PDF_w_pT0 )
	{
		buffer = new double;
		(*buffer) = ( ContributionCoeff_d+ContributionCoeff_u+ContributionCoeff_s+ContributionCoeff_c );
		ContributionCoeff_d = (ContributionCoeff_d)/(*buffer);
		ContributionCoeff_u = (ContributionCoeff_u)/(*buffer);
		ContributionCoeff_c = (ContributionCoeff_c)/(*buffer);
		ContributionCoeff_s = (ContributionCoeff_s)/(*buffer);
		delete buffer;
	}
	
	
	if( Debug_Mode )
	{
		printf( "Energy of Parton 1: %.10E\nEnergy of Parton 2: %.10E\n", p_set[0][0], p_set[1][0] );
		printf( "Four-momenta entering ME (E px py px):\n" );
		printf( "%.10E %.10E %.10E %.10E %.10E %.10E %.10E %.10E ", p_set[2][0], p_set[2][1], p_set[2][2], p_set[2][3], p_set[3][0], p_set[3][1], p_set[3][2], p_set[3][3] );
		printf( "%.10E %.10E %.10E %.10E %.10E %.10E %.10E %.10E\n", p_set[4][0], p_set[4][1], p_set[4][2], p_set[4][3], p_set[5][0], p_set[5][1], p_set[5][2], p_set[5][3] );
		printf( "Sum px=%.10E\n", (p_set[2][1]+p_set[3][1]+p_set[4][1]+p_set[5][1]) );
		printf( "Sum py=%.10E\n", (p_set[2][2]+p_set[3][2]+p_set[4][2]+p_set[5][2]) );
		printf( "Sum pz=%.10E\n", (p_set[2][3]+p_set[3][3]+p_set[4][3]+p_set[5][3]) );
		printf( "Sum E=%.10E\n", (p_set[2][0]+p_set[3][0]+p_set[4][0]+p_set[5][0]) );
	}
	
	/// Background is interesting in any case, except for the Signal Runs
	if( Test_Model[0]!='!' ) Run_MEKD_MG_MEs_BKG();
	
	/// Signal ME(s) is(are) chosen here
	if( Test_Models.size() > 0 && Test_Model[0]!='!' )
	{
		Signal_MEs.clear();
		
		for( unsigned int i=0; i<Test_Models.size(); i++)
		{
			if( Test_Models[i]=="-1" || Test_Models[i]=="ZZ" )
				{ Run_MEKD_MG_MEs_BKG(); Signal_ME=Background_ME; }
			if( Test_Models[i]=="0" || Test_Models[i]=="Custom" )
				Run_MEKD_MG_ME_Custom();
			if( Test_Models[i]=="1" || Test_Models[i]=="SMHiggs" || Test_Models[i]=="Higgs" )
				Run_MEKD_MG_ME_SMHiggs();
			if( Test_Models[i]=="2" || Test_Models[i]=="CPoddScalar" || Test_Models[i]=="CP-odd" )
				Run_MEKD_MG_ME_CPoddScalar();
			if( Test_Models[i]=="3" || Test_Models[i]=="Spin1particle" )
		 		Run_MEKD_MG_ME_Spin1();
			if( Test_Models[i]=="4" || Test_Models[i]=="Spin2particle" )
				Run_MEKD_MG_ME_Spin2();
			if( Test_Models[i]=="5" || Test_Models[i]=="CPevenScalar" || Test_Models[i]=="CP-even" )
		 		Run_MEKD_MG_ME_CPevenScalar();
				
			Signal_MEs.push_back( Signal_ME );
		}
	}
	else
	{
		if( Test_Model=="-1" || Test_Model=="ZZ" ||
			Test_Model=="!-1" || Test_Model=="!ZZ" )
			{ Run_MEKD_MG_MEs_BKG(); Signal_ME=Background_ME; }
		if( Test_Model=="0" || Test_Model=="Custom" ||
			Test_Model=="!0" || Test_Model=="!Custom" )
			Run_MEKD_MG_ME_Custom();
		if( Test_Model=="1" || Test_Model=="SMHiggs" || Test_Model=="Higgs" ||
			Test_Model=="!1" || Test_Model=="!SMHiggs" || Test_Model=="!Higgs" )
			Run_MEKD_MG_ME_SMHiggs();
		if( Test_Model=="2" || Test_Model=="CPoddScalar" || Test_Model=="CP-odd" ||
			Test_Model=="!2" || Test_Model=="!CPoddScalar" || Test_Model=="!CP-odd" )
			Run_MEKD_MG_ME_CPoddScalar();
		if( Test_Model=="3" || Test_Model=="Spin1particle" ||
			Test_Model=="!3" || Test_Model=="!Spin1particle" )
	 		Run_MEKD_MG_ME_Spin1();
		if( Test_Model=="4" || Test_Model=="Spin2particle" ||
			Test_Model=="!4" || Test_Model=="!Spin2particle" )
			Run_MEKD_MG_ME_Spin2();
		if( Test_Model=="5" || Test_Model=="CPevenScalar" || Test_Model=="CP-even" ||
			Test_Model=="!5" || Test_Model=="!CPevenScalar" || Test_Model=="!CP-even" )
	 		Run_MEKD_MG_ME_CPevenScalar();
	}
	
	
	if( Test_Model[0]!='!' ) KD = log( Signal_ME/Background_ME );
	
	return 0;
}



int MEKD_MG::Run_MEKD_MG(string Input_Model)
{
	buffer_string = Test_Model;
	Test_Model = "!";
	Test_Model += Input_Model;
	
	error_value = Run_MEKD_MG();
	
	Test_Model = buffer_string;
	return error_value;
}



///#include "MEKD_MG_RunMEs.cpp"
///////////////////////////////////
/// INCLUDED MEKD_MG_RunMEs.cpp ///
/// code follows below          ///
///                             ///
///  Part responsible for MEs   ///
///  Calculation                ///
///////////////////////////////////



int MEKD_MG::Run_MEKD_MG_ME_Custom()
{
	if( (error_value=Run_MEKD_MG_MEs_SIG_Spin0())!=0 ) return error_value;
	buffer_Custom = Signal_ME;
	if( (error_value=Run_MEKD_MG_MEs_SIG_Spin2())!=0 ) return error_value;
	Signal_ME += buffer_Custom;
	
	return 0;
}



int MEKD_MG::Run_MEKD_MG_ME_SMHiggs()
{
	if( Use_mh_eq_m4l )
	{
		Set_Of_Model_Parameters.set_block_entry( "mass", 9000006, Mass_4l );
		
		if( Use_Higgs_width ) Set_Of_Model_Parameters.set_block_entry( "decay", 9000006, static_cast<double>( MEKD_CalcHEP_Extra::Higgs_width(Mass_4l) ) );
		else Set_Of_Model_Parameters.set_block_entry( "decay", 9000006, 1 );
		
		if( Vary_signal_couplings )
		{
			Set_Of_Model_Parameters.set_block_entry( "heff", 2, 4*LmbdGG(Mass_4l) );
			Set_Of_Model_Parameters.set_block_entry( "heff", 5, hZZ_coupling );
		}
	}
	else
	{
		Set_Of_Model_Parameters.set_block_entry( "mass", 9000006, Higgs_mass );
		
		if( Use_Higgs_width ) Set_Of_Model_Parameters.set_block_entry( "decay", 9000006, Higgs_width );
		else Set_Of_Model_Parameters.set_block_entry( "decay", 9000006, 1 );
		
		if( Vary_signal_couplings )
		{
			Set_Of_Model_Parameters.set_block_entry( "heff", 2, 4*LmbdGG(Higgs_mass) );
			Set_Of_Model_Parameters.set_block_entry( "heff", 5, hZZ_coupling );
		}
	}
	
	Set_Of_Model_Parameters.set_block_entry( "heff", 1, 0 );
	Set_Of_Model_Parameters.set_block_entry( "heff", 3, 0 );
	Set_Of_Model_Parameters.set_block_entry( "heff", 4, 0 );
	Set_Of_Model_Parameters.set_block_entry( "heff", 6, 0 );
	Set_Of_Model_Parameters.set_block_entry( "heff", 7, 0 );
	Set_Of_Model_Parameters.set_block_entry( "heff", 8, 0 );
	
	return Run_MEKD_MG_MEs_SIG_Spin0();
}



int MEKD_MG::Run_MEKD_MG_ME_CPoddScalar()
{
	if( Use_mh_eq_m4l )
	{
		Set_Of_Model_Parameters.set_block_entry( "mass", 9000006, Mass_4l );
		
		if( Use_Higgs_width ) Set_Of_Model_Parameters.set_block_entry( "decay", 9000006, static_cast<double>( MEKD_CalcHEP_Extra::Higgs_width(Mass_4l) ) );
		else Set_Of_Model_Parameters.set_block_entry( "decay", 9000006, 1 );
		
		if( Vary_signal_couplings )
		{
			Set_Of_Model_Parameters.set_block_entry( "heff", 4, 4*LmbdGG(Mass_4l)*sqrt(3) );
			Set_Of_Model_Parameters.set_block_entry( "heff", 8, 2*hZZ_coupling/m_Z_temp/m_Z_temp*sqrt(3) );
		}
	}
	else
	{
		Set_Of_Model_Parameters.set_block_entry( "mass", 9000006, Higgs_mass );
		
		if( Use_Higgs_width ) Set_Of_Model_Parameters.set_block_entry( "decay", 9000006, Higgs_width );
		else Set_Of_Model_Parameters.set_block_entry( "decay", 9000006, 1 );
		
		if( Vary_signal_couplings )
		{
			Set_Of_Model_Parameters.set_block_entry( "heff", 4, 4*LmbdGG(Higgs_mass)*sqrt(3) );
			Set_Of_Model_Parameters.set_block_entry( "heff", 8, 2*hZZ_coupling/m_Z_temp/m_Z_temp*sqrt(3) );
		}
	}
	
	Set_Of_Model_Parameters.set_block_entry( "heff", 1, 0 );
	Set_Of_Model_Parameters.set_block_entry( "heff", 2, 0 );
	Set_Of_Model_Parameters.set_block_entry( "heff", 3, 0 );
	Set_Of_Model_Parameters.set_block_entry( "heff", 5, 0 );
	Set_Of_Model_Parameters.set_block_entry( "heff", 6, 0 );
	Set_Of_Model_Parameters.set_block_entry( "heff", 7, 0 );
	
	return Run_MEKD_MG_MEs_SIG_Spin0();
}



int MEKD_MG::Run_MEKD_MG_ME_CPevenScalar()
{
	if( Use_mh_eq_m4l )
	{
		Set_Of_Model_Parameters.set_block_entry( "mass", 9000006, Mass_4l );
		
		if( Use_Higgs_width ) Set_Of_Model_Parameters.set_block_entry( "decay", 9000006, static_cast<double>( MEKD_CalcHEP_Extra::Higgs_width(Mass_4l) ) );
		else Set_Of_Model_Parameters.set_block_entry( "decay", 9000006, 1 );
		
		if( Vary_signal_couplings )
		{
			Set_Of_Model_Parameters.set_block_entry( "heff", 1, 4*LmbdGG(Mass_4l)*Mass_4l*Mass_4l );
			Set_Of_Model_Parameters.set_block_entry( "heff", 2, 4*LmbdGG(Mass_4l) );
			Set_Of_Model_Parameters.set_block_entry( "heff", 3, 4*LmbdGG(Mass_4l)/Mass_4l/Mass_4l/Mass_4l/Mass_4l );
			Set_Of_Model_Parameters.set_block_entry( "heff", 5, hZZ_coupling );
			Set_Of_Model_Parameters.set_block_entry( "heff", 6, 2*hZZ_coupling/m_Z_temp/m_Z_temp );
			Set_Of_Model_Parameters.set_block_entry( "heff", 7, hZZ_coupling/m_Z_temp/m_Z_temp/m_Z_temp/m_Z_temp);
		}
	}
	else
	{
		Set_Of_Model_Parameters.set_block_entry( "mass", 9000006, Higgs_mass );
		
		if( Use_Higgs_width ) Set_Of_Model_Parameters.set_block_entry( "decay", 9000006, Higgs_width );
		else Set_Of_Model_Parameters.set_block_entry( "decay", 9000006, 1 );
		
		if( Vary_signal_couplings )
		{
			Set_Of_Model_Parameters.set_block_entry( "heff", 1, 4*LmbdGG(Higgs_mass)*Higgs_mass*Higgs_mass );
			Set_Of_Model_Parameters.set_block_entry( "heff", 2, 4*LmbdGG(Higgs_mass) );
			Set_Of_Model_Parameters.set_block_entry( "heff", 3, 4*LmbdGG(Higgs_mass)/Higgs_mass/Higgs_mass );
			Set_Of_Model_Parameters.set_block_entry( "heff", 5, hZZ_coupling );
			Set_Of_Model_Parameters.set_block_entry( "heff", 6, 2*hZZ_coupling/m_Z_temp/m_Z_temp );
			Set_Of_Model_Parameters.set_block_entry( "heff", 7, hZZ_coupling/m_Z_temp/m_Z_temp/m_Z_temp/m_Z_temp );
		}
	}
	
	Set_Of_Model_Parameters.set_block_entry( "heff", 4, 0 );
	Set_Of_Model_Parameters.set_block_entry( "heff", 8, 0 );
	
	return Run_MEKD_MG_MEs_SIG_Spin0();
}



int MEKD_MG::Run_MEKD_MG_ME_Spin1()
{
	if( Use_mh_eq_m4l )
	{
		Set_Of_Model_Parameters.set_block_entry( "mass", 300, Mass_4l );
		
		if( Use_Higgs_width ) Set_Of_Model_Parameters.set_block_entry( "decay", 300, static_cast<double>( MEKD_CalcHEP_Extra::Higgs_width(Mass_4l) ) );
		else Set_Of_Model_Parameters.set_block_entry( "decay", 300, 1 );
		
		if( Vary_signal_couplings )
		{
			Set_Of_Model_Parameters.set_block_entry( "chern", 1, Mass_4l*Mass_4l*Mass_4l*LmbdGG(Mass_4l) );
			Set_Of_Model_Parameters.set_block_entry( "chern", 3, Mass_4l*Mass_4l*0.5*hZZ_coupling );
		}
	}
	else
	{
		Set_Of_Model_Parameters.set_block_entry( "mass", 300, Higgs_mass );
		
		if( Use_Higgs_width ) Set_Of_Model_Parameters.set_block_entry( "decay", 300, Higgs_width );
		else Set_Of_Model_Parameters.set_block_entry( "decay", 300, 1 );
		
		if( Vary_signal_couplings )
		{
			Set_Of_Model_Parameters.set_block_entry( "chern", 1, Higgs_mass*Higgs_mass*Higgs_mass*LmbdGG(Higgs_mass) );
			Set_Of_Model_Parameters.set_block_entry( "chern", 3, Higgs_mass*Higgs_mass*0.5*hZZ_coupling );
		}
	}
	
	Set_Of_Model_Parameters.set_block_entry( "chern", 2, 0 );
	Set_Of_Model_Parameters.set_block_entry( "chern", 4, 0 );
	
	return Run_MEKD_MG_MEs_SIG_Spin1();
}



int MEKD_MG::Run_MEKD_MG_ME_Spin2()
{
	if( Use_mh_eq_m4l )
	{
		Set_Of_Model_Parameters.set_block_entry( "mass", 9000007, Mass_4l );
		
		if( Use_Higgs_width ) Set_Of_Model_Parameters.set_block_entry( "decay", 9000007, static_cast<double>( MEKD_CalcHEP_Extra::Higgs_width(Mass_4l) ) );
		else Set_Of_Model_Parameters.set_block_entry( "decay", 9000007, 1 );
		
		if( Vary_signal_couplings )
		{
			Set_Of_Model_Parameters.set_block_entry( "gravity", 1, 4*LmbdGG(Mass_4l) );
		}
	}
	else
	{
		Set_Of_Model_Parameters.set_block_entry( "mass", 9000007, Higgs_mass );
		
		if( Use_Higgs_width ) Set_Of_Model_Parameters.set_block_entry( "decay", 9000007, Higgs_width );
		else Set_Of_Model_Parameters.set_block_entry( "decay", 9000007, 1 );
		
		if( Vary_signal_couplings )
		{
			Set_Of_Model_Parameters.set_block_entry( "gravity", 1, 4*LmbdGG(Higgs_mass) );
		}
	}
	
	Set_Of_Model_Parameters.set_block_entry( "gravity", 2, 0 );
	Set_Of_Model_Parameters.set_block_entry( "gravity", 3, 0 );
	Set_Of_Model_Parameters.set_block_entry( "gravity", 4, 0 );
	Set_Of_Model_Parameters.set_block_entry( "gravity", 5, 0 );
	Set_Of_Model_Parameters.set_block_entry( "gravity", 6, 0 );
	Set_Of_Model_Parameters.set_block_entry( "gravity", 7, 0 );
	Set_Of_Model_Parameters.set_block_entry( "gravity", 8, 0 );
	Set_Of_Model_Parameters.set_block_entry( "gravity", 9, 0 );
	Set_Of_Model_Parameters.set_block_entry( "gravity", 10, 0 );
	Set_Of_Model_Parameters.set_block_entry( "gravity", 11, -hZZ_coupling/sqrt(2)/m_Z_temp/m_Z_temp );
	Set_Of_Model_Parameters.set_block_entry( "gravity", 12, 0 );
	Set_Of_Model_Parameters.set_block_entry( "gravity", 13, 0 );
	Set_Of_Model_Parameters.set_block_entry( "gravity", 14, 0 );
	Set_Of_Model_Parameters.set_block_entry( "gravity", 15, hZZ_coupling/sqrt(2) );
	Set_Of_Model_Parameters.set_block_entry( "gravity", 16, 0 );
	Set_Of_Model_Parameters.set_block_entry( "gravity", 17, 0 );
	Set_Of_Model_Parameters.set_block_entry( "gravity", 18, 0 );
	Set_Of_Model_Parameters.set_block_entry( "gravity", 19, 0 );
	Set_Of_Model_Parameters.set_block_entry( "gravity", 20, 0 );
	
	return Run_MEKD_MG_MEs_SIG_Spin2();
}



int MEKD_MG::Run_MEKD_MG_MEs_BKG()
{
	/// Same-flavor block. Electrons
	if( Final_state=="4e" )
	{
		/// Common mass for the same-flavor leptons
		Set_Of_Model_Parameters.set_block_entry( "mass", 13, m_e_temp );
		
		return Run_MEKD_MG_MEs_BKG_SF();
	}
	
	if( Final_state=="2e2m" || Final_state=="2e2mu" )
	{
		/// Common mass for the opposite-flavor leptons
		Set_Of_Model_Parameters.set_block_entry( "mass", 11, m_e_temp );
		Set_Of_Model_Parameters.set_block_entry( "mass", 13, m_mu_temp );
		
		return Run_MEKD_MG_MEs_BKG_OF();
	}
	
	if( Final_state=="4m" || Final_state=="4mu" )
	{
		/// Common mass for the same-flavor leptons
		Set_Of_Model_Parameters.set_block_entry( "mass", 13, m_mu_temp );
		
		return Run_MEKD_MG_MEs_BKG_SF();
	}
	
	return 1;
}



int MEKD_MG::Run_MEKD_MG_MEs_SIG_Spin0()
{
	/// Same-flavor block. Electrons
	if( Final_state=="4e" )
	{
		/// Common mass for the same-flavor leptons
		Set_Of_Model_Parameters.set_block_entry( "mass", 13, m_e_temp );
		
		return Run_MEKD_MG_MEs_SIG_Spin0_SF();
	}
	
	if( Final_state=="2e2m" || Final_state=="2e2mu" )
	{
		/// Common mass for the opposite-flavor leptons
		Set_Of_Model_Parameters.set_block_entry( "mass", 11, m_e_temp );
		Set_Of_Model_Parameters.set_block_entry( "mass", 13, m_mu_temp );
		
		return Run_MEKD_MG_MEs_SIG_Spin0_OF();
	}
	
	if( Final_state=="4m" || Final_state=="4mu" )
	{
		/// Common mass for the same-flavor leptons
		Set_Of_Model_Parameters.set_block_entry( "mass", 13, m_mu_temp );
		
		return Run_MEKD_MG_MEs_SIG_Spin0_SF();
	}
	
	return 1;
}



int MEKD_MG::Run_MEKD_MG_MEs_SIG_Spin1()
{
	/// Same-flavor block. Electrons
	if( Final_state=="4e" )
	{
		/// Common mass for the same-flavor leptons
		Set_Of_Model_Parameters.set_block_entry( "mass", 13, m_e_temp );
		
		return Run_MEKD_MG_MEs_SIG_Spin1_SF();
	}
	
	if( Final_state=="2e2m" || Final_state=="2e2mu" )
	{
		/// Common mass for the opposite-flavor leptons
		Set_Of_Model_Parameters.set_block_entry( "mass", 11, m_e_temp );
		Set_Of_Model_Parameters.set_block_entry( "mass", 13, m_mu_temp );
		
		return Run_MEKD_MG_MEs_SIG_Spin1_OF();
	}
	
	if( Final_state=="4m" || Final_state=="4mu" )
	{
		/// Common mass for the same-flavor leptons
		Set_Of_Model_Parameters.set_block_entry( "mass", 13, m_mu_temp );
		
		return Run_MEKD_MG_MEs_SIG_Spin1_SF();
	}
	
	return 1;
}



int MEKD_MG::Run_MEKD_MG_MEs_SIG_Spin2()
{
	/// Same-flavor block. Electrons
	if( Final_state=="4e" )
	{
		/// Common mass for the same-flavor leptons
		Set_Of_Model_Parameters.set_block_entry( "mass", 13, m_e_temp );
		
		return Run_MEKD_MG_MEs_SIG_Spin2_SF();
	}
	
	if( Final_state=="2e2m" || Final_state=="2e2mu" )
	{
		/// Common mass for the opposite-flavor leptons
		Set_Of_Model_Parameters.set_block_entry( "mass", 11, m_e_temp );
		Set_Of_Model_Parameters.set_block_entry( "mass", 13, m_mu_temp );
		
		return Run_MEKD_MG_MEs_SIG_Spin2_OF();
	}
	
	if( Final_state=="4m" || Final_state=="4mu" )
	{
		/// Common mass for the same-flavor leptons
		Set_Of_Model_Parameters.set_block_entry( "mass", 13, m_mu_temp );
		
		return Run_MEKD_MG_MEs_SIG_Spin2_SF();
	}
	
	return 1;
}



int MEKD_MG::Run_MEKD_MG_MEs_BKG_SF()
{
	/// Down quark block
	Set_Of_Model_Parameters.set_block_entry( "mass", 3, m_d_temp );
	ME_Background_DownType_SF.updateProc( Set_Of_Model_Parameters );
	
	p_set[0][3] = sqrt( p_set[0][0]*p_set[0][0] - m_d_temp*m_d_temp );
	p_set[1][3] = -sqrt( p_set[1][0]*p_set[1][0] - m_d_temp*m_d_temp );
	ME_Background_DownType_SF.setMomenta( p_set );
	ME_Background_DownType_SF.sigmaKin();
	buffer = const_cast<double*>( ME_Background_DownType_SF.getMatrixElements() );
	if( Use_PDF_w_pT0 )
	{
		ContributionCoeff_d = pdfreader( 1, PDFx1, Mass_4l )*pdfreader( -1, PDFx2, Mass_4l );
		Background_ME = ContributionCoeff_d*buffer[0];
		ContributionCoeff_d = pdfreader( -1, PDFx1, Mass_4l )*pdfreader( 1, PDFx2, Mass_4l );
		Background_ME += ContributionCoeff_d*buffer[1];
	}
	else Background_ME = ContributionCoeff_d*(buffer[0]+buffer[1]);
	
	
	
	/// Strange quark block
	Set_Of_Model_Parameters.set_block_entry( "mass", 3, m_s_temp );
	ME_Background_DownType_SF.updateProc( Set_Of_Model_Parameters );
	
	p_set[0][3] = sqrt( p_set[0][0]*p_set[0][0] - m_s_temp*m_s_temp );
	p_set[1][3] = -sqrt( p_set[1][0]*p_set[1][0] - m_s_temp*m_s_temp );
	ME_Background_DownType_SF.setMomenta( p_set );
	ME_Background_DownType_SF.sigmaKin();
	buffer = const_cast<double*>( ME_Background_DownType_SF.getMatrixElements() );
	if( Use_PDF_w_pT0 )
	{
		ContributionCoeff_s= pdfreader( 3, PDFx1, Mass_4l )*pdfreader( -3, PDFx2, Mass_4l );
		Background_ME += ContributionCoeff_s*buffer[0];
		ContributionCoeff_s = pdfreader( -3, PDFx1, Mass_4l )*pdfreader( 3, PDFx2, Mass_4l );
		Background_ME += ContributionCoeff_s*buffer[1];
	}
	else Background_ME += ContributionCoeff_s*(buffer[0]+buffer[1]);
	
	
	
	/// Up quark block
	Set_Of_Model_Parameters.set_block_entry( "mass", 4, m_u_temp );
	ME_Background_UpType_SF.updateProc( Set_Of_Model_Parameters );
	
	p_set[0][3] = sqrt( p_set[0][0]*p_set[0][0] - m_u_temp*m_u_temp );
	p_set[1][3] = -sqrt( p_set[1][0]*p_set[1][0] - m_u_temp*m_u_temp );
	ME_Background_UpType_SF.setMomenta( p_set );
	ME_Background_UpType_SF.sigmaKin();
	buffer = const_cast<double*>( ME_Background_UpType_SF.getMatrixElements() );
	if( Use_PDF_w_pT0 )
	{
		ContributionCoeff_u = pdfreader( 2, PDFx1, Mass_4l )*pdfreader( -2, PDFx2, Mass_4l );
		Background_ME += ContributionCoeff_u*buffer[0];
		ContributionCoeff_u = pdfreader( -2, PDFx1, Mass_4l )*pdfreader( 2, PDFx2, Mass_4l );
		Background_ME += ContributionCoeff_u*buffer[1];
	}
	else Background_ME += ContributionCoeff_u*(buffer[0]+buffer[1]);
	
	
	
	/// Charm quark block
	Set_Of_Model_Parameters.set_block_entry( "mass", 4, m_c_temp );
	ME_Background_UpType_SF.updateProc( Set_Of_Model_Parameters );
	
	p_set[0][3] = sqrt( p_set[0][0]*p_set[0][0] - m_c_temp*m_c_temp );
	p_set[1][3] = -sqrt( p_set[1][0]*p_set[1][0] - m_c_temp*m_c_temp );
	ME_Background_UpType_SF.setMomenta( p_set );
	ME_Background_UpType_SF.sigmaKin();
	buffer = const_cast<double*>( ME_Background_UpType_SF.getMatrixElements() );
	if( Use_PDF_w_pT0 )
	{
		ContributionCoeff_c = pdfreader( 4, PDFx1, Mass_4l )*pdfreader( -4, PDFx2, Mass_4l );
		Background_ME += ContributionCoeff_c*buffer[0];
		ContributionCoeff_c = pdfreader( -4, PDFx1, Mass_4l )*pdfreader( 4, PDFx2, Mass_4l );
		Background_ME += ContributionCoeff_c*buffer[1];
	}
	else Background_ME += ContributionCoeff_c*(buffer[0]+buffer[1]);
	
	
	return 0;
}



int MEKD_MG::Run_MEKD_MG_MEs_BKG_OF()
{
	/// Down quark block
	Set_Of_Model_Parameters.set_block_entry( "mass", 3, m_d_temp );
	ME_Background_DownType_OF.updateProc( Set_Of_Model_Parameters );
	
	p_set[0][3] = sqrt( p_set[0][0]*p_set[0][0] - m_d_temp*m_d_temp );
	p_set[1][3] = -sqrt( p_set[1][0]*p_set[1][0] - m_d_temp*m_d_temp );
	ME_Background_DownType_OF.setMomenta( p_set );
	ME_Background_DownType_OF.sigmaKin();
	buffer = const_cast<double*>( ME_Background_DownType_OF.getMatrixElements() );
	if( Use_PDF_w_pT0 )
	{
		ContributionCoeff_d = pdfreader( 1, PDFx1, Mass_4l )*pdfreader( -1, PDFx2, Mass_4l );
		Background_ME = ContributionCoeff_d*buffer[0];
		ContributionCoeff_d = pdfreader( -1, PDFx1, Mass_4l )*pdfreader( 1, PDFx2, Mass_4l );
		Background_ME += ContributionCoeff_d*buffer[1];
	}
	else Background_ME = ContributionCoeff_d*(buffer[0]+buffer[1]);
	
	
	
	/// Strange quark block
	Set_Of_Model_Parameters.set_block_entry( "mass", 3, m_s_temp );
	ME_Background_DownType_OF.updateProc( Set_Of_Model_Parameters );
	
	p_set[0][3] = sqrt( p_set[0][0]*p_set[0][0] - m_s_temp*m_s_temp );
	p_set[1][3] = -sqrt( p_set[1][0]*p_set[1][0] - m_s_temp*m_s_temp );
	ME_Background_DownType_OF.setMomenta( p_set );
	ME_Background_DownType_OF.sigmaKin();
	buffer = const_cast<double*>( ME_Background_DownType_OF.getMatrixElements() );
	if( Use_PDF_w_pT0 )
	{
		ContributionCoeff_s= pdfreader( 3, PDFx1, Mass_4l )*pdfreader( -3, PDFx2, Mass_4l );
		Background_ME += ContributionCoeff_s*buffer[0];
		ContributionCoeff_s = pdfreader( -3, PDFx1, Mass_4l )*pdfreader( 3, PDFx2, Mass_4l );
		Background_ME += ContributionCoeff_s*buffer[1];
	}
	else Background_ME += ContributionCoeff_s*(buffer[0]+buffer[1]);
	
	
	
	/// Up quark block
	Set_Of_Model_Parameters.set_block_entry( "mass", 4, m_u_temp );
	ME_Background_UpType_OF.updateProc( Set_Of_Model_Parameters );
	
	p_set[0][3] = sqrt( p_set[0][0]*p_set[0][0] - m_u_temp*m_u_temp );
	p_set[1][3] = -sqrt( p_set[1][0]*p_set[1][0] - m_u_temp*m_u_temp );
	ME_Background_UpType_OF.setMomenta( p_set );
	ME_Background_UpType_OF.sigmaKin();
	buffer = const_cast<double*>( ME_Background_UpType_OF.getMatrixElements() );
	if( Use_PDF_w_pT0 )
	{
		ContributionCoeff_u = pdfreader( 2, PDFx1, Mass_4l )*pdfreader( -2, PDFx2, Mass_4l );
		Background_ME += ContributionCoeff_u*buffer[0];
		ContributionCoeff_u = pdfreader( -2, PDFx1, Mass_4l )*pdfreader( 2, PDFx2, Mass_4l );
		Background_ME += ContributionCoeff_u*buffer[1];
	}
	else Background_ME += ContributionCoeff_u*(buffer[0]+buffer[1]);
	
	
	
	/// Charm quark block
	Set_Of_Model_Parameters.set_block_entry( "mass", 4, m_c_temp );
	ME_Background_UpType_OF.updateProc( Set_Of_Model_Parameters );
	
	p_set[0][3] = sqrt( p_set[0][0]*p_set[0][0] - m_c_temp*m_c_temp );
	p_set[1][3] = -sqrt( p_set[1][0]*p_set[1][0] - m_c_temp*m_c_temp );
	ME_Background_UpType_OF.setMomenta( p_set );
	ME_Background_UpType_OF.sigmaKin();
	buffer = const_cast<double*>( ME_Background_UpType_OF.getMatrixElements() );
	if( Use_PDF_w_pT0 )
	{
		ContributionCoeff_c = pdfreader( 4, PDFx1, Mass_4l )*pdfreader( -4, PDFx2, Mass_4l );
		Background_ME += ContributionCoeff_c*buffer[0];
		ContributionCoeff_c = pdfreader( -4, PDFx1, Mass_4l )*pdfreader( 4, PDFx2, Mass_4l );
		Background_ME += ContributionCoeff_c*buffer[1];
	}
	else Background_ME += ContributionCoeff_c*(buffer[0]+buffer[1]);
	
	
	return 0;
}



int MEKD_MG::Run_MEKD_MG_MEs_SIG_Spin0_SF()
{
	/// Signal block
	ME_Signal_Spin0_SF.updateProc( Set_Of_Model_Parameters );
	
	p_set[0][3] = p_set[0][0];
	p_set[1][3] = -p_set[1][0];
	ME_Signal_Spin0_SF.setMomenta( p_set );
	ME_Signal_Spin0_SF.sigmaKin();
	buffer = const_cast<double*>( ME_Signal_Spin0_SF.getMatrixElements() );
	if( Use_PDF_w_pT0 ) { Signal_ME = pdfreader( 21, PDFx1, Mass_4l )*pdfreader( 21, PDFx2, Mass_4l )*buffer[0]; }
	else Signal_ME = buffer[0];
	
	
	return 0;
}



int MEKD_MG::Run_MEKD_MG_MEs_SIG_Spin0_OF()
{
	/// Signal block
	ME_Signal_Spin0_OF.updateProc( Set_Of_Model_Parameters );
	
	p_set[0][3] = p_set[0][0];
	p_set[1][3] = -p_set[1][0];
	ME_Signal_Spin0_OF.setMomenta( p_set );
	ME_Signal_Spin0_OF.sigmaKin();
	buffer = const_cast<double*>( ME_Signal_Spin0_OF.getMatrixElements() );
	if( Use_PDF_w_pT0 ) { Signal_ME = pdfreader( 21, PDFx1, Mass_4l )*pdfreader( 21, PDFx2, Mass_4l )*buffer[0]; }
	else Signal_ME = buffer[0];
	
	
	return 0;
}



int MEKD_MG::Run_MEKD_MG_MEs_SIG_Spin1_SF()
{
	/// Signal block
	ME_Signal_Spin1_SF.updateProc( Set_Of_Model_Parameters );
	
	p_set[0][3] = p_set[0][0];
	p_set[1][3] = -p_set[1][0];
	ME_Signal_Spin1_SF.setMomenta( p_set );
	ME_Signal_Spin1_SF.sigmaKin();
	buffer = const_cast<double*>( ME_Signal_Spin1_SF.getMatrixElements() );
	if( Use_PDF_w_pT0 ) { Signal_ME = pdfreader( 21, PDFx1, Mass_4l )*pdfreader( 21, PDFx2, Mass_4l )*buffer[0]; }
	else Signal_ME = buffer[0];
	
	
	return 0;
}



int MEKD_MG::Run_MEKD_MG_MEs_SIG_Spin1_OF()
{
	/// Signal block
	ME_Signal_Spin1_OF.updateProc( Set_Of_Model_Parameters );
	
	p_set[0][3] = p_set[0][0];
	p_set[1][3] = -p_set[1][0];
	ME_Signal_Spin1_OF.setMomenta( p_set );
	ME_Signal_Spin1_OF.sigmaKin();
	buffer = const_cast<double*>( ME_Signal_Spin1_OF.getMatrixElements() );
	if( Use_PDF_w_pT0 ) { Signal_ME = pdfreader( 21, PDFx1, Mass_4l )*pdfreader( 21, PDFx2, Mass_4l )*buffer[0]; }
	else Signal_ME = buffer[0];
	
	
	return 0;
}



int MEKD_MG::Run_MEKD_MG_MEs_SIG_Spin2_SF()
{
	/// Signal block
	ME_Signal_Spin2_SF.updateProc( Set_Of_Model_Parameters );
	
	p_set[0][3] = p_set[0][0];
	p_set[1][3] = -p_set[1][0];
	ME_Signal_Spin2_SF.setMomenta( p_set );
	ME_Signal_Spin2_SF.sigmaKin();
	buffer = const_cast<double*>( ME_Signal_Spin2_SF.getMatrixElements() );
	if( Use_PDF_w_pT0 ) { Signal_ME = pdfreader( 21, PDFx1, Mass_4l )*pdfreader( 21, PDFx2, Mass_4l )*buffer[0]; }
	else Signal_ME = buffer[0];
	
	
	return 0;
}



int MEKD_MG::Run_MEKD_MG_MEs_SIG_Spin2_OF()
{
	/// Signal block
	ME_Signal_Spin2_OF.updateProc( Set_Of_Model_Parameters );
	
	p_set[0][3] = p_set[0][0];
	p_set[1][3] = -p_set[1][0];
	ME_Signal_Spin2_OF.setMomenta( p_set );
	ME_Signal_Spin2_OF.sigmaKin();
	buffer = const_cast<double*>( ME_Signal_Spin2_OF.getMatrixElements() );
	if( Use_PDF_w_pT0 ) { Signal_ME = pdfreader( 21, PDFx1, Mass_4l )*pdfreader( 21, PDFx2, Mass_4l )*buffer[0]; }
	else Signal_ME = buffer[0];
	
	
	return 0;
}



///////////////////////////////////
/// END OF MEKD_MG_RunMEs.cpp   ///
///////////////////////////////////





///#include "MEKD_MG_Sorter.cpp"
///////////////////////////////////
/// INCLUDED MEKD_MG_Sorter.cpp ///
/// code follows below          ///
///                             ///
/// Part responsible for        ///
/// momenta rearrangement       ///
///////////////////////////////////



int MEKD_MG::Arrange_Internal_pls()
{
	id_set[0]=id1; id_set[1]=id2; id_set[2]=id3; id_set[3]=id4;
	sort( id_set.begin(), id_set.end() );
	
	if( id_set[0] == -13 && id_set[1] == -11 && id_set[2] == 11 && id_set[3] == 13 )
	{
		if( id1 == 11 ) pl1_internal = pl1;
		if( id2 == 11 ) pl1_internal = pl2;
		if( id3 == 11 ) pl1_internal = pl3;
		if( id4 == 11 ) pl1_internal = pl4;
		
		if( id1 == -11 ) pl2_internal = pl1;
		if( id2 == -11 ) pl2_internal = pl2;
		if( id3 == -11 ) pl2_internal = pl3;
		if( id4 == -11 ) pl2_internal = pl4;
		
		if( id1 == 13 ) pl3_internal = pl1;
		if( id2 == 13 ) pl3_internal = pl2;
		if( id3 == 13 ) pl3_internal = pl3;
		if( id4 == 13 ) pl3_internal = pl4;
		
		if( id1 == -13 ) pl4_internal = pl1;
		if( id2 == -13 ) pl4_internal = pl2;
		if( id3 == -13 ) pl4_internal = pl3;
		if( id4 == -13 ) pl4_internal = pl4;
		Final_state = "2e2m";
		
		return 0;
	}
	
	if( id_set[0] == -13 && id_set[1] == -13 && id_set[2] == 13 && id_set[3] == 13 )
	{
		buffer_bool = false;	//first muon has beed caught
		if( id1 == 13 && !buffer_bool ) { pl1_internal=pl1; buffer_bool=true; }
		if( id2 == 13 && buffer_bool ) pl3_internal = pl2;
		if( id2 == 13 && !buffer_bool ) { pl1_internal=pl2; buffer_bool=true; }
		if( id3 == 13 && buffer_bool ) pl3_internal = pl3;
		if( id3 == 13 && !buffer_bool ) { pl1_internal=pl3; buffer_bool=true; }
		if( id4 == 13 && buffer_bool ) pl3_internal = pl4;
		
		buffer_bool = false;	//first antimuon has beed caught
		if( id1 == -13 && !buffer_bool ) { pl2_internal=pl1; buffer_bool=true; }
		if( id2 == -13 && buffer_bool ) pl4_internal = pl2;
		if( id2 == -13 && !buffer_bool ) { pl2_internal=pl2; buffer_bool=true; }
		if( id3 == -13 && buffer_bool ) pl4_internal = pl3;
		if( id3 == -13 && !buffer_bool ) { pl2_internal=pl3; buffer_bool=true; }
		if( id4 == -13 && buffer_bool ) pl4_internal = pl4;
		Final_state = "4mu";
		
		return 0;
	}
	
	if( id_set[0] == -11 && id_set[1] == -11 && id_set[2] == 11 && id_set[3] == 11 )
	{
		buffer_bool = false;	//first electron has beed caught
		if( id1 == 11 && !buffer_bool ) { pl1_internal=pl1; buffer_bool=true; }
		if( id2 == 11 && buffer_bool ) pl3_internal = pl2;
		if( id2 == 11 && !buffer_bool ) { pl1_internal=pl2; buffer_bool=true; }
		if( id3 == 11 && buffer_bool ) pl3_internal = pl3;
		if( id3 == 11 && !buffer_bool ) { pl1_internal=pl3; buffer_bool=true; }
		if( id4 == 11 && buffer_bool ) pl3_internal = pl4;
		
		buffer_bool = false;	//first positron has beed caught
		if( id1 == -11 && !buffer_bool ) { pl2_internal=pl1; buffer_bool=true; }
		if( id2 == -11 && buffer_bool ) pl4_internal = pl2;
		if( id2 == -11 && !buffer_bool ) { pl2_internal=pl2; buffer_bool=true; }
		if( id3 == -11 && buffer_bool ) pl4_internal = pl3;
		if( id3 == -11 && !buffer_bool ) { pl2_internal=pl3; buffer_bool=true; }
		if( id4 == -11 && buffer_bool ) pl4_internal = pl4;
		Final_state = "4e";
		
		return 0;
	}
	
	if( id_set[0] == 0 && id_set[1] == 0 && id_set[2] == 0 && id_set[3] == 0 )
	{
		cout << "Warning. Particle ids are not set. Assuming a proper input-lepton configuration.\n";
		pl1_internal=pl1; pl2_internal=pl2; pl3_internal=pl3; pl4_internal=pl4;
		
		return 0;
	}
	
	return 1;
}


///////////////////////////////////
/// END OF MEKD_MG_Sorter.cpp   ///
///////////////////////////////////



#endif