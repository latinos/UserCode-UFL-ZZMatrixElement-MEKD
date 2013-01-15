#ifndef MEKD_MG_h
#define MEKD_MG_h

#include <string>
#include <vector>

#include "../src/MadGraphSrc/read_slha.h"

using namespace std;


class MEKD_MG
{
public:
	/// Flags
	bool Boost_To_CM;	// for a boosted data
	bool Debug_Mode;	// Enable debugging mode
// 	bool Force_g3_running;	// unused. At some point was included for alpha_QCD
	bool Overwrite_e_and_mu_masses;	// switch for manual m_e, m_mu masses
	bool Use_Higgs_width;	//	if false, width is fixed to =1
	bool Use_mh_eq_m4l;	// Set mh to m4l for every event
	bool Use_PDF_w_pT0;	// Use PDFs in the pT=0 frame. If true, Boost_To_CM is ignored
	bool Vary_signal_couplings;	// Allow couplings to change with mass
	
	/// General parameters
	double ContributionCoeff_d;	//42	/// the value has no effect if PDF is used but the variable is always used
	double ContributionCoeff_u;	//217
	double ContributionCoeff_s;	//5
	double ContributionCoeff_c;	//3
// 	double GG;	// Assign QCD coupling, force g3 running if needed
	double Sqrt_s;	//Max energy, collision energy
	
	/// Physical parameters
	double Electron_mass;	//0.0005109989, for enabled overwriting
	double Higgs_mass;	// Works only if Use_mh_eq_m4l=false
	double Higgs_width;	// Practically not used, for future implementations
	double Muon_mass;	//0.10565837, for enabled overwriting
	double Proton_mass;	// Always used if needed
	
	/// Final-state lepton/photon information
	double *p1, *p2, *p3, *p4, *p5;
	double id1, id2, id3, id4, id5;
	
	/// String flags and file locations
	string Final_state;	// Final state, for the moment: 4e, 4mu, 2e2mu
	string Test_Model;	// -1 or ZZ; 0 or Custom; 1 or SMHiggs; 2 or CPoddScalar; 3 or CPevenScalar; 4 or Spin2particle
	vector<string> Test_Models;	// same names as for the Test_Model
	string Parameter_file;	// Location where a parameter card is stored
	string PDF_file;	// PDF/PDT table file
	
	/// Calculation results
	double Mass_4l;	//is filled after running RUN_XXXX(...)
	double Background_ME;	//may not be used if running RUN_MEKD_MG( string ) is chosen
	double Signal_ME;	//is filled after running RUN_XXXX(...)
	vector<double> Signal_MEs;	//is filled if Test_Models are set after running RUN_XXXX(...)
	double KD;	//is not filled with RUN_MEKD_MG( string )
	
	/// Functions
	void Set_Default_MEKD_MG_Parameters();
	
	int Reload_Parameters();	// reloads parameter set and updates PDF file reader
	int Run_MEKD_MG();	// main routine to evaluate matrix elements; updates "Calculation results"
	int Run_MEKD_MG(string Input_Model);	// Calculates a ME ONLY for a chosen model; ignores automatic background calculation. Updates Signal_ME
	
	/// Constructors, destructors
	MEKD_MG();
	~MEKD_MG();
	
private:
	bool Parameters_Are_Loaded, buffer_bool;
	
	int error_value;
	
	double v_expectation;	// Vacuum expectation value
	double hZZ_coupling;
	double *buffer, buffer_Custom, ml1, ml2, ml3, ml4, PDFx1, PDFx2, m_d_temp, m_u_temp, m_s_temp, m_c_temp, m_e_temp, m_mu_temp, m_Z_temp;
	double *pl1_internal, *pl2_internal, *pl3_internal, *pl4_internal, *pA1_internal;
	
	string buffer_string;
	
	vector<double> id_set;
	vector<double*> p_set;
	
	SLHAReader Set_Of_Model_Parameters;
	
	/// Internal functions ///
	int Load_Parameters();
	
	int Arrange_Internal_pls();
	
	/// Sets up particular choices
	int Run_MEKD_MG_ME_Custom();
	int Run_MEKD_MG_ME_CPevenScalar();
	int Run_MEKD_MG_ME_CPoddScalar();
	int Run_MEKD_MG_ME_SMHiggs();
	int Run_MEKD_MG_ME_Spin0PH();
	int Run_MEKD_MG_ME_Spin1();
	int Run_MEKD_MG_ME_Spin2();
	
	/// Blind-calculation functions
	int Run_MEKD_MG_MEs_BKG();
	int Run_MEKD_MG_MEs_BKG_Sub(string flavor, bool photon);
	int Run_MEKD_MG_MEs_SIG_Spin0();
	int Run_MEKD_MG_MEs_SIG_Spin0_Sub(string flavor, bool photon);
	int Run_MEKD_MG_MEs_SIG_Spin1();
	int Run_MEKD_MG_MEs_SIG_Spin1_Sub(string flavor, bool photon);
	int Run_MEKD_MG_MEs_SIG_Spin2();
	int Run_MEKD_MG_MEs_SIG_Spin2_Sub(string flavor, bool photon);
};


#endif