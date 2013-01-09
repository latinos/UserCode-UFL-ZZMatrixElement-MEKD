#ifndef MEKD_MG_Boosts_h
#define MEKD_MG_Boosts_h

#include <cmath>


void Boost_4p_2_pT0(double mass0, double *pi0, double mass1, double *pi1, double mass2, double *pi2, double mass3, double *pi3);
void Boost_4p_and_2p_2_pT0(double mass0, double *pi0, double mass1, double *pi1, double mass2, double *pi2, double mass3, double *pi3, double mass4, double *pi4, double mass5, double *pi5);
void Boost_5p_2_pT0(double mass0, double *pi0, double mass1, double *pi1, double mass2, double *pi2, double mass3, double *pi3, double mass4, double *pi4);
void Boost_5p_and_2p_2_pT0(double mass0, double *pi0, double mass1, double *pi1, double mass2, double *pi2, double mass3, double *pi3, double mass4, double *pi4, double mass5, double *pi5);
void Boost2CM(double mass0, double *pi0, double mass1, double *pi1, double mass2, double *pi2, double mass3, double *pi3);
void Boost2CM(double mass0, double *pi0, double mass1, double *pi1, double mass2, double *pi2, double mass3, double *pi3, double mass4, double *pi4);
void Boost(double *vector, double *boost);
void Boost_long(long double *vector, long double *boost);



#endif