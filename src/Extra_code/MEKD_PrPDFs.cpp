#ifndef MEKD_PrPDF_cpp
#define MEKD_PrPDF_cpp

#include <cmath>
#include <string>

using namespace std;


double PDF_SMH( double Sqrt_s, double m_H, double m4l, string Final_state );
double Pick_Greater_Nr( double x, double y);
double RooDoubleCB( double x, double mean, double width, double alpha1, double n1, double alpha2, double n2 );




double PDF_SMH( double Sqrt_s, double m_H, double m4l, string Final_state )
{
	if( Sqrt_s==8 )
	{
		if( Final_state=="4e" ) return RooDoubleCB( m4l,
				(14.8494)+(-0.405085*m_H)+(0.00422887*m_H*m_H)+(-2.16441e-05*m_H*m_H*m_H)+(5.32721e-08*m_H*m_H*m_H*m_H)+(-5.00455e-11*m_H*m_H*m_H*m_H*m_H) + m_H,
				(-15.2189)+(0.442707*m_H)+(-0.00452573*m_H*m_H)+(2.27111e-05*m_H*m_H*m_H)+(-5.39753e-08*m_H*m_H*m_H*m_H)+(4.8518e-11*m_H*m_H*m_H*m_H*m_H),
				(-22.3693)+(0.627422*m_H)+(-0.0064991*m_H*m_H)+(3.20674e-05*m_H*m_H*m_H)+(-7.51115e-08*m_H*m_H*m_H*m_H)+(6.71038e-11*m_H*m_H*m_H*m_H*m_H),
				Pick_Greater_Nr( (-197.954)+(4.6367*m_H)+(-0.0391504*m_H*m_H)+(0.000155755*m_H*m_H*m_H)+(-2.97342e-07*m_H*m_H*m_H*m_H)+(2.20176e-10*m_H*m_H*m_H*m_H*m_H), 1 ),
				(-33.5339)+(0.956795*m_H)+(-0.0100707*m_H*m_H)+(4.99974e-05*m_H*m_H*m_H)+(-1.15201e-07*m_H*m_H*m_H*m_H)+(9.94549e-11*m_H*m_H*m_H*m_H*m_H),
				20 );
		if( Final_state=="2e2m" || Final_state=="2e2mu" ) return RooDoubleCB( m4l,
				(-5.56876)+(0.155983*m_H)+(-0.00167517*m_H*m_H)+(8.37884e-06*m_H*m_H*m_H)+(-2.01225e-08*m_H*m_H*m_H*m_H)+(1.89793e-11*m_H*m_H*m_H*m_H*m_H) + m_H,
				(-23.0489)+(0.624282*m_H)+(-0.00621035*m_H*m_H)+(2.9971e-05*m_H*m_H*m_H)+(-6.8571e-08*m_H*m_H*m_H*m_H)+(5.95102e-11*m_H*m_H*m_H*m_H*m_H),
				(-14.5509)+(0.395116*m_H)+(-0.00385822*m_H*m_H)+(1.78522e-05*m_H*m_H*m_H)+(-3.88799e-08*m_H*m_H*m_H*m_H)+(3.20036e-11*m_H*m_H*m_H*m_H*m_H),
				Pick_Greater_Nr( (4.71546)+(-0.174103*m_H)+(0.00338962*m_H*m_H)+(-2.3452e-05*m_H*m_H*m_H)+(6.66662e-08*m_H*m_H*m_H*m_H)+(-6.68915e-11*m_H*m_H*m_H*m_H*m_H), 1 ),
				(-6.15641)+(0.244828*m_H)+(-0.00311039*m_H*m_H)+(1.8662e-05*m_H*m_H*m_H)+(-5.19696e-08*m_H*m_H*m_H*m_H)+(5.39877e-11*m_H*m_H*m_H*m_H*m_H),
				20 );
		if( Final_state=="4m" || Final_state=="4mu" ) return RooDoubleCB( m4l,
				(-5.71849)+(0.145626*m_H)+(-0.00138862*m_H*m_H)+(6.03825e-06*m_H*m_H*m_H)+(-1.19684e-08*m_H*m_H*m_H*m_H)+(8.75281e-12*m_H*m_H*m_H*m_H*m_H) + m_H,
				(-4.56178)+(0.123209*m_H)+(-0.00107193*m_H*m_H)+(4.5413e-06*m_H*m_H*m_H)+(-8.19429e-09*m_H*m_H*m_H*m_H)+(4.75955e-12*m_H*m_H*m_H*m_H*m_H),
				(-6.83599)+(0.196814*m_H)+(-0.00179054*m_H*m_H)+(7.59602e-06*m_H*m_H*m_H)+(-1.48893e-08*m_H*m_H*m_H*m_H)+(1.07477e-11*m_H*m_H*m_H*m_H*m_H),
				Pick_Greater_Nr( (24.2026)+(-0.5443*m_H)+(0.00517101*m_H*m_H)+(-2.3485e-05*m_H*m_H*m_H)+(5.07143e-08*m_H*m_H*m_H*m_H)+(-4.18694e-11*m_H*m_H*m_H*m_H*m_H), 1 ),
				(-26.1111)+(0.767139*m_H)+(-0.00830412*m_H*m_H)+(4.35986e-05*m_H*m_H*m_H)+(-1.10717e-07*m_H*m_H*m_H*m_H)+(1.09256e-10*m_H*m_H*m_H*m_H*m_H),
				20 );
	}
	
	return -1;
}



double Pick_Greater_Nr( double x, double y)
{ return x >= y ? x : y; }



double RooDoubleCB( double x, double mean, double width, double alpha1, double n1, double alpha2, double n2 )
{
	double t = (x-mean)/width;
	
	if( t>-alpha1 && t<alpha2 )
	{
		return exp( -0.5*t*t );
	}
	else if( t<-alpha1 )
	{
		double A1 = pow( n1/fabs(alpha1), n1 )*exp( -alpha1*alpha1/2 );
		double B1 = n1/fabs( alpha1 )-fabs( alpha1 );
		
		return A1*pow( B1-t, -n1 );
	}
	else if( t>alpha2 )
	{
		double A2 = pow( n2/fabs(alpha2), n2 )*exp( -alpha2*alpha2/2 );
		double B2 = n2/fabs(alpha2)-fabs( alpha2 );
		
		return A2*pow( B2+t,-n2 );
	}
	
//	cout << "ERROR evaluating range..." << endl;
	return 99;
}

#endif