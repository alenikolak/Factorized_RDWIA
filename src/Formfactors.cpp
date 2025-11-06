#include "Constants.h"

//Constants for Kelly parametrization
const double rmup=2.79285;
const double rmun=-1.91304;

const double As=1.70;
const double Bs=3.30;

const double a0=1.0;
const double a1Gep=-0.24;
const double b1Gep=10.98;
const double b2Gep=12.82;
const double b3Gep=21.97;

const double a1Gmp=0.12;
const double b1Gmp=10.97;
const double b2Gmp=18.86;
const double b3Gmp=6.55;

const double a1Gmn=2.33;
const double b1Gmn=14.72;
const double b2Gmn=24.20;
const double b3Gmn=84.1;

const double MV2 = 0.710649*1E6; //MV=0.843 GeV -> in MeV^2

//Constants for axial dipole FF:
const double MA2 = (1.0*1.0)*1E6; //MeV2


void Fp_kelly(double Q2, double &F1, double &F2)
{
    double tau = Q2/(4.*MN2);

    double GEp = (a0 + a1Gep*tau)/(1. + b1Gep*tau + b2Gep*pow(tau,2) + b3Gep*pow(tau,3));
    double GMp = rmup*(a0 + a1Gmp*tau)/(1. + b1Gmp*tau + b2Gmp*pow(tau,2) + b3Gmp*pow(tau,3));

    
    F1 = (GEp+tau*GMp)/(1.+tau);
    F2 = (GMp-GEp)/(1.+tau);

}


void Fn_kelly(double Q2, double &F1, double &F2)
{

    	double tau = Q2/(4.*MN2);

    	double DipV = 1./pow( 1+Q2/MV2, 2 ); 
	double GEn = As*tau/(1.+Bs*tau)*DipV;
	double GMn = rmun*(a0 + a1Gmn*tau)/(1. + b1Gmn*tau + b2Gmn*pow(tau,2) + b3Gmn*pow(tau,3));

    	F1 = (GEn+tau*GMn)/(1.+tau);
    	F2 = (GMn-GEn)/(1.+tau);
	
}


void FIV_kelly(double Q2, double &F1, double &F2)
{
	//isovector form factors
	double F1p, F2p, F1n, F2n;
	Fp_kelly(Q2,F1p,F2p);
	Fn_kelly(Q2,F1n,F2n);
	
	F1 = F1p - F1n;
	F2 = F2p - F2n;

}

void F1_only(double Q2, double &F1, double &F2)
{
	Fp_kelly(Q2, F1, F2);
	F2=0;
}

void F2_only(double Q2, double &F1, double &F2)
{
	Fp_kelly(Q2, F1, F2);
	F1=0;
}

void F2_one(double Q2, double &F1, double &F2)
{
	F1=0;
	F2=1.;
}


void FA_dipole(double Q2, double &GA, double &GP)
{

    double DipA = 1./pow( 1+Q2/MA2, 2 ); 

    GA = gA * DipA ;
    GP = GA*4*MN2/(Q2 + Mpi2);
}

void FA_only(double Q2, double &GA, double &GP)
{
	FA_dipole(Q2, GA, GP);
	GP=0.;

}


void GP_only(double Q2, double &GA, double &GP)
{
	FA_dipole(Q2, GA, GP);
	GA=0.;
}

void Fzero(double Q2, double &f, double &g)
{
	f=0.;
	g=0.;
}
