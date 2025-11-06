#ifndef CONSTANTS_H
#define CONSTANTS_H


/*
The "Constants.h" header file, when included, provides some useful constant objects that will be used throughout the code.
*/


#include <complex>
using namespace std;

// The complex number unit
const complex<double> I = complex<double> (0, 1);

// Some real constants
const double MN = 938.918695;  // average nucleon mass
const double Mp = 938.918695;
const double MN2 = 881568.315821;//MN^2

const double Mpi = 138.0389867; // average pion mass (2*Mpi_chrgd+Mpi_ntrl)/3
const double Mpi2 = 19054.761849; //Mpi^2

const double muonmass = 105.658369;  // muon mass (MeV)
const double electronmass = 0.511;   // electron mass

const double Pi = 3.141592654;  // \pi
const double hbc = 197.3270;  // hbar*c (MeV*fm)

// // // // // // // // // // // // // // // 
const double x2Pi = 2*Pi;
// // // // // // // // // // // // // // // 

// Electroweak constants
const double M_W = 80385.;  // W-boson mass (MeV)
const double M_Z = 91187.6; // Z-boson mass (MeV)
const double Cabibbo = 0.974;  // cosine Cabibbo angle
const double G_Fermi = 1.16637e-11;  // Fermi constant (MeV^{-2})
const double Alpha = 0.007297353;  // fine structure constant
const double gA = 1.26; // axial charge 

const int down = 0; //spin -1/2
const int up = 1;   //spin 1/2

const int proton = 1;
const int neutron = 2;



#endif
