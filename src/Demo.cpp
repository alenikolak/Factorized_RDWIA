#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <cstdarg>

// // // // // 
#include <chrono>
// // // // //

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Constants.h"
#include "Formfactors.h"
#include "S_tables.h"
#include "Leptonic.h"
#include "Hadronic.h"
#include "Crosssections.h"

using namespace std;

int main(int argc, char* argv[])
{

	//We set form factors, proton form factors from Kelly
	FormFactors FF;
	FF.F_vector = Fp_kelly;
	FF.F_axial = Fzero; //No axial FF

	//Reading the S-state tables
	S_level S_low;
	S_low.read_S("./Overlaps/EDRMF/C12/S_lev1_N1"); //Low energy table with TN steps of 1 MeV TN < 15 MeV

	S_level S_high;
	S_high.read_S("./Overlaps/EDRMF/C12/highTS_lev1_N1"); //Higher energies with TN steps of 5 MeV, TN >= 15 MeV


	//And real also the P-state
	S_level P_low;
	P_low.read_S("./Overlaps/EDRMF/C12/S_lev2_N1");

	S_level P_high;
	P_high.read_S("./Overlaps/EDRMF/C12/highTS_lev2_N1");


	double w = 100.;
	double q = 380.;

	cosdep_Responses(w,q, &S_high, &FF);


	return 0;
	//The reading takes most of the time

	//The functions 'cross_section' in Crosssections.cpp give the exclusive cross section, but they are not so usefull usually

	//Let's instead do the inclusive response
	q=380.;
	for (double w = 0 ; w < q ; w+=1)
	{
		double Responses_S1[5];
		Inclusive_Responses(w, q, &S_low, &FF, Responses_S1); //S state, low enery table

		double Responses_S2[5];
		Inclusive_Responses(w, q, &S_high, &FF, Responses_S2); //S state high energy table

		double Responses_P1[5];
		Inclusive_Responses(w, q, &P_low, &FF, Responses_P1); //P state low energy table

		double Responses_P2[5];
		Inclusive_Responses(w, q, &P_high, &FF, Responses_P2); //high energy table

		double f_p = 3.3/4.; //Change the occupation numbers of the states, why not ?
		double f_s = 1.8/2.;

		f_p = 1;
		f_s = 1;
	
		cout << w << "  ";
		for (int iR = 0 ; iR < 5 ; iR++)
		{
			cout << f_p*(Responses_P1[iR] +  Responses_P2[iR]) + f_s*(Responses_S1[iR] + Responses_S2[iR]) << "  "; //Take the sum of the responses from low- and high-T, their physical region doesn't overlap
		}
		cout << endl;
	}

	cout << endl << endl;

	//We can do also the cross section
	double costheta_l = cos(36.*Pi/180.);
	double E_i = 680.;

	//Now we need a lepton tensor
	LeptonTensor Lep;
	
	Lep.set_elec(); // This sets the masses to electronmass, helicity to 0 and the correct prefactor


	//To compare with other calculations, we ignore nuclear recoil. 
	//Easily done by making Mass of residual system very big
	S_high.MB = 1e10;
	S_low.MB = 1e10;
	P_high.MB = 1e10;
	P_low.MB = 1e10;
	
	double CS_i[5]; //Stores the partial cross sections, not used now

	double units = (hbc*hbc)*10*1e6; //All outputs are in units MeV^n. This makes it nb/MeV
	for (double w = 0 ; w < 250 ; w+=2.)
	{
		double E_f = E_i - w;
		double CS1 = CS_inclusive(E_i, E_f, costheta_l, &Lep, &FF, &S_low, CS_i);
		double CS2 = CS_inclusive(E_i, E_f, costheta_l, &Lep, &FF, &S_high, CS_i);
		double CS3 = CS_inclusive(E_i, E_f, costheta_l, &Lep, &FF, &P_low, CS_i);
		double CS4 = CS_inclusive(E_i, E_f, costheta_l, &Lep, &FF, &P_high, CS_i);
		cout << w << "  " << (CS1 + CS2)*units << "  " << (CS3 + CS4)*units <<  endl;  //S and P-shell contributions

	}


	return 0;

}

