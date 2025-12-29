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

	//Getting the level information from stdin
	//
	std::string nucleus;
	int lev;
	int nuc;

	cout << "#Give nucleus (\"C12\" or \"O16\")" << endl;
	cin >> nucleus;
	cout << "#Give level ( int from 1 - 6) )" << endl;
	cin >> lev;
	cout << "#Give final nucleon( 1 (proton) or 2 (neutron) )" << endl;
	cin >> nuc;
	cout << "#INFO: " << nucleus << " " << lev << " " << nuc << endl;
	std::string fnm = "./Overlaps/EDRMF/" + nucleus + "/S_lev" + to_string(lev) + "_N"+to_string(nuc);
	cout << "# Reading " << fnm << endl;
	S_level S_low;
	S_low.read_S(fnm); //Low energy table with TN steps of 1 MeV TN < 15 MeV


	fnm = "./Overlaps/EDRMF/" + nucleus + "/highTS_lev" + to_string(lev) + "_N"+to_string(nuc);
	cout << "# Reading " << fnm << endl;
	S_level S_high;
	S_high.read_S(fnm); //high energy table with TN steps of 5 MeV TN >= 15 MeV


	//We set form factors, full vector + axial
	FormFactors FF;
	FF.F_vector = FIV_kelly; //Isovector form factors from kelly
	FF.F_axial = FA_dipole;  


	//Output cosine dependence of all responses
	//If these match, everything should match
	
	//For low energy response:
	double E_kappa = S_low.E_kappa;
	double w = E_kappa + 10; //Low energy table, not interpolated in energy: T_N = w - E_kappa = 10

	int n_q = 10;
	double q_step = 15.;
	double q_min = w + 10.;
	for (int iq = 0; iq < n_q ; iq++)
	{
		double q = q_min + iq*q_step; 	
		cosdep_Responses(w,q, &S_low, &FF); //Ignores recoil
	}
	std::cout << endl << endl;


	double q = 100.;
	for (double w = E_kappa + 1.22 ; w < E_kappa + 14; w += 1) //Interpolating in T_N
	{
		cosdep_Responses(w,q, &S_low, &FF);
	}

	std::cout << endl << endl;

	//For high-TN response:
	w = E_kappa + 100.;
	n_q = 10;
	q_step = 35;
	q_min = w + 10.;
	for (int iq = 0; iq < n_q ; iq++)
	{
		double q = q_min + iq*q_step; 	
		cosdep_Responses(w,q, &S_high, &FF); //Ignores recoil
	}

	std::cout << endl << endl;
	
	q = 400.;
	for (double w = E_kappa + 50.22 ; w < q ; w += 20.) //With the smaller tables there will be a piece zero above T_N = 260, but ok
	{
		cosdep_Responses(w,q, &S_high, &FF); //Ignores recoil
	}
	

	std::cout << endl << endl;

	//Output of integrated inclusive responses as well
	for (double q = 100. ; q < 600 ; q+=100)
	{

		for (double w = 0 ; w < q ; w+=1)
		{
			double Responses_S1[5];
			Inclusive_Responses(w, q, &S_low, &FF, Responses_S1); // low enery table

			double Responses_S2[5];
			Inclusive_Responses(w, q, &S_high, &FF, Responses_S2); // high energy table

//			double f_p = 3.3/4.; //Change the occupation numbers of the states, why not ?
//			double f_s = 1.8/2.;
		
			cout << w << "  " << q << "  ";
			for (int iR = 0 ; iR < 5 ; iR++)
			{
				cout << 1*(Responses_S1[iR] +  Responses_S2[iR]) << "  "; //Take the sum of the responses from low- and high-T, their physical region doesn't overlap
			}

			//Add also the responses where we ignore nuclear recoil:
			Inclusive_Responses_th(w, q, &S_low, &FF, Responses_S1); // low enery table

			Inclusive_Responses_th(w, q, &S_high, &FF, Responses_S2); // high energy table
			for (int iR = 0 ; iR < 5 ; iR++)
			{
				cout << 1*(Responses_S1[iR] +  Responses_S2[iR]) << "  "; 
			}

			cout << endl;
		}
		std::cout << endl << endl;

	}

	return 0;
}

