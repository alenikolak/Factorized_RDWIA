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

	cout << "This 'demo' code explains the top-level function calls and inputs to compute the nucleon knockout cross section from overlap matrices" << endl;
	cout << "The code is designed to be read, the output otherwise will not make much sense" << endl; 
	cout << "Now reading tables of overlaps from ./Overlaps/EDRMF/C12/S_lev1_N1 and ./Overlaps/EDRMF/C12/highTS_lev1_N1" << endl;
	cout << "This takes most of the runtime" << endl;
	//Reading the S-state tables
	S_level S_low;
	S_low.read_S("./Overlaps/EDRMF/C12/S_lev1_N1"); //Low energy table with TN steps of 1 MeV TN < 15 MeV

	S_level S_high;
	S_high.read_S("./Overlaps/EDRMF/C12/highTS_lev1_N1"); //Higher energies with TN steps of 5 MeV, TN >= 15 MeV

	cout << endl << endl;


	//We set form factors, e.g. proton form factors from Kelly
	FormFactors FF;
	FF.F_vector = Fp_kelly;
	FF.F_axial = Fzero; //No axial FF

	//And define a lepton tensor
	LeptonTensor Lep;
	Lep.set_elec(); // This sets the masses to electronmass, helicity to 0 and sets the prefactor F_X(Q^2)


	cout << "The function 'cross_section_cos' takes as input lepton and nucleon kinematics, the overlap matrix, lepton tensor, and form factors." << endl;
	cout << "These are combined as explained in section III to produce the cross section" << endl;

	//For given kinematics:
	double E_i = 1000.;  //initial lepton energy
	double E_f = 850.;  //Final-state lepton energy
	double costheta_l = 0.8; //Angle between k_i and k_f
	
	double costheta_N = 0.995; //Scattering angle of nucleon with respect to q = k_i - k_f
	double phi_N = 0; //Azimuth angle

	//There are in principle up to two possible solutions for the outgoing nucleon momentum
	double k_N_solutions[2];
	double CS_k[2];

	//The function call
	int N_sol = cross_section_cos(E_i, E_f, costheta_l , costheta_N , phi_N, &Lep, &FF, &S_high, CS_k, k_N_solutions);
	// returns the number of solutions and stores the cross sections and associated nucleon energies in CS and k_N_solutions
	// All differential cross sections in the code are in units MeV^{n} with n the appropriate power.

	cout << "As explained in appendix B, two solutions are sometimes possible for the nucleon momentum given lepton kinematics and a nucleon scattering angle" << endl;
	cout << "This is only possible in a narrow region of phase space q_{two_sol}(w) < q < q_max(w)" << endl;
	string fnm_out = "demo_two_solutions.out";
	ofstream out;
	out.open(fnm_out);
	cout << "This is illustrated for fixed w, the results are written to " << fnm_out << endl; 


	//As explained in Appendix B of the paper, two solutions are only possible in a small kinematic region where q_twosol(w) < q < q_max(w). 
	//No solutions are possible for q > q_max(w) 
	//We illustrate this in an ineficient way using the cross section function
	//For fixed energy transfer:
	double w = 120.;

	//And nuclear masses
	double M_A = S_high.MA; //Nucleus mass
	double M_B = S_high.MB; //Final-state nucleus mass

	//The bounds are
	double q_two_sol = sqrt( pow(w + M_A - MN,2) - M_B*M_B);
	double q_max = sqrt( (w+M_A)*(w+M_A) - (MN + M_B)*(MN+M_B)  );

	E_i = q_max*4; //This E_i is arbitrary, but set large enough so that the lepton kinematics are physical for any q chosen below.

	out << "# q (MeV) ,  Number solutions for k_N (costheta_N = 0.995) , Number solutions for k_N (costheta_N = 0.85) , Maximum number of solutions for any angle" << endl;
	out << "#For fixed energy transfer w = E_i - E_f = " << w << endl;
	out << "#q_{two_sol} =  " << q_two_sol << "  q_max = " << q_max << endl;
	for (double q = w + 10 ; q < q_max + 20 ; q+=2)
	{
		E_f = E_i - w;
		costheta_l = (E_i*E_i + E_f*E_f - q*q)/(2.*E_i*E_f);

		int N_sol_max = 1; 
		if (q > q_two_sol){N_sol_max = 2;} //In this region up to two solutions are possible for some angles

		if (abs(costheta_l) > 1){N_sol_max = 0;} //invalid lepton kinematics no solutions, should never happen here
		if (q > q_max){N_sol_max=0;} //Beyond this bound there are no solutions for any angle

		int N_sol =  cross_section_cos(E_i, E_f, costheta_l , costheta_N , phi_N, &Lep, &FF, &S_high, CS_k, k_N_solutions); //Number solutions for cos\theta_N = 0.995
		int N_sol_cossmaller =  cross_section_cos(E_i, E_f, costheta_l , costheta_N*0.85 , phi_N, &Lep, &FF, &S_high, CS_k, k_N_solutions); //Number of solutions for smaller cos\theta_N

		out << q << "  " << N_sol << "  " << N_sol_cossmaller << "  " << N_sol_max << endl;

	}		
	out.close();

	cout << endl << endl;
	cout << "To avoid these two solutions one can use missing momentum p_m as input instead of the scattering angle. " << endl;
	cout << "This is implemented in function 'cross_section'" << endl;
	cout << "This is still inefficient because the angular dependence of the cross section factorizes as in Eqs. (22 - 24)" << endl;
	cout << "For integrated cross sections or efficient event generation this factorization should be exploited, see e.g. [arxiv:2011.05269]" << endl;
	cout << endl << endl;

	//The two solutions can be avoided by using missing momentum as the independent variable instead of scattering angle
	
	//For given kinematics:
	E_i = 1000.;  //initial lepton energy
	E_f = 850.;  //Final-state lepton energy
	costheta_l = 0.8; //Angle between k_i and k_f
	
	double p_miss = 150.; //Missing momentum
	phi_N = 0; //Azimuth angle

	//The function call
	double CS = cross_section(E_i, E_f, costheta_l, p_miss, phi_N, &Lep, &FF, &S_high);
	//Returns a single cross section, since now p_miss determines uniquely the kinetic energy of the residual system
	
	//Of course, this function is not very efficient because the phi_N-dependence of the cross section factorizes
	//Using the notation of [arxiv:1807.11281] we define the structure functions A,B,C,D,E which are independent of phi_N
	double ABCDE[5];
	
	//Computed in 
	int ERR = Get_ABCDE(E_i, E_f, costheta_l, p_miss, &Lep, &FF, &S_high, ABCDE);
	//In an identical way as done for the cross section
	// The cross section is then  CS =  A + B*cos(phi) + C*cos(2phi) + D*sin(phi) + E*sin(2phi)
	// The definition of the functions can be read from Eqs. (22 - 24) of our paper.
	// This factorization may be used to more efficiently generate events in the same way as in Ref. [2011.05269]
	
	fnm_out = "phi_electron.out";
	out.open(fnm_out);

	cout << "The phi-dependence for (unpolarized) electron scattering is written to " << fnm_out << endl;
	cout << "In this case there is no odd phi-dependence in the cross section as explained in section IV" << endl;

	out << "#Phi_N (rad),  CS(phi),  A , B*cos(phi) , C*cos(2phi), D*sin(phi), E*sin(2phi) " << endl; 
	for (double phi = 0 ; phi < 2*Pi ; phi+=0.05)
	{
		CS = cross_section(E_i, E_f, costheta_l, p_miss, phi, &Lep, &FF, &S_high);
		out << phi << "  " << CS << "  ";
		out << ABCDE[0] << "  " ;
		out << ABCDE[1]*cos(phi) << "  " << ABCDE[2]*cos(2*phi) << "  ";
		out << ABCDE[3]*sin(phi) << "  " << ABCDE[4]*sin(2*phi) << endl; //For unpolarized electron scattering these are zero
	}

	out.close();
	
	cout << endl << endl;

	
	fnm_out = "phi_neutrino.out";
	out.open(fnm_out);
	cout << "We can easily change the lepton tensor and form factors to treat the neutrino scattering case written to " << fnm_out << endl;
	cout << "There now is odd-phi dependence." << endl;
	cout << "This dependence does not arise in the PWIA, as explained in Section IV.C" << endl;
	cout << "See also the discussion in appendix C" << endl;

	//To treat neutrino scattering we simply change the form factors:
	FF.F_vector = FIV_kelly; //Isovector form factors from Kelly
	FF.F_axial = FA_dipole; //Dipole form factors

	//And change the lepton tensor
	Lep.set_CC_numu(); // This sets initial mass to zero, final mass to muon, helicity to -1 and sets the prefactor F_X(Q^2)

	//We can now compute ABCDE again
	ERR = Get_ABCDE(E_i, E_f, costheta_l, p_miss, &Lep, &FF, &S_high, ABCDE);

	//and output the cross section	
	out << "#Phi_N (rad),  CS(phi),  A , B*cos(phi) , C*cos(2phi), D*sin(phi), E*sin(2phi) " << endl; 
	for (double phi = 0 ; phi < 2*Pi ; phi+=0.05)
	{
		CS = cross_section(E_i, E_f, costheta_l, p_miss, phi, &Lep, &FF, &S_high);
		out << phi << "  " << CS << "  ";
		out << ABCDE[0] << "  " ;
		out << ABCDE[1]*cos(phi) << "  " << ABCDE[2]*cos(2*phi) << "  ";
		out << ABCDE[3]*sin(phi) << "  "; //In the PWIA these do not arise
		out << ABCDE[4]*sin(2*phi) << endl; //In the PWIA these do not arise
		//The odd phi-dependence , in particular the sin(2phi) term, is uniquely due to nucleon distortion, see discussion in appendix C.
	}

	out.close();

	cout << endl << endl;
	fnm_out = "CS_inclusive_nu.out";
	out.open(fnm_out);
	cout << "Finally, the inclusive cross section is simply obtained from the inclusive responses (Eq. 22 integrated over p_m and phi) and written to " << fnm_out << endl;


//
//	//For the inclusive neutrino cross section at familiar kinematics:
	costheta_l = cos(36.*Pi/180.);
	E_i = 680.;
	E_f = 100.;

	double CS_partial[5]; //Stores the partial cross sections : CC,LL,CL, T, Tp contributions


	//We simply call
	CS = CS_inclusive(E_i, E_f, costheta_l, &Lep, &FF, &S_high, CS_partial);
	//Which integrates the inclusive responses and contracts with the lepton tensor

	out << "# w (MeV) , CS " << endl;
	for (double w = 0 ; w < 270 ; w+=2.)
	{
		double E_f = E_i - w;
		double CS1 = CS_inclusive(E_i, E_f, costheta_l, &Lep, &FF, &S_low, CS_partial); //Cross section from low-energy table
		double CS2 = CS_inclusive(E_i, E_f, costheta_l, &Lep, &FF, &S_high, CS_partial); //From high energy table
		out << w  << "  " << CS1 + CS2 << endl; //The physical region of the tables doesn't overlap so we can just add the cross sections	
	}

	out.close();


	fnm_out = "CS_inclusive_e.out";
	out.open(fnm_out);
	cout << "Again, one can easily change the lepton tensor and form factors to obtain the electron cross section written to " << fnm_out << endl;

	//For electron scattering on a proton:
	FF.F_vector = Fp_kelly; //Isovector form factors from Kelly
	FF.F_axial = Fzero; //No axial

	//And change the lepton tensor
	Lep.set_elec();

	out << "# w (MeV) , CS " << endl;
	for (double w = 0 ; w < 270 ; w+=2.)
	{
		double E_f = E_i - w;
		double CS1 = CS_inclusive(E_i, E_f, costheta_l, &Lep, &FF, &S_low, CS_partial); //Cross section from low-energy table
		double CS2 = CS_inclusive(E_i, E_f, costheta_l, &Lep, &FF, &S_high, CS_partial); //From high energy table
		out << w  << "  " << CS1 + CS2 << endl; //The physical region of the tables doesn't overlap so we can just add the cross sections	
	}

	out.close();

	cout << endl << endl;
	cout << "If one adds the p-shell contributions to the folder ./Overlaps/EDRMF/C12/ one can reproduce Fig. 5 from [arxiv:1904.10696] by uncommenting the last block of code" << endl;
	//Commented out below between the stars, will work if you put the p-shell overlaps in the correct folder.

/***
	//We read the P-shell tables
	cout << "Reading again ... " << endl;
	S_level P_low;
	P_low.read_S("./Overlaps/EDRMF/C12/S_lev2_N1"); //Low energy table with TN steps of 1 MeV TN < 15 MeV
	S_level P_high;
	P_high.read_S("./Overlaps/EDRMF/C12/highTS_lev2_N1"); //Higher energies with TN steps of 5 MeV, TN >= 15 MeV

	//We will approximate all outgoing nucleons as protons (N=1) above. 
	//In this way we make the mistake of including the coulomb potential for what should be a neutron state. 
	//In principle one can use the neutron tables and do everything correctly, effect is small.
	

	//To compare with [arxiv:1904.10696] we ignore nuclear recoil. 
	//Easily done by making Mass of residual system very big
	S_high.MB = 1e10;
	S_low.MB = 1e10;
	P_high.MB = 1e10;
	P_low.MB = 1e10;

	
	double units = (hbc*hbc)*10*1e6; //All outputs are in units MeV^n. This makes it nb/MeV

	//Also get a form factor struct for the neutrons:
	FormFactors FFn;
	FFn.F_vector = Fn_kelly;
	FFn.F_axial = Fzero; //No axial FF
	

	fnm_out = "CS_inclusive_e_pn.out";
	//The result is panel b of Fig. 5 in [arxiv:1904.10696]
	out.open(fnm_out);
	for (double w = 0 ; w < 250 ; w+=2.)
	{
		double E_f = E_i - w;
		double CS1 = CS_inclusive(E_i, E_f, costheta_l, &Lep, &FF, &S_low, CS_partial);
		double CS2 = CS_inclusive(E_i, E_f, costheta_l, &Lep, &FF, &S_high, CS_partial);
		double CS3 = CS_inclusive(E_i, E_f, costheta_l, &Lep, &FF, &P_low, CS_partial);
		double CS4 = CS_inclusive(E_i, E_f, costheta_l, &Lep, &FF, &P_high, CS_partial);

		//One could change the occupation number to take into account SRC of course:
		//f_p = 3.3/4.
		//f_s = 1.8/2.  from [arix:2306.10823] based on Benhar SF.
		//but this is not done here

		double CS_protons = (CS1+CS2+CS3+CS4)*units;

		//Now evaualte the neutron contribution:
		CS1 = CS_inclusive(E_i, E_f, costheta_l, &Lep, &FFn, &S_low, CS_partial);
		CS2 = CS_inclusive(E_i, E_f, costheta_l, &Lep, &FFn, &S_high, CS_partial);
		CS3 = CS_inclusive(E_i, E_f, costheta_l, &Lep, &FFn, &P_low, CS_partial);
		CS4 = CS_inclusive(E_i, E_f, costheta_l, &Lep, &FFn, &P_high, CS_partial);

		double CS_neutrons = (CS1+CS2+CS3+CS4)*units;
		out << w << "  " << CS_protons << "  " << CS_neutrons << endl;  

	}

	out.close();

***/

	return 0;

}

