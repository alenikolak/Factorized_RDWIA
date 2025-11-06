#ifndef S_TABLES_H
#define S_TABLES_H

#include <complex>
#include <string>
#include <memory>

using namespace std;

struct S_level
{	
	//Grids to store the matrices
	const static int m_max = 2; //Enough for the j=1/2 and j=3/2 cases, we only store the m_j > 0 case
	const static int n_pm_max = 100; 
	const static int n_cos_max = 30;
	const static int n_TN_max = 100;
	
	complex<double> (*S_kappa)[n_TN_max][n_pm_max][n_cos_max][2][m_max] = new complex<double>[16][n_TN_max][n_pm_max][n_cos_max][2][m_max];
	//The 16 components of the matrix are:
	//scalar, pseudoscalar, V^\mu, A^\mu, T^{\mu\nu} with nu > mu

	double TN_grid[n_TN_max];
	double cos_grid[n_cos_max];
	double pm_grid[n_pm_max];

	int N_pm = 0;
	int N_TN = 0;
	int N_cos = 0;

	//Information on the level
	double E_kappa = 0; //
	double MA = 0;
	double MB = 0;
	int two_J;
	int nMJ;

	//Output, overkill to store all components, but clearer
	complex<double> sc[2][m_max*2];
	complex<double> ps[2][m_max*2];
	complex<double> V[4][2][m_max*2];
	complex<double> A[4][2][m_max*2];
	complex<double> T[4][4][2][m_max*2];

	
	int eval(double TN, double pm, double cos);
	int read_S(string fnm);


	void copy_all(int spin, int imj, int iT, int ip, int icos, complex<double> s_el[16][2][2][2]);

	~S_level() {
		delete[] S_kappa;
	}
	
};

int read_S(string fnm);

#endif
