#include "S_tables.h"
#include "Constants.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <cstdarg>

int ERR_f(int line, string expected, string found)
{

	cerr << "File error at line : " << line << "  expected : " << expected << "  but actually found : " << found << endl;

	return -1;

}


int ERR_f(int line, int expected, int found)
{

	cerr << "File error at line : " << line << "  expected : " << expected << "  but actually found : " << found << endl;

	return -1;

}


int ERR_f(int line, double expected, double found)
{

	cerr << "File error at line : " << line << "  expected : " << expected << "  but actually found : " << found << endl;

	return -1;

}

int S_level::read_S(string fnm)
{
	ifstream ing;
	ing.open(fnm);
	
	int line=1;
	ing >> this->two_J >> this->E_kappa >> this->MA >>this-> N_TN >> this->N_pm >> this->N_cos;
	ing.ignore(1000, '\n');

	if (this->two_J%2 == 0){return ERR_f(line, "odd value for 2J", to_string(this->two_J));}
	this->nMJ = (this->two_J+1);
	
	//Mass of the residual system
	this->MB = this->MA - MN + this->E_kappa; //Only place where a 'constant' enters

	int nMJ_infile = this->nMJ/2; //Only half of the m_j states need to be stored in the file
	int block_length = this->N_pm*this->N_cos*nMJ_infile*2; //For every TN there is a grid with this size

	//Chech here also if the sizes are too big maybe

	//Read the grids
	line++;
	int n;
	ing >> n;
	ing.ignore(1000, '\n');
	if (n != this->N_TN){ return ERR_f(line, N_TN, n); }

	line++;
	for (int iT = 0 ; iT < n ; iT++)
	{	
		int i;
		double T;
		ing >> i >> T;
		ing.ignore(1000, '\n');
		if (i != iT){  return ERR_f(line, iT, i); }
		this->TN_grid[iT] = T;

		line++;
	}


	ing >> n;
	ing.ignore(1000, '\n');
	if (n != this->N_pm){ return ERR_f(line, N_pm, n); }

	line++;
	for (int iT = 0 ; iT < n ; iT++)
	{	
		int i;
		double p;
		ing >> i >> p;
		ing.ignore(1000, '\n');
		if (i != iT){  return ERR_f(line, iT, i); }
		this->pm_grid[iT] = p;

		line++;
	}


	ing >> n;
	ing.ignore(1000, '\n');
	if (n != this->N_cos){ return ERR_f(line, this->N_cos, n); }

	line++;
	for (int iT = 0 ; iT < n ; iT++)
	{	
		int i;
		double cos;
		ing >> i >> cos;
		ing.ignore(1000, '\n');
		if (i != iT){  return ERR_f(line, iT, i); }
		this->cos_grid[iT] = cos;

		line++;
	}

	//Now we can read the actual data
	for (int iT = 0 ; iT < this->N_TN ; iT++)
	{
		//Every block should start with the TN value
		int i;
		double T;
		ing >> i >> T;
		if (iT != i){return ERR_f(line,iT,i);}
		if (T != this->TN_grid[iT]){ return ERR_f(line,this->TN_grid[iT],T); }
		ing.ignore(1000, '\n');
		line++;

		for (int i = 0 ; i < block_length ; i++)
		{
			int ip, ic, im, is;
			ing >> ip >> ic >> im >> is; //Ok
			for ( int mu = 0 ; mu < 16 ; mu++)
			{	
				double R, I;
				ing >> R >> I; 
				this->S_kappa[mu][iT][ip][ic][is][im] = complex<double>(R,I);
			}
			ing.ignore(1000, '\n');
			line++;
		}

	}
	ing.close();
	return 0;
}

int index_min_const(double x, double *vec_x, int N_x)
{
	//Get the index such that vec_x[i] <= x < vec_x[i+1]
	//Assumption is that the grid has a constant step
	if (x > vec_x[N_x-1]){return -1;}

	if (x < vec_x[0]){ return -1;} //This might leave a hole at low p_m, could include pm=0 case in tables

	double step = vec_x[4] - vec_x[3];

	int i = static_cast <int> (floor( ( x- vec_x[0])/step));	

	return i;
}

complex<double> interpolate_3d(complex<double> S[2][2][2], double del_1, double del_2, double del_3)
{
	complex<double> S_2d[2][2];
	complex<double> S_1d[2];

	S_2d[0][0] = S[0][0][0]*(1.-del_1) + S[1][0][0]*del_1;
	S_2d[0][1] = S[0][0][1]*(1.-del_1) + S[1][0][1]*del_1;
	S_2d[1][0] = S[0][1][0]*(1.-del_1) + S[1][1][0]*del_1;
	S_2d[1][1] = S[0][1][1]*(1.-del_1) + S[1][1][1]*del_1;

	S_1d[0] = S_2d[0][0]*(1.-del_2) + S_2d[1][0]*del_2;
	S_1d[1] = S_2d[0][1]*(1.-del_2) + S_2d[1][1]*del_2;

	return S_1d[0]*(1.-del_3) + S_1d[1]*del_3;
}

void S_level::copy_all(int spin, int imj, int iT, int ip, int icos, complex<double> s_el[16][2][2][2])
{
	//Probably unnecesary but clearer using this
	for (int mu = 0 ; mu < 16 ; mu ++)
	{
		for (int itt = 0 ; itt <2 ; itt++)
		{
			for (int ipp = 0 ; ipp <2 ; ipp++)
			{
				for (int icc = 0 ; icc <2 ; icc++)
				{
					s_el[mu][itt][ipp][icc] = this->S_kappa[mu][iT+itt][ip + ipp][icos+icc][spin][imj];
				}
			}
		}
	}
}


int S_level::eval(double TN, double pm, double cos)
{

	//evaluate by 3d linear interpolation	

	int iT = index_min_const(TN, this->TN_grid, N_TN);
	if (iT < 0){return -1;}
	int ip = index_min_const(pm, this->pm_grid, N_pm);
	if (ip < 0){return -1;}
	int ic = index_min_const(cos, this->cos_grid, N_cos);
	if (ic < 0){return -1;}


	double del_T = this->TN_grid[iT+1] - this->TN_grid[iT];
	double del_p = this->pm_grid[ip+1] - this->pm_grid[ip];
	double del_c = this->cos_grid[ic+1] - this->cos_grid[ic];

	double T_L = (TN - TN_grid[iT])/del_T;
	double p_L = (pm - pm_grid[ip])/del_p;
	double c_L = (cos - cos_grid[ic])/del_c;

	complex<double> S_int_grid[16][2][2][2];

	//number of m-states
	int nM = this->two_J+1;
	for (int im = 0 ; im < nM/2 ; im++ ) //Only loop over half the m states, the others are determined by a flip of angular momenta
	{
		for (int spin = 0 ; spin < 2 ; spin++)
		{
			int spin_flip = (spin+1)%2;
			int m_flip = im + nM/2;

			//Put all the relevant data to be interpolated in S_int_grid
			this->copy_all(spin, im, iT, ip, ic, S_int_grid);
			
			//scalar
			this->sc[spin][im] = interpolate_3d(S_int_grid[0], T_L, p_L, c_L);
			this->sc[spin_flip][m_flip] = this->sc[spin][im];


			//Pseudoscalar
			this->ps[spin][im] = interpolate_3d(S_int_grid[1], T_L, p_L, c_L);
			this->ps[spin_flip][m_flip] = -this->ps[spin][im];  //We choose to flip the axial, unless there is a \gamma^2 in the tensor. One could do the opposite, overall sign doesn't matter


			//counter for the tensor one:
			int iTensor = 0;
			for (int mu = 0 ; mu < 4 ; mu++)
			{
				//Vector & axialvector
				this->V[mu][spin][im] = interpolate_3d(S_int_grid[2+mu], T_L, p_L, c_L);
				this->A[mu][spin][im] = interpolate_3d(S_int_grid[6+mu], T_L, p_L, c_L);

				this->A[mu][spin_flip][m_flip] = -A[mu][spin][im];
				this->V[mu][spin_flip][m_flip] =  V[mu][spin][im];

				if (mu == 2)
				{
					this->A[mu][spin_flip][m_flip] = A[mu][spin][im];
					this->V[mu][spin_flip][m_flip] = -V[mu][spin][im];
				}

				//Tensor:

				this->T[mu][mu][spin][im] = 0.;
				this->T[mu][mu][spin_flip][m_flip] = 0.; //Just to be sure

				for (int nu = mu+1 ; nu < 4 ; nu++)
				{
					this->T[mu][nu][spin][im] = interpolate_3d(S_int_grid[10 + iTensor], T_L, p_L, c_L);
					iTensor++;
				
					//Setting the flipped one	
					this->T[mu][nu][spin_flip][m_flip] = this->T[mu][nu][spin][im];
					if (nu == 2 || mu == 2){
						this->T[mu][nu][spin_flip][m_flip] = -this->T[mu][nu][spin][im];
					}

					this->T[nu][mu][spin][im] = -this->T[mu][nu][spin][im]; //Antisymmetric tensor, setting the lower diagonal elements
					this->T[nu][mu][spin_flip][m_flip] = -this->T[mu][nu][spin_flip][m_flip];  								
				}
				
	
			}
			
		}//spin
	
	}//m

	return 0;

}
