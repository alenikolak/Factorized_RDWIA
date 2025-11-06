#ifndef KINEMATICS_H
#define KINEMATICS_H


int Kproton(double MA, double MB, double w, double q, double costheta_qN, double kN_sols[])
{
	//Proton momentum including nuclear recoil given the scattering angle

  //Here one cannot set MB large to ignore recoil
  //To be sure nothing bad is going on we check:
  if (MB > MA){return 0;}
	
  //We solve the equation
  // \sqrt(c) = w + E_a = sqrt(p^2 + 2*p*a_1 + a_0) + sqrt(p^2 + M2) for p

  // with a_o = M_B^2 + (q)^2 
  // a_1 = - qcostheta_qN 

  double c = pow(w + MA,2);
  double a_0 = pow(MB,2) + pow(q,2);
  double a_1 = -1.*q*costheta_qN;

  //This is 2nd order poly in p:
  //#p^0 * [ c^2 - 2*c*(M^2 +a_0) + (a_0 - M^2)^2 ] 
  //#p^1 * 4a_1[ a_0 - M^2 - c ]
  //#p^2 * 4 ( a_1 - c)
  //

  double A = 4*(a_1*a_1 - c); //Negative ~ -4*M_A^2 , unless q is absolutely huge
  double B = 4*a_1*(a_0 - MN2 - c); //usually positive
  double C = c*c - 2*c*(a_0 + MN2) + pow(a_0 - MN2,2);


  double det2 = B*B - 4.*A*C;

  if (det2 < 0){ return 0; } //No solutions

  double det = sqrt(det2);



  double k_plus = -1.*(B+det)/(2.*A);

  if (k_plus < 0){ return 0; } //No solutions: k_pi_min will be smaller unless A is positive which shouldn't happen

  //Check E_conservation:
  double kN = k_plus;
  double kN2 = kN*kN;

  if( fabs( w + MA  - sqrt(kN2 + MN2) - sqrt( MB*MB + q*q -2.*kN*q*costheta_qN + kN2)) > 0.001 ){return 0;} //For testing, should always be ok

  kN_sols[0] = kN;
  int nr_solutions = 1;

  double k_pi_min = (-1.*B+det)/(2.*A);
 
  if (k_pi_min < 0){return nr_solutions;} 
  kN = k_pi_min;
  kN2 = kN*kN;

  if( fabs( w + MA  - sqrt(kN2 + MN2) - sqrt( MB*MB + q*q -2.*kN*q*costheta_qN + kN2)) > 0.001 ){return nr_solutions;}

  kN_sols[1] = kN;
  nr_solutions = 2;

  return nr_solutions;

}




#endif
