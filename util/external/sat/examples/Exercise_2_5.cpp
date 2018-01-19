//------------------------------------------------------------------------------
//
// Exercise_2_5.cpp
// 
// Purpose: 
//
//   Satellite Orbits - Models, Methods, and Applications
//   Exercise 2-5: Sunsynchronous repeat orbits
//
// Notes:
//
//   This software is protected by national and international copyright. 
//   Any unauthorized use, reproduction or modificaton is unlawful and 
//   will be prosecuted. Commercial and non-private application of the 
//   software in any form is strictly prohibited unless otherwise granted
//   by the authors.
//
//   The code is provided without any warranty; without even the implied 
//   warranty of merchantibility or fitness for a particular purpose.
//
// Last modified:
//
//   2000/03/04  OMO  Final version (1st edition)
//   2005/04/14  OMO  Final version (2nd reprint)
//
// (c) 1999-2005  O. Montenbruck, E. Gill
//
//------------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <iomanip>

#include "GNU_iomanip.h"

#include "SAT_Const.h"

using namespace std;


// x mod y
double Modulo(double x, double y) {
  return y*((x/y)-floor(x/y));
};


//------------------------------------------------------------------------------
//
// Main program
//
//------------------------------------------------------------------------------

int main() {

  const int    K   =  3;                              // Number of cycles
  const int    N   = 43;                              // Number of orbits
  const double T_N = (double) K / (double) N;         // Draconic period [d]
  const double J_2 = 1.08263e-3;                      // Oblateness coefficient

  const double Omega_dot = +0.985647240*Rad;          // Nodal rate [rad/d]
  
  double omega_dot = 0.0;                             // Arg. of latitude rate
  double Delta_n   = 0.0;                             // Perturb. of mean motion
  double n_0       = 0.0;
  double a         = 0.0; 
  double h         = 0.0;
  double i         = 0.0;

  // Header 

  cout << "Exercise 2-5: Sunsynchronous repeat orbits" << endl << endl;

  cout << "  Draconic period          T_N = " 
       << fixed << setprecision(7) << setw(12) << T_N << " d" << endl 
       << "  Rate                 2pi/T_N = "
       << setprecision(5) << setw(12) << 360.0/T_N << " deg/d" << endl
       << "  Node offset  Delta lam_Omega = "
       << setprecision(5) << setw(12) << -K*360.0/N << " deg/Orbit" << endl
       << endl;

  // Iteration
  
  for (int iterat=0; iterat<=2; iterat++) {
  
    // Secular rates of argument of perigee and mean anomaly [rad/d]
    if (iterat==0) {
      omega_dot = 0.0; 
      Delta_n   = 0.0; 
    }
    else {
      omega_dot = -0.75*n_0*J_2*pow(R_Earth/a,2)*(1.0-5*pow(cos(i),2)); 
      Delta_n   = -0.75*n_0*J_2*pow(R_Earth/a,2)*(1.0-3*pow(cos(i),2));
    };
    
    // Mean motion, semimajor axis and altitude
    n_0 = pi2/T_N - Delta_n - omega_dot;                // [rad/d]
    a   = pow( GM_Earth/pow(n_0/86400.0,2), 1.0/3.0 );  // [m]
    h   = a - R_Earth;                                  // [m]
  
    // Inclination [rad}
    i   = acos(-2.0*Omega_dot/(3.0*n_0*J_2)*pow(a/R_Earth,2));
                                                      
    cout << "  Iteration " << iterat << endl << endl;
    cout << "  Arg. perigee rate  omega_dot = " 
         << setprecision(5) << setw(12) << Deg*omega_dot << " deg/d" << endl 
         << "  Perturb. mean motion   n-n_0 = " 
         << setprecision(5) << setw(12) << Deg*Delta_n << " deg/d" << endl  
         << "  Mean motion              n_0 = " 
         << setprecision(5) << setw(12) << Deg*n_0 << " deg/d" << endl 
         << "  Semimajor axis             a = " 
         << setprecision(3) << setw(12) << a/1000.0 << " km" << endl 
         << "  Altitude                   h = " 
         << setw(12) << h/1000.0 << " km" << endl 
         << "  Inclination                i = " 
         << setw(12) << Deg*i << " deg" << endl 
         << endl;
  
  };

  // Ascending nodes

  cout << "  Greenwich longitude of ascending node" << endl << endl;
  cout << "      Day 1          Day 2          Day 3    " << endl;
  cout << "  Orbit   [deg]  Orbit   [deg]  Orbit   [deg] " << endl;
  cout << fixed << setprecision(2);

  for (int I=0; I<15; I++) {
    for (int J=0; J<3; J++) {
      cout << setw(5) << 15*J+I << setw(10) 
           << Modulo(-(15*J+I)*K*360.0/N+180.0,360.0)-180.0;
    };
    cout << endl;
  };

  return 0;
  
}
