//------------------------------------------------------------------------------
//
// Exercise_2_1.cpp
// 
// Purpose: 
//
//   Satellite Orbits - Models, Methods, and Applications
//   Exercise 2-1: Orbit raising using Hohmann transfer
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


//------------------------------------------------------------------------------
//
// Main program
//
//------------------------------------------------------------------------------

int main() {

  // Constants

  const double h_1 = 750.0e3;     // Initial altitude [m]
  const double h_2 = 775.0e3;     // Final   altitude [m]

  const double r_1 = R_Earth+h_1; // Initial radius [m]
  const double r_2 = R_Earth+h_2; // Final   radius [m]

  // Variables

  double v_1,v_2,a_t,v_p,v_a;

  // Circular velocities

  v_1 = sqrt(GM_Earth/r_1);  // [m/s]
  v_2 = sqrt(GM_Earth/r_2);  // [m/s]

  // Transfer orbit 

  a_t = 0.5*(r_1+r_2);                // [m]
  v_p = sqrt(GM_Earth*r_2/(a_t*r_1)); // [m/s]
  v_a = sqrt(GM_Earth*r_1/(a_t*r_2)); // [m/s]

  // Header

  cout << "Exercise 2-1: Orbit raising using Hohmann transfer" << endl << endl;
  
  // Results

  cout << fixed << setprecision(2) 
    << " Initial altitude     h_1    " << setw(10) << h_1/1000 << " km" << endl
    << " Final   altitude     h_2    " << setw(10) << h_2/1000 << " km" << endl     
    << " Circular velocity    v_1    " << setw(10) << v_1 << " m/s" << endl
    << " Circular velocity    v_2    " << setw(10) << v_2 << " m/s" << endl
    << " Difference                  " << setw(10) << v_1-v_2 << " m/s" << endl
    << endl;

  cout << fixed << setprecision(2) 
    << " Transfer orbit sma   a_t    " << setw(10) << a_t/1000 << " km" << endl
    << " Pericenter velocity  v_p    " << setw(10) << v_p << " m/s" << endl
    << " Apocenter  velocity  v_a    " << setw(10) << v_a << " m/s" << endl
    << " Difference           v_p-v_1" << setw(10) << v_p-v_1 << " m/s" << endl
    << " Difference           v_2-v_a" << setw(10) << v_2-v_a << " m/s" << endl
    << " Total velocity diff.        " << setw(10) << v_2-v_a+v_p-v_1 << " m/s" 
    << endl << endl;

  return 0;

}
