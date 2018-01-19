//------------------------------------------------------------------------------
//
// Exercise_2_2.cpp
// 
// Purpose: 
//
//   Satellite Orbits - Models, Methods, and Applications
//   Exercise 2-2: Solution of Kepler's equation
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
#include "SAT_Kepler.h"

using namespace std;


//------------------------------------------------------------------------------
//
// Main program
//
//------------------------------------------------------------------------------

int main() {


  // Variables

  const double  MeanAnom[] = {4.0,50.0};   // [deg]
  int     i;  
  double  M,E,E_ref,e;


  for (int iCase=0; iCase<=1; iCase++) {
  
    // Test Case

    M     = MeanAnom[iCase]*Rad;
    e     = 0.72; 
    E_ref = EccAnom(M,e);

    // Header

    cout << "Exercise 2-2: Solution of Kepler's equation" << endl << endl;
  
    cout << fixed << setprecision(11) 
         << "  M =" << setw(14) << M << endl
         << "  e =" << setw(14) << e << endl
         << "  E =" << setw(14) << E_ref << endl;

    // Newton's iteration

    cout << endl
         << "  a) Newton's iteration" << endl << endl
         << "  i         E         Accuracy  sin/cos" << endl;

    E = M;
    i = 0;
    do {
      i++;
      E = E - (E-e*sin(E)-M)/(1.0-e*cos(E)); 
      cout << setw(3) << i 
           << fixed      << setprecision(11) << setw(16) << E 
           << scientific << setprecision(2)  << setw(11) << fabs(E-E_ref)
           << setw(6) << 2*i
           << endl;
    } while ( fabs(E-E_ref)> 1.0e-10 );

    // Fixed point iteration

    cout << endl
         << "  b) Fixed point iteration" << endl << endl
         << "  i         E         Accuracy  sin/cos" << endl;

    E = M;
    i = 0;
    do {
      i++;
      E = M + e*sin(E); 
      cout << setw(3) << i 
           << fixed      << setprecision(11) << setw(16) << E 
           << scientific << setprecision(2)  << setw(11) << fabs(E-E_ref)
           << setw(6) << i
           << endl;
    } while ( fabs(E-E_ref)> 1.0e-10 );

    cout << endl << endl;   

  };

  return 0;
  
}
