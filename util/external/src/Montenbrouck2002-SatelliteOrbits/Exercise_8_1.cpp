//------------------------------------------------------------------------------
//
// Exercise_8_1.cpp
// 
// Purpose: 
//
//   Satellite Orbits - Models, Methods, and Applications 
//   Exercise 8-1: Least-squares fit using Givens rotations
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

#include <iostream>
#include <iomanip>

#include "GNU_iomanip.h"

#include "SAT_Const.h"
#include "SAT_Filter.h"
#include "SAT_VecMat.h"

using namespace std;


//------------------------------------------------------------------------------
//
// Main program
//
//------------------------------------------------------------------------------

int main() {

  // Constants

  const double t[7] = { 0.04, 0.32, 0.51, 0.73, 1.03, 1.42, 1.60 };
  const double z[7] = { 2.63, 1.18, 1.16, 1.54, 2.65, 5.41, 7.67 };

  // Variables

  int     i;
  double  b;
  Vector  a(3),c(3),d(3);
  Matrix  R(3,3);
  LSQ     PolyFit(3);


  // Header 

  cout << "Exercise 8-1: Least-squares fit using Givens rotation" 
       << endl << endl;


  // Accumulation of data equations

  for (i=0;i<7;i++) {
    
    // Data equation

    a(0) = 1.0;
    a(1) = t[i];
    a(2) = t[i]*t[i];
    b    = z[i];
    
    cout << "Observation " << noshowpos << setw(1) << i << endl << endl
         << fixed << setprecision(4) << showpos
         << "  a = " << setw(8) << a << "    "
         << "  b = " << setw(8)  << b << endl << endl;
    
    // Process data equation

    PolyFit.Accumulate(a,b);
  
    // Square-root information matrix and transformed data 

    R = PolyFit.SRIM();
    d = PolyFit.Data();

    cout << "      " << setw(8) << R.Row(0)  << "    "
         << "      " << setw(8)  << d(0) << endl;
    cout << "  R = " << setw(8) << R.Row(1)  << "    "
         << "  d = " << setw(8)  << d(1) << endl;
    cout << "      " << setw(8) << R.Row(2)  << "    "
         << "      " << setw(8)  << d(2) << endl << endl;

  }
  
  // Solution of least squares system

  PolyFit.Solve(c);

  cout << endl
       << "Adjusted polynomial coefficients" << endl << endl;
  for (i=0;i<3;i++) 
    cout << "  c(" << noshowpos << setw(1) << i << ") = " 
         << setprecision(6) << setw(10) << c(i) << endl;


  return 0;
  
}
