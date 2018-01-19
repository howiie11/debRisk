//------------------------------------------------------------------------------
//
// Exercise_4_1.cpp
// 
// Purpose: 
//
//   Satellite Orbits - Models, Methods, and Applications
//   Exercise 4-1: Runge-Kutta 4th-order Integration
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

#include "SAT_DE.h"
#include "SAT_Kepler.h"
#include "SAT_VecMat.h"

using namespace std;



//------------------------------------------------------------------------------
//
// f_Kep6D
//
// Purpose:
// 
//   Computes the derivative of the state vector for the normalized (GM=1)
//   Kepler's problem in three dimensions
//
// Note:
//
//   pAux is expected to point to an integer variable that will be incremented
//   by one on each call of Deriv
//
//------------------------------------------------------------------------------

void f_Kep6D ( double t, const Vector& y, Vector& yp, void* pAux )
{
  
  // Pointer to auxiliary integer variable used as function call counter
  int* pCalls = static_cast<int*>(pAux);
  
  // State vector derivative
  Vector r = y.slice(0,2);
  Vector v = y.slice(3,5);
  yp = Stack ( v, -r/(pow(Norm(r),3)) );
  
  // Increment function call count
  (*pCalls)++;

};


//------------------------------------------------------------------------------
//
// Main program
//
//------------------------------------------------------------------------------

int main() {

  // Constants
  
  const double  GM    = 1.0;                   // Gravitational coefficient
  const double  e     = 0.1;                   // Eccentricity
  const double  t_end = 20.0;                  // End time
  const Vector  Kep(1.0,e,0.0,0.0,0.0,0.0);    // (a,e,i,Omega,omega,M)
  const Vector  y_ref = State(GM,Kep,t_end);   // Reference solution

  const int     Steps[] = { 50, 100, 250, 500, 750, 1000, 1500, 2000 };
  
  // Variables
  
  int     nCalls;                              // Function call count
  int     iCase;
  double  t,h;                                 // Time and step size
  Vector  y(6);                                // State vector
  
  RK4     Orbit(f_Kep6D,6,&nCalls);            // Object for integrating the
                                               // differential equation
                                               // defined by f_Kep6D using the
                                               // 4th-order Runge-Kutta method
  // Header 

  cout << "Exercise 4-1: Runge-Kutta 4th-order integration" << endl << endl;
  cout << "  Problem D1 (e=0.1)" << endl << endl;
  cout << "  N_fnc   Accuracy   Digits " << endl;
    
  // Loop over test cases

  for (iCase=0; iCase<8; iCase++) {
  
    // Step size
    h = t_end/Steps[iCase];

    // Initial values
    t = 0.0;
    y = Vector( 1.0-e,0.0,0.0, 0.0,sqrt((1+e)/(1-e)),0.0 );
    nCalls = 0;

    // Integration from t=t to t=t_end
    for (int i=1; i<=Steps[iCase]; i++)  
      Orbit.Step( t, y, h );
 
    // Output
    cout << fixed  << setw(6) << nCalls
         << scientific << setprecision(3) << setw(13)
         << Norm(y-y_ref) 
         << fixed << setprecision(2) << setw(7)
         << -log10(Norm(y-y_ref)) << endl;
  
  };

  return 0;
  
}
