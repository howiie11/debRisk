//------------------------------------------------------------------------------
//
// Exercise_4_3.cpp
// 
// Purpose: 
//
//   Satellite Orbits - Models, Methods and Applications
//   Exercise 4-3: Step size control of DE multistep method
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

// Record for passing global data between f_Kep6D_ and the calling program 

struct AuxDataRecord {
  int     n_step;
  double  t;
};

//------------------------------------------------------------------------------
//
// f_Kep6D_
//
// Purpose:
// 
//   Computes the derivative of the state vector for the normalized (GM=1)
//   Kepler's problem in three dimensions
//
// Note:
//
//   pAux is expected to point to a variable of type AuxDataRecord, which is
//   used to communicate with the other program sections and to hold data 
//   between subsequent calls of this function
//
//------------------------------------------------------------------------------

void f_Kep6D_ ( double t, const Vector& y, Vector& yp, void* pAux )
{
   
  // State vector derivative

  Vector r = y.slice(0,2);
  Vector v = y.slice(3,5);
  yp = Stack ( v, -r/(pow(Norm(r),3)) );
  
  // Pointer to auxiliary data record
  AuxDataRecord* p = static_cast<AuxDataRecord*>(pAux);

  // Write current time, step size and radius; store time for next step
  if (t-(*p).t>1.0e-10) {
    (*p).n_step++;
    cout << setw(5) << (*p).n_step
         << fixed << setprecision(6) << setw(12) << t
         << fixed << setprecision(6) << setw(12) << t-(*p).t
         << fixed << setprecision(3) << setw(10) << Norm(r)
         << endl;
    (*p).t = t;
  };
  
};


//------------------------------------------------------------------------------
//
// Main program
//
//------------------------------------------------------------------------------

int main() {

  // Constants
  
  const double  GM    = 1.0;                   // Gravitational coefficient
  const double  e     = 0.9;                   // Eccentricity
  const double  t_end = 20.0;                  // End time
  const Vector  Kep(1.0,e,0.0,0.0,0.0,0.0);    // (a,e,i,Omega,omega,M)
  const Vector  y_ref = State(GM,Kep,t_end);   // Reference solution

  // Variables
  
  double         t;                            // Time
  double         relerr,abserr;                // Accuracy requirements
  Vector         y(6);                         // State vector
  AuxDataRecord  Aux;                          // Auxiliary data
  DE             Orbit(f_Kep6D_,6,&Aux);       // Object for integrating the
                                               // differential equation
                                               // defined by f_Kep6D_
  // Header 

  cout << "Exercise 4-3: Step size control of DE multistep method" 
       << endl << endl;
  cout << "  Step       t           h          r " << endl;
    
  // Initial values
  t = 0.0;
  y = Vector( 1.0-e,0.0,0.0, 0.0,sqrt((1+e)/(1-e)),0.0 );
  
  Aux.n_step = 0;
  Aux.t      = 0;

  relerr = 0.0;
  abserr = 1.0e-8;

  // Integration from t=t to t=t_end
  Orbit.Init (t,relerr,abserr);
  Orbit.Integ (t_end,y);

  return 0;

}
