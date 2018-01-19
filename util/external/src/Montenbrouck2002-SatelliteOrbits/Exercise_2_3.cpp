//------------------------------------------------------------------------------
//
// Exercise_2_3.cpp
// 
// Purpose: 
//
//   Satellite Orbits - Models, Methods, and Applications
//   Exercise 2-3: Osculating elements
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
#include "SAT_Kepler.h"
#include "SAT_RefSys.h"
#include "SAT_Time.h"
#include "SAT_VecMat.h"

using namespace std;


//------------------------------------------------------------------------------
//
// Main program
//
//------------------------------------------------------------------------------

int main() {

  // Position and velocity 

  const Vector  r(+10000.0e3,+40000.0e3, -5000.0e3);  // [m]
  const Vector  v(  -1.500e3,  +1.000e3,  -0.100e3);  // [m/s]

  // Variables

  int     i;
  Vector  y(6), Kep(6); 

  // Orbital elements

  y = Stack(r,v);

  Kep = Elements ( GM_Earth, y );

  // Output

  cout << "Exercise 2-3: Osculating elements" << endl << endl;

  cout << "State vector:" << endl << endl;
  cout << "  Position       " << fixed << setprecision(3);
  for (i=0; i<3; i++) { cout << setw(12) << r(i)/1000.0; }; 
  cout << "  [km]" << endl;
  cout << "  Velocity       " << setprecision(6);
  for (i=0; i<3; i++) { cout << setw(12) << v(i)/1000.0; };
  cout << "  [km/s]" << endl;

  cout << endl;
  
  cout << "Orbital elements:" << endl << endl
       << setprecision(3)
       << "  Semimajor axis   " << setw(10) << Kep(0)/1000.0 << " km" << endl
       << setprecision(7)
       << "  Eccentricity     " << setw(10) << Kep(1)<< endl
       << setprecision(3)
       << "  Inclination      " << setw(10) << Kep(2)*Deg << " deg"<< endl
       << "  RA ascend. node  " << setw(10) << Kep(3)*Deg << " deg"<< endl
       << "  Arg. of perigee  " << setw(10) << Kep(4)*Deg << " deg"<< endl
       << "  Mean anomaly     " << setw(10) << Kep(5)*Deg << " deg"<< endl
       << endl;

  return 0;

}
