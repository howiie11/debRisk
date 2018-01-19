//------------------------------------------------------------------------------
//
// Exercise_2_6.cpp
// 
// Purpose: 
//
//   Satellite Orbits - Models, Methods, and Applications
//   Exercise 2-6: Initial orbit determination
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

  // Ground station

  const Vector   R_Sta(+1344.143e3,+6068.601e3,+1429.311e3);  // Position vector
  const Geodetic Sta(R_Sta,R_Earth,f_Earth);             // Geodetic coordinates

  // Observations

  struct ObsType {
    double Mjd_UTC;
    double Azim, Elev, Range;
  };

  const ObsType  Obs[2] = {  { Mjd ( 1999,04,02, 00,30,00.0 ),
                               132.67*Rad, 32.44*Rad, 16945.450e3 },
                             { Mjd ( 1999,04,02, 03,00,00.0 ),
                               123.08*Rad, 50.06*Rad, 37350.340e3 } };

  // Variables

  int       i,j;
  double    Az,El,d;
  Vector    s(3), Kep(6); 
  Matrix    E(3,3), U(3,3);
  Vector    r[2];

  // Transformation to local tangent coordinates

  E = Sta.LTC_Matrix();                                   

  // Convert observations

  for (i=0; i<2; i++) {
    // Earth rotation
    U = R_z(GMST(Obs[i].Mjd_UTC));
    // Topocentric position vector
    Az = Obs[i].Azim;  El = Obs[i].Elev;  d = Obs[i].Range;
    s = d*VecPolar(pi/2-Az,El); 
    // Inertial position vector
    r[i] = Transp(U)*(Transp(E)*s + R_Sta);                
  }
                                    
  // Orbital elements

  Kep = Elements ( GM_Earth, Obs[0].Mjd_UTC, Obs[1].Mjd_UTC, r[0],r[1] );


  // Output

  cout << "Exercise 2-6: Initial orbit determination" << endl << endl;

  cout << "Inertial positions:" << endl << endl
       << setw(36) << "[km]" << setw(12) << "[km]" << setw(12) << "[km]"
       << endl;
  for (i=0; i<2; i++) { 
    cout << "  " << Date(Obs[i].Mjd_UTC) << fixed << setprecision(3);
    for (j=0; j<3; j++) { cout << setw(12) << r[i](j)/1000.0; };
    cout << endl;
  };  
  cout << endl;
  
  cout << "Orbital elements:" << endl << endl
       << "  Epoch (1st obs.)  " << Date(Obs[0].Mjd_UTC)  << endl
       << fixed << setprecision(3)
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
