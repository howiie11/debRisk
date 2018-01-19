//------------------------------------------------------------------------------
//
// Exercise_2_4.cpp
// 
// Purpose: 
//
//   Satellite Orbits - Models, Methods, and Applications
//   Exercise 2-4: Topocentric satellite motion
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

  const double   lon_Sta = 11.0*Rad;           // [rad]
  const double   lat_Sta = 48.0*Rad;           // [rad]
  const double   alt_h   =  0.0e3;             // [m]
  
  const Geodetic Sta(lon_Sta,lat_Sta,alt_h);   // Geodetic coordinates
  

  // Spacecraft orbit

  const double   Mjd_Epoch = Mjd(1997,01,01);  // Epoch

  const double   a     = 960.0e3 + R_Earth;    // Semimajor axis [m]
  const double   e     =   0.0;                // Eccentricity
  const double   i     =  97.0*Rad;            // Inclination [rad]
  const double   Omega = 130.7*Rad;            // RA ascend. node [rad]
  const double   omega =   0.0*Rad;            // Argument of latitude [rad]
  const double   M0    =   0.0*Rad;            // Mean anomaly at epoch [rad]

  const Vector   Kep(a,e,i,Omega,omega,M0);    // Keplerian elements

  // Variables

  double    Mjd_UTC, dt;
  double    Azim, Elev, Dist;
  Vector    R_Sta(3), s(3);
  Vector    r(3); 
  Matrix    E(3,3), U(3,3);


  // Station

  R_Sta = Sta.Position(R_Earth,f_Earth);        // Geocentric position vector
  E     = Sta.LTC_Matrix();                     // Transformation to 
                                                // local tangent coordinates
  // Header

  cout << "Exercise 2-4: Topocentric satellite motion" << endl << endl
       << "   Date         UTC           Az         El      Dist" << endl
       << "yyyy/mm/dd  hh:mm:ss.sss     [deg]     [deg]     [km]" << endl;
  
  // Orbit

  for (int Minute=6; Minute<=24; Minute++) {
  
    Mjd_UTC = Mjd_Epoch + Minute/1440.0;        // Time

    dt = (Mjd_UTC-Mjd_Epoch)*86400.0;           // Time since epoch [s] 

    r = State(GM_Earth,Kep,dt).slice(0,2);      // Inertial position vector
    
    U = R_z(GMST(Mjd_UTC));                     // Earth rotation 
    s = E * ( U*r - R_Sta );                    // Topocentric position vector

    Dist = Norm(s);                             // Distance
    AzEl ( s, Azim, Elev );                     // Azimuth, Elevation

    cout << Date(Mjd_UTC) 
         << fixed << setprecision(1) 
         << setw(10) << Azim*Deg << setw(10) << Elev*Deg
         << setw(10) << Dist/1000.0 << endl;

  };

  return 0;

}
