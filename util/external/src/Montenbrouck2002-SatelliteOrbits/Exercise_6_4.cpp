//------------------------------------------------------------------------------
//
// Exercise_6_4.cpp
// 
// Purpose: 
//
//   Satellite Orbits - Models, Methods, and Applications 
//   Exercise 6-4: Tropospheric Refraction
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
//   2000/03/04  EGO  Final version (1st edition)
//   2005/04/14  OMO  Final version (2nd reprint)
//
// (c) 1999-2005  O. Montenbruck, E. Gill
//
//------------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <iomanip>

#include "SAT_Const.h"
#include "SAT_Kepler.h"
#include "SAT_RefSys.h"
#include "SAT_Time.h"
#include "SAT_VecMat.h"

#include "GNU_iomanip.h"

using namespace std;


//------------------------------------------------------------------------------
//
// Main program
//
//------------------------------------------------------------------------------

int main() {

  // Ground station

  const double   lon_Sta = 11.0*Rad;             // [rad]
  const double   lat_Sta = 48.0*Rad;             // [rad]
  const double   alt_h   =  0.0e3;               // [m]
  
  const Geodetic Sta(lon_Sta,lat_Sta,alt_h);     // Geodetic coordinates
  

  // Fixed media data at ground site

  const double   T0 =  273.2;                    // Temperature at 0 deg C [K]
  const double   pa = 1024.0;                    // Partial pressure of dry air [mb]
  const double   fh =    0.7;                    // Relative humidity 

  // Spacecraft orbit

  const double   Mjd_Epoch = Mjd(1997,01,01);    // Epoch

  const double   a     = 42164.0e3;              // Semimajor axis [m]
  const double   e     =   0.000296;             // Eccentricity
  const double   i     =  0.05*Rad;              // Inclination [rad]
  const double   Omega = 150.7*Rad;              // RA ascend. node [rad]
  const double   omega =   0.0*Rad;              // Argument of latitude [rad]
  const double   M0    =   0.0*Rad;              // Mean anomaly at epoch [rad]

  const Vector   Kep(a,e,i,Omega,omega,M0);      // Keplerian elements


  // Variables

  int       Hour;                                
  double    Mjd_UTC, dt;                         
  double    Azim, Elev, Elev0=0, Dist;
  double    Ns, eh, T, TC;
  Vector    dElev(2);                            
  Vector    R_Sta(3);                            
  Vector    r(3), s(3);                                
  Matrix    U(3,3), E(3,3);                       

  double    Tv[2]   = { 303.0, 283.0 };          // Temperature [K] 
  
  // Station

  R_Sta = Sta.Position(R_Earth,f_Earth);         // Geocentric position vector
  E     = Sta.LTC_Matrix();                      // Transformation to 
                                                 // local tangent coordinates


  // Header

  cout << "Exercise 6-4: Tropospheric Refraction" << endl << endl;
  

  // Orbit 

  for (Hour=0; Hour<=8; Hour++) {
           
    Mjd_UTC = Mjd_Epoch + 3.0*Hour/24.0;         // Modified Julian Date [UTC]

    dt = (Mjd_UTC-Mjd_Epoch)*86400.0;            // Time since epoch [s] 

    r = State(GM_Earth,Kep,dt).slice(0,2);       // Inertial position vector
    
    U = R_z(GMST(Mjd_UTC));                      // Earth rotation 
    s = E * ( U*r - R_Sta );                     // Topocentric position vector

    Dist = Norm(s);                              // Distance
    AzEl ( s, Azim, Elev );                      // Azimuth, Elevation

    if (Hour==0) {
      Elev0 = Elev;                              // Store initial elevation
      cout << "E0 [deg] " << fixed << setprecision(3) << setw(10) << Elev0*Deg << endl << endl
       << "   Date         UTC          E-E0      dE_t1     dE_t2 " << endl
       << "yyyy/mm/dd  hh:mm:ss.sss     [deg]     [deg]     [deg]" << endl;
    };
      
    for (int Ti=0; Ti<=1; Ti++) {                // Evaluate at 2 temperatures
      T  = Tv[Ti];                               // Map to scalar 
      TC = T-T0;                                 // Temperature [C]
      eh = 6.10*fh*exp(17.15*TC/(234.7+TC));     // Partial water pressure
      Ns = 77.64*pa/T + 3.734e5*eh/(T*T);        // Refractivity
      dElev(Ti) = Ns*1.0e-6/tan(Elev);           // Tropospheric refraction
    };


    cout << Date(Mjd_UTC) 
         << fixed << setprecision(3) 
         << setw(10) << (Elev-Elev0)*Deg << setw(10) << dElev(0)*Deg
         << setw(10) << dElev(1)*Deg << endl;

  };

  return 0;

}
