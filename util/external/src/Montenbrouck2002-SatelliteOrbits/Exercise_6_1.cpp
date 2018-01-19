//------------------------------------------------------------------------------
//
// Exercise_6_1.cpp
// 
// Purpose: 
//
//   Satellite Orbits - Models, Methods, and Applications 
//   Exercise 6-1: Light Time Iteration
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

  const double   lon_Sta = 11.0*Rad;             // [rad]
  const double   lat_Sta = 48.0*Rad;             // [rad]
  const double   alt_h   =  0.0e3;               // [m]
  
  const Geodetic Sta(lon_Sta,lat_Sta,alt_h);     // Geodetic coordinates
  

  // Spacecraft orbit

  const double   Mjd_Epoch = Mjd(1997,01,01);    // Epoch

  const double   a     = 960.0e3 + R_Earth;      // Semimajor axis [m]
  const double   e     =   0.0;                  // Eccentricity
  const double   i     =  97.0*Rad;              // Inclination [rad]
  const double   Omega = 130.7*Rad;              // RA ascend. node [rad]
  const double   omega =   0.0*Rad;              // Argument of latitude [rad]
  const double   M0    =   0.0*Rad;              // Mean anomaly at epoch [rad]

  const Vector   Kep(a,e,i,Omega,omega,M0);      // Keplerian elements


  // Light time iteration 

  const int      I_max = 2;                      // Maxim. number of iterations


  // Variables

  int       Iteration, Step;                     // Loop counters
  double    Mjd_UTC, t;                          // Time
  double    rho,range;                           // Range 1-way/2-way
  double    tau_up,tau_down;                     // Upleg/downleg light time
  Vector    R_Sta(3);                            // Earth-fixed station position
  Vector    r_Sta(3);                            // Inertial station position
  Vector    r(3);                                // Inertial satellite position 
  Vector    rho_up(I_max+1), rho_down(I_max+1);  // Upleg/downleg range
  Matrix    U(3,3);                              // Earth rotation matrix


  // Station

  R_Sta = Sta.Position(R_Earth,f_Earth);         // Geocentric position vector


  // Header

  cout << "Exercise 6-1: Light time iteration" << endl << endl
       << "   Date         UTC        Distance   " << 
          "Down It 1   It 2   Up It 1    Range" << endl
       << "yyyy/mm/dd  hh:mm:ss.sss      [m]     " << 
          "    [m]     [mm]      [m]      [m] " << endl;
  
  // Orbit 

  for (Step=0; Step<=6; Step++) {
  
    // Ground-received time

    t = 360.0 + 180.0*Step;                      // Time since epoch [s]
    
    Mjd_UTC = Mjd_Epoch + t/86400.0;             // Modified Julian Date [UTC]

    U = R_z(GMST(Mjd_UTC));                      // Earth rotation matrix

    r_Sta = Transp(U)*R_Sta;                     // Inertial station position

    // Light time iteration for downleg satellite -> station

    tau_down = 0.0;
    for (Iteration=0;Iteration<=I_max;Iteration++) {
      r = State(GM_Earth,Kep,t-tau_down).slice(0,2);  // Spacecraft position
      rho = Norm(r-r_Sta);                            // Downleg range
      tau_down = rho/c_light;                         // Downleg light time
      rho_down(Iteration) = rho;
    };

    // Light time iteration for upleg station -> satellite 

    tau_up = 0.0;
    for (Iteration=0;Iteration<=I_max;Iteration++) {
      U = R_z(GMST(Mjd_UTC-(tau_down+tau_up)/86400.0));
      r_Sta = Transp(U)*R_Sta;                        // Inertial station pos.
      rho = Norm(r-r_Sta);                            // at ground transmit time
      tau_up = rho/c_light;                           // Upleg light time
      rho_up(Iteration) = rho;
    };

    // Two-way range 

    range = 0.5 * ( rho_down(I_max) + rho_up(I_max) ); 

    cout << Date(Mjd_UTC) 
         << fixed << setprecision(1) 
         << setw(12) << rho_down(0)                   // Geometric range
         << setw( 9) << rho_down(1)-rho_down(0)   
         << setw( 9) <<(rho_down(2)-rho_down(1))*1000.0 
         << setw( 9) << rho_up(1)-rho_up(0)   
         << setw(12) << range << endl;

  };

  return 0;

}
