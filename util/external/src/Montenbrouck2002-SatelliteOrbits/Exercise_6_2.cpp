//------------------------------------------------------------------------------
//
// Exercise_6_2.cpp
// 
// Purpose: 
//
//   Satellite Orbits - Models, Methods, and Applications 
//   Exercise 6-2: Range Rate Modelling
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
  
  // Earth rotation

  const Vector   omega_vec(0.0,0.0,omega_Earth); // Earth rotation vector

  // Spacecraft orbit

  const double   Mjd_Epoch = Mjd(1997,01,01);    // Epoch

  const double   a     = 960.0e3 + R_Earth;      // Semimajor axis [m]
  const double   e     =   0.0;                  // Eccentricity
  const double   i     =  97.0*Rad;              // Inclination [rad]
  const double   Omega = 130.7*Rad;              // RA ascend. node [rad]
  const double   omega =   0.0*Rad;              // Argument of latitude [rad]
  const double   M0    =   0.0*Rad;              // Mean anomaly at epoch [rad]

  const Vector   Kep(a,e,i,Omega,omega,M0);      // Keplerian elements


  // Radar Modelling

  const int      I_max = 3;                      // Maximum light time iterations
  const double   Count = 1.0;                    // Doppler count time [s]


  // Variables

  int       Iteration, Step;                     // Loop counters
  double    Mjd_UTC, t;                          // Time
  double    rho;                                 // Range 1-way
  double    range1,range0;                       // Range 2-way at end, begin of count
  double    range_rate;                          // Range rate
  double    Doppler;                             // Instantaneous Doppler
  double    tau_up,tau_down;                     // Upleg/downleg light time
  double    rho_up, rho_down;                    // Upleg/downleg range
  Vector    R_Sta(3);                            // Earth-fixed station position
  Vector    r_Sta(3);                            // Inertial station position
  Vector    r(3);                                // Inertial satellite position
  Vector    x(3),v(3);                           // Earth-fixed satellite position, velocity 
  Vector    u(3);                                // Unit vector satellite station
  Matrix    U(3,3);                              // Earth rotation matrix


  // Station

  R_Sta = Sta.Position(R_Earth,f_Earth);         // Geocentric position vector


  // Header

  cout << "Exercise 6-2: Range Rate Modelling" << endl << endl
       << "   Date         UTC       Range Rate     Doppler  Difference" << endl
       << "yyyy/mm/dd  hh:mm:ss.sss      [m/s]       [m/s]       [m/s] " << endl;
  
  // Orbit 

  for (Step=0; Step<=6; Step++) {
  
    // Ground-received time

    t = 360.0 + 180.0*Step;                      // Time since epoch [s]
    
    Mjd_UTC = Mjd_Epoch + t/86400.0;             // Modified Julian Date [UTC]

    U = R_z(GMST(Mjd_UTC));                      // Earth rotation matrix

    r_Sta = Transp(U)*R_Sta;                     // Inertial station position

    // Light time iteration at count interval end for downleg satellite -> station

    tau_down = 0.0;
    for (Iteration=0;Iteration<=I_max;Iteration++) {
      r = State(GM_Earth,Kep,t-tau_down).slice(0,2);  // Spacecraft position
      rho = Norm(r-r_Sta);                            // Downleg range
      tau_down = rho/c_light;                         // Downleg light time
      rho_down = rho;
    };

    // Light time iteration at count interval end for upleg station -> satellite 

    tau_up = 0.0;
    for (Iteration=0;Iteration<=I_max;Iteration++) {
      U = R_z(GMST(Mjd_UTC-(tau_down+tau_up)/86400.0));
      r_Sta = Transp(U)*R_Sta;                        // Inertial station pos.
      rho = Norm(r-r_Sta);                            // at ground transmit time
      tau_up = rho/c_light;                           // Upleg light time
      rho_up = rho;
    };

    // Two-way range at end of count interval 

    range1 = 0.5 * ( rho_down + rho_up ); 

    // Station position at begin of count interval

    U = R_z(GMST(Mjd_UTC-Count/86400.0));        // Earth rotation matrix

    r_Sta = Transp(U)*R_Sta;                     // Inertial station position

    // Light time iteration at count interval begin for downleg satellite -> station

    tau_down = 0.0;
    for (Iteration=0;Iteration<=I_max;Iteration++) {
      r = State(GM_Earth,Kep,t-tau_down-Count).slice(0,2);  // Spacecraft position
      rho = Norm(r-r_Sta);                            // Downleg range
      tau_down = rho/c_light;                         // Downleg light time
      rho_down = rho;
    };

    // Light time iteration at count interval begin for upleg station -> satellite 

    tau_up = 0.0;
    for (Iteration=0;Iteration<=I_max;Iteration++) {
      U = R_z(GMST(Mjd_UTC-(tau_down+tau_up+Count)/86400.0));
      r_Sta = Transp(U)*R_Sta;                        // Inertial station pos.
      rho = Norm(r-r_Sta);                            // at ground transmit time
      tau_up = rho/c_light;                           // Upleg light time
      rho_up = rho;
    };

    // Two-way range at begin of count interval 

    range0 = 0.5 * ( rho_down + rho_up ); 

    // Two-way average range rate

    range_rate = (range1 - range0)/Count;


    // Instantaneous Doppler modelling at mid of count interval 
 
    U = R_z(GMST(Mjd_UTC-(Count/2.0)/86400));            // Earth rotation matrix
    x = U*State(GM_Earth,Kep,t-(Count/2.0)).slice(0,2);  // Spacecraft position
    v = U*State(GM_Earth,Kep,t-(Count/2.0)).slice(3,5)   // Spacecraft velocity
        - Cross(omega_vec,x);
    u = (x-R_Sta)/Norm(x-R_Sta);                         // Unit vector s/c-station 

    Doppler = Dot(v,u);                                  // Instantaneous Doppler


    // Output

    cout << Date(Mjd_UTC) 
         << fixed << setprecision(3) 
         << setw(12) << range_rate 
         << setw(12) << Doppler 
         << setw(12) << range_rate - Doppler << endl;

  };

  return 0;

}
