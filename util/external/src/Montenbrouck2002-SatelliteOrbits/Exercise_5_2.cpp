//------------------------------------------------------------------------------
//
// Exercise_5_2.cpp
// 
// Purpose: 
//
//   Satellite Orbits - Models, Methods, and Applications
//   Exercise 5-2: Velocity in the Earth-fixed frame
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


  // Variables
  
  int    i;                  // Loop counter
  double MJD_GPS,MJD_TT;     // Modified Julian Date (GPS,TT)
  double MJD_UTC,MJD_UT1;    // Modified Julian Date (UTC,UT1)
  Matrix P(3,3),N(3,3);      // Precession/nutation matrix 
  Matrix Theta(3,3);         // Sidereal Time matrix 
  Matrix S(3,3),dTheta(3,3); // and derivative
  Matrix Pi(3,3);            // Polar motion matrix 
  Matrix U(3,3),dU(3,3);     // ICRS to ITRS transformation and derivative
  Vector r_WGS(3),v_WGS(3);  // Position/velocity in the Earth-fixed frame
  Vector r(3),v(3);          // Position/velocity in the ICRS 
  Vector y(6),Kep(6);        // Satte vector and Keplerian elements
  

  // Header 

  cout << "Exercise 5-2: Velocity in the Earth-fixed frame"
       << endl << endl;
    

  // Earth Orientation Parameters (UT1-UTC[s],UTC-TAI[s], x["], y["])
  // (from IERS Bulletin B #135 and C #16; valid for 1999/03/04 0:00 UTC)

  IERS::Set ( 0.649232, -32.0, 0.06740, 0.24173 ); 

  // Date

  MJD_GPS = Mjd ( 1999,03,04, 0,0,0.0 );

  MJD_UTC = MJD_GPS - IERS::GPS_UTC(MJD_GPS)/86400.0;
  MJD_UT1 = MJD_UTC + IERS::UT1_UTC(MJD_UTC)/86400.0;
  MJD_TT  = MJD_UTC + IERS::TT_UTC(MJD_UTC)/86400.0; 

  // Earth-fixed state vector of GPS satellite #PRN15
  // (from NIMA ephemeris nim09994.eph; WGS84(G873) system)
  
  r_WGS = Vector(19440.953805e+3,16881.609273e+3, -6777.115092e+3); // [m]
  v_WGS = Vector(-8111.827456e-1,-2573.799137e-1,-30689.508125e-1); // [m/s]


  // ICRS to ITRS transformation matrix and derivative

  P      = PrecMatrix(MJD_J2000,MJD_TT);    // IAU 1976 Precession
  N      = NutMatrix(MJD_TT);               // IAU 1980 Nutation
  Theta  = GHAMatrix(MJD_UT1);              // Earth rotation
  Pi     = PoleMatrix(MJD_UTC);             // Polar motion
  
  S(0,1) = 1.0; S(1,0) = -1.0;              // Derivative of Earth rotation 
  dTheta = omega_Earth*S*Theta;             // matrix [1/s]
                               
  U      = Pi*Theta*N*P;                    // ICRS to ITRS transformation
  dU     = Pi*dTheta*N*P;                   // Derivative [1/s]
 
  // Transformation from WGS to ICRS

  r = Transp(U)*r_WGS;
  v = Transp(U)*v_WGS + Transp(dU)*r_WGS;

  // Orbital elements

  y   = Stack(r,v);
  Kep = Elements ( GM_Earth, y );


  // Output
  
  cout << "Date" << endl << endl
       << " " << Date(MJD_GPS) << " GPS" << endl
       << " " << Date(MJD_UTC) << " UTC" << endl
       << " " << Date(MJD_UT1) << " UT1" << endl
       << " " << Date(MJD_TT ) << " TT " << endl << endl;
       
  cout << "WGS84 (G873) State vector:" << endl << endl;
  cout << "  Position       " << fixed << setprecision(3);
  for (i=0; i<3; i++) { cout << setw(12) << r_WGS(i)/1000.0; }; 
  cout << "  [km]" << endl;
  cout << "  Velocity       " << setprecision(6);
  for (i=0; i<3; i++) { cout << setw(12) << v_WGS(i)/1000.0; };
  cout << "  [km/s]" << endl << endl;

  cout << "ICRS-ITRS transformation" << endl << endl
       << setprecision(8) << fixed << showpos << setw(12) << U << endl;

  cout << "Derivative of ICRS-ITRS transformation [10^(-4)/s]" << endl << endl
       << setprecision(8) << fixed << showpos << setw(12) << dU*1.0e4 << endl;

  cout << "ICRS State vector:" << endl << endl;
  cout << "  Position       " << fixed << setprecision(3);
  for (i=0; i<3; i++) { cout << setw(12) << r(i)/1000.0; }; 
  cout << "  [km]" << endl;
  cout << "  Velocity       " << setprecision(6);
  for (i=0; i<3; i++) { cout << setw(12) << v(i)/1000.0; };
  cout << "  [km/s]" << endl << endl;
  
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
