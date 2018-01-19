//------------------------------------------------------------------------------
//
// Exercise_5_1.cpp
// 
// Purpose: 
//
//   Satellite Orbits - Models, Methods, and Applications
//   Exercise 5-1: Transformation from celestial to terrestrial reference system
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
  
  double MJD_UTC;    // Modified Julian Date (UTC)
  double MJD_UT1;    // Modified Julian Date (UTC)
  double MJD_TT;     // Modified Julian Date (TT)
  Matrix P(3,3);     // Precession matrix (ICRS -> mean-of-date)
  Matrix N(3,3);     // Nutation matrix (mean-of-date -> true-of-date)
  Matrix Theta(3,3); // Sidereal Time matrix (tod -> pseudo-Earth-fixed)
  Matrix Pi(3,3);    // Polar motion matrix (pseudo-Earth-fixed -> ITRS)


  // Header 

  cout << "Exercise 5-1: Transformation from celestial "
       << "to terrestrial reference system" 
       << endl << endl;
    

  // Earth Orientation Parameters (UT1-UTC[s],UTC-TAI[s], x["], y["])
  // (from IERS Bulletin B #135 and C #16; valid for 1999/03/04 0:00 UTC)

  IERS::Set ( 0.649232, -32.0, 0.06740, 0.24173 ); 


  // Date

  MJD_UTC = Mjd ( 1999,03,04, 0,0,0.0 );
  MJD_UT1 = MJD_UTC + IERS::UT1_UTC(MJD_UTC)/86400.0;
  MJD_TT  = MJD_UTC + IERS::TT_UTC(MJD_UTC)/86400.0;
  
  // IAU 1976 Precession
  // (ICRF to mean equator and equinox of date)

  P = PrecMatrix(MJD_J2000,MJD_TT);

  // IAU 1980 Nutation
  // (Transformation to the true equator and equinox)

  N = NutMatrix(MJD_TT);

  // Apparent Sidereal Time
  // Rotation about the Celestial Ephemeris Pole

  Theta = GHAMatrix(MJD_UT1);   // Note: here we evaluate the equation of the
                                // equinoxes with the MJD_UT1 time argument 
                                // (instead of MJD_TT)
                                
  // Polar motion
  // (Transformation from the CEP to the IRP of the ITRS)

  Pi = PoleMatrix(MJD_UTC);     // Note: the time argument of polar motion series
                                // is not rigorously defined, but any differences
                                // are negligible

  // Output
  
  cout << "Date" << endl
       << " " << Date(MJD_UTC) << " UTC" << endl
       << " " << Date(MJD_UT1) << " UT1" << endl
       << " " << Date(MJD_TT ) << " TT " << endl << endl << endl;

  cout << "IAU 1976 Precession matrix (ICRS to tod)" << endl
       << setprecision(8) << fixed << showpos << setw(12) 
       << P << endl;

  cout << "IAU 1980 Nutation matrix (tod to mod)" << endl
       << setprecision(8) << fixed << showpos << setw(12) 
       << N << endl;

  cout << "Earth Rotation matrix" << endl
       << setprecision(8) << fixed << showpos << setw(12) 
       << Theta << endl;

  cout << "Polar motion matrix" << endl
       << setprecision(8) << fixed << showpos << setw(12) 
       << Pi << endl << endl;
       
  cout << "ICRS-ITRS transformation" << endl
       << setprecision(8) << fixed << showpos << setw(12) 
       << Pi*Theta*N*P;

  return 0;

}
