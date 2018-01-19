//------------------------------------------------------------------------------
//
// Exercise_5_3.cpp
// 
// Purpose: 
//
//   Satellite Orbits - Models, Methods, and Applications
//   Exercise 5-3: Geodetic coordinates
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
#include "SAT_VecMat.h"

using namespace std;


//------------------------------------------------------------------------------
//
// Main program
//
//------------------------------------------------------------------------------

int main() {

  // WGS84 geoid parameters
  
  const double R_WGS84     =   6378.137e3;      // Radius Earth [m]
  const double f_WGS84     = 1.0/298.257223563; // Flattening   

  // Variables
  
  Vector    R_Sta(3);   // Cartesian station coordinates
  Geodetic  Sta;        // Station coordinates

  // Header 

  cout << "Exercise 5-3: Geodetic coordinates" << endl << endl;
    
  // Coordinates of NIMA GPS station at Diego Garcia (WGS84(G873); epoch 1997.0)

  R_Sta = Vector ( 1917032.190, 6029782.349, -801376.113 );   // [m]

  // Geodetic coordinates

  Sta = Geodetic(R_Sta,R_WGS84,f_WGS84);
  
  // Output
  
  cout << "Cartesian station coordinates (WGS84) [m]" << endl << endl
       << setprecision(3) << fixed << showpos << setw(13) 
       << Sta.Position(R_WGS84,f_WGS84) << endl 
       << endl;

  cout << "Geodetic station coordinates (WGS84)" << endl << endl
       << " longitude " << setprecision(8)
                        << setw(12) << Deg*Sta.lon << " deg" << endl
       << " latitude  " << setw(12) << Deg*Sta.lat << " deg" << endl
       << " height    " << setprecision(3)
                        << setw(12) << Sta.h   << " m" << endl;

  return 0;
  
}
