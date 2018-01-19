//------------------------------------------------------------------------------
//
// Exercise_8_2.cpp
// 
// Purpose: 
//
//   Satellite Orbits - Models, Methods, and Applications 
//   Exercise 8-2: Least-squares orbit determination
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
#include <cmath>

#include "GNU_iomanip.h"

#include "SAT_Const.h"
#include "SAT_Kepler.h"
#include "SAT_Filter.h"
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

  // Constants

  const int     N_obs = 6;
  const double  Step  = 1200.0;

  const Vector  Null3D(0.0,0.0,0.0);
  
  const double  sigma_range = 10.0;        // [m]
  const double  sigma_angle = 0.01*Rad;    // [rad] (=36")

  const char*   Label[6] = { "x  [m]  ", "y  [m]  ", "z  [m]  ",
                             "vx [m/s]", "vy [m/s]", "vz [m/s]"  };
  
  // Variables

  int     i,iterat;
  double  Mjd0,t,MjdUTC,Theta;
  double  Azim,Elev,Dist;
  Vector  Y0_ref(6),Y0_apr(6),Y0(6),Y(6),r(3),R(3),s(3);
  Vector  dAds(3),dEds(3),dDds(3);
  Vector  dAdY0(6),dEdY0(6),dDdY0(6);
  Matrix  dYdY0(6,6),U(3,3),E(3,3);
  LSQ     OrbEst(6);
  Vector  dY0(6),SigY0(6);

  struct ObsType {
    double Mjd_UTC;
    double Azim, Elev, Dist;
  };

  ObsType Obs[N_obs];


  // Ground station

  R = Vector(+1344.0e3,+6069.0e3,1429.0e3);    // [m]
  E = Geodetic(R).LTC_Matrix();

  // Header 

  cout << "Exercise 8-2: Least-squares orbit determination" << endl << endl;


  // Generation of artificial observations from given epoch state 

  Mjd0 = Mjd(1995,03,30, 00,00,00.0);       // Epoch (UTC)

  Y0_ref = Vector(-6345.000e3, -3723.000e3,  -580.000e3,     // [m]
                  +2.169000e3, -9.266000e3, -1.079000e3 );   // [m/s]
  Y0 = Y0_ref;

  cout << "Measurements" << endl << endl
       << "     Date          UTC       Az[deg]   El[deg]   Range[km]" << endl;

  for (i=0; i<N_obs; i++) {

    // Time increment and propagation
    t      = (i+1)*Step;                    // Time since epoch [s]
    MjdUTC = Mjd0 + t/86400.0;              // Modified Julian Date
    TwoBody ( GM_Earth,Y0_ref,t, Y,dYdY0 ); // State vector
    
    // Topocentric coordinates
    Theta  = GMST(MjdUTC);                  // Earth rotation
    U = R_z(Theta);        
    r = Y.slice(0,2);
    s = E*(U*r-R);                          // Topocentric position [m]
    AzEl(s,Azim,Elev); Dist=Norm(s);        // Azimuth, Elevation, Range    
    
    // Observation record
    Obs[i].Mjd_UTC = MjdUTC;
    Obs[i].Azim = Azim;
    Obs[i].Elev = Elev;
    Obs[i].Dist = Dist;

    // Output
    cout << "  " << Date(MjdUTC) << fixed << setprecision(3)
         << setw(10) << Deg*Azim << setw(10) << Deg*Elev 
         << setw(12) << Dist/1000.0 << endl;

  };

  cout << endl;


  //
  // Orbit determination
  //

  Mjd0   = Mjd(1995,03,30, 00,00,00.0);       // Epoch (UTC)

  Y0_apr = Y0_ref + Vector(+10.0e3,-5.0e3,+1.0e3,-1.0,+3.0,-0.5);
  Y0     = Y0_apr;

  // Iteration
  
  for (iterat=1; iterat<=3; iterat++) {

    OrbEst.Init();

    cout << "Iteration Nr. " << iterat << endl << endl
       << "  Residuals:" << endl << endl
       << "     Date          UTC       Az[deg]   El[deg]  Range[m]" << endl;

    for (i=0; i<N_obs; i++) {

      // Time increment and propagation
      MjdUTC = Obs[i].Mjd_UTC;                // Modified Julian Date
      t      = (MjdUTC-Mjd0)*86400.0;         // Time since epoch [s]
      TwoBody ( GM_Earth,Y0,t, Y,dYdY0 );     // State vector
    
      // Topocentric coordinates
      Theta  = GMST(MjdUTC);                  // Earth rotation
      U = R_z(Theta);        
      r = Y.slice(0,2);
      s = E*(U*r-R);                          // Topocentric position [m]
    
      // Observations and partials
      AzEl(s,Azim,Elev,dAds,dEds);            // Azimuth, Elevation
      Dist=Norm(s); dDds=s/Dist;              // Range
    
      dAdY0 = Stack(dAds*E*U,Null3D)*dYdY0;
      dEdY0 = Stack(dEds*E*U,Null3D)*dYdY0;
      dDdY0 = Stack(dDds*E*U,Null3D)*dYdY0;
    
      // Accumulate least-squares system

      OrbEst.Accumulate(dAdY0,(Obs[i].Azim-Azim),sigma_angle/cos(Elev));
      OrbEst.Accumulate(dEdY0,(Obs[i].Elev-Elev),sigma_angle);
      OrbEst.Accumulate(dDdY0,(Obs[i].Dist-Dist),sigma_range);

      // Output
      cout << "  " << Date(MjdUTC) << fixed << setprecision(3)
           << setw(10) << Deg*(Obs[i].Azim-Azim) 
           << setw(10) << Deg*(Obs[i].Elev-Elev) 
           << setprecision(1)
           << setw(10) << Obs[i].Dist-Dist << endl;

    };

    // Solve least-squares system

    OrbEst.Solve(dY0); 
    SigY0 = OrbEst.StdDev(); 

    cout << endl << "  Correction:" << endl << endl
         << "  Pos" << setprecision(1) << setw(10) << dY0.slice(0,2) 
         << "  m  " << endl
         << "  Vel" << setprecision(4) << setw(10) << dY0.slice(3,5) 
         << "  m/s" << endl << endl;  

    // Correct epoch state

    Y0 = Y0 + dY0;

  };

  // Summary

  cout << "Summary:" << endl
       << "             a priori   correction      final        sigma" 
       << endl;
  for (i=0;i<6;i++) {
    cout << "  " << Label[i];
    if (i<3) cout << setprecision(1); else cout << setprecision(4);
    cout << setw(12) << Y0_apr(i)
         << setw(11) << Y0(i)-Y0_apr(i) 
         << setw(14) << Y0(i)
         << setw(11) << SigY0(i)
         << endl;
  }

  return 0;
  
}
