//------------------------------------------------------------------------------
//
// Exercise_8_3.cpp
// 
// Purpose: 
//
//   Satellite Orbits - Models, Methods, and Applications 
//   Exercise 8-3: Orbit Determination using Extended Kalman Filter
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

  const int     N_obs =   6;
  const double  Step  =  1200.0;

  const Vector  Null3D(0.0,0.0,0.0);
  
  const double  sigma_range = 10.0;        // [m]
  const double  sigma_angle = 0.01*Rad;    // [rad] (0.01 deg = 36")

  
  // Variables

  int     i;
  double  Mjd0,t,t_old,MjdUTC,Theta;
  double  Azim,Elev,Dist;
  Vector  Y0_true(6),Y_true(6),Y(6),Y_old(6);
  Vector  ErrY(6),SigY(6);
  Vector  r(3),R(3),s(3);
  Vector  dAds(3),dEds(3),dDds(3);
  Vector  dAdY(6),dEdY(6),dDdY(6);
  Matrix  U(3,3),E(3,3);
  Matrix  Phi(6,6),Phi_true(6,6), P(6,6);
  EKF     Filter(6);

  struct ObsType {
    double Mjd_UTC;
    double Azim, Elev, Dist;
  };

  ObsType Obs[N_obs];


  // Ground station

  R = Vector(+1344.0e3,+6069.0e3,1429.0e3);    // [m] Bangalore
  E = Geodetic(R).LTC_Matrix();

  // Header 

  cout << "Exercise 8-3: Sequential orbit determination" << endl << endl;


  // Generation of artificial observations from given epoch state 

  Mjd0 = Mjd(1995,03,30, 00,00,00.0);       // Epoch (UTC)

  Y0_true = Vector(-6345.000e3, -3723.000e3,  -580.000e3,     // [m]
                   +2.169000e3, -9.266000e3, -1.079000e3 );   // [m/s]

  cout << "Measurements" << endl << endl
       << "     Date          UTC       Az[deg]   El[deg]   Range[km]" << endl;

  for (i=0; i<N_obs; i++) {

    // Time increment and propagation
    t      = (i+1)*Step;                    // Time since epoch [s]
    MjdUTC = Mjd0 + t/86400.0;              // Modified Julian Date
    TwoBody ( GM_Earth,Y0_true,t, Y,Phi );  // State vector
    
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

  cout << "State errors" << endl << endl
       << "                                     Pos[m]             Vel[m/s]    " 
       << endl
       << "     Date          UTC    Upd.   Error     Sigma     Error     Sigma"
       << endl;

  // Initialization
  
  Mjd0 = Mjd(1995,03,30, 00,00,00.0);        // Epoch (UTC)
  
  t = 0.0;
  
  Y = Y0_true + Vector(+10.0e3,-5.0e3,+1.0e3,-1.0,+3.0,-0.5);
  
  P = 0.0;
  for (i=0; i<3; i++) P(i,i)=1.0e8;
  for (i=3; i<6; i++) P(i,i)=1.0e2;

  Filter.Init(t,Y,P);

  // Measurement loop

  for (i=0; i<N_obs; i++) {

    // Previous step
    t_old = Filter.Time();
    Y_old = Filter.State();

    // Propagation to measurement epoch
    MjdUTC = Obs[i].Mjd_UTC;                  // Modified Julian Date
    t      = (MjdUTC-Mjd0)*86400.0;           // Time since epoch [s]
    TwoBody ( GM_Earth,Y_old,t-t_old, Y,Phi); // State vector
    Theta  = GMST(MjdUTC);                    // Earth rotation
    U = R_z(Theta);        
    
    // Time update
    Filter.TimeUpdate(t,Y,Phi);

    // Truth orbit
    TwoBody( GM_Earth,Y0_true,t, Y_true,Phi_true );
    
    // State error and standard deviation
    ErrY = Filter.State()-Y_true;
    SigY = Filter.StdDev();
    cout << Date(MjdUTC) << "  t  " << fixed 
         << setprecision(1) << setw(10) << Norm(ErrY.slice(0,2))
         << setprecision(1) << setw(10) << Norm(SigY.slice(0,2))
         << setprecision(4) << setw(10) << Norm(ErrY.slice(3,5))
         << setprecision(4) << setw(10) << Norm(SigY.slice(3,5))
         << endl;

    // Azimuth and partials
    r = Filter.State().slice(0,2);
    s = E*(U*r-R);                            // Topocentric position [m]
    AzEl(s,Azim,Elev,dAds,dEds);              // Azimuth, Elevation
    dAdY = Stack(dAds*E*U,Null3D);

    // Measurement update
    Filter.MeasUpdate ( Obs[i].Azim, Azim, sigma_angle/cos(Elev), dAdY );
    ErrY = Filter.State()-Y_true;
    SigY = Filter.StdDev();
    cout << setw(29) << "  Az " << fixed 
         << setprecision(1) << setw(10) << Norm(ErrY.slice(0,2))
         << setprecision(1) << setw(10) << Norm(SigY.slice(0,2))
         << setprecision(4) << setw(10) << Norm(ErrY.slice(3,5))
         << setprecision(4) << setw(10) << Norm(SigY.slice(3,5))
         << endl;

    // Elevation and partials
    r = Filter.State().slice(0,2);
    s = E*(U*r-R);                            // Topocentric position [m]
    AzEl(s,Azim,Elev,dAds,dEds);              // Azimuth, Elevation
    dEdY = Stack(dEds*E*U,Null3D);
    
    // Measurement update
    Filter.MeasUpdate ( Obs[i].Elev, Elev, sigma_angle, dEdY );
    ErrY = Filter.State()-Y_true;
    SigY = Filter.StdDev();
    cout << setw(29) << "  El " << fixed 
         << setprecision(1) << setw(10) << Norm(ErrY.slice(0,2))
         << setprecision(1) << setw(10) << Norm(SigY.slice(0,2))
         << setprecision(4) << setw(10) << Norm(ErrY.slice(3,5))
         << setprecision(4) << setw(10) << Norm(SigY.slice(3,5))
         << endl;
    
    // Range and partials
    r = Filter.State().slice(0,2);
    s = E*(U*r-R);                            // Topocentric position [m]
    Dist=Norm(s); dDds=s/Dist;                // Range
    dDdY = Stack(dDds*E*U,Null3D);

    // Measurement update
    Filter.MeasUpdate ( Obs[i].Dist, Dist, sigma_range, dDdY );
    ErrY = Filter.State()-Y_true;
    SigY = Filter.StdDev();
    cout << setw(29) << "  rho" << fixed 
         << setprecision(1) << setw(10) << Norm(ErrY.slice(0,2))
         << setprecision(1) << setw(10) << Norm(SigY.slice(0,2))
         << setprecision(4) << setw(10) << Norm(ErrY.slice(3,5))
         << setprecision(4) << setw(10) << Norm(SigY.slice(3,5))
         << endl;
    
  };

  return 0;
  
}
