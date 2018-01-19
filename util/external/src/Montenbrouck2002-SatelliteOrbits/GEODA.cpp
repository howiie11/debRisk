//------------------------------------------------------------------------------
//
// GEODA.cpp
// 
// Purpose: 
//
//   Satellite Orbits - Models, Methods, and Applications 
//   Geostationary satellite Orbit Determination error Analysis
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
#include <fstream>
#include <cmath>

#include "SAT_Const.h"
#include "SAT_Filter.h"
#include "SAT_RefSys.h"
#include "SAT_Time.h"
#include "SAT_VecMat.h"

#include "GNU_iomanip.h"

using namespace std;


//------------------------------------------------------------------------------
//
// Global types and data
//
//------------------------------------------------------------------------------

struct Meas {
  double   sig_noise;    // Measurement accuracy
  double   sig_bias;     // Bias uncertainty
  double   step;         // Interval
  bool     est_bias;     // Flag for bias estimation
  int      index;        // Bias parameter index
};

struct Station {
  Geodetic Geod;         // Geodetic coordinates
  Meas     ang;          // Angle measurement parameters
  Meas     rng;          // Range measurement parameters
};


//------------------------------------------------------------------------------
//
// Phi
//
// Purpose:
//
//   State transition matrix 
//
//------------------------------------------------------------------------------

Matrix Phi(double n, double t) {

  const double c = cos(n*t);
  const double s = sin(n*t);

  Matrix T(6,6);

  T(0,0) = 4-3*c;      T(0,1) = 0;    T(0,3) = s/n;        T(0,4) = 2*(1-c)/n;
  T(1,0) = 6*(s-n*t);  T(1,1) = 1;    T(1,3) = 2*(c-1)/n;  T(1,4) = 4*s/n-3*t;
  T(2,2) = c;          T(2,5) = s/n;
 
  T(3,0) = 3*n*s;      T(3,1) = 0;    T(3,3) = c;          T(3,4) = 2*s;
  T(4,0) = 6*n*(c-1);  T(4,1) = 0;    T(4,3) = -2*s;       T(4,4) = 4*c-3;
  T(5,2) = -n*s;       T(5,5) = c;

  return T;

};


//------------------------------------------------------------------------------
//
// GetSetup
//
// Purpose:
// 
//   Reads GEODA setup parameters from file
//
//------------------------------------------------------------------------------

void GetSetup ( ifstream& inp,
                double& lon, Station& Sta1, Station& Sta2, 
                int& n_est, int& n_con,
                double& T_track, double& T_pred ) {

  // Variables

  int       k;
  double    lam,phi,h;

  // Input

  inp.ignore(50); inp >> lon;  inp.ignore(81,'\n');
  
  inp.ignore(50); inp >> lam >> phi >> h;  inp.ignore(81,'\n'); 
  Sta1.Geod=Geodetic(Rad*(lam-lon),Rad*phi,h);
  
  inp.ignore(50); inp >> Sta1.ang.sig_noise >> Sta1.ang.sig_bias 
                      >> Sta1.ang.step >> Sta1.ang.est_bias; 
  inp.ignore(81,'\n'); 
  
  inp.ignore(50); inp >> Sta1.rng.sig_noise >> Sta1.rng.sig_bias 
                      >> Sta1.rng.step >> Sta1.rng.est_bias; 
  inp.ignore(81,'\n'); 
  
  inp.ignore(50); inp >> lam >> phi >> h;  inp.ignore(81,'\n'); 
  Sta2.Geod=Geodetic(Rad*(lam-lon),Rad*phi,h);

  inp.ignore(50); inp >> Sta2.ang.sig_noise >> Sta2.ang.sig_bias 
                      >> Sta2.ang.step >> Sta2.ang.est_bias; 
  inp.ignore(81,'\n'); 
  
  inp.ignore(50); inp >> Sta2.rng.sig_noise >> Sta2.rng.sig_bias 
                      >> Sta2.rng.step >> Sta2.rng.est_bias; 
  inp.ignore(81,'\n'); 

  inp.ignore(50); inp >> T_track;   inp.ignore(81,'\n'); 
  inp.ignore(50); inp >> T_pred;    inp.ignore(81,'\n'); 

  // Change units

  lon*=Rad; 
  Sta1.ang.sig_noise*=Rad; Sta1.ang.sig_bias*=Rad; 
  Sta2.ang.sig_noise*=Rad; Sta2.ang.sig_bias*=Rad;
  
  // Define bias parameter indizes 
  // (estimation parameters first, consider parameters last)

  n_est =  6;   // Estimation parameter count (default: state vector only)
  n_con =  0;   // Consider parameter count
  k     =  0;   // Bias parameter index

  if (  Sta1.ang.est_bias ) { Sta1.ang.index=k; k+=2; n_est+=2; };
  if (  Sta1.rng.est_bias ) { Sta1.rng.index=k; k+=1; n_est+=1; };
  if (  Sta2.ang.est_bias ) { Sta2.ang.index=k; k+=2; n_est+=2; };
  if (  Sta2.rng.est_bias ) { Sta2.rng.index=k; k+=1; n_est+=1; };

  if ( !Sta1.ang.est_bias ) { Sta1.ang.index=k; k+=2; n_con+=2; };
  if ( !Sta1.rng.est_bias ) { Sta1.rng.index=k; k+=1; n_con+=1; };
  if ( !Sta2.ang.est_bias ) { Sta2.ang.index=k; k+=2; n_con+=2; };
  if ( !Sta2.rng.est_bias ) { Sta2.rng.index=k; k+=1; n_con+=1; };


}


//------------------------------------------------------------------------------
//
// Main program
//
//------------------------------------------------------------------------------

int main(int argc, char* argv[]) {

  // Constants
    
  const int     n_orb  = 6;              // Number of dynamical parameters
  const int     n_bias = 6;              // Number of bias parameters
  const int     n_dim  = n_orb + n_bias; // Dimension of least squares system

  const double  a = 42164000.0;          // Geostationary radius [m]
  const double  n = 2*pi*1.0027/24.0;    // Mean angular velocity [rad/h]

  const Vector  dady(4,0,0,0,2/n,0);     // Partials of semi-major axis
                                         // w.r.t. state vector

  // Variables

  ifstream inp;                          // Input file

  double   lon;                          // Subsatellite longitude
  Station  Sta1, Sta2;                   // Station specific parameters
  int      n_est;                        // Number of estimation parameters
  int      n_con;                        // Number of consider parameters
  double   T_track, T_pred;              // Tracking and prediction interval

  int      i_Sta;                        // Counter
  Station  Sta;                          // Station parameters
  double   t;                            // Time [h]
  double   Az,El,D;                      // Modelled observations
  Vector   r(3);                         // Position vector
  Vector   dAds(3),dEds(3);              // Partials wrt. topocentric 
  Vector   dAdr(3),dEdr(3),dDdr(3);      // and equatorial coordinates
  Vector   dhdy(n_orb);                  // Partials wrt. epoch state 
  Vector   dhdb(n_bias);                 // Partials wrt. bias parameters
  Vector   Var_b(n_bias);                // Bias parameter variance
  Matrix   E(3,3);                       // Transf. to topocentric coordinates
  Matrix   P_y0(n_orb,n_orb);            // State covariance matrix
  Matrix   P_y(n_orb,n_orb);             // State covariance matrix
  Matrix   T(n_orb,n_orb);               // State transition matrix
  Matrix   U(n_orb,n_orb);               // Covariance transformation
  double   P_a;                          // Semi-major axis covariance
  

  // Control parameters from setup file

  if (argc>1) inp.open(argv[1]); else inp.open("GEODA.inp"); 
  if (!inp) { 
    cerr << "ERROR: Could not open GEODA setup file" << endl;
    exit(1);
  };
  
  GetSetup ( inp, lon, Sta1,Sta2, n_est,n_con, T_track,T_pred );
  
  inp.close();
  
  // Dynamic variables

  LSQ    Lsq(n_dim);                     // Least squares system

  Matrix R_xx(n_est,n_est);              // Square root information matrix
  Matrix R_xx_inv(n_est,n_est);          // Inverse of SRIM
  Matrix P_x(n_est,n_est);               // Estimation parameter covariance
  Matrix R_xc(n_est,n_dim-n_est);        // Transformed consider param. partials
  Matrix C(n_dim-n_est,n_dim-n_est);     // Consider parameter covariance
  Matrix Pc_x(n_est,n_est);              // Consider covariance matrix

  // Bias parameter variance

  Var_b(Sta1.ang.index  ) = pow(Sta1.ang.sig_bias,2);  // Azimuth bias station 1
  Var_b(Sta1.ang.index+1) = pow(Sta1.ang.sig_bias,2);  // Elevation bias station 1
  Var_b(Sta1.rng.index  ) = pow(Sta1.rng.sig_bias,2);  // Range bias station 1
  Var_b(Sta2.ang.index  ) = pow(Sta2.ang.sig_bias,2);  // Azimuth bias station 2
  Var_b(Sta2.ang.index+1) = pow(Sta2.ang.sig_bias,2);  // Elevation bias station 2
  Var_b(Sta2.rng.index  ) = pow(Sta2.rng.sig_bias,2);  // Range bias station 2

  // Consider parameter covariance
  
  C = Diag(Var_b.slice(n_est-n_orb,n_bias-1)); 

  // Partial derivatives of inertial w.r.t. rotating state vector

  U = Id(6); U(4,0) = +n; U(3,1) = -n;
  
  // Header
  
  cout << "GEO Orbit Determination Error Analysis" << endl << endl
       << endl << endl
       << fixed << setprecision(3)
       << "  Setup parameters:" << endl << endl 
       << "  Subsatellite longitude [deg]                   :"
       << setw(8) << Deg*lon << endl
       << "  Station 1 (lon [deg], lat [deg], alt [m])      :"
       << setw(8) << Deg*(Sta1.Geod.lon+lon) << setw(9) << Deg*Sta1.Geod.lat 
       << setw(9) << Sta1.Geod.h << endl
       << "  Angles (noise & bias [deg], step [h], est.bias):"
       << setw(8) << Deg*Sta1.ang.sig_noise << setw(9) << Deg*Sta1.ang.sig_bias
       << setw(9) << Sta1.ang.step << setw(3) << Sta1.ang.est_bias << endl
       << "  Range (noise & bias [deg], step [h], est.bias) :" 
       << setw(8) << Sta1.rng.sig_noise << setw(9) << Sta1.rng.sig_bias
       << setw(9) << Sta1.rng.step << setw(3) << Sta1.rng.est_bias << endl
       << "  Station 2 (lon [deg], lat [deg], alt [m])      :"
       << setw(8) << Deg*(Sta2.Geod.lon+lon) << setw(9) << Deg*Sta2.Geod.lat 
       << setw(9) << Sta2.Geod.h << endl
       << "  Angles (noise & bias [deg], step [h], est.bias):"
       << setw(8) << Deg*Sta2.ang.sig_noise << setw(9) << Deg*Sta2.ang.sig_bias
       << setw(9) << Sta2.ang.step << setw(3) << Sta2.ang.est_bias << endl
       << "  Range (noise & bias [deg], step [h], est.bias) :" 
       << setw(8) << Sta2.rng.sig_noise << setw(9) << Sta2.rng.sig_bias
       << setw(9) << Sta2.rng.step << setw(3) << Sta2.rng.est_bias << endl
       << "  Tracking interval [h]                          :" 
       << setw(8) << T_track << endl
       << "  Prediction interval [h]                        :"
       << setw(8) << T_pred << endl
       << endl;

  
  // Initialization of least squares system
  
  Lsq.Init();

  // Process partials of measurements from both stations

  for (i_Sta=1;i_Sta<=2; i_Sta++) {
   
    // Select station

    if (i_Sta==1) Sta=Sta1; else Sta=Sta2; 

    // Azimuth, elevation, range and partials

    E = Sta.Geod.LTC_Matrix();
    r = Vector(a,0.0,0.0) - Sta.Geod.Position();
    AzEl(E*r,Az,El,dAds,dEds); D=Norm(r); 
    dAdr=dAds*E; dEdr=dEds*E; dDdr=r/D;

    // Angle measurements
    t = 0.0;
    if (Sta.ang.step>0.0) for(;;) {
      // Partials of azimuth measurement
      dhdy = dAdr*Phi(n,t).slice(0,2,0,5);       // Partials w.r.t. epoch state
      dhdb = 0.0;  dhdb(Sta.ang.index) = 1.0;    // Partials w.r.t. biases
      // Add partials to least squares system
      Lsq.Accumulate ( Stack(dhdy,dhdb), 0.0, Sta.ang.sig_noise  );
      // Partials of elevation measurement
      dhdy = dEdr*Phi(n,t).slice(0,2,0,5);       // Partials w.r.t. epoch state
      dhdb = 0.0;  dhdb(Sta.ang.index+1) = 1.0;  // Partials w.r.t. biases
      // Add partials to least squares system
      Lsq.Accumulate ( Stack(dhdy,dhdb), 0.0, Sta.ang.sig_noise );
      // Increment time of measurement and exit if past end of data arc
      t += Sta.ang.step;
      if (t>T_track) break;
    };

    // Range measurements
    t = 0.0;
    if (Sta.rng.step>0.0) for(;;) {
      // Partials of range measurement
      dhdy = dDdr*Phi(n,t).slice(0,2,0,5);       // Partials w.r.t. epoch state
      dhdb = 0.0; dhdb(Sta.rng.index) = 1.0;     // Partials w.r.t. biases
      Lsq.Accumulate ( Stack(dhdy,dhdb), 0.0, Sta.rng.sig_noise );
      // Increment time of measurements and exit if past end of data arc
      t += Sta.rng.step;
      if (t>T_track) break;
    };

  };


  // Decomposition of square root information matrix

  R_xx = Lsq.SRIM().slice(0,n_est-1,0    ,n_est-1); // Square root inf. matrix
                                                    // of estimated parameters
  R_xc = Lsq.SRIM().slice(0,n_est-1,n_est,n_dim-1); // Transformed partials  
                                                    // wrt. consider parameters

  // Computed covariance

  InvUpper(R_xx,R_xx_inv);                      // Inverse of SRIM
  P_x  = R_xx_inv*Transp(R_xx_inv);             // Estim. parameter covariance
  P_y0 = P_x.slice(0,5,0,5);                    // Epoch state covariance
  P_a  = Dot(dady*P_y0,dady);                   // Semi-major axis covariance
  P_y  = U*P_y0*Transp(U);                      // Inertial state covariance

  // Output

  cout << "Epoch state standard deviation (unmodelled effects ignored)" 
       << endl << endl;
  cout << "  Position (xyz)  " 
       << fixed << setprecision(2) << setw(10) 
       << P_y.Diag().Sqrt().slice(0,2) << " [m]" << endl
       << "  Velocity (xyz)  " 
       << fixed << setprecision(5) << setw(10) 
       << P_y.Diag().Sqrt().slice(3,5)/3600.0 << " [m/s]" << endl
       <<  endl;
  cout << "  Semi-major axis "
       << fixed << setprecision(2) << setw(10) 
       << sqrt(P_a) << " [m]" << endl << endl;

  // Consider covariance

  Pc_x = P_x + R_xx_inv*R_xc*C*Transp(R_xc)*Transp(R_xx_inv);
  P_y0 = Pc_x.slice(0,5,0,5);                   // Epoch state covariance
  P_a  = Dot(dady*P_y0,dady);                   // Semi-major axis covariance
  P_y  = U*P_y0*Transp(U);                      // Inertial state covariance

  // Output

  cout << "Epoch state standard deviation (including unmodelled effects)" 
       << endl << endl;
  cout << "  Position (xyz)  " 
       << fixed << setprecision(2) << setw(10) 
       << P_y.Diag().Sqrt().slice(0,2) << " [m]" << endl
       << "  Velocity (xyz)  " 
       << fixed << setprecision(5) << setw(10) 
       << P_y.Diag().Sqrt().slice(3,5)/3600.0 << " [m/s]" << endl
       << endl;
  cout << "  Semi-major axis "
       << fixed << setprecision(2) << setw(10) 
       << sqrt(P_a) << " [m]" << endl << endl;

/*
  // Debug print 
  cout << "Effect of unestimated parameters" << endl << endl;
  for (int i=0; i<n_con; i++) {
    cout << fixed << setprecision(3) << setw(10) 
         << -((R_xx_inv*R_xc).Col(i)*sqrt(Var_b(n_est-n_orb+i))).slice(0,2) << endl;
  };
  cout << endl;
*/

  // Position/velocity standard deviation 

  cout << "Position/velocity standard deviation (including unmodelled effects)" 
       << endl << endl
       << "                   Position                      Velocity         " 
       << endl
       << "  Time    radial     tang     normal    radial     tang     normal" 
       << endl
       << "   [h]      [m]       [m]       [m]      [m/s]     [m/s]     [m/s]" 
       << endl;

  t = 0.0;

  for(;;) {
    // Transition matrix
    T = Phi(n,t);
    // Covariance 
    P_y = U*T*P_y0*Transp(T)*Transp(U);
    // Output
    cout << fixed 
         << setprecision(1) << setw(6)  << t 
         << setprecision(1) << setw(10) << P_y.Diag().Sqrt().slice(0,2)
         << setprecision(4) << setw(10) << P_y.Diag().Sqrt().slice(3,5)/3600.0 
         << endl;
    // Increment time and exit if past end of prediction arc
    t += 3.0;
    if (t>T_pred) break;
  };

  return 0;

};
