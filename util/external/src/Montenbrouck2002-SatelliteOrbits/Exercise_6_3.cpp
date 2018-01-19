//------------------------------------------------------------------------------
//
// Exercise_6_3.cpp
// 
// Purpose: 
//
//   Satellite Orbits - Models, Methods, and Applications 
//   Exercise 6-3: User Clock Error from GPS Pseudorange.
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
#include "SAT_Time.h"
#include "SAT_VecMat.h"

using namespace std;


//------------------------------------------------------------------------------
//
// Intpol
// 
// Purpose:
//
//   Interpolation using Nevilles's algorithm
//
// Input/output:
//
//   n         Number of points
//   x         Abscissa values x_i
//   y         Ordinate values y_i=y(x_i)
//   x0        Interpolation point
//   <return>  Interpolated value y(x0)
//
// References:
//
//   Schwarz H.R.; Numerische Mathematik; B. G. Teubner, Stuttgart (1988).
// 
//------------------------------------------------------------------------------

double Intpol ( int n, const double *x, const double *y, double x0 )
{
  int     i,k;                    // Loop counters
  double* p = new double [n];     // Interpolation tableau
  double  y0;
  // Neville interpolation
  for (i=0;i<n;i++) p[i]=y[i];    // Copy of ordinate values
  for (k=1;k<n;k++)
    for (i=n-1;i>=k;i--)
      p[i] += (x0-x[i])*(p[i]-p[i-1])/(x[i]-x[i-k]);
  y0 = p[n-1];
  delete []p;       // Clean-up 
  return y0;        // Result
}


//------------------------------------------------------------------------------
//
// Main program
//
//------------------------------------------------------------------------------

int main() {

  // Constants
    
  const int    n_eph =  13;                       // Number of ephemeris points
  const int    n_obs =   6;                       // Number of observations
  const int    n_it  =   3;                       // Light time iterations steps

  const double L1    = 154.0;                     // Frequency multiplier L1
  const double L2    = 120.0;                     // Frequency multiplier L2
  const double f_rel = 1.0/(1.0-(L2*L2)/(L1*L1)); // Ionosphere multiplier

  // Goldstone station coordinates WGS-84 [m]

  const Vector r_Sta ( -2353614.1280, -4641385.4470, 3676976.5010 ); 

  // Variables

  int    i, i_it;                                 // Loop counters
  double tau;                                     // Light time [s]
  double t_snd,t_rcv;                             // Send and receive time [h]
  double dt_User;                                 // User clock offset [s]
  double dt_Sat;                                  // Satellite clock offset [s]
  double mean = 0.0;                              // Mean clock error [m]
  double range;                                   // Observed range [m]
  double rho;                                     // Computed range [m]
  Vector r_GPS;                                   // GPS pos., Earth-fixed [m]
  double res[n_obs];                              // Pseudorange residuals [m]


  // SP3 PRN1 position coordinates (x,y,z) [m] and clock error [s]
  // 1998/02/19 08:00:00.0 - 11:00:00.0, at 15 m intervals [GPS time]

  double t[n_eph]   = {  8.00, 8.25, 8.50, 8.75, 9.00, 9.25, 9.50, 9.75, 
                         10.00, 10.25, 10.50, 10.75, 11.00 };
                       
  double x[n_eph]   = { -15504291.797,-15284290.679,-14871711.829,-14242843.546,
                        -13380818.523,-12276418.004,-10928585.710, -9344633.744, 
                         -7540134.384, -5538503.062, -3370289.205, -1072201.838,
                          1314093.678 };
  
  double y[n_eph]   = { -21530763.883,-21684703.684,-21600510.259,-21306712.708,
                        -20837175.663,-20229688.085,-19524421.024,-18762314.034,
                        -17983451.817,-17225491.970,-16522202.377,-15902162.018,
                        -15387672.739 };   

  double z[n_eph]   = {  -1271498.273,  1573435.406,  4391350.089,  7133948.741,  
                          9754366.309, 12207953.668, 14453015.617, 16451492.281, 
                         18169574.686, 19578246.580, 20653745.961, 21377940.941, 
                         21738615.794 };             
 
  double clk[n_eph] = {  40.018233e-6, 40.097295e-6, 40.028697e-6, 40.154941e-6, 
                         40.193626e-6, 40.039288e-6, 40.012677e-6, 39.883106e-6, 
                         40.181357e-6, 40.328261e-6, 40.039533e-6, 40.052642e-6, 
                         40.025493e-6  };


  // Goldstone pseudo range observations P2, C1 [m] from PRN1
  // 1998/02/19 08:30:00.0 - 11:00:00.0, at 30 m intervals [GPS time]

  double    t_obs[n_obs]    = { 8.5, 9.0, 9.5, 10.0, 10.5, 11.0 };
  
  double    range_P2[n_obs] = { 21096577.475, 20519964.850, 20282706.954, 
                                20375838.496, 20751678.769, 21340055.129 };

  double    range_C1[n_obs] = { 21096579.501, 20519966.875, 20282709.233,
                                20375840.613, 20751680.997, 21340057.362 };


  // Process pseudorange measurements
  
  for (i=0;i<=n_obs-1;i++) {

    dt_Sat = Intpol(n_eph,t,clk,t_obs[i]);        // Satellite clock offset [s]
    
    range  = range_P2[i]                          // Pseudrorange corrected for
             - (range_P2[i]-range_C1[i])*f_rel    // - ionosphere
             + c_light*dt_Sat;                    // - satellite clock offset
             
    // Light time iteration for downleg GPS -> station

    dt_User = 0.0;
    tau     = 0.0;
    for (i_it=1;i_it<=n_it;i_it++) {
      t_rcv = t_obs[i] - dt_User/3600.0;          // Receive time (GPS)
      t_snd = t_rcv - tau/3600.0;                 // Transmit time (GPS) [h]
      r_GPS = Vector ( Intpol(n_eph,t,x,t_snd),   // Position of GPS satellite
                       Intpol(n_eph,t,y,t_snd),   // at t_snd in Earth-fixed
                       Intpol(n_eph,t,z,t_snd) ); // system at t_snd and in 
      r_GPS = R_z(omega_Earth*tau) * r_GPS;       // Earth-fixed frame at t_obs
      rho   = Norm(r_GPS-r_Sta);                  // Range
      tau   = rho/c_light;                        // Light time
      res[i]  = range - rho;                      // Residual
      dt_User = res[i]/c_light;                   // Approx. user clock error 
    };

    mean += res[i]; 

  };
  
  mean = mean/n_obs;                              // Scaling of mean value


  // Output

  cout << "Exercise 6-3: User Clock Error from GPS Pseudorange" << endl 
       << endl
       << "                               Error    Err-Mean" << endl
       << "                                [m]        [m]  " << endl;

  for (i=0;i<=n_obs-1;i++) {
    cout << " " << Date(Mjd(1998,2,19)+t_obs[i]/24.0) 
         << fixed << setprecision(2) << setw(12) << res[i]
         << setw(10) << res[i] - mean << endl; 
    };
  cout << endl 
       << " Mean value " << setw(25) << mean << endl; 

  return 0;
  
}
