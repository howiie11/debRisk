//------------------------------------------------------------------------------
//
// RTOD.cpp
// 
// Purpose: 
//
//   Satellite Orbits - Models, Methods, and Applications 
//   Real Time Orbit Determination based on GPS navigation data 
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
#include "SAT_DE.h"
#include "SAT_Kepler.h"
#include "SAT_Filter.h"
#include "SAT_Force.h"
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

// GPS data record

struct ObsType {
  double  Mjd_GPS;    // GPS time
  Vector  r_obs;      // Measured position vector (WGS84, [m])
  Vector  v_obs;      // Measured velocity vector (WGS84, [m/s])
  Vector  r_ref;      // True position vector (WGS84, [m])
  Vector  v_ref;      // True velocity vector (WGS84, [m/s])
};

// Record for passing global data between Deriv and the calling program 

struct AuxParam {
  double  Mjd0_GPS;   // Reference epoch (GPS time)
  int     n_grav;     // Gravity model degree and order
};



//------------------------------------------------------------------------------
//
// Accel
//
// Purpose:
//
//   Computes the acceleration of an Earth orbiting satellite due to 
//   the Earth's harmonic gravity field
//
// Input/Output:
//
//   Mjd_GPS     GPS Time (Modified Julian Date)
//   r           Satellite position vector in the pseudo-true-of-date system
//   n           Degree and order of gravity field
//   <return>    Acceleration (a=d^2r/dt^2) in the pseudo-true-of-date system
//
//------------------------------------------------------------------------------

Vector Accel ( double Mjd_GPS, const Vector& r, int n )
{

  // Variables 

  double  Mjd_UT1;
  Vector  a(3);
  Matrix  U(3,3);

  // Acceleration due to harmonic gravity field

  Mjd_UT1 = Mjd_GPS + (IERS::UT1_UTC(Mjd_GPS)-IERS::GPS_UTC(Mjd_GPS))/86400.0;
  
  U = R_z(GMST(Mjd_UT1));

  a = AccelHarmonic ( r,U, GM_Earth, Grav.R_ref,Grav.CS, n,n );

  // Acceleration
  
  return a;

}


//------------------------------------------------------------------------------
//
// Deriv
//
// Purpose:
// 
//   Computes the derivative of the state vector 
//
// Note:
//
//   pAux is expected to point to a variable of type AuxDataRecord, which is
//   used to communicate with the other program sections and to hold data 
//   between subsequent calls of this function
//
//------------------------------------------------------------------------------

void Deriv ( double t, const Vector& y, Vector& yp, void* pAux )
{

  // Pointer to auxiliary data record
  
  AuxParam* p = static_cast<AuxParam*>(pAux);

  // Time

  double  Mjd_GPS = (*p).Mjd0_GPS + t/86400.0;

  // State vector derivative

  Vector r = y.slice(0,2);
  Vector v = y.slice(3,5);
  
  yp = Stack ( v, Accel(Mjd_GPS,r,(*p).n_grav) );
   
};


//------------------------------------------------------------------------------
//
// GetSetup
//
// Purpose:
// 
//   Reads RTOD setup parameters from file
//
//------------------------------------------------------------------------------

void GetSetup ( ifstream& inp,
                int& n_grav, double& Step, 
                double& sig_xyz,double& sig_pos, double& sig_vel, 
                double& w_pos, double& w_vel, double& EditLevel ) {

  // Input

  inp.ignore(30); inp >> n_grav;    inp.ignore(133,'\n'); 
  inp.ignore(30); inp >> Step;      inp.ignore(133,'\n'); 
  inp.ignore(30); inp >> sig_xyz;   inp.ignore(133,'\n'); 
  inp.ignore(30); inp >> sig_pos;   inp.ignore(133,'\n'); 
  inp.ignore(30); inp >> sig_vel;   inp.ignore(133,'\n'); 
  inp.ignore(30); inp >> w_pos;     inp.ignore(133,'\n'); 
  inp.ignore(30); inp >> w_vel;     inp.ignore(133,'\n'); 
  inp.ignore(30); inp >> EditLevel; inp.ignore(133,'\n'); 

}


//------------------------------------------------------------------------------
//
// GetObs
//
// Purpose:
// 
//   Read observation and reference orbit from input file
//
//------------------------------------------------------------------------------

void GetObs(ifstream& inp, ObsType& Obs) {

  // Variables

  int      Y,M,D,h,m;
  double   s;
  Vector   r_obs(3),v_obs(3),r_ref(3),v_ref(3);

  // Input

  inp >> Y; inp.ignore(1);
  inp >> M; inp.ignore(1);
  inp >> D; inp.ignore(1);
  inp >> h; inp.ignore(1);
  inp >> m; inp.ignore(1);
  inp >> s 
      >> r_obs(0) >> r_obs(1) >> r_obs(2)
      >> v_obs(0) >> v_obs(1) >> v_obs(2)
      >> r_ref(0) >> r_ref(1) >> r_ref(2)
      >> v_ref(0) >> v_ref(1) >> v_ref(2);
  inp.ignore(255,'\n'); 

  // Data copy

  Obs.Mjd_GPS = Mjd(Y,M,D,h,m,s);
  Obs.r_obs   = r_obs;
  Obs.v_obs   = v_obs;
  Obs.r_ref   = r_ref;
  Obs.v_ref   = v_ref;

}


//------------------------------------------------------------------------------
//
// Main program
//
//------------------------------------------------------------------------------

int main(int argc, char* argv[]) {

  // Constants
    
  const Vector  dxdY(1.0,0.0,0.0,0.0,0.0,0.0);  // Partials
  const Vector  dydY(0.0,1.0,0.0,0.0,0.0,0.0);  
  const Vector  dzdY(0.0,0.0,1.0,0.0,0.0,0.0);  

  // Variables

  bool      Reject;                       // Flag for rejected data
  int       i;                            // Loop counter
  double    Step;                         // Nominal RK4 step size [s]
  double    sigma_xyz;                    // Measurement std. dev. [m]
  double    sigma_pos;                    // A priori std. dev. [m]
  double    sigma_vel;                    // A priori std. dev. [m/s]
  double    w_pos;                        // State noise std. dev. [m]
  double    w_vel;                        // State noise std. dev. [m/s]
  double    Mjd0,Mjd_UT1;                 // Time 
  double    t,t_old;
  double    EditLevel;                    // Data edit level
  double    a_obs, a_est, a_ref;          // Semi-major axis
  double    Sig_a;                        // Standard deviation
  Vector    Y(6),Y_old(6),r(3),v(3);      // State vector
  Vector    r_est(3),v_est(3);            // State vector
  Vector    SigY(6);                      // Standard deviation
  Vector    dadY(6);                      // Partials of sma. wrt. state 
  Matrix    U(3,3),dU(3,3),S(3,3);        // Sidereal Time matrix and derivative
  Matrix    Phi(6,6), P(6,6), Q(6,6);     // State transition, covariance and
                                          // process  noise matrix
  EKF       Filter(6);                    // Extended Kalman filter object
  AuxParam  Aux;                          // Auxiliary data
  double    h;                            // Step size 
  RK4       Orbit(Deriv,6,&Aux);          // Object for integrating the
                                          // equation of motion
  ifstream  inp;                          // Input file 
  ObsType   Obs;                          // Observations



  // Initialize UT1-UTC and UTC-TAI time difference

  IERS::Set ( -0.05,-30.00, 0.0, 0.0 ); 

  // Control parameters from setup file

  if (argc>1) inp.open(argv[1]); else inp.open("RTOD.inp"); 
  if (!inp) { 
    cerr << "ERROR: Could not open RTOD setup file" << endl;
    exit(1);
  };
  
  GetSetup ( inp, Aux.n_grav,Step, sigma_xyz,sigma_pos,sigma_vel, 
             w_pos,w_vel, EditLevel );
  
  inp.close();

  // Open GPS data input file, skip header and read first observation
  
  if (argc>2) inp.open(argv[2]); else inp.open("RTOD.dat"); 
  if (!inp) { 
    cerr << "ERROR: Could not open RTOD data file" << endl;
    exit(1);
  };
  
  inp.ignore(255,'\n'); 
  GetObs (inp, Obs);

  // Reference epoch (GPS time)

  Mjd0 = Obs.Mjd_GPS;
  Aux.Mjd0_GPS = Obs.Mjd_GPS;
  t = 0.0;

  // Transformation from WGS to pseudo-true-of-date system

  Mjd_UT1 = Mjd0 + ( IERS::UT1_UTC(Mjd0)
                      - IERS::GPS_UTC(Mjd0) )/86400.0;

  S(0,1) = 1.0; S(1,0) = -1.0; 
  U  = R_z(GMST(Mjd_UT1));            // Earth rotation matrix
  dU = omega_Earth*S*U;               // and derivative [1/s]
 
  // Pseudo true of date position vector

  r = Transp(U)*Obs.r_obs;
  v = Transp(U)*Obs.v_obs + Transp(dU)*Obs.r_obs;
  
  Y = Stack(r,v);

  // A priori covariance

  SigY = Vector ( sigma_pos,sigma_pos,sigma_pos,    
                  sigma_vel,sigma_vel,sigma_vel );
  P = Diag(SigY);                    
  P = P*P;
  
  // Initialization

  Filter.Init(t,Y,P);

  // Header
  
  cout << "GPS/MET Real Time Orbit Determination" << endl << endl
       << " Initial Epoch:        " << "   " << Date(Obs.Mjd_GPS) << endl
       << fixed 
       << " Inertial position:    " << setw(15) << setprecision(3) << r << endl
       << " Inertial velocity:    " << setw(15) << setprecision(6) << v << endl
       << endl
       << " Gravity field order:  " << Aux.n_grav << endl
       << setprecision(1) 
       << " Step size:            " << Step << " s" << endl
       << " Measurement sigma:    " << sigma_xyz << " m" << endl
       << " A priori pos. sigma:  " << sigma_pos << " m" << endl
       << setprecision(4)
       << " A priori vel. sigma:  " << sigma_vel << " m/s" << endl
       << " State noise position: " << setprecision(1)<< w_pos << " m" << endl
       << " State noise velocity: " << setprecision(4)<< w_vel <<" m/s"<< endl
       << " Edit level (sigma):   " << setprecision(1) << EditLevel << endl
       << endl << endl
       << "   Time      x(WGS)     y(WGS)     z(WGS)" 
       << "       Position [m]        Velocity [m/s]    Semi-major axis [m] "
       << endl 
       << "    [s]        [m]        [m]        [m] "
       << "    Meas. Sigma  Estim. Meas.  Sigma  Estim.  Meas. Sigma  Estim."
       << endl;

  // Measurement loop

  for (i=1;i<=2000;i++) {

    // Previous step
    t_old = Filter.Time();
    Y_old = Filter.State();

    // Next observation
    GetObs (inp, Obs);
    if (inp.fail()) break;
    
    // Propagation to measurement epoch
    t = (Obs.Mjd_GPS-Mjd0)*86400.0;             // Time since epoch [s]
    TwoBody ( GM_Earth, Y_old,t-t_old, Y,Phi);  // State vector

    // Process Noise
    SigY = Vector(w_pos,w_pos,w_pos,w_vel,w_vel,w_vel);
    Q = Diag(SigY);
    Q = Q*Q;

    // Integration to time of measurement  
    Y = Y_old;
    while (t_old<t) {            
      h = (t-t_old);             
      if (h>Step) h=Step;        
      Orbit.Step(t_old,Y,h);
    };
         
    // Time update
    Filter.TimeUpdate(t,Y,Phi,Q);

    // Earth rotation matrix and derivative
    Mjd_UT1 = Obs.Mjd_GPS + ( IERS::UT1_UTC(Obs.Mjd_GPS)
                              - IERS::GPS_UTC(Obs.Mjd_GPS) )/86400.0;
    U  = R_z(GMST(Mjd_UT1));            
    dU = omega_Earth*S*U;   

    // Transformation of measurements to inertial system
    r = Transp(U)*Obs.r_obs;

    // Data editing limit 
    Reject = ( Norm(r-Y.slice(0,2)) > EditLevel*sigma_xyz );
    
    // Measurement updates
    if ( !Reject ) {
      Filter.MeasUpdate ( r(0), Filter.State()(0), sigma_xyz, dxdY );
      Filter.MeasUpdate ( r(1), Filter.State()(1), sigma_xyz, dydY );
      Filter.MeasUpdate ( r(2), Filter.State()(2), sigma_xyz, dzdY );
    }
    
    // Estimated state vector in rotating, Earth-fixed system

    r_est = U*Filter.State().slice(0,2);
    v_est = U*Filter.State().slice(3,5) + dU*Filter.State().slice(0,2);

    // Semi-major axis

    r = Transp(U)*Obs.r_obs;
    v = Transp(U)*Obs.v_obs + Transp(dU)*Obs.r_obs;
    a_obs = 1.0/(2.0/Norm(r)-Dot(v,v)/GM_Earth);

    r = Transp(U)*Obs.r_ref;
    v = Transp(U)*Obs.v_ref + Transp(dU)*Obs.r_ref;
    a_ref = 1.0/(2.0/Norm(r)-Dot(v,v)/GM_Earth);

    r = Filter.State().slice(0,2);
    v = Filter.State().slice(3,5);
    a_est = 1.0/(2.0/Norm(r)-Dot(v,v)/GM_Earth);

    dadY = Stack ( (2.0*pow(a_est,2)/pow(Norm(r),3))*r,
                   (2.0*pow(a_est,2)/GM_Earth      )*v );

    Sig_a = sqrt ( Dot(dadY*Filter.Cov(),dadY) );

    // Output

    cout << setw(8) << setprecision(1) << t << " "     // Time
         << setw(11) << r_est                          // Estimated position
         << setw(7) << setprecision(1)
         << Norm(Obs.r_obs-Obs.r_ref)                  // Measurement errors
         << (Reject? "*":" ")                          // Data editing flag
         << setw(6) << setprecision(1)
         << Norm(Filter.StdDev().slice(0,2))           // Standard deviation
         << setw(7) << setprecision(1)
         << Norm(r_est-Obs.r_ref)                      // Position error
         << setw(7) << setprecision(3)
         << Norm(Obs.v_obs-Obs.v_ref)                  // Measurement errors
         << setw(7) << setprecision(3)
         << Norm(Filter.StdDev().slice(3,5))           // Standard deviation
         << setw(7) << setprecision(3)
         << Norm(v_est-Obs.v_ref)                      // Velocity error
         << setw(7) << setprecision(1)
         << fabs(a_obs-a_ref)                          // Measurement errors
         << setw(7) << setprecision(1)
         << Sig_a                                      // Standard deviation
         << setw(7) << setprecision(1)
         << fabs(a_est-a_ref)                          // Semi-major axis error
         << endl;

  }

  // Close data file

  inp.close();

  return 0;

}
