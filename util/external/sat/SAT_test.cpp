//------------------------------------------------------------------------------
//
// Exercise_3_4.cpp
// 
// Purpose: 
//
//   Satellite Orbits - Models, Methods, and Applications
//   Exercise 3-4: Orbit Perturbations
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
#include "SAT_DE.h"
#include "SAT_Kepler.h"
#include "SAT_Force.h"
#include "SAT_RefSys.h"
#include "SAT_Time.h"
#include "SAT_VecMat.h"

using namespace std;

//------------------------------------------------------------------------------
//
// Global types and data
//
//------------------------------------------------------------------------------


// Record for passing global data between Deriv and the calling program 

struct AuxParam {
  double  Mjd0_TT;
  double  Area,mass,CR,CD;
  int     n,m;
  bool    Sun,Moon,SRad,Drag;
};



//------------------------------------------------------------------------------
//
// Accel
//
// Purpose:
//
//   Computes the acceleration of an Earth orbiting satellite due to 
//    - the Earth's harmonic gravity field, 
//    - the gravitational perturbations of the Sun and Moon
//    - the solar radiation pressure and
//    - the atmospheric drag
//
// Input/Output:
//
//   Mjd_TT      Terrestrial Time (Modified Julian Date)
//   r           Satellite position vector in the ICRF/EME2000 system
//   v           Satellite velocity vector in the ICRF/EME2000 system
//   Area        Cross-section 
//   mass        Spacecraft mass
//   CR          Radiation pressure coefficient
//   CD          Drag coefficient
//   <return>    Acceleration (a=d^2r/dt^2) in the ICRF/EME2000 system
//
//------------------------------------------------------------------------------

Vector Accel ( double Mjd_TT, const Vector& r, const Vector& v, 
               double Area, double mass, double CR, double CD,
               int n, int m, 
               bool FlagSun, bool FlagMoon, bool FlagSRad, bool FlagDrag )
{

  double Mjd_UT1;
  Vector a(3), r_Sun(3), r_Moon(3);
  Matrix T(3,3), E(3,3);

  // Acceleration due to harmonic gravity field

  Mjd_UT1 = Mjd_TT + (IERS::UT1_UTC(Mjd_TT)-IERS::TT_UTC(Mjd_TT))/86400.0;

  T = NutMatrix(Mjd_TT) * PrecMatrix(MJD_J2000,Mjd_TT);
  E = GHAMatrix(Mjd_UT1) * T;

  a = AccelHarmonic ( r,E, Grav.GM,Grav.R_ref,Grav.CS, n,m );

  // Luni-solar perturbations 

  r_Sun  = Sun(Mjd_TT);
  r_Moon = Moon(Mjd_TT);

  if (FlagSun)  a += AccelPointMass ( r, r_Sun,  GM_Sun  ); 
  if (FlagMoon) a += AccelPointMass ( r, r_Moon, GM_Moon ); 

  // Solar radiation pressure

  if (FlagSRad) a += AccelSolrad ( r, r_Sun, Area, mass, CR, P_Sol, AU );

  // Atmospheric drag

  if (FlagDrag) a += AccelDrag ( Mjd_TT, r, v, T, Area, mass, CD );

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

  double  Mjd_TT = (*p).Mjd0_TT + t/86400.0;

  // State vector components

  Vector r = y.slice(0,2);
  Vector v = y.slice(3,5);

  // Acceleration 

  Vector a(3);

  a = Accel ( Mjd_TT,r,v, (*p).Area, (*p).mass, (*p).CR, (*p).CD,
              (*p).n, (*p).m, (*p).Sun, (*p).Moon, (*p).SRad, (*p).Drag );

  // State vector derivative
  
  yp = Stack ( v, a );
   
};


//------------------------------------------------------------------------------
//
// Ephemeris computation
//
//------------------------------------------------------------------------------

void Ephemeris ( const Vector& Y0, int N_Step, double Step, AuxParam p, 
                 Vector Eph[] )

{

  int       i;
  double    t,t_end;
  double    relerr,abserr;        // Accuracy requirements
  DE        Orb(Deriv,6,&p);      // Object for integrating the eq. of motion
  Vector    Y(6);

  relerr = 1.0e-13;
  abserr = 1.0e-6;
  t      = 0.0;
  Y      = Y0;
  Orb.Init(t,relerr,abserr);
  for (i=0;i<=N_Step;i++) {
    t_end = Step*i;  
    Orb.Integ ( t_end, Y);
    Eph[i] = Y;
  };


}

//------------------------------------------------------------------------------
//
// Maximum computation
//
//------------------------------------------------------------------------------

template<class T> const T& Max (const T& a, const T& b)
{
    return (a<b) ? b:a;
}


//------------------------------------------------------------------------------
//
// Main program
//
//------------------------------------------------------------------------------

int main() {

  // Constants

  const int N_Step = 720;   // Maximum number of steps
  
  // Variables

  int       i;
  int       N_Step1;
  int       N_Step2;
  double    Step;
  double    Mjd0_UTC;
  Vector    Y0(6);
  Vector    Kep(6);
  AuxParam  Aux_ref,Aux;              // Auxiliary parameters
  double    Max_J20,Max_J22,Max_J44,Max_J1010;
  double    Max_Sun,Max_Moon,Max_SRad,Max_Drag;
  Vector    Eph_Ref  [N_Step+1];
  Vector    Eph_J20  [N_Step+1];
  Vector    Eph_J22  [N_Step+1];
  Vector    Eph_J44  [N_Step+1];
  Vector    Eph_J1010[N_Step+1];
  Vector    Eph_Sun  [N_Step+1];
  Vector    Eph_Moon [N_Step+1];
  Vector    Eph_SRad [N_Step+1];
  Vector    Eph_Drag [N_Step+1];
  

  // Initialize UT1-UTC and UTC-TAI time difference

  IERS::Set ( -0.05,-30.00, 0.0, 0.0 ); 

  
  // Epoch state (remote sensing satellite)

  Mjd0_UTC = Mjd(1999,03,01,00,00,0.0);
  Kep = Vector ( 7178.0e3, 0.0010, 98.57*Rad, 0.0, 0.0, 0.0 ); 
  Y0 = State ( GM_Earth, Kep, 0.0 );

  // Model parameters

  Aux_ref.Mjd0_TT = Mjd0_UTC + IERS::TT_UTC(Mjd0_UTC)/86400.0;
  Aux_ref.Area    = 5.0;     // [m^2]  Remote sensing satellite
  Aux_ref.mass    = 1000.0;  // [kg]
  Aux_ref.CR      = 1.3;     
  Aux_ref.CD      = 2.3;
  Aux_ref.n       = 20;
  Aux_ref.m       = 20;
  Aux_ref.Sun     = true;
  Aux_ref.Moon    = true;
  Aux_ref.SRad    = true;
  Aux_ref.Drag    = true;

  // Reference orbit

  Step    =  120.0; // [s]
  N_Step1 =     50; // 100 mins
  N_Step2 =    720; // 1 day

  Aux = Aux_ref; 
  Ephemeris ( Y0, N_Step2, Step, Aux, Eph_Ref );
  
  // J2,0 perturbations
  Aux.n = 2; Aux.m = 0;
  Ephemeris ( Y0, N_Step2, Step, Aux, Eph_J20 );

  // J2,2 perturbations
  Aux.n = 2; Aux.m = 2;
  Ephemeris ( Y0, N_Step2, Step, Aux, Eph_J22 );

  // J4,4 perturbations
  Aux.n = 4; Aux.m = 4;
  Ephemeris ( Y0, N_Step2, Step, Aux, Eph_J44 );

  // J10,10 perturbations
  Aux.n = 10; Aux.m = 10;
  Ephemeris ( Y0, N_Step2, Step, Aux, Eph_J1010 );
  Aux.n = 20; Aux.m = 20;

  // Solar perturbations
  Aux.Sun = false;
  Ephemeris ( Y0, N_Step2, Step, Aux, Eph_Sun );
  Aux.Sun = true;

  // Lunar perturbations
  Aux.Moon = false;
  Ephemeris ( Y0, N_Step2, Step, Aux, Eph_Moon );
  Aux.Moon = true;

  // Solar radiation pressure perturbations
  Aux.SRad = false;
  Ephemeris ( Y0, N_Step2, Step, Aux, Eph_SRad );
  Aux.SRad = true;

  // Drag perturbations
  Aux.Drag = false;
  Ephemeris ( Y0, N_Step2, Step, Aux, Eph_Drag );
  Aux.Drag = true;

  // Find maximum over N_Step1 steps

  Max_J20=Max_J22=Max_J44=Max_J1010=Max_Sun=Max_Moon=Max_SRad=Max_Drag = 0.0;
  for (i=0;i<=N_Step1;i++) {
    Max_J20   = Max(Norm (Eph_J20  [i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_J20);
    Max_J22   = Max(Norm (Eph_J22  [i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_J22);
    Max_J44   = Max(Norm (Eph_J44  [i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_J44);
    Max_J1010 = Max(Norm (Eph_J1010[i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_J1010);
    Max_Sun   = Max(Norm (Eph_Sun  [i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_Sun);
    Max_Moon  = Max(Norm (Eph_Moon [i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_Moon);
    Max_SRad  = Max(Norm (Eph_SRad [i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_SRad);
    Max_Drag  = Max(Norm (Eph_Drag [i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_Drag);
  };

  // Output
  
  cout << "Exercise 3-4: Orbit Perturbations " << endl << endl;

  cout << "Remote sensing satellite: " << endl << endl;
  
  cout << "  Maximum position errors within "
       << fixed << setprecision(1)
       << setw(5) << N_Step1*Step/60.0 
       << " min propagation interval " << endl << endl;
  cout << "    J2,0    J2,2    J4,4  J10,10" 
       << "     Sun    Moon  SolRad    Drag" << endl
       << "     [m]     [m]     [m]     [m]"
       << "     [m]     [m]     [m]     [m]" << endl;
  cout << setw(8) << Max_J20   << setw(8) << Max_J22  << setw(8) << Max_J44
       << setw(8) << Max_J1010 << setw(8) << Max_Sun  << setw(8) << Max_Moon
       << setw(8) << Max_SRad  << setw(8) << Max_Drag << endl << endl;

  // Find maximum over N_Step2 steps

  for (i=N_Step1+1;i<=N_Step2;i++) {
    Max_J20   = Max(Norm (Eph_J20  [i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_J20);
    Max_J22   = Max(Norm (Eph_J22  [i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_J22);
    Max_J44   = Max(Norm (Eph_J44  [i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_J44);
    Max_J1010 = Max(Norm (Eph_J1010[i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_J1010);
    Max_Sun   = Max(Norm (Eph_Sun  [i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_Sun);
    Max_Moon  = Max(Norm (Eph_Moon [i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_Moon);
    Max_SRad  = Max(Norm (Eph_SRad [i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_SRad);
    Max_Drag  = Max(Norm (Eph_Drag [i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_Drag);
  };

  // Output

  cout << "  Maximum position errors within "
       << setw(5) << N_Step2*Step/60.0 
       << " min propagation interval " << endl << endl;
  cout << "    J2,0    J2,2    J4,4  J10,10" 
       << "     Sun    Moon  SolRad    Drag" << endl
       << "     [m]     [m]     [m]     [m]"
       << "     [m]     [m]     [m]     [m]" << endl;
  cout << setw(8) << Max_J20   << setw(8) << Max_J22  << setw(8) << Max_J44
       << setw(8) << Max_J1010 << setw(8) << Max_Sun  << setw(8) << Max_Moon
       << setw(8) << Max_SRad  << setw(8) << Max_Drag << endl << endl;



  // Epoch state (geostationary satellite)

  Kep = Vector (42166.0e3, 0.0004,  0.02*Rad, 0.0, 0.0, 0.0 );
  Y0 = State ( GM_Earth, Kep, 0.0 );

  // Model parameters

  Aux_ref.Area    = 10.0;    // [m^2]  Geostationary satellite

  // Reference orbit

  Step    = 1200.0; // [s]
  N_Step1 =     72; // 1 day
  N_Step2 =    144; // 2 days 

  Aux = Aux_ref; 
  Ephemeris ( Y0, N_Step2, Step, Aux, Eph_Ref );
  
  // J2,0 perturbations
  Aux.n = 2; Aux.m = 0;
  Ephemeris ( Y0, N_Step2, Step, Aux, Eph_J20 );

  // J2,2 perturbations
  Aux.n = 2; Aux.m = 2;
  Ephemeris ( Y0, N_Step2, Step, Aux, Eph_J22 );

  // J4,4 perturbations
  Aux.n = 4; Aux.m = 4;
  Ephemeris ( Y0, N_Step2, Step, Aux, Eph_J44 );

  // J10,10 perturbations
  Aux.n = 10; Aux.m = 10;
  Ephemeris ( Y0, N_Step2, Step, Aux, Eph_J1010 );
  Aux.n = 20; Aux.m = 20;

  // Solar perturbations
  Aux.Sun = false;
  Ephemeris ( Y0, N_Step2, Step, Aux, Eph_Sun );
  Aux.Sun = true;

  // Lunar perturbations
  Aux.Moon = false;
  Ephemeris ( Y0, N_Step2, Step, Aux, Eph_Moon );
  Aux.Moon = true;

  // Solar radiation pressure perturbations
  Aux.SRad = false;
  Ephemeris ( Y0, N_Step2, Step, Aux, Eph_SRad );
  Aux.SRad = true;

  // Drag perturbations
  Aux.Drag = false;
  Ephemeris ( Y0, N_Step2, Step, Aux, Eph_Drag );
  Aux.Drag = true;

  // Find maximum over N_Step1 steps

  Max_J20=Max_J22=Max_J44=Max_J1010=Max_Sun=Max_Moon=Max_SRad=Max_Drag = 0.0;
  for (i=0;i<=N_Step1;i++) {
    Max_J20   = Max(Norm (Eph_J20  [i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_J20);
    Max_J22   = Max(Norm (Eph_J22  [i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_J22);
    Max_J44   = Max(Norm (Eph_J44  [i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_J44);
    Max_J1010 = Max(Norm (Eph_J1010[i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_J1010);
    Max_Sun   = Max(Norm (Eph_Sun  [i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_Sun);
    Max_Moon  = Max(Norm (Eph_Moon [i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_Moon);
    Max_SRad  = Max(Norm (Eph_SRad [i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_SRad);
    Max_Drag  = Max(Norm (Eph_Drag [i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_Drag);
  };

  // Output
  
  cout << "Geostationary satellite: " << endl << endl;
  
  cout << "  Maximum position errors within "
       << fixed << setprecision(1)
       << setw(5) << N_Step1*Step/60.0 
       << " min propagation interval " << endl << endl;
  cout << "    J2,0    J2,2    J4,4  J10,10" 
       << "     Sun    Moon  SolRad    Drag" << endl
       << "     [m]     [m]     [m]     [m]"
       << "     [m]     [m]     [m]     [m]" << endl;
  cout << setw(8) << Max_J20   << setw(8) << Max_J22  << setw(8) << Max_J44
       << setw(8) << Max_J1010 << setw(8) << Max_Sun  << setw(8) << Max_Moon
       << setw(8) << Max_SRad  << setw(8) << Max_Drag << endl << endl;

  // Find maximum over N_Step2 steps

  for (i=N_Step1+1;i<=N_Step2;i++) {
    Max_J20   = Max(Norm (Eph_J20  [i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_J20);
    Max_J22   = Max(Norm (Eph_J22  [i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_J22);
    Max_J44   = Max(Norm (Eph_J44  [i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_J44);
    Max_J1010 = Max(Norm (Eph_J1010[i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_J1010);
    Max_Sun   = Max(Norm (Eph_Sun  [i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_Sun);
    Max_Moon  = Max(Norm (Eph_Moon [i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_Moon);
    Max_SRad  = Max(Norm (Eph_SRad [i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_SRad);
    Max_Drag  = Max(Norm (Eph_Drag [i].slice(0,2)-Eph_Ref[i].slice(0,2)), Max_Drag);
  };

  // Output

  cout << "  Maximum position errors within "
       << setw(5) << N_Step2*Step/60.0 
       << " min propagation interval " << endl << endl;
  cout << "    J2,0    J2,2    J4,4  J10,10" 
       << "     Sun    Moon  SolRad    Drag" << endl
       << "     [m]     [m]     [m]     [m]"
       << "     [m]     [m]     [m]     [m]" << endl;
  cout << setw(8) << Max_J20   << setw(8) << Max_J22  << setw(8) << Max_J44
       << setw(8) << Max_J1010 << setw(8) << Max_Sun  << setw(8) << Max_Moon
       << setw(8) << Max_SRad  << setw(8) << Max_Drag << endl << endl;

  return 0;

}
