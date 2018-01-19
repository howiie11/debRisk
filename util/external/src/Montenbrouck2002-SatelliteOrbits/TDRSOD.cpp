//------------------------------------------------------------------------------
//
// TDRSOD.cpp
// 
// Purpose: 
//
//   Satellite Orbits - Models, Methods, and Applications 
//   Orbit Determination from Tracking and Data Relay Satellite measurements
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
#include <vector>

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

// TDRS tracking data record

struct ObsType {
  double  Mjd_UTC;    // Ground received time
  int     StaNo;      // Station number
  int     TdrsNo;     // TDRS number
  double  Range_4w;   // 4-way range measurement
};

// Dummy relational operators required for use with STL template class vector<>

bool operator==(ObsType a, ObsType b) { return true; };
bool operator< (ObsType a, ObsType b) { return true; };


// Record for passing global data between the differential equation function 
// and external program units

struct AuxParam {
  double  Mjd0_TT;       // Reference epoch
  double  Area,mass;     // Spacecraft cross-section [m^2] and mass [kg]
  double  CR,CD;         // Radiation pressure and drag coefficient
};


// Station parameters

struct StaParam {
  int       StaNo;       // Station number
  Geodetic  Loc;         // Station location
};


//------------------------------------------------------------------------------
//
// TrjData class (specification and implementation)
//
// Purpose:
//
//   The TrjDate class provides a framework for the integration of spacecraft
//   trajectories including both the state and the variational equations. 
//   It handles the required initialization steps, the packing and unpacking
//   of the state vector, state transition matrix and sensitivity matrix and
//   the mapping between the external time scale (Modified Julian Date TT)
//   and the internal time scale (seconds since reference epoch)
//
//------------------------------------------------------------------------------

class TrjData {
  
  public:
  
    // Constructor
    
    TrjData(): N_var(0),p(0) {};         // Default constructor
    
    TrjData ( DEfunct    Deriv,          // Differential equation
              DEfunct    VarEqn,         // Variational equations
              int        n_var,          // Dimension of variat. eqns.
              AuxParam*  pAux   )        // Pointer to auxiliary data
    {
      // Store dimension and aux. data pointer
      N_var = n_var;                
      p     = pAux;
      // Allocate objects
      Orb.Define(Deriv ,6      ,pAux);   // State vector integration
      Var.Define(VarEqn,7*n_var,pAux);   // Variational equations integration
      y     = Vector(6);                 // State vector
      Y     = Vector(7*n_var);           // State and partials
      PhiS_ = Matrix(6,n_var);           // transition and sensitivity matrix
    };
    
    // Definition
    void Define ( DEfunct    Deriv,      // Differential equation
                  DEfunct    VarEqn,     // Variational equations
                  int        n_var,      // Dimension of variat. eqns.
                  AuxParam*  pAux   )    // Pointer to auxiliary data
    {
      // Store dimension and aux. data pointer
      N_var = n_var;                
      p     = pAux;
      // Allocate objects
      Orb.Define(Deriv ,6      ,pAux);   // State vector integration
      Var.Define(VarEqn,7*n_var,pAux);   // Variational equations integration
      y     = Vector(6);                 // State vector
      Y     = Vector(7*n_var);           // State and partials
      PhiS_ = Matrix(6,n_var);           // Transition and sensitivity matrix
    };
  
    // Initialization
    
    void Init (double Mjd0_TT, Vector y0) 
    {
      double    t = ( Mjd0_TT - (*p).Mjd0_TT )*86400.0;
      // Define initial epoch relative to reference time 
      // and accuracy requirements
      Orb.Init(t,1.0e-13,1.0e-6);
      Var.Init(t,1.0e-13,1.0e-6);
      // Copy epoch state vector
      y = y0;
      // Create combined vector from epoch state, epoch transition matrix (=1)
      // and epoch sensitivity matrix (=0)
      for (int i=0; i<=5; i++) {
        Y(i) = y0(i);
        for (int j=0; j<N_var; j++)  Y(6*(j+1)+i) = ( i==j ? 1 : 0 ); 
      };
    };
  
    // Internal integration to given time 
    
    void Integ(double Mjd_TT) {
      double    t = ( Mjd_TT - (*p).Mjd0_TT )*86400.0;
      Orb.Integ ( t, y );
      Var.Integ ( t, Y );
    };

    // State interpolation
    
    Vector State(double Mjd_TT) {
      double    t = ( Mjd_TT - (*p).Mjd0_TT )*86400.0;
      Orb.Intrp ( t, y );
      return y;
    };
    
    // State transition and sensitivity matrix interpolation

    Matrix PhiS (double Mjd_TT) {
      double    t = ( Mjd_TT - (*p).Mjd0_TT )*86400.0;
      Var.Intrp ( t, Y );
      for (int j=0; j<N_var; j++) 
        PhiS_.SetCol(j,Y.slice(6*(j+1),6*(j+1)+5));
      return PhiS_;
    };
  
  private:

    // Elements

    int       N_var;      // Dimension of variational equations
    DE        Orb;        // Object for integrating the state vector
    DE        Var;        // Object for integrating the variational equations
    AuxParam* p;          // Pointer to auxiliary parameters
    Vector    y;          // Auxiliary state vector
    Vector    Y;          // Auxiliary state vector and partials
    Matrix    PhiS_;      // State transition and sensitivity matrix

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
               double Area, double mass, double CR, double CD  )
{

  // Constants 
  
  const int n = 20;  // Degree and order of gravity field (<=Grav.n_max)
  const int m = 20;
  
  // Variables 

  double  Mjd_UT1;
  Vector  a(3), r_Sun(3), r_Moon(3);
  Matrix  T(3,3), U(3,3);

  // Acceleration due to harmonic gravity field

  Mjd_UT1 = Mjd_TT + (IERS::UT1_UTC(Mjd_TT)-IERS::TT_UTC(Mjd_TT))/86400.0;

  T = NutMatrix(Mjd_TT) * PrecMatrix(MJD_J2000,Mjd_TT);
  U = GHAMatrix(Mjd_UT1) * T;

  a = AccelHarmonic ( r,U, Grav.GM,Grav.R_ref,Grav.CS, n,m );

  // Luni-solar perturbations 

  r_Sun  = Sun(Mjd_TT);
  r_Moon = Moon(Mjd_TT);

  a += AccelPointMass ( r, r_Sun,  GM_Sun  ); 
  a += AccelPointMass ( r, r_Moon, GM_Moon ); 

  // Solar radiation pressure

  a += Illumination ( r, r_Sun )
       * AccelSolrad ( r, r_Sun, Area, mass, CR, P_Sol, AU );

  // Atmospheric drag

  a += AccelDrag ( Mjd_TT, r, v, T, Area, mass, CD );

  // Acceleration
  
  return a;

}


//------------------------------------------------------------------------------
//
// Accel
//
// Purpose:
//
//   Computes the acceleration and gradient of an Earth orbiting satellite
//   due to the Earth's gravity field (low degree and order)
//
// Input/Output:
//
//   Mjd_TT      Terrestrial Time (Modified Julian Date)
//   r           Satellite position vector in the ICRF/EME2000 system
//   v           Satellite velocity vector in the ICRF/EME2000 system
//   a           Acceleration (a=d^2r/dt^2) in the ICRF/EME2000 system
//   Area        Cross-section 
//   mass        Spacecraft mass
//   G           Gradient (G=da/dr) in the ICRF/EME2000 system
//   dadCD       Partials of acceleration w.r.t. to drag coefficient
//   dadCR       Partials of acceleration w.r.t. to solrad coefficient
//
//------------------------------------------------------------------------------

void Accel ( double Mjd_TT, const Vector& r, const Vector& v,
             double Area, double mass, 
             Vector& a, Matrix& G, Vector& dadCD, Vector& dadCR )
{

  // Constants

  const int    n = 2;     // Degree and order of gravity field (<=Grav.n_max)
  const int    m = 0;
  const double d = 1.0;   // Position increment [m]

  // Variables 

  double  Mjd_UT1;
  Vector  da(3),dr(3),r_Sun(3);
  Matrix  T(3,3), U(3,3);

  // Acceleration due to the Earth's gravity field

  Mjd_UT1 = Mjd_TT + (IERS::UT1_UTC(Mjd_TT)-IERS::TT_UTC(Mjd_TT))/86400.0;

  T = NutMatrix(Mjd_TT) * PrecMatrix(MJD_J2000,Mjd_TT);
  U = GHAMatrix(Mjd_UT1) * T;

  a = AccelHarmonic ( r,U, Grav.GM,Grav.R_ref,Grav.CS, n,m );

  // Gradient

  for (int i=0; i<=2; i++) {
    // Set offset in i-th component of the position vector
    dr = 0.0; dr(i) = d;
    // Acceleration difference
    da = AccelHarmonic ( r+dr,U, GM_Earth, Grav.R_ref,Grav.CS, n,m ) -  a;
    // Derivative with respect to i-th axis
    G.SetCol(i,da/d);
  }

  // Drag coefficient partials

  dadCD = AccelDrag ( Mjd_TT, r, v, T, Area, mass, 1.0 );

  // Radiation pressure coefficient partials
  
  r_Sun = Sun(Mjd_TT);
  dadCR = Illumination ( r, r_Sun )
          * AccelSolrad ( r, r_Sun, Area, mass, 1.0, P_Sol, AU );

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

  a = Accel ( Mjd_TT,r,v, (*p).Area, (*p).mass, (*p).CR, (*p).CD );

  // State vector derivative
  
  yp = Stack ( v, a );
   
};


//------------------------------------------------------------------------------
//
// VarEqn
//
// Purpose:
// 
//   Computes the variational equations, i.e. the derivative of the state vector
//   and the state transition matrix
//
// Input/Output:
//
//   t           Time since epoch (*pAux).Mjd_0 in [s]
//   yPhiS       (6+36+6)-dim vector comprising the state vector (y), the
//               state transition matrix (Phi) and the sensitivity matrix
//               in column wise storage order
//   yPhiSp      Derivative of yPhiS
//   pAux        Pointer; pAux is expected to point to a variable of type
//               AuxDataRecord, which is used to communicate with the other
//               program sections and to hold data between subsequent calls
//               of this function
//
//------------------------------------------------------------------------------

void VarEqn ( double t, const Vector& yPhiS, Vector& yPhiSp, void* pAux )
{

  // Variables

  AuxParam* p;             // Auxiliary data pointer
  int       i,j;           // Loop counters
  double    Mjd_TT;        // Modified Julian Date
  Vector    r(3), v(3);    // Position, velocity
  Vector    a(3);          // Acceleration
  Vector    dadCD(3);      // Drag coefficient partials
  Vector    dadCR(3);      // Radiation pressure coefficient partials
  Matrix    G(3,3);        // Gradient of acceleration
  Matrix    Phi(6,6);      // State transition matrix 
  Matrix    S(6,2);        // Sensitivity matrix
  Matrix    Phip(6,6);     // Time derivative of state transition matrix
  Matrix    Sp(6,2);       // Time derivative of sensitivity matrix
  Matrix    dfdy(6,6);     // 
  Matrix    dfdp(6,2);     // 

  // Pointer to auxiliary data record
  
  p = static_cast<AuxParam*>(pAux);

  // Time

  Mjd_TT = (*p).Mjd0_TT + t/86400.0;

  // State vector components

  r = yPhiS.slice(0,2);
  v = yPhiS.slice(3,5);

  // State transition matrix

  for (j=0; j<=5; j++) Phi.SetCol(j,yPhiS.slice(6*(j+1),6*(j+1)+5));
  
  // Sensitivity matrix 

  for (j=6; j<=7; j++) S.SetCol(j-6,yPhiS.slice(6*(j+1),6*(j+1)+5));
  
  // Acceleration and gradient

  Accel ( Mjd_TT,r,v, (*p).Area,(*p).mass, a, G,dadCD,dadCR );

  // Time derivative of state transition matrix

  for (i=0; i<=2; i++) {
    for (j=0; j<=2; j++) {
      dfdy(i  ,j) = 0.0;                  // dv/dr(i,j)
      dfdy(i+3,j) = G(i,j);               // da/dr(i,j)
      dfdy(i  ,j+3) = ( i==j ? 1 : 0 );   // dv/dv(i,j)
      dfdy(i+3,j+3) = 0.0;                // da/dv(i,j)
    };
  };

  Phip = dfdy*Phi;
  
  // Time derivative of sensitivity matrix

  for (i=0; i<=2; i++) {
    dfdp(i  ,0) = 0.0;                    // dv/dCD(i)
    dfdp(i+3,0) = dadCD(i);               // da/dCD(i)
    dfdp(i  ,1) = 0.0;                    // dv/dCR(i)
    dfdp(i+3,1) = dadCR(i);               // da/dCR(i)
  };

  Sp = dfdy*S + dfdp;

  // Derivative of combined state vector and state transition matrix
  
  for (i=0; i<=2; i++) {
    yPhiSp(i)   = v(i);                    // dr/dt(i)
    yPhiSp(i+3) = a(i);                    // dv/dt(i)
  };
  for (i=0; i<=5; i++) 
    for (j=0; j<=5; j++)
      yPhiSp(6*(j+1)+i  ) = Phip(i,j);     // dPhi/dt(i,j)
  for (i=0; i<=5; i++) 
    for (j=6; j<=7; j++)
      yPhiSp(6*(j+1)+i  ) = Sp(i,j-6);     // dS/dt(i,j)

};


//------------------------------------------------------------------------------
//
// Range_4W
//
// Purpose:
//
//   Computes the TDRS 4-way range
//
// Input/output:
//
//   Mjd_UTC     Ground-received time of measurement
//   Trj_User    User satellite trajectory integration object
//   Trj_TDRS    Relay satellite trajectory integration object
//   Sta         Ground station
//   rho         4-way range
//   drho_dUser  Partials w.r.t. user satellite epoch state and parameters
//   drho_dTDRS  Partials w.r.t. TDRS satellite epoch state and parameters
//
//------------------------------------------------------------------------------

void Range_4W ( double MJD_UTC, TrjData& Trj_User, TrjData& Trj_TDRS, 
                const Geodetic& Sta, 
                double& rho, Vector& drho_dUser, Vector& drho_dTDRS )

{
  // Constants

  const int i_max = 2;             // Order of light-time correction

  // Variables

  int     i;                       // Counter
  double  MJD_TT;                  // Terrestrial Time
  double  UT1_TT;                  // UT1-TT time difference [d]
  double  tau1,tau2,tau3,tau4;     // Light-times
  double  rho1,rho2,rho3,rho4;     // Legs
  Vector  R_Sta(3), r_Sta(3);      // Earth-fixed and inertial station position
  Vector  r_User(3),r_TDRS(3);     // User and relay satellite position
  Vector  r1(3),r2(3),r3(3),r4(3); // Legs
  Matrix  PN(3,3);                 // Precession/nutation matrix
  
  // Ground received time (Terrestrial Time)

  MJD_TT = MJD_UTC + IERS::TT_UTC (MJD_UTC)/86400.0;
  UT1_TT = ( IERS::UT1_UTC(MJD_UTC)- IERS::TT_UTC(MJD_UTC) ) / 86400.0;

  // Precession, nutation
  
  PN = Transp ( NutMatrix(MJD_TT) * PrecMatrix(MJD_J2000,MJD_TT) );

  // Station (pseudo Earth-fixed)

  R_Sta = Transp(PoleMatrix(MJD_UTC)) * Sta.Position();

  // Downleg TRDS -> station
  
  tau1 = 0.0;
  r_Sta = PN * Transp(GHAMatrix(MJD_TT+UT1_TT)) * R_Sta;
  for (i=0;i<=i_max;i++) {
    r_TDRS = Trj_TDRS.State(MJD_TT-tau1).slice(0,2);
    r1     = r_TDRS-r_Sta;
    rho1   = Norm(r1);
    tau1   = (rho1/c_light)/86400.0;
  };

  // Downleg User -> TDRS

  tau2 = 0.0;
  for (i=0;i<=i_max;i++) {
    r_User = Trj_User.State(MJD_TT-tau1-tau2).slice(0,2);
    r2     = r_User-r_TDRS;
    rho2   = Norm(r2);
    tau2   = (rho2/c_light)/86400.0;
  };

  // Upleg TDRS -> User

  tau3 = 0.0;
  for (i=0;i<=i_max;i++) {
    r_TDRS = Trj_TDRS.State(MJD_TT-tau1-tau2-tau3).slice(0,2);
    r3     = r_TDRS-r_User;
    rho3   = Norm(r3);
    tau3   = (rho3/c_light)/86400.0;
  };

  // Upleg station -> TDRS

  tau4 = 0.0;
  for (i=0;i<=i_max;i++) {
    r_Sta  = PN * Transp(GHAMatrix(MJD_TT-tau1-tau2-tau3-tau4+UT1_TT)) * R_Sta;
    r4     = r_Sta-r_TDRS;
    rho4   = Norm(r4);
    tau4   = (rho3/c_light)/86400.0;
  };

  // Range

  rho = 0.5 * ( rho1 + rho2 + rho3 + rho4 );

  // Range partials w.r.t. parameters of user satellite 

  drho_dUser =  0.5 * ( r2/rho2 - r3/rho3 ) 
                    * Trj_User.PhiS(MJD_TT-tau1-tau2).slice(0,2,0,7);
  drho_dTDRS =  0.5 * ( r1/rho1 - r2/rho2 ) 
                    * Trj_TDRS.PhiS(MJD_TT-tau1).slice(0,2,0,7)
              + 0.5 * ( r3/rho3 - r4/rho4 )
                    * Trj_TDRS.PhiS(MJD_TT-tau1-tau2-tau3).slice(0,2,0,7);

};


//------------------------------------------------------------------------------
//
// GetSetup
//
// Purpose:
// 
//   Reads TDRSOD setup parameters from file
//
//------------------------------------------------------------------------------

void GetSetup ( ifstream& inp,
                double& Mjd0_TT, int& n_iterat,
                int& n_sat, int TdrsNo[], Vector y0[], Vector Sigy0[], 
                Vector p[], Vector Sigp[], AuxParam Aux[], 
                int& n_Sta, StaParam Sta[]
               ) {

  // Variables

  int       i,i_sat,i_Sta;
  int       Y,M,D,h,m;
  double    s, Mjd0_UTC;
  double    CD, CR;
  double    UT1_UTC, UTC_TAI, x_pole, y_pole;
  Vector    R(3);

  // TDRS satellite numbers

  inp.ignore(25); inp >> n_sat;    // Number of TDRS satellites
  n_sat += 1;                      // Total number of satellites
  for (i_sat=1;i_sat<n_sat;i_sat++) 
    inp >> TdrsNo[i_sat];
  TdrsNo[0] = 0;
  inp.ignore(81,'\n');

  // Station numbers

  inp.ignore(25); inp >> n_Sta; 
  for (i_Sta=0;i_Sta<n_Sta;i_Sta++) inp >> Sta[i_Sta].StaNo;
  inp.ignore(81,'\n');

  // Number of iterations

  inp.ignore(25); inp >> n_iterat;  inp.ignore(81,'\n');

  // Earth rotation parameters

  inp.ignore(25); inp >> UT1_UTC >> UTC_TAI;  inp.ignore(81,'\n'); 
  inp.ignore(25); inp >> x_pole >> y_pole;  inp.ignore(81,'\n'); 
  
  IERS::Set ( UT1_UTC, UTC_TAI, x_pole, y_pole );

  // Epoch 

  inp.ignore(25);
  inp >> Y; inp.ignore(1);  inp >> M; inp.ignore(1);  inp >> D; 
  inp >> h; inp.ignore(1);  inp >> m; inp.ignore(1);  inp >> s ;
  inp.ignore(81,'\n');
   
  Mjd0_UTC = Mjd(Y,M,D,h,m,s);                            // UTC  
  Mjd0_TT  = Mjd0_UTC + IERS::TT_UTC(Mjd0_UTC)/86400.0;   // TT

  for (i_sat=0;i_sat<n_sat;i_sat++)  Aux[i_sat].Mjd0_TT = Mjd0_TT;

  // Epoch state vector, s/c parameters and a priori standard deviations
  // for user s/c and TDRS satellites

  for (i_sat=0;i_sat<n_sat;i_sat++) {
  
    // Instantiation of vector objects

    y0[i_sat]      = Vector(6);                     // Epoch state (xyz,vxyz)
    Sigy0[i_sat]   = Vector(6);                     // Standard deviation
    p[i_sat]       = Vector(2);                     // Parameters (CD,CR)
    Sigp[i_sat]    = Vector(2);                     // Standard deviation

    // Epoch state and standard deviation

    for (i=0;i<6;i++) {
      inp.ignore(25); inp >> y0[i_sat](i) >> Sigy0[i_sat](i);  
      inp.ignore(81,'\n'); 
    };
    
    // Mass, area, CD, CR, sig(CD), sig(CR)

    inp.ignore(25); inp >> Aux[i_sat].mass >> Aux[i_sat].Area;  
    inp.ignore(81,'\n');
    inp.ignore(25); inp >> CD >> Sigp[i_sat](0);  inp.ignore(81,'\n'); 
    inp.ignore(25); inp >> CR >> Sigp[i_sat](1);  inp.ignore(81,'\n'); 
    p[i_sat](0) = CD;  Aux[i_sat].CD = CD;
    p[i_sat](1) = CR;  Aux[i_sat].CR = CR;

  };

  // Station location

  for (i_Sta=0;i_Sta<n_Sta;i_Sta++) {
    inp.ignore(25); inp >> R(0) >> R(1) >> R(2);  inp.ignore(81,'\n'); 
    Sta[i_Sta].Loc = Geodetic(R);
  }

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
  inp >> s ;
  
  inp >> Obs.StaNo >> Obs.TdrsNo >> Obs.Range_4w;
  inp.ignore(255,'\n'); 

  // Copy data

  Obs.Mjd_UTC = Mjd(Y,M,D,h,m,s);   // UTC
  Obs.Range_4w *=1000.0;            // [m]

}


//------------------------------------------------------------------------------
//
// Main program
//
//------------------------------------------------------------------------------

int main(int argc, char* argv[]) {

  
  // Constants
  
  const int    N_sat    =  3;           // Max. numb of satellites (user+TDRS)
  const int    N_Sta    =  2;           // Max. numb of stations 

  const double sigma_range = 10.0;      // TDRS measurement accuracy [m]

  const int    Digits[8] =              // Number of decimals for 
               { 2,2,2,5,5,5, 4,4 };    // printout of estimation parameters

  const char*  Label[8] =           // Estimation parameter names 
               { "x  [m]  ", "y  [m]  ", "z  [m]  ",
                 "vx [m/s]", "vy [m/s]", "vz [m/s]",
                 "CD      ", "CR      " };

  // Variables

  ifstream  inp;                         // Input file
  int       i,k,i_sat,i_Sta,iterat;      // Loop counters
  int       n_sat;                       // Actual num. of sats. (user+TDRS)
  int       n_Sta;                       // Actual num. of ground stations
  int       n_iterat;                    // Number of iterations
  int       n_obs;                       // Number of observations
  int       TdrsNo[N_sat];               // TDRS s/c number 
  int       I_TDRS;                      // Relay satellite index (1,..,n_sat)
  int       I_Sta;                       // Station index
  int       offset;                      // Index offset
  double    Mjd0_TT;                     // Epoch (TT)
  double    Mjd_UTC, Mjd_TT;             // Time (UTC,TT)
  Vector    y0[N_sat];                   // Epoch states of user and relay sats
  Vector    Sigy0[N_sat];                // Standard deviation
  Vector    p[N_sat];                    // S/C parameters of user/real sats
  Vector    Sigp[N_sat];                 // Standard deviation
  AuxParam  Aux[N_sat];                  // Auxil. params. of user and relay sat
  Vector    drho_dx[N_sat];              // Partials w.r.t. state and params
  TrjData   Trj[N_sat];                  // Trajectory of user and relay sat
  StaParam  Sta[N_Sta];                  // Ground stations
  double    rho;                         // Range
  double    res, Sum_res, Sum_ressqr;    // Residual, sum, sum of squares
  double    Mean, RMS;                   // Residuals statistics
  ObsType   tmp;                         // Auxiliary variable
  vector<ObsType>  Obs;                  // Observations list
  
  
  // Read setup parameters from file
   
  if (argc>1) inp.open(argv[1]); else inp.open("TDRSOD.inp"); 
  if (!inp) { 
    cerr << "ERROR: Could not open TDRSOD setup file" << endl;
    exit(1);
  };
  
  GetSetup ( inp, 
             Mjd0_TT, n_iterat,
             n_sat, TdrsNo, y0, Sigy0, p, Sigp, Aux,
             n_Sta, Sta );
  
  inp.close();

  // Allocation of dynamic variables and objects

  const int n_est = 8*n_sat;             // Number of estimation parameters

  LSQ       OrbEst(n_est);               // Least squares system
  Vector    X_apr(n_est),X(n_est);       // Estimation parameter vector
  Vector    dzdX(n_est);                 // Partials of modelled observations
  Vector    dX(n_est),SigX(n_est);       // Correction and standard deviation
  Matrix    P_apr(n_est,n_est);          // A priori covariance

  // Instantiation of integrator objects

  for (i_sat=0;i_sat<n_sat;i_sat++) {
    Trj[i_sat].Define(Deriv,VarEqn,8,&Aux[i_sat]);
  };


  // A priori estimation parameter vector and covariance

  for (i_sat=0;i_sat<n_sat;i_sat++) {
    offset = 8*i_sat;
    for (i=0;i<6;i++) X_apr(offset+i  ) = y0[i_sat](i);
    for (i=0;i<2;i++) X_apr(offset+6+i) = p[i_sat](i);
    for (i=0;i<6;i++) SigX(offset+i  )  = Sigy0[i_sat](i);
    for (i=0;i<2;i++) SigX(offset+6+i)  = Sigp[i_sat](i);
  };

  P_apr = Diag(SigX)*Diag(SigX);

  X = X_apr;


  // Read all observations from input file to observations vector
  
  if (argc>2) inp.open(argv[2]); else inp.open("TDRSOD.dat"); 
  if (!inp) { 
    cerr << "ERROR: Could not open TDRSOD tracking data file" << endl;
    exit(1);
  };
  
  inp.ignore(255,'\n');       // Skip header

  for (;;) {
    GetObs (inp, tmp);
    if (inp.fail()) break;    // Terminate at end of file
    Obs.push_back(tmp);
  }

  inp.close();
  

  // Header
  
  cout << "TDRS Orbit Determination" << endl 
       << endl
       << Obs.size() << " measurements read from tracking data file" << endl
       << endl;


  // Iteration loop

  for (iterat=1; iterat<=n_iterat; iterat++) {

    // Iteration header

    cout << "Iteration " << iterat << endl << endl
         << "     Date         UTC        Sta TDRS"
         << "   obs  [m]    comp [m]    o-c [m]"
         << endl;

    // Initialize least squares system and measurement statistics

    OrbEst.Init(X_apr-X,P_apr);

    n_obs      = 1;
    Sum_res    = 0.0;
    Sum_ressqr = 0.0;

    // Initialize trajectory objects with epoch and epoch state

    for (i_sat=0;i_sat<n_sat;i_sat++)  
      Trj[i_sat].Init(Mjd0_TT,y0[i_sat]);

    // Process measurements

    for (i=0;i<(int)Obs.size();i++) {
    
      // Time step
      Mjd_UTC = Obs[i].Mjd_UTC;
      Mjd_TT  = Mjd_UTC + IERS::TT_UTC(Mjd_UTC)/86400.0;
    
      // Integrate user and relay satellite trajectories to ground-received time
      // of current measurement
      for (i_sat=0;i_sat<n_sat;i_sat++)  
        Trj[i_sat].Integ(Mjd_TT);
    
      // Select TDRS satellite and ground station
      I_TDRS = -1;
      for (i_sat=1;i_sat<n_sat;i_sat++) 
        if (TdrsNo[i_sat]==Obs[i].TdrsNo) {
          I_TDRS = i_sat; break;
        };
      if (I_TDRS<0) continue;   // Unknown TDRS satellite; skip observation
      
      // Select ground station
      I_Sta = -1;
      for (i_Sta=0;i_Sta<n_Sta;i_Sta++) 
        if (Sta[i_Sta].StaNo==Obs[i].StaNo) {
          I_Sta = i_Sta; break;
        };
      if (I_Sta<0) continue;    // Unknown station; skip this observation

      // Compute 4-way range and partials using interpolated trajectories
      // of user and relay satellite
      Range_4W ( Mjd_UTC, Trj[0], Trj[I_TDRS], Sta[I_Sta].Loc, 
                 rho, drho_dx[0], drho_dx[I_TDRS] );

      // Residual (observed - computed) 
      res = Obs[i].Range_4w - rho;

      // Partials of modelled observations

      dzdX   = 0;
      offset = 8*I_TDRS;
      for (k=0;k<8;k++){
        dzdX(k)        = drho_dx[0](k);       // Partials w.r.t. user satellite
        dzdX(offset+k) = drho_dx[I_TDRS](k);  // Partials w.r.t. relay satellite
      }

      // Statistics
      n_obs      += 1;
      Sum_res    += res;
      Sum_ressqr += res*res;

      // Process data equation
      OrbEst.Accumulate ( dzdX, res, sigma_range );

      // Output
      cout << "  " << Date(Mjd_UTC)
           << fixed
           << setw(6) << Obs[i].StaNo
           << setw(4) << Obs[i].TdrsNo
           << setprecision(4) << setw(13) << Obs[i].Range_4w/1000.0
           << setprecision(4) << setw(12) << (rho)/1000.0
           << setprecision(2) << setw( 9) << (Obs[i].Range_4w-rho)
           << endl;

    };

    // Solve least-squares system

    OrbEst.Solve(dX); 
    SigX = OrbEst.StdDev(); 
    X += dX;

    // Update epoch states and parameters

    for (i_sat=0;i_sat<n_sat;i_sat++) {
      offset = 8*i_sat;
      y0[i_sat]     = X.slice(offset+0,offset+5);
      Aux[i_sat].CD = X(offset+6);
      Aux[i_sat].CR = X(offset+7);
    };

    // Residuals statistics

    Mean = Sum_res / n_obs;
    RMS  = sqrt ( Sum_ressqr / n_obs  -  Mean*Mean );

    // Print results

    cout << setw(70) << "_______" << endl
         << setw(55) << " " << "Mean" << setprecision(2) << setw(11) << Mean
         << " m" << endl
         << setw(55) << " " << "RMS " << setprecision(2) << setw(11) << RMS
         << " m" << endl
         << endl << endl
         << "Results of iteration " << iterat << endl << endl
         << "  Parameter            " 
         << "a priori   tot.corr. last corr.       final       sigma"
         << endl;
    for (i=0;i<n_est;i++) {
      cout << "  " << Label[i%8];
      if (i/8==0) cout << " User  ";
      else        cout << " TDRS-" << TdrsNo[i/8];
      cout << setprecision(Digits[i%8])
           << setw(14) << X_apr(i)
           << setw(11) << X(i)-X_apr(i) 
           << setw(11) << dX(i) 
           << setw(14) << X(i)
           << setw(11) << SigX(i)
           << endl;
    };
    cout << endl << endl;

  }

  return 0;

}
