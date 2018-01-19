//------------------------------------------------------------------------------
//
// Exercise_7_1.cpp
// 
// Purpose: 
//
//   Satellite Orbits - Models, Methods, and Applications 
//   Exercise 7-1: State transition matrix
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
#include "SAT_VecMat.h"

using namespace std;

//------------------------------------------------------------------------------
//
// Global types and data
//
//------------------------------------------------------------------------------

// Record for passing global data between Deriv and the calling program 

struct AuxParam {
  double  Mjd_0;     // Reference epoch
  int     n_a,m_a;   // Degree and order of gravity field for acceleration
  int     n_G,m_G;   // Degree and order of gravity field for gradient
};


//------------------------------------------------------------------------------
//
// Accel
//
// Purpose:
//
//   Computes the acceleration of an Earth orbiting satellite due to 
//   the Earth's harmonic gravity field up to degree and order 10
//
// Input/Output:
//
//   Mjd_UT      Modified Julian Date (Universal Time)
//   r           Satellite position vector in the true-of-date system
//   n,m         Gravity model degree and order 
//   <return>    Acceleration (a=d^2r/dt^2) in the true-of-date system
//
//------------------------------------------------------------------------------

Vector Accel ( double Mjd_UT, const Vector& r, int n, int m )
{

  // Variables 

  Matrix  U(3,3);

  // Earth rotation matrix 

  U = R_z(GMST(Mjd_UT));

  // Acceleration due to harmonic gravity field

  return  AccelHarmonic ( r,U, GM_Earth, Grav.R_ref,Grav.CS, n,m );

}


//------------------------------------------------------------------------------
//
// Gradient
//
// Purpose:
//
//   Computes the gradient of the Earth's harmonic gravity field 
//
// Input/Output:
//
//   Mjd_UT      Modified Julian Date (Universal Time)
//   r           Satellite position vector in the true-of-date system
//   n,m         Gravity model degree and order 
//   <return>    Gradient (G=da/dr) in the true-of-date system
//
//------------------------------------------------------------------------------

Matrix Gradient ( double Mjd_UT, const Vector& r, int n, int m )
{

  // Constants

  const double d = 1.0;   // Position increment [m]
  
  // Variables 

  int     i;
  Vector  a(3),da(3),dr(3);
  Matrix  U(3,3),G(3,3);

  // Earth rotation matrix 

  U = R_z(GMST(Mjd_UT));

  // Acceleration
  
  a = AccelHarmonic ( r,U, GM_Earth, Grav.R_ref,Grav.CS, n,m );

  // Gradient

  for (i=0; i<=2; i++) {
    // Set offset in i-th component of the position vector
    dr = 0.0; dr(i) = d;
    // Acceleration difference
    da = AccelHarmonic ( r+dr/2,U, GM_Earth, Grav.R_ref,Grav.CS, n,m ) -   
         AccelHarmonic ( r-dr/2,U, GM_Earth, Grav.R_ref,Grav.CS, n,m );
//  da = AccelHarmonic ( r+dr,U, GM_Earth, Grav.R_ref,Grav.CS, n,m ) -  a;
    // Derivative with respect to i-th axis
    G.SetCol(i,da/d);
  }

  return G;

}


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
//   yPhi        (6+36)-dim vector comprising the state vector (y) and the
//               state transition matrix (Phi) in column wise storage order
//   yPhip       Derivative of yPhi
//   pAux        Pointer; pAux is expected to point to a variable of type
//               AuxDataRecord, which is used to communicate with the other
//               program sections and to hold data between subsequent calls
//               of this function
//
//------------------------------------------------------------------------------

void VarEqn ( double t, const Vector& yPhi, Vector& yPhip, void* pAux )
{

  // Variables

  AuxParam* p;             // Auxiliary data pointer
  int       i,j;           // Loop counters
  double    Mjd;           // Modified Julian Date
  Vector    r(3), v(3);    // Position, velocity
  Vector    a(3);          // Acceleration
  Matrix    G(3,3);        // Gradient of acceleration
  Matrix    Phi(6,6);      // State transition matrix
  Matrix    Phip(6,6);     // Time derivative of state transition matrix
  Matrix    dfdy(6,6);     // 
      
  // Pointer to auxiliary data record
  
  p = static_cast<AuxParam*>(pAux);

  // Time

  Mjd = (*p).Mjd_0 + t/86400.0;

  // State vector components

  r = yPhi.slice(0,2);
  v = yPhi.slice(3,5);

  // State transition matrix

  for (j=0; j<=5; j++) Phi.SetCol(j,yPhi.slice(6*(j+1),6*(j+1)+5));
  
  // Acceleration and gradient

  a = Accel    ( Mjd,r, (*p).n_a, (*p).m_a );
  G = Gradient ( Mjd,r, (*p).n_G, (*p).m_G );

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
  
  // Derivative of combined state vector and state transition matrix
  
  for (i=0; i<=2; i++) {
    yPhip(i)   = v(i);                    // dr/dt(i)
    yPhip(i+3) = a(i);                    // dv/dt(i)
  };
  for (i=0; i<=5; i++) 
    for (j=0; j<=5; j++)
      yPhip(6*(j+1)+i  ) = Phip(i,j);     // dPhi/dt(i,j)

};


//------------------------------------------------------------------------------
//
// Max: finds the largest element of a matrix
//
//------------------------------------------------------------------------------

double Max(const Matrix& M) {

  const int n = M.size1();
  const int m = M.size2();
  double max = 0.0;

  for (int i=0; i<n; i++) 
    for (int j=0;j<m; j++)
      if (fabs(M(i,j))>max) max=fabs(M(i,j));

  return max;
}

//------------------------------------------------------------------------------
//
// InvSymp: 
//
// Purpose:
//
//   Inverts a symplectic matrix
//
// Input/output:
//
//   A         Sympletctic (even-dimensional square) matrix 
//   <return>  Inverse of A
//
//                      ( +(A_22)^T  -(A_12)^T )            ( A_11  A_12 )
//             A^(-1) = (                      )  with  A = (            )
//                      ( -(A_21)^T  +(A_11)^T )            ( A_21  A_22 )
//
//------------------------------------------------------------------------------

Matrix InvSymp (const Matrix& A) {

  const int N = A.size1();
  const int n = N/2;
  Matrix Inv(N,N);
  
  if (A.size2()!=N || 2*n!=N ) {
    cout << "ERROR: Invalid shape in InvSymp(Matrix)";
    exit(1);
  };

  for (int i=0; i<n; i++) {
    for (int j=0;j<n; j++) {
      Inv(i  ,j  ) = +A(j+n,i+n);    
      Inv(i  ,j+n) = -A(j  ,i+n);
      Inv(i+n,j  ) = -A(j+n,i  );
      Inv(i+n,j+n) = +A(j  ,i  );
    };
  };

  return Inv;

}


//------------------------------------------------------------------------------
//
// Main program
//
//------------------------------------------------------------------------------

int main() {

  // Constants
  
  const double  relerr = 1.0e-13;                     // Relative and absolute
  const double  abserr = 1.0e-6;                      // accuracy requirement

  const double  a   = R_Earth+650.0e3;                // Semi-major axis [m]
  const double  e   = 0.001;                          // Eccentricity
  const double  inc = Rad*51.0;                       // Inclination [rad]
  
  const double  n   = sqrt(GM_Earth/(a*a*a));         // Mean motion

  const Matrix Scale=Diag(Vector(1,1,1,1/n,1/n,1/n)); // Velocity
  const Matrix scale=Diag(Vector(1,1,1,  n,  n,  n)); // normalization

  const double t_end = 86400.0; 
  const double step  =   300.0;

  // Variables

  int       i,j;                     // Loop counters
  double    Mjd0;                    // Epoch (Modified Julian Date)
  double    t;                       // Integration time [s]
  double    max, max_Kep;            // Matrix error norm
  AuxParam  Aux1,Aux2;               // Auxiliary parameter records
  DE        Int1(VarEqn,42,&Aux1);   // Object for integrating the variat. eqns.
  DE        Int2(VarEqn,42,&Aux2);   // 
  Vector    y0(6), y(6);             // State vector
  Vector    yPhi1(42), yPhi2(42);    // State and transition matrix vector
  Matrix    Phi1(6,6), Phi2(6,6);    // State transition matrix
  Matrix    Phi_Kep(6,6);            // State transition matrix (Keplerian)
  
  
  // Epoch state and transition matrix

  Mjd0 = MJD_J2000;
  
  y0 = State ( GM_Earth, Vector(a,e,inc,0.0,0.0,0.0) );

  for (i=0; i<=5; i++) {
    yPhi1(i) = y0(i);
    for (j=0; j<=5; j++)  yPhi1(6*(j+1)+i) = ( i==j ? 1 : 0 ); 
  };

  yPhi2 = yPhi1;


  // Model parameters for reference solution (full 10x10 gravity model) and
  // simplified variational equations

  Aux1.Mjd_0 = Mjd0;  // Reference epoch
  Aux1.n_a   = 10;    // Degree of gravity field for trajectory computation
  Aux1.m_a   = 10;    // Order of gravity field for trajectory computation
  Aux1.n_G   = 10;    // Degree of gravity field for variational equations
  Aux1.m_G   = 10;    // Order of gravity field for variational equations

  Aux2.Mjd_0 = Mjd0;  // Reference epoch
  Aux2.n_a   =  2;    // Degree of gravity field for trajectory computation
  Aux2.m_a   =  0;    // Order of gravity field for trajectory computation
  Aux2.n_G   =  2;    // Degree of gravity field for variational equations
  Aux2.m_G   =  0;    // Order of gravity field for variational equations

  // Initialization
    
  t = 0.0;
  Int1.Init ( t, relerr, abserr );
  Int2.Init ( t, relerr, abserr );
  
  // Header
  
  cout << "Exercise 7-1: State transition matrix" << endl
       << endl 
       << "  Time [s]  Error (J2)  Error(Kep)" << endl;
 
  // Steps
  
  while (t<t_end) {

    // New output time

    t += step;

    // Integration

    Int1.Integ ( t, yPhi1 ); // Reference
    Int2.Integ ( t, yPhi2 ); // Simplified
  
    // Extract and normalize state transition matrices

    for (j=0; j<=5; j++) Phi1.SetCol(j,yPhi1.slice(6*(j+1),6*(j+1)+5));
    for (j=0; j<=5; j++) Phi2.SetCol(j,yPhi2.slice(6*(j+1),6*(j+1)+5));

    Phi1 = Scale*Phi1*scale;
    Phi2 = Scale*Phi2*scale;

    // Keplerian state transition matrix (normalized)

    TwoBody ( GM_Earth, y0,t, y,Phi_Kep); 
    Phi_Kep = Scale*Phi_Kep*scale;

    // Matrix error norms

    max     = Max ( Id(6)-Phi1*InvSymp(Phi2   ) );
    max_Kep = Max ( Id(6)-Phi1*InvSymp(Phi_Kep) );
  
    // Output 
    
    cout << fixed << setprecision(2) 
         << setw(10) << t
         << setw(10) << max << setw(12) << max_Kep << endl;

  };

  // Output

  cout << endl
       << "Time since epoch [s]" << endl
       << endl
       << fixed << setprecision(2) << setw(12) << t << endl
       << endl;
  
  cout << "State transition matrix (reference)" << endl 
       << endl
       << fixed << setprecision(5) << setw(12) << Phi1 << endl;

  cout << "J2 State transition matrix" << endl 
       << endl
       << fixed << setprecision(5) << setw(12) << Phi2 << endl;

  cout << "Keplerian State transition matrix" << endl 
       << endl
       << fixed << setprecision(5) << setw(12) << Phi_Kep << endl;

  return 0;

}
