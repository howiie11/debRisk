//------------------------------------------------------------------------------
//
// Exercise_4_2.cpp
// 
// Purpose: 
//
//   Satellite Orbits - Models, Methods, and Applications
//   Exercise 4-2: Gauss-Jackson 4th-order predictor 
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

#include <cmath>
#include <iostream>
#include <iomanip>

#include "GNU_iomanip.h"

#include "SAT_Const.h"
#include "SAT_Kepler.h"
#include "SAT_VecMat.h"

using namespace std;


//------------------------------------------------------------------------------
//
// GJ4P class 
//
//------------------------------------------------------------------------------

// Function prototype for second order differential equations
// void f (double t, const Vector& r, const Vector& v, Vector& a)

typedef void (*GJ4Pfunct)(
  double        t,     // Independent variable
  const Vector& r,     // Position vector 
  const Vector& v,     // Velocity vector  r'=v
  Vector&       a,     // Acceleration     r''=a=f(t,r,v)
  void*         pAux   // Pointer to auxiliary data used within f
);


// Specification

class GJ4P
{
  public:

    // Constructor
    GJ4P (
      GJ4Pfunct   f_,        // Differential equation
      int         n_eqn_,    // Dimension
      void*       pAux_      // Pointer to auxiliary data
      )
    : n_eqn(n_eqn_), f(f_), pAux(pAux_)
    {};
    
    // Initialization step
    void Init (
      double        t_0,     // Initial value of independent variable
      const Vector& r_0,     // Initial value r_0=r(t_0)
      const Vector& v_0,     // Initial value v_0=dr/dt(t_0)
      double        h_       // Step size
    );

    // Integration step
    void Step (         
      double&  t,            // Independent variable; updated by t+h
      Vector&  r,            // Value of r(t); updated by r(t+h)
      Vector&  v             // Value of v(t)=dr/dt(t); updated by v(t+h)
    );

  private:

    // 4th order Runge-Kutta step
    void RK4 (         
      double&  t,            // Independent variable; updated by t+h
      Vector&  r,            // Value of r(t); updated by r(t+h)
      Vector&  v,            // Value of v(t)=dr/dt(t); updated by v(t+h)
      double   h             // Step size
    );

    // Elements
    int         n_eqn;       // Dimension
    GJ4Pfunct   f;           // Differential equation
    double      h;           // Step size
    void*       pAux;        // Pointer to auxiliary data requird by f
    Vector      S2,S1;       // First and second sum of acceleration
    Vector      D[4];        // Backward differences of acceleration at t
    Vector      d[4];        // Backward differences of acceleration at t+h
    Vector      r_p,v_p;     // Predictor

};


//
// 4th order Runge-Kutta step for 2nd order differential equation
//

void GJ4P::RK4 (double&  t, Vector&  r, Vector&  v, double h )
{
  Vector v_1, v_2, v_3, v_4;
  Vector a_1, a_2, a_3, a_4;

  v_1 = v;              f( t      , r            , v_1, a_1, pAux ); 
  v_2 = v+(h/2.0)*a_1;  f( t+h/2.0, r+(h/2.0)*v_1, v_2, a_2, pAux );
  v_3 = v+(h/2.0)*a_2;  f( t+h/2.0, r+(h/2.0)*v_2, v_3, a_3, pAux ); 
  v_4 = v+h*a_3;        f( t+h    , r+h*v_3      , v_4, a_4, pAux );
  
  t = t + h;
  r = r + (h/6.0)*( v_1 + 2.0*v_2 + 2.0*v_3 + v_4 );
  v = v + (h/6.0)*( a_1 + 2.0*a_2 + 2.0*a_3 + a_4 );

};


//
// Initialization of backwards differences from initial conditions
//

void GJ4P::Init(double t_0, const Vector& r_0, const Vector& v_0, double h_)
{
  // Order of method
  
  const int m = 4;

  // Coefficients gamma/delta of 1st/2nd order Moulton/Cowell corrector method

  const double gc[m+1] = {+1.0, -1/2.0, -1/12.0, -1/24.0, -19/720.0 };
  const double dc[m+2] = {+1.0,   -1.0, +1/12.0,     0.0,  -1/240.0, -1/240.0 };

  int       i,j;
  double    t = t_0;
  Vector    r = r_0;
  Vector    v = v_0;

  // Save step size  

  h = h_;     

  // Create table of accelerations at past times t-3h, t-2h, and t-h using
  // RK4 steps

  f(t,r,v,D[0],pAux);     // D[i]=a(t-ih)
  for (i=1;i<=m-1;i++) {
    RK4(t,r,v,-h);  f(t,r,v,D[i],pAux);   
  };

  // Compute backwards differences
  
  for (i=1;i<=m-1;i++) 
    for (j=m-1;j>=i;j--) D[j] = D[j-1]-D[j];

  // Initialize backwards sums using 4th order GJ corrector

  S1 = v_0/h;              for (i=1;i<=m  ;i++) S1 -= gc[i]*D[i-1];
  S2 = r_0/(h*h)-dc[1]*S1; for (i=2;i<=m+1;i++) S2 -= dc[i]*D[i-2];

};


//
// Step from t to t+h
//

void GJ4P::Step (double& t, Vector& r, Vector& v) 
{
  // Order of method
  
  const int m = 4;  

  // Coefficients gamma/delta of 1st/2nd order Bashforth/Stoermr predictor
  
  const double gp[m+1] = {+1.0, +1/2.0, +5/12.0,  +3/8.0, +251/720.0 };
  const double dp[m+2] = {+1.0,    0.0, +1/12.0, +1/12.0,  +19/240.0,  +3/40.0 };

  int i;
  
  // 4th order predictor

  r_p = dp[0]*S2; for(i=2;i<=m+1;i++) r_p += dp[i]*D[i-2]; r_p = (h*h)*r_p;
  v_p = gp[0]*S1; for(i=1;i<=m  ;i++) v_p += gp[i]*D[i-1]; v_p =     h*v_p;

  // Update backwards difference table

  f ( t+h, r_p,v_p, d[0], pAux );               // Acceleration at t+h
  for (i=1;i<=m-1;i++) d[i]=d[i-1]-D[i-1];      // New differences at t+h
  for (i=0;i<=m-1;i++) D[i]=d[i];               // Update differences 
  S1 += d[0];  S2 += S1;                        // Update sums

  // Update independent variable and solution

  t = t + h;
  r = r_p;
  v = v_p;

};


//------------------------------------------------------------------------------
//
// f_Kep3D
//
// Purpose:
// 
//   Computes the second time derivative of the position vector for the 
//   normalized (GM=1) Kepler's problem in three dimensions
//
// Note:
//
//   pAux is expected to point to an integer variable that will be incremented
//   by one on each call of f_Kep3D
//
//------------------------------------------------------------------------------

void f_Kep3D ( double t, const Vector& r, const Vector& v, 
               Vector& a, void* pAux )
{
  
  // Pointer to auxiliary integer variable used as function call counter
  int* pCalls = static_cast<int*>(pAux);
  
  // 2nd order derivative d^2(r)/dt^2 of the position vector
  a = -r/(pow(Norm(r),3));
  
  // Increment function call count
  (*pCalls)++;

};


//------------------------------------------------------------------------------
//
// Main program
//
//------------------------------------------------------------------------------

int main() {

  // Constants
  
  const double  GM    = 1.0;                   // Gravitational coefficient
  const double  e     = 0.1;                   // Eccentricity
  const double  t_end = 20.0;                  // End time
  const Vector  Kep(1.0,e,0.0,0.0,0.0,0.0);    // (a,e,i,Omega,omega,M)
  const Vector  y_ref = State(GM,Kep,t_end);   // Reference solution

  const int     Steps[] = { 100, 300, 600, 1000, 1500, 2000, 3000, 4000 };
  
  // Variables
  
  int     nCalls;                              // Function call count
  int     iCase;
  double  t,h;                                 // Time and step size
  Vector  y(6);                                // State vector
  Vector  r(3);
  Vector  v(3);

  GJ4P    Orbit(f_Kep3D,3,&nCalls);            // Object for integrating the
                                               // 2nd order diff. equation
                                               // defined by f_Kep3D using the
                                               // 4th-order GJ predictor
  
  // Header 

  cout << "Exercise 4-2: Gauss-Jackson 4th-order predictor" 
       << endl << endl;
  cout << "  Problem D1 (e=0.1)" << endl << endl;
  cout << "  N_fnc   Accuracy   Digits " << endl;
    
  // Loop over test cases

  for (iCase=0; iCase<8; iCase++) {
  
    // Step size
    h = t_end/Steps[iCase];

    // Initial values
    t = 0.0;
    r = Vector( 1.0-e, 0.0, 0.0 );
    v = Vector( 0.0, sqrt((1+e)/(1-e)), 0.0 );
    nCalls = 0;

    // Integration from t=t to t=t_end
    Orbit.Init ( t, r, v, h );
    for (int i=1; i<=Steps[iCase]; i++)  
      Orbit.Step( t, r, v );
    y = Stack(r,v);

    // Output
    cout << fixed  << setw(6) << nCalls
         << scientific << setprecision(3) << setw(13)
         << Norm(y-y_ref) 
         << fixed << setprecision(2) << setw(7)
         << -log10(Norm(y-y_ref)) << endl;
  
  };

  return 0;

}
