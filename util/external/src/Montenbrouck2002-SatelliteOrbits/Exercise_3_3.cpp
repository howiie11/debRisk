//------------------------------------------------------------------------------
//
// Exercise_3_3.cpp
// 
// Purpose: 
//
//   Satellite Orbits - Models, Methods, and Applications
//   Exercise 3-3: Accelerations
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

#include <cmath>
#include <iostream>
#include <iomanip>

#include "GNU_iomanip.h"

#include "SAT_Const.h"
#include "SAT_Force.h"
#include "SAT_VecMat.h"

using namespace std;

// Pegasus function prototype: double f(double x);
typedef double (*PegasusFunct) (double x);


//------------------------------------------------------------------------------
// A_Diff: Function computes the difference of acceleration due to SRP and DRG
//   h         Height [m]
//   <return>  Difference of acceleration
//------------------------------------------------------------------------------
double A_Diff (double h)
{
  double CD     =     2.3;      // Spacecraft parameters
  double CR     =     1.3;
  double Mjd_TT = 51269.0;      // State epoch

  Vector r(3);
  double dens;

  r = Vector (1.0, 0.0, 0.0 ) * (Grav.R_ref + h);
  dens = Density_HP (Mjd_TT,r);

  return (0.5*CD*dens*Grav.GM/(Grav.R_ref + h)) - (P_Sol*CR);
}


//------------------------------------------------------------------------------
//
// Pegasus: Root finder using the Pegasus method
//
// Input:
//
//   PegasusFunct  Pointer to the function to be examined
//
//   LowerBound    Lower bound of search interval
//   UpperBound    Upper bound of search interval
//   Accuracy      Desired accuracy for the root
//
// Output:
//
//   Root          Root found (valid only if Success is true)
//   Success       Flag indicating success of the routine
//
// References:                                                               
//
//   Dowell M., Jarratt P., 'A modified Regula Falsi Method for Computing    
//     the root of an equation', BIT 11, p.168-174 (1971).                   
//   Dowell M., Jarratt P., 'The "PEGASUS Method for Computing the root      
//     of an equation', BIT 12, p.503-508 (1972).                            
//   Engeln-Muellges G., Reutter F., 'Formelsammlung zur Numerischen           
//     Mathematik mit FORTRAN77-Programmen', Bibliogr. Institut,             
//     Zuerich (1986).                                                       
//
// Notes:
//
//   Pegasus assumes that the root to be found is bracketed in the interval
//   [LowerBound, UpperBound]. Ordinates for these abscissae must therefore
//   have different signs.
//
//------------------------------------------------------------------------------
void Pegasus ( PegasusFunct f,
               double LowerBound, double UpperBound, double  Accuracy,
               double& Root, bool& Success )
{
  //
  // Constants
  //
  const int MaxIterat = 30; 

  //
  // Variables
  //
  double x1 = LowerBound; double f1 = f(x1);
  double x2 = UpperBound; double f2 = f(x2);
  double x3 = 0.0;        double f3 = 0.0;
  
  int Iterat = 0; 


  // Initialization
  Success = false;
  Root    = x1;

  
  // Iteration
  if ( f1 * f2 < 0.0 )
    do 
    {
      // Approximation of the root by interpolation
      x3 = x2 - f2/( (f2-f1)/(x2-x1) ); f3 = f(x3);

      // Replace (x1,f2) and (x2,f2) by new values, such that
      // the root is again within the interval [x1,x2]
      if ( f3 * f2 <= 0.0 ) {
        // Root in [x2,x3]
        x1 = x2; f1 = f2; // Replace (x1,f1) by (x2,f2)
        x2 = x3; f2 = f3; // Replace (x2,f2) by (x3,f3)
      }
      else {
        // Root in [x1,x3]
        f1 = f1 * f2/(f2+f3); // Replace (x1,f1) by (x1,f1')
        x2 = x3; f2 = f3;     // Replace (x2,f2) by (x3,f3)
      }

      if (fabs(f1) < fabs(f2))
        Root = x1;
      else
        Root = x2;

      Success = (fabs(x2-x1) <= Accuracy);
      Iterat++;
    }
    while ( !Success && (Iterat<MaxIterat) );
}


//------------------------------------------------------------------------------
//
// Main program
//
//------------------------------------------------------------------------------

int main() {
 
  const double r_Moon = 384400.0e+03; // Geocentric Moon distance [m]
  const double J20_norm =  4.841e-04; // Normalized geopotential coefficients
  const double J22_norm =  2.812e-06;

  const double h1  =  150.0e+3;       // Start of interval [m]
  const double h2  = 2000.0e+3;       // Stop of interval [m]
  const double eps =  100.0;          // Accuracy [m]  

  bool Success = false;

  double r5,r,h;


  cout << "Exercise 3-3: Accelerations " << endl << endl;

  // Balance of accelerations due to drag and solar radiation pressure

  Pegasus(A_Diff, h1, h2, eps, h, Success);
  
  cout << " a_DRG > a_SRP  for " << setprecision(0) << fixed << setw(6) 
       << h/1000.0 << " km" << endl;

  // Balance of accelerations due to J22 and Moon

  r5 = 3.0/2.0 * Grav.GM/GM_Moon * pow(Grav.R_ref,2) * pow(r_Moon,3) * J22_norm;
  r = pow(r5,0.2);
  
  cout << " a_J22 > a_Moon for " << setprecision(0) << fixed << setw(6) 
       << (r - Grav.R_ref)/1000.0 << " km" << endl;

  // Balance of accelerations due to J22 and Sun

  r5 = 3.0/2.0 * Grav.GM/GM_Sun * pow(Grav.R_ref,2) * pow(AU,3) * J22_norm;
  r = pow(r5,0.2);

  cout << " a_J22 > a_Sun  for " << setprecision(0) << fixed << setw(6) 
       << (r - Grav.R_ref)/1000.0 << " km" << endl;

  // Balance of accelerations due to J20 and Moon

  r5 = 3.0/2.0 * Grav.GM/GM_Moon * pow(Grav.R_ref,2) * pow(r_Moon,3) * J20_norm;
  r = pow(r5,0.2);
  
  cout << " a_J20 > a_Moon for " << setprecision(0) << fixed << setw(6) 
       << (r - Grav.R_ref)/1000.0 << " km" << endl;

  // Balance of accelerations due to J20 and Sun

  r5 = 3.0/2.0 * Grav.GM/GM_Sun * pow(Grav.R_ref,2) * pow(AU,3) * J20_norm;
  r = pow(r5,0.2);

  cout << " a_J20 > a_Sun  for " << setprecision(0) << fixed << setw(6) 
       << (r - Grav.R_ref)/1000.0 << " km" << endl;

  return 0;

}
