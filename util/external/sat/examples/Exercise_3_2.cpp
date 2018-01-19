//------------------------------------------------------------------------------
//
// Exercise_3_2.cpp
// 
// Purpose: 
//
//   Satellite Orbits - Models, Methods, and Applications
//   Exercise 3-2: Lunar ephemerides
//
// Note:
//
//   This program rigorously evaluates the heights were different acceleration
//   types will balance for a specific spacecraft state vector and parameters. 
//   The results will thus deviate from the estimates given in the book, that
//   are based on analytical estimates.
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
#include "SAT_Force.h"
#include "SAT_Time.h"
#include "SAT_VecMat.h"

using namespace std;


//------------------------------------------------------------------------------
//
// Cheb3D class (specification and implementation)
//
//------------------------------------------------------------------------------

// Chebyshev approximation of 3-dimensional vectors

class Cheb3D {
  public:
    // Constructor
    Cheb3D ( int n, double ta, double tb, 
             const Vector& cx, const Vector& cy, const Vector& cz )
      : N(n), Ta(ta), Tb(tb)
    {
      Cx=cx; Cy=cy; Cz=cz;
    };
    // Evaluation 
    Vector  Value (double t) const;
  private:
    // Elements
    int     N;       // Number of coefficients
    double  Ta;      // Begin interval
    double  Tb;      // End interval
    Vector  Cx;      // Coefficients of Chebyshev polyomial (x-coordinate)
    Vector  Cy;      // Coefficients of Chebyshev polyomial (y-coordinate)
    Vector  Cz;      // Coefficients of Chebyshev polyomial (z-coordinate)
};


//
//  Evaluation of the Chebyshev approximation of a three-dimensional vector
//

Vector Cheb3D::Value (double t) const
{
  
  // Variables
  int     i;
  Vector  f1(3), f2(3), old_f1(3);
  double  tau;

  // Check validity
  if ( (t<Ta) || (Tb<t) ) {
    cerr << "ERROR: Time out of range in Cheb3D::Value" << endl;
    exit(1);
  };

  // Clenshaw algorithm
  tau = (2.0*t-Ta-Tb)/(Tb-Ta);  

  f1 = 0.0;
  f2 = 0.0;  

  for (i=N-1; i>=1; i--) {
    old_f1 = f1; 
    f1 = 2.0*tau*f1-f2+Vector(Cx(i),Cy(i),Cz(i));  
    f2 = old_f1;
  };
  
  return tau*f1-f2+Vector(Cx(0),Cy(0),Cz(0));

}


//------------------------------------------------------------------------------
//
// Main program
//
//------------------------------------------------------------------------------

int main() {

  // Constants
  
  const int    N_Step =   8;
  const double Step   = 0.5; // [d]


  // Sample Chebyshev coefficients from JPL DE405 for geocentric lunar 
  // coordinates in the ITRF(EME2000) frame valid from 2006/03/14 TDB 
  // to 2006/03/18 TDB 

  const int    N_coeff = 13;          // Number of coefficients 
  const double t1 = 53808.0;          // MJD at start of interval 
  const double t2 = 53812.0;          // MJD at end   of interval 
  const double Cx_moon[N_coeff] = {   // Coefficients for x-coordinate
 -0.383089044877016277e+06, 0.218158411754834669e+05, 0.179067292901463843e+05,
 -0.836928063411765777e+02,-0.628266733052023696e+02,-0.459274434235101225e+00,
  0.491167202819885532e-01, 0.770804039287614762e-03,-0.125935992206166816e-03,
  0.500271026610991370e-05, 0.107044869185752331e-05, 0.172472464343636242e-08,
 -0.269667589576924680e-08                                                   };
  const double Cy_moon[N_coeff] = {   // Coefficients for y-coordinate
 -0.379891721705081436e+05,-0.143611643157166138e+06, 0.187126702787245881e+04, 
  0.112734362473135207e+04, 0.932891213817359177e+00,-0.191932684130578513e+01, 
 -0.266517663331897990e-01, 0.104558913448630337e-02,-0.359077689123857890e-04, 
 -0.123405162037249834e-04, 0.180479239596339495e-06, 0.525522632333670539e-07, 
  0.543313967008773005e-09                                                   }; 
  const double Cz_moon[N_coeff] = {   // Coefficients for z-coordinate
 -0.178496690739133737e+05,-0.788257550331743259e+05, 0.880684692614081882e+03, 
  0.618395886330471512e+03, 0.103331218594995988e+01,-0.104949867328178592e+01,
 -0.150337371962561087e-01, 0.569056416308259317e-03,-0.186297523286550968e-04,
 -0.680012420653791955e-05, 0.902057208454410917e-07, 0.287891446432139173e-07, 
  0.319822827699973363e-09                                                   };

  // Chebyshev approximation of lunar coordinates

  const Cheb3D MoonCheb ( N_coeff, t1,t2,                 // Order and interval
                          Vector(&Cx_moon[0],N_coeff),    // Coefficients
                          Vector(&Cy_moon[0],N_coeff),      
                          Vector(&Cz_moon[0],N_coeff)  );     
  
  // Variables

  int       i;
  double    Mjd0, Mjd_TT;
  Vector    r(3);

  // Epoch

  Mjd0 = Mjd(2006,03,14,00,00,0.0);

  // Output
  
  cout << "Exercise 3-2: Lunar Ephemerides " << endl << endl;

  cout << " Moon position from low precision analytical theory" << endl;
  cout << " Date [TT]                 " << " Position [km] " << endl;

  for (i=0;i<=N_Step;i++) {
    Mjd_TT = Mjd0 + i*Step;
    r = Moon (Mjd_TT)/1000.0;
    cout << " " << Date(Mjd_TT) 
         << fixed << setprecision(3) << setw(14) << r << endl; 
  };

  cout << endl << " Moon position from DE405" << endl;
  cout << " Date [TT]                 " << " Position [km] " << endl;

  for (i=0;i<=N_Step;i++) {
    Mjd_TT = Mjd0 + i*Step;
    r = MoonCheb.Value (Mjd_TT);
    cout << " " << Date(Mjd_TT) 
         << fixed << setprecision(3) << setw(14) << r << endl; 
  };

  return 0;

}
