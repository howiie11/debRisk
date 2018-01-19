//------------------------------------------------------------------------------
//
// SAT_Time.h
// 
// Purpose:
//
//    Time and date computation
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
//   2005/04/14  OMO  Final version (2nd reprint)
//
// (c) 1999-2005  O. Montenbruck, E. Gill
//
//------------------------------------------------------------------------------

#ifndef INC_SAT_TIME_H
#define INC_SAT_TIME_H

#include <iostream>


//------------------------------------------------------------------------------
//
// Mjd
//
// Purpose:
//
//   Modified Julian Date from calendar date and time
//
// Input/output:
//
//   Year      Calendar date components
//   Month
//   Day
//   Hour      Time components (optional)
//   Min
//   Sec
//   <return>  Modified Julian Date
//
//------------------------------------------------------------------------------

double Mjd ( int Year,   int Month, int Day, 
             int Hour=0, int Min=0, double Sec=0.0 );


//------------------------------------------------------------------------------
//
// CalDat
//
// Purpose:
//
//   Calendar date and time from Modified Julian Date
//
// Input/output:
//
//   Mjd       Modified Julian Date
//   Year      Calendar date components
//   Month
//   Day
//   Hour      Time components
//   Min
//   Sec
//
//------------------------------------------------------------------------------

void CalDat ( double Mjd, 
              int& Year, int& Month, int& Day,
              int& Hour, int& Min, double& Sec );



//------------------------------------------------------------------------------
//
// Date (class definition)
//
// Purpose:
//
//   Auxiliary class for date and time output
//
//------------------------------------------------------------------------------

class Date
{
  public:
    Date(double Mjd);                                             // Constructor    
    friend std::ostream& operator<< (std::ostream& os, const Date& D); // Output
  private:
    double mjd;
};


#endif  // include blocker
