//------------------------------------------------------------------------------
//
// GNU_iomanip.h
//
// Purpose:    
//
//    Temporaray implementation of ostream manipulators from the 
//    C++ Standard Library, which are not contained in <iomanip> 
//    as provided with GNU C++.
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

#ifndef INC_GNU_IOMANIP_H
#define INC_GNU_IOMANIP_H


#if defined __GNUC__ && defined __GNUC_MINOR__
# define _GNUC_BEFORE(maj,min) ((__GNUC__<<16)+__GNUC_MINOR__<((maj)<<16)+(min))
#else
# define _GNUC_BEFORE(maj,min) 0
#endif


#if _GNUC_BEFORE(3,0)  

   // ... code required up to gcc 2.99 

#  include <iomanip>
#  include <iostream>

   namespace {
     ostream& left (ostream& os) {
       os.setf(ios::left ,ios::adjustfield); return os;
     };
     ostream& right(ostream& os) {
       os.setf(ios::right,ios::adjustfield); return os;
     };
     ostream& fixed(ostream& os) {
       os.setf(ios::fixed,ios::floatfield); return os;
     };
     ostream& scientific(ostream& os) {
       os.setf(ios::scientific,ios::floatfield);  return os;
     };
     ostream& showpos  (ostream& os) {
       os.setf(ios::showpos); return os;
     };
     ostream& noshowpos(ostream& os) {
       os.unsetf(ios::showpos); return os;
     };
   }

#endif  // version check


#endif  // include blocker
