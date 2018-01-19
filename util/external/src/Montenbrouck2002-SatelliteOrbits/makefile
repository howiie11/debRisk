#-------------------------------------------------------------------------------
#
# makefile
#
# Purpose:
#
#   Linux make file for C++ library and programs of
#   "O. Montenbruck, E. Gill; Satellite Orbits - Models, Methods, and 
#   Applications; Springer Verlag (2000)."
#
# Notes:
# 
#   To build the executable programs, copy all source and header files 
#   (*.cpp,*.h) as well as this makefile into the working directory:
#      >tar -xvf /cdrom/SAT/Linux/satsrc.tar
#      >cp /cdrom/SAT/Linux/makefile .
#   Then start the make tool:
#      >make
#
# Last modified:
#
#   2000/03/04  OMO  Final version (1st edition)
#   2005/04/14  OMO  Final version (2nd reprint)
#
#-------------------------------------------------------------------------------

# Compilation options
CXXFLAGS += -Wall -O3 -I.

# The library is:
LIB=libSAT.a

# and is made from the following object files
LIBOBJS= \
  SAT_Force.o SAT_Kepler.o \
  SAT_RefSys.o SAT_Time.o SAT_VecMat.o

# The application files themselves are:
APPLICS= \
  Exercise_2_1 Exercise_2_2 Exercise_2_3 Exercise_2_4 \
  Exercise_2_5 Exercise_2_6 \
  Exercise_3_1 Exercise_3_2 Exercise_3_3 Exercise_3_4 \
  Exercise_4_1 Exercise_4_2 Exercise_4_3 \
  Exercise_5_1 Exercise_5_2 Exercise_5_3 \
  Exercise_6_1 Exercise_6_2 Exercise_6_3 Exercise_6_4 \
  Exercise_7_1 \
  Exercise_8_1 Exercise_8_2 Exercise_8_3 \
  GEODA RTOD TDRSOD

# Our main target is to make all the applications
all: $(APPLICS)

# The library is made from the object files by using 'ar'
libSAT.a: $(LIBOBJS)
	ar rc $(LIB) $(LIBOBJS)

# Each of the application files, in addition to being dependant on its
# own source file, needs the library
$(APPLICS): $(LIB)

