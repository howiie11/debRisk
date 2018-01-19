include ../../../compiler.in 

CC=g++
CXXFLAGS += -Wall -O3 -I. -w

LIBDIR=../../lib/$(ARCH)/
INCDIR=../../include

LIB=libSAT.a
LIBOBJS=SAT_Force.o SAT_RefSys.o SAT_VecMat.o SAT_DE.o SAT_Kepler.o SAT_Time.o
HFILES=$(shell ls *.h)

$(LIB): $(LIBOBJS)
	ar rc $(LIB) $(LIBOBJS)

clean:
	@-rm -rf *~ &> /dev/null
	@-rm -rf $(LIBOBJS) $(LIB) *.out &> /dev/null

%.opp:%.cpp $(LIB) $(HFILES)
	$(CC) $(CXXFLAGS) -c -o $@ $<

%.out:%.opp 
	$(CC) $< $(LIB) -o $@

test:SAT_test.out

all:$(LIB) test

install:
	cp *.h $(INCDIR)
	cp *.a $(LIBDIR)
