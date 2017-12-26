include ../../../compiler.in 

CFLAGS = -g -w -DINLINE
MFLAGS = -lm
CC = g++ -fpermissive $(CFLAGS)

LIBDIR=../../lib/$(ARCH)/
INCDIR=../../include

LIB=libnrlmsise.a
LIBOBJS=nrlmsise-00.o nrlmsise-00_test.o nrlmsise-00_data.o

$(LIB): $(LIBOBJS)
	ar rc $(LIB) $(LIBOBJS)

%.o:%.c nrlmsise-00.h
	$(CC) -c -o $@ $<

%.opp:%.cpp $(LIB) nrlmsise-00.h
	$(CC) -c -o $@ $<

%.out:%.opp 
	$(CC) $< $(LIB) -o $@

clean:
	@-rm -rf *~ &> /dev/null
	@-rm -rf $(LIBOBJS) $(LIB) *.out &> /dev/null

test:nrlmsise-00_test.out

all:$(LIB) test

install:
	@-cp *.h $(INCDIR)
	@-cp *.a $(LIBDIR)

