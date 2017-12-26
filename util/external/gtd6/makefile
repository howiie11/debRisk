include ../../../compiler.in 
CFLAGS = -g -w -DINLINE
MFLAGS = -lm
CC = gfortran -fpermissive $(CFLAGS)

LIBDIR=../../lib/$(ARCH)/
INCDIR=../../include

LIB=libgd6.a
LIBOBJS=gtd6.o

$(LIB): $(LIBOBJS)
	ar rc $(LIB) $(LIBOBJS)

%.o:%.f
	$(CC) -c $(CFLAGS) $^

clean:
	@-rm -rf *~ &> /dev/null
	@-rm -rf $(LIBOBJS) $(LIB) &> /dev/null

test:

all:$(LIB) test

install:
	@-cp *.h $(INCDIR)
	@-cp *.a $(LIBDIR)


