CFLAGS = -Wall -g -DINLINE
MFLAGS = -lm
#CC = gcc $(CFLAGS)
CC = g++ -fpermissive $(CFLAGS)

nrlmsise-test:nrlmsise-00.o nrlmsise-00_test.o nrlmsise-00_data.o gtd6.o
	$(CC) -o nrlmsise-test nrlmsise-00.o nrlmsise-00_test.o nrlmsise-00_data.o $(MFLAGS)

nrlmsise-00.o:nrlmsise-00.c nrlmsise-00.h
	$(CC) -c nrlmsise-00.c

nrlmsise-00_test.o:nrlmsise-00_test.c nrlmsise-00.h
	$(CC) -c  nrlmsise-00_test.c

nrlmsise-00_data.o:nrlmsise-00_data.c nrlmsise-00.h
	$(CC) -c  nrlmsise-00_data.c

gtd6.o:gtd6.f
	gfortran -c $^

clean:
	rm -rf  nrlmsise-test nrlmsise-00.o nrlmsise-00_test.o nrlmsise-00_data.o gtd6.o
