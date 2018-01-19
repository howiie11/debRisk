include compiler.in

UTIL=util
LIBDIR=$(UTIL)/lib/$(ARCH)

CPP=g++
CFLAGS=-I. -w -I$(UTIL)/include -I$(UTIL)/include/gsl -I$(UTIL)/kernels
LFLAGS=-lm $(LIBDIR)/cspice.a $(LIBDIR)/csupport.a $(LIBDIR)/libgsl.a $(LIBDIR)/libgslcblas.a \
	   $(LIBDIR)/libSAT.a $(LIBDIR)/libnrlmsise.a $(LIBDIR)/libgd6.a $(LIBDIR)/libconfig++.a 

show:
	@echo "It works!"


install:
	make unpack
	make -C util/external all install

cleancrap:
	@echo "Cleaning crap..."
	@find . -name "*~" -exec rm -rf {} \;
	@find . -name "#*#" -exec rm -rf {} \;
	@find . -name "#*" -exec rm -rf {} \;

cleanexe:
	@echo "Cleaning executable and log..."
	@rm -rf *.pyc *.out *.exe *.log *.opp

clean:cleancrap cleanexe
	@echo "Cleaning everything..."
	@rm -rf *.png *.dat

%.out:%.opp 
	$(CPP) $< $(LFLAGS) -o $@
	./$@

%.opp:%.cpp *.hpp
	$(CPP) -c $(CFLAGS) $< -o $@

external:
	$(MAKE) -C $(EXTDIR)

commit:
	@echo "Commiting..."
	@-git commit -am "Commit"
	@-git push origin master

pull:
	@-git reset --hard HEAD
	@-git pull

pack:
	@echo "Packing data..."
	@bash .store/pack.sh pack

unpack:
	@echo "Unpacking data..."
	@bash .store/pack.sh unpack
