include compiler.in

show:
	@echo "It works!"

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
