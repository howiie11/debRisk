cleancrap:
	@-find . -name "*~" -exec rm {} \;
	@-find . -name "#*#" -exec rm {} \;

%:
	for dir in $(shell find . -type d);do if [ $$dir != "." ];then make -C $$dir $@;fi;done
