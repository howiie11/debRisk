include ../../compiler.in

PACKAGES=$(shell find . -type d)

show:
	echo $(PACKAGES)

%:
	for dir in $(PACKAGES);do if [ $$dir != "." ];then make -C $$dir $@;fi;done
