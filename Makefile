.PHONY: all compile doc clean clean-doc clean-compile

all: compile doc

compile:
	./Allwmake

doc:
	doxygen doc/Doxygen/Doxyfile 2>&1 | tee doc/Doxygen/doxy.log

clean: clean-compile clean-doc

clean-compile:
	./Allwclean

clean-doc:
	rm -fr doc/output
