#Makefile

.PHONY: default test docs clean

default:
	cd src; $(MAKE)

clean:
	rm -f *.a *.exe
	cd src/; $(MAKE) clean
	cd docs/; $(MAKE) clean

docs:
	cd docs; $(MAKE)

test:
	cd test; $(MAKE)
