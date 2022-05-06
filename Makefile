#Makefile

.PHONY: default test docs clean

default:
	cd src; $(MAKE)

clean:
	cd src/; $(MAKE) clean
	cd docs/; $(MAKE) clean
	cd tests/; $(MAKE) clean

docs:
	cd docs; $(MAKE)

test:
	cd tests; $(MAKE)
