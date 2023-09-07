# Makefile for Sphinx documentation
# Variables for this project
DOXY_SOURCE = ../doxygen
DOXY_BUILDDIR = ../doxygen/_build
BREATHE_PROJECT_NAME = LibRPA

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?= -Dbreathe_projects.$(BREATHE_PROJECT_NAME)=$(DOXY_BUILDDIR)/xml
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

# Targets
.PHONY: html help clean

html: $(DOXY_BUILDDIR)/xml
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS)

$(DOXY_BUILDDIR)/xml:
	cd $(DOXY_SOURCE) && $(MAKE)

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS)

clean:
	rm -rf $(BUILDDIR)
