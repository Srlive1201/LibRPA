DOXYFILE_TEMPLATE = Doxyfile.in
DOXYFILE = Doxyfile
DOXY_OUTPUT_DIRECTORY = librpa
HTML_OUTPUT = html
DOT_PATH = $(shell which dot)

DOXYGEN_AWESOME_VERSION = 2.3.3
DOXYGEN_AWESOME_DIRECTORY = doxygen-awesome-css

.PHONY: default clean

# CAIRO_DEBUG_PDF=1 is a workaround to fix the
# "Failed to extract bounding box from generated diagram file" 
# error before doxygen 1.9.7.
# For details please see https://github.com/doxygen/doxygen/issues/9319
default: $(DOXYFILE) $(DOXYGEN_AWESOME_DIRECTORY)/doxygen-awesome.css
	CAIRO_DEBUG_PDF="1" doxygen $<

$(DOXYFILE): $(DOXYFILE_TEMPLATE) Makefile
	sed -e "s|\(OUTPUT_DIRECTORY *=\) *.*|\1 $(DOXY_OUTPUT_DIRECTORY)|g" \
		-e "s|\(INPUT *=\) *.*|\1 README.md ../../src ../../driver|g" \
		-e "s|\(GENERATE_TAGFILE *=\) *.*|\1 $(DOXY_OUTPUT_DIRECTORY)/$(HTML_OUTPUT)/tagfile.xml|g" \
		-e "s|\(DOT_PATH *=\) *.*|\1 $(DOT_PATH)|g" $< > $@


# https://github.com/boschglobal/doxysphinx/blob/main/external/get_doxygen_awesome.sh
$(DOXYGEN_AWESOME_DIRECTORY)/doxygen-awesome.css:
	mkdir -p $(DOXYGEN_AWESOME_DIRECTORY)
	curl -sSL https://github.com/jothepro/doxygen-awesome-css/archive/refs/tags/v$(DOXYGEN_AWESOME_VERSION).tar.gz \
		| tar -xzvf - --strip-components=1 -C $(DOXYGEN_AWESOME_DIRECTORY) \
		doxygen-awesome-css-$(DOXYGEN_AWESOME_VERSION)/doxygen-awesome.css doxygen-awesome-css-$(DOXYGEN_AWESOME_VERSION)/LICENSE

clean:
	rm -rf $(DOXYFILE)

veryclean: clean
	rm -rf $(DOXY_OUTPUT_DIRECTORY)
