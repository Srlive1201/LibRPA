DOXYFILE_TEMPLATE = Doxyfile.in
DOXYFILE = Doxyfile
DOXY_OUTPUT_DIRECTORY = _build
DOT_PATH = \/usr\/bin

.PHONY: default extract clean

# CAIRO_DEBUG_PDF=1 is a workaround to fix the
# "Failed to extract bounding box from generated diagram file" 
# error before doxygen 1.9.7.
# For details please see https://github.com/doxygen/doxygen/issues/9319
default: $(DOXYFILE)
	CAIRO_DEBUG_PDF="1" doxygen $<

$(DOXYFILE): $(DOXYFILE_TEMPLATE) Makefile
	sed -e "s/\(OUTPUT_DIRECTORY *=\) *.*/\1 $(DOXY_OUTPUT_DIRECTORY)/g" \
		-e "s|\(INPUT *=\) *.*|\1 ../../src ../../driver|g" \
		-e "s/\(DOT_PATH *=\) *.*/\1 $(DOT_PATH)/g" $< > $@

clean:
	rm -rf $(DOXYFILE) $(DOXY_OUTPUT_DIRECTORY)
