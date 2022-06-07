DOXYFILE_TEMPLATE = Doxyfile.in
DOXYFILE = Doxyfile
DOXY_OUTPUT_DIRECTORY = _build
DOT_PATH = \/usr\/bin

.PHONY: default extract clean

default: $(DOXYFILE)
	doxygen $<

$(DOXYFILE): $(DOXYFILE_TEMPLATE) Makefile
	# $(shell sed -e "s/\(OUTPUT_DIRECTORY *=\) *.*/\1 $(DOXY_OUTPUT_DIRECTORY)/g" $< > $@)
	# $(shell sed "s/\(INPUT *=\) *.*/\1 $(PWD)\/..\/src/g" $< > $@)
	# $(shell sed "s/\(DOT_PATH *=\) *.*/\1 $(shell which dot)/g" $< > $@)
	sed -e "s/\(OUTPUT_DIRECTORY *=\) *.*/\1 $(DOXY_OUTPUT_DIRECTORY)/g" \
		-e "s|\(INPUT *=\) *.*|\1 ../src|g" \
		-e "s/\(DOT_PATH *=\) *.*/\1 $(DOT_PATH)/g" $< > $@

clean:
	rm -rf $(DOXYFILE) $(DOXY_OUTPUT_DIRECTORY)
