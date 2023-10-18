# Create man pages and html pages from the md files
# requires groff and pandoc

VERSION_MAJOR != awk -F' ' '{if($$2=="DW_VERSION_MAJOR") print $$3}' ../src/dw_version.h | sed 's/"//g'
VERSION_MINOR != awk -F' ' '{if($$2=="DW_VERSION_MINOR") print $$3}' ../src/dw_version.h | sed 's/"//g'
VERSION_PATCH != awk -F' ' '{if($$2=="DW_VERSION_PATCH") print $$3}' ../src/dw_version.h | sed 's/"//g'
VERSION="$(VERSION_MAJOR).$(VERSION_MINOR).$(VERSION_PATCH)"


all: dw.txt dw.1 dw_bw.1 dw_bw.txt dw_dots.txt dw_dots.1

%.txt: %.1
	man ./$< | col -b > $@

%.1: %.md
	pandoc $< -s -t man -o $@
	echo "VERSION = $(VERSION)"
	# on MacOS, use:
	# sed -i '' "s/ABCDE/$(VERSION)/g" $@
	sed -i "s/ABCDE/$(VERSION)/g" $@
