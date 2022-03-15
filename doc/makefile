# Create man pages and html pages from the md files
# requires groff and pandoc

VERSION != awk -F' ' '{if($$2=="deconwolf_version") print $$3}' ../src/dw_version.h | sed 's/"//g'

all: dw.txt dw_bw.txt

%.txt: %.1
	man ./$< | col -b > $@

%.1: %.md
	pandoc $< -s -t man -o $@
	echo "VERSION = $(VERSION)"
	# on MacOS, use:
	# sed -i '' "s/ABCDE/$(VERSION)/g" $@
	sed -i "s/ABCDE/$(VERSION)/g" $@
