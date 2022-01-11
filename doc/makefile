# Create man pages and html pages from the md files
# requires groff and pandoc

all: dw.1 dw_bw.1
	man ./deconwolf.1 | col -b > dw.1.txt
	man ./dw_bw.1 | col -b > dw_bw.1.txt
	# groff -m man -T xhtml dw.1 > dw.html
	# groff -m man -T xhtml dw_bw.1 > dw_bw.html

dw.1: dw.1.md
	pandoc dw.1.md -s -t man -o dw.1

dw_bw.1: dw_bw.1.md
	pandoc dw_bw.1.md -s -t man -o dw_bw.1
