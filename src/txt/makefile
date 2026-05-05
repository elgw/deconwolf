sources = $(wildcard *.txt)
targets = $(patsubst %.txt,%.txt.h,$(sources))

%.txt.h: %.txt
	$(info $< -> $@)
	xxd -i "$<" > "$@"

build: $(targets)
	$(info done)
