# Normal build:
# make -B
# To build with debug flags:
# make debug=1 -B

DEBUG ?= 0
ifeq (DEBUG, 1)
    CFLAGS =-g3 -gdwarf2 -DDEBUG
else
    CFLAGS=-DNDEBUG -O3 -flto -march=native
endif

CC = gcc $(CFLAGS)

EXECUTABLE = bin/deconwolf
OBJECTS = fim.o tiling.o fft.o fim_tiff.o dw.o deconwolf.o
LIBRARIES = -lm -lfftw3f -lfftw3f_threads -ltiff
SRCDIR = src/

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) -o $@ $^ $(LIBRARIES)

%.o: $(SRCDIR)%.c
	$(CC) -c $<

clean:
	rm -f $(EXECUTABLE) $(OBJECTS)
