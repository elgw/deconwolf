cc=gcc
CFLAGS=-lm -Wno-unknown-pragmas -Wall -g
CFLAGS2=-lm -Wno-unknown-pragmas -Wall -O3 -flto -march=native -ftree-vectorize

cflags+=-fopenmp
cflags+=`pkg-config cairo --cflags`
ldflags+=`pkg-config cairo --libs`
ldflags+=-lm -ltiff -lfftw3f -lfftw3f_threads

all: tiling_ut fft_ut fim_ut

tiling_ut:
	$(cc) tiling_ut.c tiling.c fim_tiff.c fim.c $(CFLAGS) -ltiff -lfftw3f -o tiling_ut
fft_ut:
	$(cc) fft_ut.c fft.c fim.c $(CFLAGS) -lm -lfftw3f -lfftw3f_threads -o fft_ut

dw_render: dw_render.c dw_render.h fim.c ftab.c fim_tiff.c
	$(cc) $(cflags) -DSTANDALONE dw_render.c fim.c fim_tiff.c fft.c ftab.c $(ldflags) -o dw_render

# Test the float image module
fim_ut: fim_ut.c fim.c fim.h
	$(cc) -fopenmp fim_ut.c fft.c fim.c fim_tiff.c $(CFLAGS) -lm -o fim_ut -lfftw3f -lfftw3f_threads -ltiff

ftab_ut: ftab.c
	$(cc) ftab_ut.c ftab.c $(CFLAGS) -o ftab_ut

# Read a tif file and write it back as foo.tif
fim_tiff_ut:
	$(cc)  -Wall fim_tiff.c fim.c fft.c -g -fno-inline -Dunittest -O3 -march=native -ltiff -lm -lfftw3f -o fim_tiff_ut

dw_tiff_max:
	$(cc) fim.c fim_tiff.c deconwolf_tif_max.c $(CFLAGS2) -lm -ltiff -lfftw3f -lfftw3f_threads -o ../bin/dw_tiff_max

tiff_from_raw:
	$(cc) tiff_from_raw.c fim_tiff.c $(CFLAGS2) -lm -ltiff -o ../bin/dw_tiff_from_raw

li:
	$(cc) -DLI_TEST li.c  -o li_ut -lm `pkg-config --cflags --libs gsl`
