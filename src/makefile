cc=gcc
CFLAGS=-lm -Wno-unknown-pragmas -Wall -g
CFLAGS2=-lm -Wno-unknown-pragmas -Wall -O3 -flto -march=native -ftree-vectorize

cflags+=-fopenmp -Wall -Wextra -g
cflags+=`pkg-config cairo libpng --cflags`
ldflags+=`pkg-config cairo libpng --libs`
ldflags+=-lm -ltiff -lfftw3f -lfftw3f_threads -lfftw3 -lfftw3_threads

gsl_ldflags=`pkg-config --libs gsl`
gsl_cflags=`pkg-config --cflags  gsl`

all: tiling_ut fft_ut fim_ut

tiling_ut:
	$(cc) tiling_ut.c tiling.c fim_tiff.c fim.c $(CFLAGS) -ltiff -lfftw3f -o tiling_ut
fft_ut:
	$(cc) fft_ut.c fft.c fim.c $(CFLAGS) -lm -lfftw3f -lfftw3f_threads -o fft_ut

dw_render_files=dw_render.c dw_render.h fim.c ftab.c fim_tiff.c dw_util.c fft.c
dw_render: $(dw_render_files)
	$(cc) $(cflags) -DSTANDALONE $(dw_render_files) $(ldflags) -o dw_render

dw_nuclei: dw_nuclei.c fim.c ftab.c fim_tiff.c prf_forest.o
	$(cc) $(cflags) -DSTANDALONE dw_nuclei.c fim.c dw_util.c fim_tiff.c fft.c ftab.c prf_forest.o prf_tree.o prf_util.o qsort.o $(ldflags) -o dw_nuclei

dw_dots_files=dw_dots.c fim.c ftab.c fim_tiff.c fwhm.c dw_util.c fft.c
dots_ldflags=$(ldflags) $(gsl_ldflags)
dots_cflags=$(cflags) $(gsl_cflags)
dw_dots: $(dw_dots_files)
	$(cc) $(dots_cflags) -DSTANDALONE $(dw_dots_files) $(dots_ldflags) -o dw_dots

psf_files = dw_psf.c fim.c ftab.c fim_tiff.c fwhm.c dw_util.c fft.c
psf_ldflags=$(ldflags) $(gsl_ldflags)
psf_cflags=$(cflags) $(gsl_cflags)
dw_psf: $(psf_files)
	$(cc) $(psf_cflags) -DSTANDALONE $(psf_files) $(psf_ldflags) -o dw_psf

prf_forest.o: random_forest/prf_forest.c
	$(cc) -c $(cflags) random_forest/prf*.c random_forest/qsort.c

# Test the float image module
fim_ut: fim_ut.c fim.c fim.h
	$(cc) -fopenmp fim_ut.c fft.c fim.c fim_tiff.c ftab.c $(CFLAGS) -lm -o fim_ut -lfftw3f -lfftw3f_threads -ltiff

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
