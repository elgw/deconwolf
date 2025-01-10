CC=gcc -std=gnu99

CFLAGS=-Wall -Wextra

ASAN?=0
ifeq ($(ASAN),1)
CFLAGS+=-g3 -fanalyzer -fno-inline -fstack-protector-all -fsanitize=address -fno-omit-frame-pointer
endif


DEBUG?=0
ifeq ($(DEBUG),1)
CFLAGS+=-g3 -fno-inline -fstack-protector-all -fno-omit-frame-pointer
endif

ifeq ($(DEBUG),0)
ifeq ($(ASAN),0)
CFLAGS+=-O3
LDFLAGS+=-flto
endif
endif

FORCE: ;

LDFLAGS+=-lm

# lib Cairo
CFLAGS+=`pkg-config cairo libpng --cflags`
LDFLAGS+=`pkg-config cairo libpng --libs`

# lib GSL
CFLAGS+=`pkg-config --cflags gsl`
LDFLAGS+=`pkg-config --libs gsl`

# Tiff
CFLAGS+=`pkg-config libtiff-4 --cflags`
LDFLAGS+=`pkg-config libtiff-4 --libs`

# FFTW3
CFLAGS += `pkg-config fftw3 fftw3f --cflags`
LDFLAGS += `pkg-config fftw3 fftw3f --libs` -lfftw3_omp -lfftw3f_omp

# OpenMP
CFLAGS+=-fopenmp

all: tiling_ut fft_ut fim_ut

sparse_preprocess_files=sparse_preprocess_cli.c fim_tiff.o fim.o dw_util.o ftab.o fft.o sparse_preprocess.c
sparse_preprocess_cli: $(sparse_preprocess_files)
	$(CC) $(CFLAGS) -DSTANDALONE $(sparse_preprocess_files) $(LDFLAGS) -o sparse_preprocess_cli

tiling_files=fft.o tiling_ut.o tiling.o fim_tiff.o fim.o ftab.o dw_util.o
tiling_ut: $(tiling_files)
	$(CC) $(CFLAGS) $(tiling_files) $(LDFLAGS) -o tiling_ut

fft_ut_files=fft_ut.o fft.o fim.o ftab.o fim_tiff.o dw_util.o
fft_ut: $(fft_ut_files)
	$(CC) $(fft_ut_files) $(CFLAGS) $(LDFLAGS) -o fft_ut

dw_render_files=dw_render.c fim.o ftab.o fim_tiff.o dw_util.o fft.o
dw_render: $(dw_render_files)
	$(CC) $(CFLAGS) -DSTANDALONE $(dw_render_files) $(LDFLAGS) -o dw_render

## Note: prf_* and qsort comes from pixel_random_forest repo

dw_nuclei_files=dw_nuclei.c fim.o ftab.o fim_tiff.o dw_util.o \
trafo/libtrafo.a fft.o dw_png.o

dw_nuclei: $(dw_nuclei_files)
	$(CC) $(CFLAGS) -DSTANDALONE $(dw_nuclei_files) $(LDFLAGS) -o dw_nuclei

dw_png_ut_files = dw_png.c dw_png_ut.c
dw_png_ut: $(dw_png_ut_files)
	$(CC) $(CFLAGS) $(dw_png_ut_files) $(LDFLAGS) -o dw_png_ut

dw_dots_files=dw_dots.c fim.o ftab.o fim_tiff.o  dw_util.o fft.o gmlfit.o
dw_dots: $(dw_dots_files)
	$(CC) $(CFLAGS) -DSTANDALONE $(dw_dots_files) $(LDFLAGS) -o dw_dots

psf_files = dw_psf.c fim.o ftab.o fim_tiff.o  dw_util.o fft.o
dw_psf: $(psf_files)
	$(CC) $(CFLAGS) -DSTANDALONE $(psf_files) $(LDFLAGS) -o dw_psf

spsf_files = dw_psf_sted.c fim.o ftab.o fim_tiff.o  dw_util.o fft.o
dw_psf_sted: $(spsf_files)
	$(CC) $(CFLAGS) -DSTANDALONE $(spsf_files) $(LDFLAGS) -o dw_psf_sted

# Test the float image module
fim_ut_files=fim_ut.o fim.o dw_util.o ftab.o fft.o fim_tiff.o
fim_ut: $(fim_ut_files)
	$(CC)  $(CFLAGS) $(fim_ut_files) $(LDFLAGS) -o fim_ut

ftab_ut: ftab.c
	$(CC) ftab_ut.c ftab.c $(CFLAGS) -o ftab_ut

fim_tiff_ut_files=fim_tiff.c fim.o fft.o dw_util.o ftab.o
fim_tiff_ut: $(fim_tiff_ut_files)
	$(CC) $(CFLAGS) -Dunittest $(fim_tiff_ut_files) $(LDFLAGS) -o fim_tiff_ut

dw_tiff_max_files=fim.o fim_tiff.o dw_maxproj.c ftab.o fft.o dw_util.o
dw_tiff_max: $(dw_tiff_max_files)
	$(CC) -DSTANDALONE $(CFLAGS) $(dw_tiff_max_files) $(LDFLAGS) -o dw_tiff_max

tiff_from_raw_files=tiff_from_raw.c fim.o fim_tiff.o ftab.o fft.o dw_util.o
tiff_from_raw:
	$(CC) $(CFLAGS) $(tiff_from_raw_files) $(LDFLAGS) -o tiff_from_raw

li:
	$(CC) -DLI_TEST li.c  -o li_ut -lm `pkg-config --CFLAGS --libs gsl`

%.o: %.c
	$(CC) $(CFLAGS) -c $<

trafo/libtrafo.a: FORCE
	make -C trafo libtrafo.a
