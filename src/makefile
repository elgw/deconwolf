cc=gcc
CFLAGS=-lm -Wno-unknown-pragmas -Wall -g

all: tiling_ut fft_ut fim_ut

tiling_ut:
	$(cc) tiling_ut.c tiling.c fim_tiff.c fim.c $(CFLAGS) -ltiff -lfftw3f -o tiling_ut 
fft_ut:
	$(cc) fft_ut.c fft.c fim.c $(CFLAGS) -lm -lfftw3f -lfftw3f_threads -o fft_ut

fim_ut:
	$(cc) fim_ut.c fim.c $(CFLAGS) -lm -o fim_ut -lfftw3f


