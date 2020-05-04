cc=gcc
cflags=-Wall -g
#cflags=-Wall -O3 -flto

deconwolf: tiffio fft
	$(cc) src/deconwolf.c $(cflags) fft.o tiffio.o -lm -lfftw3f -ltiff -lfftw3f_threads  -o bin/deconwolf

fft:
	$(cc) -c src/fft.c $(cflags) -lm -lfftw3f -lfftw3f_threads -o fft.o

tiffio:
	$(cc) -c src/tiffio.c $(cflags) -L/usr/lib/x86_64-linux-gnu/ -ltiff -lm -o tiffio.o
