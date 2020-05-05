cc=gcc
cflags_dbg=-Wall -g
cflags=-Wall -O3 -flto -DNDEBUG

deconwolf: tiffio fft
	$(cc) src/deconwolf.c $(cflags) fft.o tiffio.o -lm -lfftw3f -ltiff -lfftw3f_threads  -o bin/deconwolf

debug: tiffio_dbg fft_dbg
	$(cc) src/deconwolf.c $(cflags_dbg) fft_dbg.o tiffio_dbg.o -lm -lfftw3f -ltiff -lfftw3f_threads  -o bin/deconwolf

fft:
	$(cc) -c src/fft.c $(cflags) -lm -lfftw3f -lfftw3f_threads -o fft.o
fft_dbg:
	$(cc) -c src/fft.c $(cflags) -lm -lfftw3f -lfftw3f_threads -o fft_dbg.o

tiffio:
	$(cc) -c src/tiffio.c $(cflags) -L/usr/lib/x86_64-linux-gnu/ -ltiff -lm -o tiffio.o:

tiffio_dbg:
	$(cc) -c src/tiffio.c $(cflags) -L/usr/lib/x86_64-linux-gnu/ -ltiff -lm -o tiffio_dbg.o
