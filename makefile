cc=gcc
cflags_dbg=-Wall -g
cflags=-Wall -O3 -flto -DNDEBUG

deconwolf: tiffio fft tiling
	$(cc) src/deconwolf.c $(cflags) fft.o tiffio.o tiling.o -lm -lfftw3f -ltiff -lfftw3f_threads  -o bin/deconwolf

debug: tiffio_dbg fft_dbg tiling_dbg
	$(cc) src/deconwolf.c $(cflags_dbg) fft_dbg.o tiffio_dbg.o tiling_dbg.o -lm -lfftw3f -ltiff -lfftw3f_threads  -o bin/deconwolf

tiling:
	$(cc) -c src/tiling.c $(cflags) -lm -o tiling.o

tiling_dbg:
	$(cc) -c src/tiling.c $(clags_dbg) -lm -o tiling_dbg.o

fft:
	$(cc) -c src/fft.c $(cflags) -lm -lfftw3f -lfftw3f_threads -o fft.o
fft_dbg:
	$(cc) -c src/fft.c $(cflags_dbg) -lm -lfftw3f -lfftw3f_threads -o fft_dbg.o

tiffio:
	$(cc) -c src/tiffio.c $(cflags) -L/usr/lib/x86_64-linux-gnu/ -ltiff -lm -o tiffio.o

tiffio_dbg:
	$(cc) -c src/tiffio.c $(cflags_dbg) -L/usr/lib/x86_64-linux-gnu/ -ltiff -lm -o tiffio_dbg.o
