cc=gcc

CC_VERSION = "$(shell gcc --version | head -n 1)"
GIT_VERSION = "$(shell git log --pretty=format:'%aD:%H' -n 1)"

f1 = -DCC_VERSION=\"$(CC_VERSION)\"
f2 = $(f1) -DGIT_VERSION=\"$(GIT_VERSION)\"

cflags_dbg=-Wall -g $(f2)
cflags=-Wall -O3 -flto -DNDEBUG -march=native $(f2)

deconwolf: fft tiling fim fim_tiff
	$(cc) src/deconwolf.c $(cflags) fft.o tiling.o fim.o fim_tiff.o -lm -lfftw3f -ltiff -lfftw3f_threads  -o bin/deconwolf

debug: fim_tiff_dbg fft_dbg tiling_dbg fim_dbg
	$(cc) src/deconwolf.c $(cflags_dbg) fft_dbg.o fim_tiff_dbg.o tiling_dbg.o fim_dbg.o -lm -lfftw3f -ltiff -lfftw3f_threads  -o bin/deconwolf

fim:
	$(cc) -c src/fim.c $(cflags) -o fim.o

fim_dbg:
	$(cc) -c src/fim.c $(clags_dbg) -o fim_dbg.o

tiling:
	$(cc) -c src/tiling.c $(cflags) -o tiling.o

tiling_dbg:
	$(cc) -c src/tiling.c $(clags_dbg) -o tiling_dbg.o

fft:
	$(cc) -c src/fft.c $(cflags) -o fft.o
fft_dbg:
	$(cc) -c src/fft.c $(cflags_dbg) -o fft_dbg.o

fim_tiff:
	$(cc) -c src/fim_tiff.c $(cflags) -L/usr/lib/x86_64-linux-gnu/ -o fim_tiff.o

fim_tiff_dbg:
	$(cc) -c src/fim_tiff.c $(cflags_dbg) -L/usr/lib/x86_64-linux-gnu/ -o fim_tif_dbg.o
