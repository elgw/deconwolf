# Normal build:
# make -B
# To build with debug flags and no OpenMP
# make DEBUG=1 OMP=0 -B
# For normal build
# make -B

UNAME_S := $(shell uname -s)
$(info Host type: $(UNAME_S))

dw = bin/dw
dwbw = bin/dw_bw

CFLAGS = -Wall -Wextra -std=gnu99

CC_VERSION = "$(shell gcc --version | head -n 1)"
GIT_VERSION = "$(shell git log --pretty=format:'%aD:%H' -n 1)"

CFLAGS += -DCC_VERSION=\"$(CC_VERSION)\"
CFLAGS += -DGIT_VERSION=\"$(GIT_VERSION)\"

# Since we don't check for these errors, see man math_error
# we can disable this feature
CFLAGS += -fno-math-errno

DESTDIR?=/usr/local/bin
DEBUG?=0
ifeq ($(DEBUG),1)
    CFLAGS += -g3 -DDEBUG
else
   #CFLAGS +=  -g -O2 -ftree-vectorize -Wno-unknown-pragmas -flto
   CFLAGS +=  -O3 -Wno-unknown-pragmas -flto -DNDEBUG
   #-fno-math-errno no relevant performance gain
endif

dw_LIBRARIES =  -lm -ltiff
dwtm_LIBRARIES =  -lm -ltiff
dwbw_LIBRARIES = -lm -ltiff -lpthread

### GSL
CFLAGS += `gsl-config --cflags`
dwbw_LIBRARIES += `gsl-config --libs`

### FFT Backend
FFTW3=1
MKL?=0
CUDA?=0

ifeq ($(MKL), 1)
dw=bin/dw-mkl
FFTW3=0
CUDA=0
endif

ifeq ($(CUDA), 1)
dw=bin/dw-cuda
FFTW3=0
MKL=0
endif

ifeq ($(MKL),1)
$(info FFTW backend: MKL)
CFLAGS += -DMKL `pkg-config mkl-static-lp64-iomp --cflags`
dw_LIBRARIES += `pkg-config mkl-static-lp64-iomp --cflags --libs`
dwbw_LIBRARIES += `pkg-config mkl-static-lp64-iomp --cflags --libs`
endif

ifeq ($(FFTW3), 1)
$(info FFTW backend: FFTW3)
CFLAGS += `pkg-config fftw3 fftw3f --cflags`
dw_LIBRARIES += `pkg-config fftw3 fftw3f --libs`
dwbw_LIBRARIES += `pkg-config fftw3 fftw3f --libs`
ifneq ($(WINDOWS),1)
dw_LIBRARIES += -lfftw3f_threads
dwbw_LIBRARIES += -lfftw3f_threads
endif
endif

ifeq ($(CUDA),1)
$(info FFTW backend: CUDA)
$(info FFT BACKEND: CUDA)
CUDA_DIR=/usr/local/cuda/
CFLAGS += -I$(CUDA_DIR)/include/ -L$(CUDA_DIR)/lib64/ -DCUDA
dw_LIBRARIES +=  -lcufftw
dwbw_LIBRARIES += -lcufftw
endif

## OpenMP
OMP?=1
ifeq ($(OMP), 1)
ifeq ($(UNAME_S), Darwin)
dw_LIBRARIES += -Xpreprocessor -lomp
else
CFLAGS += -fopenmp
endif
else
CFLAGS += -fno-openmp
endif

ifneq ($(UNAME_S),Darwin)
    CFLAGS+=-march=native -mtune=native
endif

ifeq ($(WINDOWS),1)
	CFLAGS += -DWINDOWS
endif

## MANPATH
MANPATH=/usr/share/man/man1/
ifeq ($(UNAME_S),Darwin)
	MANPATH=/usr/local/share/man/man1
endif


CC = cc $(CFLAGS)
SRCDIR = src/

dw_OBJECTS = fim.o tiling.o fft.o fim_tiff.o dw.o deconwolf.o deconwolf_tif_max.o method_eve.o method_identity.o method_rl.o
dwbw_OBJECTS = fim.o fim_tiff.o dw_bwpsf.o bw_gsl.o lanczos.o li.o

# dwtm = bin/dw_tiffmax
# dwtm_OBJECTS = fim.o fim_tiff.o deconwolf_tif_max.o

all: $(dw) $(dwtm) $(dwbw)
# all: $(dw) $(dwtm) $(dwbw)


$(dw): $(dw_OBJECTS)
	$(CC) -o $@ $^ $(dw_LIBRARIES)

$(dwbw): $(dwbw_OBJECTS)
	$(CC) -o $@ $^ $(dwbw_LIBRARIES)

%.o: $(SRCDIR)%.c
	$(CC) -c $<

clean:
	rm -f $(dw) $(dw_OBJECTS)
	rm -f $(dwbw) $(dwbw_OBJECTS)

install:
	# Binaries
	cp bin/dw_bw $(DESTDIR)/dw_bw
	cp bin/dw $(DESTDIR)/dw
	# Man pages
	cp doc/dw.1 $(MANPATH)/dw.1
	cp doc/dw.1 $(MANPATH)/deconwolf.1
	cp doc/dw_bw.1 $(MANPATH)/dw_bw.1


uninstall:
	rm $(DESTDIR)/dw
	rm $(DESTDIR)/dw_bw
	rm $(DESTDIR)/dw_batch
	rm $(MANPATH)/dw.1
	rm $(MANPATH)/deconwolf.1
	rm $(MANPATH)/dw_bw.1
