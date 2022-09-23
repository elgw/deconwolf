# Normal build:
# make -B
# To build with debug flags and no OpenMP
# make DEBUG=1 OMP=0 -B
# For normal build
# make -B

# To enable the method shbcl use
# make OPENCL=1

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
    CFLAGS += -g3 -O1 -DNDEBUG -fno-inline
else
    #CFLAGS +=  -g -O2 -ftree-vectorize -Wno-unknown-pragmas -flto
    CFLAGS += -O3 -Wno-unknown-pragmas -flto -DNDEBUG
    #-fno-math-errno no relevant performance gain
endif

dw_LIBRARIES =  -lm -ltiff
dwtm_LIBRARIES =  -lm -ltiff
dwbw_LIBRARIES = -lm -ltiff -lpthread

### GSL
CFLAGS += `gsl-config --cflags`
dwbw_LIBRARIES += `gsl-config --libs`
dw_LIBRARIES += `gsl-config --libs`

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
$(info OMP enabled)
ifeq ($(UNAME_S), Darwin)
dw_LIBRARIES += -Xpreprocessor -lomp
else
CFLAGS += -fopenmp
endif
else
CFLAGS += -fno-openmp
$(info OMP disabled)
endif

## OpenCL
OPENCL?=0
ifeq ($(OPENCL), 1)
$(info OpenCL enabled)
CFLAGS+=-DOPENCL
CFLAGS+=-I/opt/nvidia/hpc_sdk/Linux_x86_64/21.3/cuda/11.2/targets/x86_64-linux/include/
LDFLAGS=-L/opt/nvidia/hpc_sdk/Linux_x86_64/21.3/cuda/11.2/targets/x86_64-linux/lib/
dw_LIBRARIES+=-lclFFT -lOpenCL
dw_OBJECTS+=method_shb_cl.o cl_util.o
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

dw_OBJECTS += fim.o tiling.o fft.o fim_tiff.o dw.o deconwolf.o dw_maxproj.o dw_util.o method_eve.o method_identity.o method_rl.o method_ave.o method_shb.o dw_imshift.o fft.o dw_nuclei.o dw_dots.o fwhm.o ftab.o dw_psf.o
dwbw_OBJECTS = fim.o fim_tiff.o dw_bwpsf.o bw_gsl.o lanczos.o li.o fft.o dw_util.o

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

kernels:
	# cl_complex_square
	cp src/kernels/cl_complex_square.c cl_complex_square
	xxd -i cl_complex_square > src/kernels/cl_complex_square.h
	rm cl_complex_square
	# cl_complex_mul
	cp src/kernels/cl_complex_mul.c cl_complex_mul
	xxd -i cl_complex_mul > src/kernels/cl_complex_mul.h
	rm cl_complex_mul
	# cl_complex_mul_inplace
	cp src/kernels/cl_complex_mul_inplace.c cl_complex_mul_inplace
	xxd -i cl_complex_mul_inplace > src/kernels/cl_complex_mul_inplace.h
	rm cl_complex_mul_inplace
	# cl_complex_mul_conj
	cp src/kernels/cl_complex_mul_conj.c cl_complex_mul_conj
	xxd -i cl_complex_mul_conj > src/kernels/cl_complex_mul_conj.h
	rm cl_complex_mul_conj
	# cl_complex_mul_conj_inplace
	cp src/kernels/cl_complex_mul_conj_inplace.c cl_complex_mul_conj_inplace
	xxd -i cl_complex_mul_conj_inplace > src/kernels/cl_complex_mul_conj_inplace.h
	rm cl_complex_mul_conj_inplace


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
