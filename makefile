# Normal build:
# make -B
# To build with debug flags and no OpenMP
# make DEBUG=1 OMP=0 -B
# For normal build
# make -B

# To enable the methods shbcl1, shbcl2 use
# make -B OPENCL=1
#
# To build with clang, specify CC=clang

CC?=gcc
# CC=clang # also requires the package libomp-14-dev

UNAME_S := $(shell uname -s)
$(info Host type: $(UNAME_S))

dw = bin/dw
dwbw = bin/dw_bw

CFLAGS = -Wall -Wextra -std=gnu99

CC_VERSION = "$(shell $(CC) --version | head -n 1)"
GIT_VERSION = "$(shell git log --pretty=format:'%aD:%H' -n 1)"

CFLAGS += -DCC_VERSION=\"$(CC_VERSION)\"
CFLAGS += -DGIT_VERSION=\"$(GIT_VERSION)\"

# Since we don't check for these errors, see man math_error
# we can disable this feature
CFLAGS += -fno-math-errno

DESTDIR?=/usr/local/bin
DEBUG?=0
ifeq ($(DEBUG),1)
    # TODO: only compiles with flto on, why?
    # TODO: does not compile without -O1 or above
    CFLAGS += -O1 -g3 -Wno-unknown-pragmas -flto
else
    CFLAGS += -O3 -Wno-unknown-pragmas -flto
    # -O2 -ftree-vectorize and -O3 give about the same performance
    # -DNDEBUG will not make it faster so don't use it
    # -fno-math-errno no relevant performance gain
endif

dw_LIBRARIES =  -lm -ltiff
dwtm_LIBRARIES =  -lm -ltiff
dwbw_LIBRARIES = -lm -ltiff -lpthread

### GSL
CFLAGS += `gsl-config --cflags`
MOSTLYSTATIC?=0
ifeq ($(MOSTLYSTATIC), 0)
    dwbw_LIBRARIES += `gsl-config --libs`
    dw_LIBRARIES += `gsl-config --libs`
else
   dwbw_LIBRARIES += -l:libgsl.a -l:libgslcblas.a
   dw_LIBRARIES += -l:libgsl.a -l:libgslcblas.a
endif

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
ifeq ($(MOSTLYSTATIC), 0)
dw_LIBRARIES += `pkg-config fftw3 fftw3f --libs`
dwbw_LIBRARIES += `pkg-config fftw3 fftw3f --libs`
else
dw_LIBRARIES += -l:libfftw3.a -l:libfftw3f.a
dwbw_LIBRARIES += -l:libfftw3.a -l:libfftw3f.a
endif
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
CFLAGS+=-Xpreprocessor -fopenmp
dw_LIBRARIES += -lomp
else
CFLAGS += -fopenmp
ifeq ($(CC),clang)
CFLAGS+=-fopenmp=libiomp5
endif
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
ifeq ($(UNAME_S), Darwin)
dw_LIBRARIES+=-framework OpenCL
else
CFLAGS+=-I/opt/nvidia/hpc_sdk/Linux_x86_64/21.3/cuda/11.2/targets/x86_64-linux/include/
LDFLAGS=-L/opt/nvidia/hpc_sdk/Linux_x86_64/21.3/cuda/11.2/targets/x86_64-linux/lib/
dw_LIBRARIES+=-lOpenCL
endif
dw_OBJECTS+=method_shb_cl.o method_shb_cl2.o cl_util.o
endif

# clFFT
ifeq ($(OPENCL), 1)
dw_LIBRARIES+=-lclFFT
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

# Extra protection
CC+=-g3

CCF = $(CC) $(CFLAGS)
SRCDIR = src/

dw_OBJECTS += fim.o tiling.o fft.o fim_tiff.o dw.o deconwolf.o dw_maxproj.o dw_util.o method_eve.o method_identity.o method_rl.o method_ave.o method_shb.o dw_imshift.o fft.o dw_nuclei.o dw_dots.o fwhm.o ftab.o dw_psf.o dw_tiff_merge.o dw_psf_sted.o

dwbw_OBJECTS = fim.o fim_tiff.o dw_bwpsf.o bw_gsl.o lanczos.o li.o fft.o dw_util.o

# dwtm = bin/dw_tiffmax
# dwtm_OBJECTS = fim.o fim_tiff.o deconwolf_tif_max.o

all: $(dw) $(dwtm) $(dwbw)
# all: $(dw) $(dwtm) $(dwbw)


$(dw): $(dw_OBJECTS)
	$(CCF) -o $@ $^ $(dw_LIBRARIES)

$(dwbw): $(dwbw_OBJECTS)
	$(CCF) -o $@ $^ $(dwbw_LIBRARIES)

%.o: $(SRCDIR)%.c
	$(CCF) -c $<

clean:
	rm -f *.o

# the kernels are mostly included by cl_util.c some in method_shb_cl*
# TODO: this is a silly list
kernels:
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
	# for iDiv
	cp src/kernels/cl_idiv_kernel.c cl_idiv_kernel
	xxd -i cl_idiv_kernel > src/kernels/cl_idiv_kernel.h
	rm cl_idiv_kernel
	# y = im/y
	cp src/kernels/cl_update_y_kernel.c cl_update_y_kernel
	xxd -i cl_update_y_kernel > src/kernels/cl_update_y_kernel.h
	rm cl_update_y_kernel
	#
	xxd -i src/kernels/cl_real_mul_inplace.c > src/kernels/cl_real_mul_inplace.h
	xxd -i src/kernels/cl_positivity.c > src/kernels/cl_positivity.h
	xxd -i src/kernels/cl_shb_update.c > src/kernels/cl_shb_update.h
	xxd -i src/kernels/cl_preprocess_image.c > src/kernels/cl_preprocess_image.h

install:
	# Binaries
	cp bin/dw_bw $(DESTDIR)/dw_bw
	cp bin/dw $(DESTDIR)/dw
	# Man pages
	cp doc/dw.1 $(MANPATH)/dw.1
	cp doc/dw.1 $(MANPATH)/deconwolf.1
	cp doc/dw_bw.1 $(MANPATH)/dw_bw.1


uninstall:
	rm -f $(DESTDIR)/dw
	rm -f $(DESTDIR)/dw_bw
	rm -f $(DESTDIR)/dw_batch
	rm -f $(MANPATH)/dw.1
	rm -f $(MANPATH)/deconwolf.1
	rm -f $(MANPATH)/dw_bw.1
