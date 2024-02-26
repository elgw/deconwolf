## Normal build:
# make -B

## To build with debug flags and no OpenMP
# make DEBUG=1 OMP=0 DEBUG=1

## With GPU acceleration
# make kernels
# make -B VKFFT=1
#

CC=gcc -std=gnu11
# CC=clang # also requires the package libomp-14-dev

UNAME_S := $(shell uname -s)
$(info Host type: $(UNAME_S))

dw = bin/dw
dwbw = bin/dw_bw

CFLAGS = -Wall -Wextra

CC_VERSION = "$(shell $(CC) --version | head -n 1)"
GIT_VERSION = "$(shell git log --pretty=format:'%aD:%H' -n 1)"

CFLAGS += -DCC_VERSION=\"$(CC_VERSION)\"
CFLAGS += -DGIT_VERSION=\"$(GIT_VERSION)\"

DESTDIR?=/usr/local/bin
DEBUG?=0

#
# Optimization flags
#

# Notes:
# -O2 -ftree-vectorize and -O3 give about the same performance
# -DNDEBUG turns off some self-tests
# -fno-math-errno gives no relevant performance gain

ifeq ($(DEBUG),1)
    CFLAGS += -O0 -Wno-unknown-pragmas -fanalyzer -g3
else
    CFLAGS += -O3 -flto -DNDEBUG
endif

#
# Check if the linker can perform link time optimizations
#

define check_ld_flag
$(shell $(LD) $(1) -e 0 /dev/null 2>/dev/null && echo $(1))
endef

$(info --testing -flto and -flto=auto)
LD_HAS_FLTO=$(call check_ld_flag,-flto)
LD_HAS_FLTO_AUTO=$(call check_ld_flag,-flto=auto)

ifneq ($(LD_HAS_FLTO_AUTO),)
$(info Enabling -flto=auto)
CFLAGS+=-flto=auto
else
ifneq ($(LD_HAS_FLTO),)
$(info Enabling -flto)
CFLAGS+=-flto
else
$(info Neither -flto nor -flto=auto available)
endif
endif

#
# Set up for linking
#

dw_LIBRARIES=
dwbw_LIBRARIES=

#
# Math library
#

$(info -- Math (-lm) library?)
LD_HAS_LM=$(call check_ld_flag,-lm)
ifneq ($(LD_HAS_LM),)
$(info found. Will use -lm)
dw_libraries+=-lm
else
$(info NOT FOUND)
endif

#
# Tiff library
#

$(info -- Looking for libtiff)
TIFFLIB=

ifeq ($(TIFFLIB),)
LIBTIFF6 = $(shell pkg-config libtiff-6 --exists ; echo $$?)
ifeq ($(LIBTIFF6),0)
TIFFLIB=libtiff-6
endif
endif


ifeq ($(TIFFLIB),)
LIBTIFF5 = $(shell pkg-config libtiff-5 --exists ; echo $$?)
ifeq ($(LIBTIFF4),0)
TIFFLIB=libtiff-5
endif
endif

ifeq ($(TIFFLIB),)
LIBTIFF4 = $(shell pkg-config libtiff-4 --exists ; echo $$?)
ifeq ($(LIBTIFF4),0)
	TIFFLIB=libtiff-4
endif
endif

ifeq ($(TIFFLIB),)
LIBTIFF = $(shell pkg-config libtiff --exists ; echo $$?)
ifeq ($(LIBTIFF4),0)
TIFFLIB=libtiff
endif
endif


ifneq ($(TIFFLIB),)
$(info found $(TIFFLIB))
CFLAGS+=$(shell pkg-config ${TIFFLIB} --cflags)
dw_LIBRARIES+=$(shell pkg-config ${TIFFLIB} --libs)
dwbw_LIBRARIES+=$(shell pkg-config ${TIFFLIB} --libs)
endif

ifeq ($(TIFFLIB),)
$(error tiff library not found)
endif

##
## GSL
##

$(info -- Looking for GSL)
HAS_GSL= $(shell gsl-config --version)
ifneq ($(HAS_GSL),)
GSL_VERSION = $(shell gsl-config --version )
$(info found GSL ${GSL_VERSION})
else
$(error Could not find GSL)
endif

GSLFLAGS=`gsl-config --cflags`
CFLAGS += $(GSLFLAGS)
MOSTLYSTATIC?=0
ifeq ($(MOSTLYSTATIC), 0)
    dwbw_LIBRARIES += `gsl-config --libs`
    dw_LIBRARIES += `gsl-config --libs`
else
   dwbw_LIBRARIES += -l:libgsl.a -l:libgslcblas.a
   dw_LIBRARIES += -l:libgsl.a -l:libgslcblas.a
endif

##
## FFT Backend
##
FFTW3=1
MKL?=0

ifeq ($(MKL), 1)
dw=bin/dw-mkl
FFTW3=0
endif

ifeq ($(MKL),1)
$(info FFTW backend: MKL)
CFLAGS += -DMKL `pkg-config mkl-static-lp64-iomp --cflags`
dw_LIBRARIES += `pkg-config mkl-static-lp64-iomp --cflags --libs`
dwbw_LIBRARIES += `pkg-config mkl-static-lp64-iomp --cflags --libs`
endif

ifeq ($(FFTW3), 1)
$(info -- Looking for FFTW3)
FFTW_EXISTS = $(shell pkg-config libtiff-4 --exists ; echo $$?)
ifeq ($(FFTW_EXISTS), 0)
$(info found)
else
$(error Could not find FFTW3)
endif

CFLAGS += `pkg-config fftw3 fftw3f --cflags`
ifeq ($(MOSTLYSTATIC), 0)
dw_LIBRARIES += `pkg-config fftw3 fftw3f --libs`
dwbw_LIBRARIES += `pkg-config fftw3 fftw3f --libs`
else
dw_LIBRARIES += -l:libfftw3.a -l:libfftw3f.a
dwbw_LIBRARIES += -l:libfftw3.a -l:libfftw3f.a
endif
ifneq ($(WINDOWS),1)
# _omp is faster than _threads on my machines
dw_LIBRARIES += -lfftw3f_omp
dwbw_LIBRARIES += -lfftw3f_omp
#dw_LIBRARIES += -lfftw3f_threads
#dwbw_LIBRARIES += -lfftw3f_threads
endif
endif

##
## OpenMP
##

OMP?=1
ifeq ($(OMP), 1)
$(info --Enabling OpenMP (OMP))
ifeq ($(UNAME_S), Darwin)
CFLAGS+=-Xpreprocessor -fopenmp
dw_LIBRARIES += -lomp
dwbw_LIBRARIES += -lomp
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


##
## GPU VKFFT + OpenCL
##


VKFFT?=0
ifeq ($(VKFFT), 1)
$(info -- Including VkFFT)
CFLAGS+=-DVKFFT_BACKEND=3
CFLAGS+=-Isrc/VkFFT/vkFFT/
OPENCL=1
endif

#
# OpenCL
#

OPENCL?=0
ifeq ($(OPENCL), 1)
$(info -- OpenCL enabled)
CFLAGS+=-DOPENCL
ifeq ($(UNAME_S), Darwin)
dw_LIBRARIES+=-framework OpenCL
else
$(info -- Looking for OpenCL)
OPENCL_EXISTS = $(shell pkg-config OpenCL --exists ; echo $$?)
ifeq ($(OPENCL_EXISTS),0)
$(info found)
CFLAGS+=`pkg-config OpenCL --cflags`
dw_LIBRARIES+=`pkg-config OpenCL --libs`
else
$(error ERROR: Could not find OpenCL)
endif
endif
dw_OBJECTS+=method_shb_cl.o method_shb_cl2.o cl_util.o
endif


###########################
# Platform specific extras
###########################

ifneq ($(UNAME_S),Darwin)
    CFLAGS+=
endif

ifeq ($(WINDOWS),1)
	CFLAGS += -DWINDOWS
endif

## MANPATH
MANPATH=/usr/share/man/man1/
ifeq ($(UNAME_S),Darwin)
	MANPATH=/usr/local/share/man/man1
endif

SRCDIR = src/

dw_OBJECTS += fim.o \
tiling.o \
fft.o \
fim_tiff.o \
dw.o deconwolf.o \
dw_maxproj.o \
dw_util.o \
method_identity.o \
method_rl.o \
method_shb.o \
dw_imshift.o \
fft.o \
dw_dots.o \
fwhm.o \
ftab.o \
dw_psf.o \
dw_tiff_merge.o \
dw_psf_sted.o
#dw_nuclei.o

dwbw_OBJECTS = fim.o \
fim_tiff.o \
dw_bwpsf.o \
bw_gsl.o \
lanczos.o \
li.o fft.o \
dw_util.o \
ftab.o

$(info Everything looks ok)

all: $(dw) $(dwbw)

$(dw): $(dw_OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(dw_LIBRARIES)

$(dwbw): $(dwbw_OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(dwbw_LIBRARIES)

%.o: $(SRCDIR)%.c
	$(CC) $(CFLAGS) -c $<

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
