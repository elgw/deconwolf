# Note: Please use cmake to build deconwolf
#       the makefile is only for testing
#
# -- Normal build:
# make -B
#
# -- With GPU acceleration
# make kernels
# make -B VKFFT=1
#
# To build with debug flags and no OpenMP
# make DEBUG=1 OMP=0 DEBUG=1
#
# There is also ...
# -D_FORTIFY_SOURCE=2
# -D_FORTIFY_SOURCE=3

DESTDIR?=/usr/local/bin
DEBUG?=0

UNAME_S := $(shell uname -s)
$(info Host type: $(UNAME_S))

dw = bin/dw
dwbw = bin/dw_bw

CFLAGS += -std=gnu11 -Wall -Wextra

CFLAGS+=-Isrc/kdtree/include

FORTIFY?=0
ifeq ($(FORTIFY),1)
CFLAGS+=-D_FORTIFY_SOURCE=3
endif

SANITIZE?=0
ifeq ($(SANITIZE),1)
CFLAGS+=-fsanitize=address -static-libasan
endif

ANALYZE?=0
ifeq ($(ANALYZE),1)
CFLAGS+=-fanalyzer
endif

# FFT Backend, pick __one__
FFTW3=1
MKL?=0

# GPU acceleration with OpenCL/VkFFT
VKFFT?=1

## MANPATH
MANPATH=/usr/share/man/man1/
ifeq ($(UNAME_S),Darwin)
MANPATH=/usr/local/share/man/man1
endif

PKGCONF=pkg-config

#
# Bake some information
#

CC_VERSION = "$(shell $(CC) --version | head -n 1)"

$(info -- Checking git commit version)
ifneq ($(wildcard .git/.*),)
GIT_VERSION = "$(shell git log --pretty=format:'%aD:%H' -n 1)"
else
$(info Building otside of the git repository)
GIT_VERSION= "unknown (.git not available)"
endif
$(info found $(GIT_VERSION))

CFLAGS += -DCC_VERSION=\"$(CC_VERSION)\"
CFLAGS += -DGIT_VERSION=\"$(GIT_VERSION)\"

#
# Tools available?
#

$(info -- Looking for compiler ($(CC)))
ifeq (, $(shell which $(CC)))
	$(error Can not find $(CC))
endif

$(info -- Looking pkg-config)
ifeq (, $(shell which pkg-config))
ifeq (, $(shell which pkgconf))
$(error Can not find pkg-config or pkgconf)
else
PKGCONF=pkgconf
endif
endif

# Only needed for developers that edit the kernels
$(info -- Looking for xxd)
ifeq (, $(shell which xxd))
$(info xxd not found)
endif


#
# Optimization flags
#

# Notes:
# -O2 -ftree-vectorize and -O3 give about the same performance
# -DNDEBUG turns off some self-tests
# -fno-math-errno gives no relevant performance gain

ifeq ($(DEBUG),1)
    CFLAGS += -O0 -Wno-unknown-pragmas -g3
else
    CFLAGS += -O3 -DNDEBUG
endif

# possibly good to compile with -fanalyzer now and then,
# however without VkFFT ...

#
# Check if the linker can perform link time optimizations
#
ifeq ($(DEBUG),0)
define check_ld_flag
$(shell $(LD) $(1) -e 0 /dev/null 2>/dev/null && echo $(1))
endef

$(info -- testing -flto and -flto=auto)
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
endif

#
# Set up for linking
#

dw_LIBRARIES=
dwbw_LIBRARIES=

#
# kd tree
#
dw_LIBRARIES+=-Lsrc/kdtree/ -lkdtree

FORCE: ;

src/kdtree/libkdtree.a: FORCE
	$(MAKE) -C $(@D) libkdtree.a


#
# Math library
#

$(info -- Math (-lm) library?)
LD_HAS_LM=$(call check_ld_flag,-lm)
ifneq ($(LD_HAS_LM),)
$(info will use -lm)
dw_libraries+=-lm
else
$(info Not available, might be built in)
endif

#
# Tiff library
#

$(info -- Looking for libtiff)
TIFFLIB=

ifeq ($(TIFFLIB),)
LIBTIFF6 = $(shell $(PKGCONF) libtiff-6 --exists ; echo $$?)
ifeq ($(LIBTIFF6),0)
TIFFLIB=libtiff-6
endif
endif


ifeq ($(TIFFLIB),)
LIBTIFF5 = $(shell $(PKGCONF) libtiff-5 --exists ; echo $$?)
ifeq ($(LIBTIFF4),0)
TIFFLIB=libtiff-5
endif
endif

ifeq ($(TIFFLIB),)
LIBTIFF4 = $(shell $(PKGCONF) libtiff-4 --exists ; echo $$?)
ifeq ($(LIBTIFF4),0)
	TIFFLIB=libtiff-4
endif
endif

ifeq ($(TIFFLIB),)
LIBTIFF = $(shell $(PKGCONF) libtiff --exists ; echo $$?)
ifeq ($(LIBTIFF4),0)
TIFFLIB=libtiff
endif
endif


ifneq ($(TIFFLIB),)
$(info found $(TIFFLIB))
CFLAGS+=$(shell $(PKGCONF) ${TIFFLIB} --cflags)
dw_LIBRARIES+=$(shell $(PKGCONF) ${TIFFLIB} --libs)
dwbw_LIBRARIES+=$(shell $(PKGCONF) ${TIFFLIB} --libs)
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

# There can be only one
ifeq ($(FFTW3), 1)
MKL=0
endif

ifeq ($(MKL), 1)
FFTW3=0
endif


ifeq ($(MKL), 1)
$(info FFTW backend: MKL)
CFLAGS += -DMKL $(shell $(PKGCONF) mkl-static-lp64-iomp --cflags )
dw_LIBRARIES += $(shell $(PKGCONF) mkl-static-lp64-iomp --cflags --libs)
dwbw_LIBRARIES += $(shell $(PKGCONF) mkl-static-lp64-iomp --cflags --libs)
endif

ifeq ($(FFTW3), 1)
$(info -- Looking for FFTW3)
FFTW_EXISTS = $(shell $(PKGCONF) fftw3 --exists ; echo $$?)
ifeq ($(FFTW_EXISTS), 0)

else
$(error Could not find FFTW3)
endif

CFLAGS += $(shell $(PKGCONF) fftw3 fftw3f --cflags)
ifeq ($(MOSTLYSTATIC), 0)
dw_LIBRARIES += $(shell $(PKGCONF) fftw3 fftw3f --libs)
dwbw_LIBRARIES += $(shell $(PKGCONF) fftw3 fftw3f --libs)
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

# OpenMP is required. On windows, compile
# with clang

$(info -- Enabling OpenMP (OMP))
ifeq ($(UNAME_S), Darwin)
CFLAGS+=-Xpreprocessor -fopenmp
dw_LIBRARIES += -lomp
dwbw_LIBRARIES += -lomp
else
ifeq ($(CC),clang)
# Ubuntu 24: requires also
# sudo apt-get install libomp-dev libclang-rt-dev
CFLAGS+=-fopenmp=libomp
else
CFLAGS += -fopenmp
endif
endif



##
## GPU VKFFT + OpenCL
##



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
OPENCL_EXISTS = $(shell $(PKGCONF) OpenCL --exists ; echo $$?)
ifeq ($(OPENCL_EXISTS),0)

CFLAGS+=$(shell $(PKGCONF) OpenCL --cflags)
dw_LIBRARIES+=$(shell $(PKGCONF) OpenCL --libs)
else
$(error ERROR: Could not find OpenCL)
endif
endif
dw_OBJECTS+=method_shb_cl.o method_shb_cl2.o cl_util.o
endif


##
## Platform specific extras
##


ifneq ($(UNAME_S),Darwin)
    CFLAGS+=
endif

ifeq ($(WINDOWS),1)
	CFLAGS += -DWINDOWS
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
ftab.o \
dw_psf.o \
dw_tiff_merge.o \
dw_psf_sted.o \
sparse_preprocess.o \
sparse_preprocess_cli.o \
gmlfit.o \
dw_align_dots.o
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
