# Build deconwolf
#
# Usage:
# meson setup builddir
# cd builddir
# meson compile
#
# Inspired by
# https://gitlab.surrey.ac.uk/ri0005/binary_c/-/blob/455a3598c506cc157b45b90016ce896645d6ca0d/meson.build
#
# TODO:
# - OpenCL
# - -DNDEBUG -DVKFFT -DOpenCL
# 

project('deconwolf', 'c',
                 default_options : ['c_std=gnu11'])

subdir('src/')
subdir('src/random_forest/')
subdir('doc/')

# basic system requirements
system_requirements = {
                            'meson':'--version',
                            'ninja':'--version',
                            'gsl-config':'--version',
                        }
                        
foreach cmd, version_args : system_requirements
    _program = find_program(cmd, required:true)
    if _program.found()
           _program_path = _program.path()
           if version_args != ''
               _program_version = run_command(_program_path,version_args).stdout().strip()
               message(cmd + ' version "' + _program_version + '"')
           endif
    else
        message(cmd + 'not found : it is a requirement')
        exit
    endif
endforeach

########################################################
# C-flags and ld-flags

dw_cflags = []
dw_bw_cflags = []
dw_deps = []
dw_bw_deps = []


thread_dep = dependency('threads')
cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required : false)


dw_deps += [m_dep, thread_dep]
dw_bw_deps += [m_dep, thread_dep]


gsl_dep = dependency(
    'gsl',
    required : true)

dw_bw_deps += gsl_dep
dw_deps += gsl_dep

# pkg-config --list-all | grep tiff
tiff_dep = dependency(
    'libtiff-4',
    required : true)
    
dw_deps += tiff_dep
dw_bw_deps += tiff_dep

fftw3_dep = dependency(
    'fftw3',
    required : true)
dw_deps += fftw3_dep
dw_bw_deps += fftw3_dep

fftw3f_dep = dependency(
    'fftw3f',
    required : true)
dw_deps += fftw3f_dep
dw_bw_deps += fftw3f_dep

dw_ldflags = ['-lfftw3f_threads']
dw_bw_ldflags = ['-lfftw3f_threads']

compiler = meson.get_compiler('c')

if compiler.has_argument('-mtune=native')
   dw_cflags += ['-mtune=native']
endif

if compiler.has_argument('-march=native')
   dw_cflags += ['-march=native']
endif

if compiler.has_argument('-fno-math-errno')
   dw_cflags += ['-fno-math-errno']
endif


if cc.links('int main () { return 0; }', name: '-flto', args: '-flto')
   dw_ldflags += ['-flto']
   dw_bw_ldflags += ['-flto']
endif

omp_deps = dependency('openmp', required : false)
if omp_deps.found()
   dw_deps += omp_deps
   dw_bw_deps += omp_deps
endif



os = host_machine.system()
if os == 'darwin'
   # for OpenMP
   # installed with brew install libomp
   dw_cflags += ['-Xpreprocessor', '-fopenmp']
   dw_ldflags += ['-lomp']
endif

_program = find_program('git', required:false)
if _program.found()
   git_version = run_command(_program, ['log', '--pretty=format:\'%aD:%H\'', '-n', '1']).stdout().strip()
   message('git_version = ' + git_version)
   dw_cflags += '-DGIT_VERSION="' + git_version + '"'
else
      dw_cflags += '-DGIT_VERSION=UNKNOWN'
endif

executable('dw', dw_src,
                    dependencies: dw_deps,
                    install: true,
                    c_args: [dw_cflags],
                    link_args: [dw_ldflags])
                    
executable('dw_bw', dw_bw_src,
                        dependencies : dw_bw_deps,
                        install: true,
                        c_args: [dw_bw_cflags],
                        link_args: [dw_bw_ldflags])

