project('deconwolf', 'c',
        default_options : ['c_std=gnu99', 'optimization=3', 
        'warning_level=1', 'b_ndebug=true', 'b_lto=true'])

cc = meson.get_compiler('c')

src = ['src/deconwolf_cli.c', 'src/deconwolf.c', 'src/fim_tiff.c', 'src/tiling.c', 'src/fft.c', 'src/fim.c']

fftw3f_dep = dependency('fftw3f')
fftw3f_threads_dep = cc.find_library('fftw3f_threads')
tiff_dep = dependency('libtiff-4')
m_dep = cc.find_library('m', required : false)
dep = [fftw3f_dep, fftw3f_threads_dep, tiff_dep, m_dep]

executable('deconwolf', src, dependencies : dep, install : true)

