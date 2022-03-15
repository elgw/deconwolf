# import deconwolf
# d = deconwolf.dw
# d.run(image_file = 'asdasf', psf_file = 'asdfasdf', ...)

import subprocess
import os

class dw:

    def run(image_file='', psf_file='', iters=50, prefix=''):
        assert(os.path.isfile(image_file))
        assert(os.path.isfile(psf_file))
        assert(iter > 0)

        cmd = ['dw', '--iter', str(iters), image_file, psf_file];
        if len(prefix) > 0:
            cmd.append('--prefix')
            cmd.append(str(iter))

        print('Will run ' + ' '.join(cmd))
        #subprocess.run(['dw', '--iter', str(iter), image_file, psf_file])
