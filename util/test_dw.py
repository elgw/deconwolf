main_help_msg = """
Tests expected behavior of a deconwolf binary

Usage:

Always specify what binary to use, for example:

$ python test_dw.py --dw ../build/dw

Notes:

- Most modules have their own unit tests which can be built
within the src/ folder.

"""
import os, sys, argparse, subprocess, tempfile, shutil
import numpy as np

def errol():
    print("🔴 ", end="")

def check_return_status_when_output_exists(dw):
    cmd = [dw, 'temp/input.tif', 'temp/psf.tif']

    print("")
    print(f"Checking the return value when the output file exists")
    print(f"$ {' '.join(cmd)}")

    if not os.path.isdir('temp'):
        os.mkdir('temp')
    open('temp/input.tif', 'a').close()
    open('temp/psf.tif', 'a').close()
    open('temp/dw_input.tif', 'a').close()

    status = subprocess.run(cmd, capture_output=True)
    if(status.returncode != 0):
        errol()
        print(f"dw should return status 0 if output file already exists")
        print(f"CMD: {' '.join(cmd)}")

def check_help_section(command):
    """Check that there is some consistency in how the help sections
    are written """

    print("")
    print(f"Checking the help section for")
    print(f"$ {' '.join(command)}")

    res = subprocess.run(command, capture_output=True)
    if res.returncode != 0:
        errol()
        print(f"Non-zero return code: {res.returncode}")
        return

    assert(len(res.stderr) == 0)
    for n, line in enumerate(res.stdout.splitlines()):

        line = line.decode("utf-8")
        line = line.replace('\t', '        ')
        # Not too wide output
        if len(line) >80:
            # tab to 8 whitespaces
            errol()
            print(f"line {n+1} is {len(line)} chars, >80")
            print(line)
        # Not to many leading white spaces
        if len(line) > 8:
            if line[0:9] == '         ':
                errol()
                print(f"Too many initial white spaces on line {n+1}")
                print(line)
        # If the the first characters are '--' then there should be no
        # space before that
        test = line.strip()
        if len(test) > 3:
            #breakpoint()
            if(test[0:1] == '-'):
                if(line[0:3] != '  -'):
                    errol()
                    print(f" Not 2 ws before initial '-' on line {n+1}")
                    print(line)

def test_tif_npy(image):
    """ Convert
    image/tif1 -> npy1 -> tif2 -> npy2
    and finally compare npy1 and npy2
    """
    print("")
    print("Checking tif2npy and npy2tif")
    tempdir0 = tempfile.TemporaryDirectory(prefix='test_dw_')
    tempdir = tempdir0.name
    # print(f"tempdir={tempdir}")
    shutil.copy(args.image, tempdir)
    # Check that tif2npy and npy2tif works
    imfile = os.path.basename(args.image)
    tiffile = tempdir + os.path.sep + imfile
    tiffile2 = tempdir + os.path.sep + "resaved_" + imfile
    npyfile = tiffile + ".npy"
    npyfile2 = tiffile2 + ".npy"
    command = [dw, 'tif2npy', tiffile, npyfile]
    res = subprocess.run(command, capture_output=True)
    if res.returncode != 0:
        errol()
        print(f"Non-zero return code: {res.returncode}")

    command = [dw, 'npy2tif', npyfile, tiffile2]
    res = subprocess.run(command, capture_output=True)
    if res.returncode != 0:
        errol()
        print(f"Non-zero return code: {res.returncode}")

    command = [dw, 'tif2npy', tiffile2, npyfile2]
    res = subprocess.run(command, capture_output=True)
    if res.returncode != 0:
        errol()
        print(f"Non-zero return code: {res.returncode}")

    A = np.load(npyfile).astype('f8')
    B = np.load(npyfile2).astype('f8')
    assert(A.shape == B.shape)
    assert(abs(np.mean(A) - np.mean(B)) < 0.1)
    c = np.corrcoef(A.flatten(), B.flatten())
    assert(c[0][1] > 0.99)
    tempdir0.cleanup()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=main_help_msg)
    parser.add_argument('--dw',
                        help='path to the dw binary',
                        required=True)
    parser.add_argument('--image',
                        help='path to a small test image',
                        required=False)

    args = parser.parse_args()

    if not args.image:
        print('No --image given, skipping some tests')

    dw = args.dw

    # Check the consistency of the help sections

    check_help_section([dw, '--help'])
    subcmds = ['maxproj', 'merge', 'dots', 'psf',
               'psf-STED', 'nuclei', 'background',
               'tif2npy', 'npy2tif', 'imshift', 'align-dots'];
    for subcmd in subcmds:
        check_help_section([dw, subcmd, '--help'])

    # Test conversion between tif and npy
    if args.image:
        test_tif_npy(args.image)

    #
    # Expected behavior
    #

    check_return_status_when_output_exists(dw)

    # Check if CPU and GPU implementations give similar results

    # Check that tiling works

    # at least that it does not crash

    # ... and whatmore ...
