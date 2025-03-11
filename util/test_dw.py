"""
Some tests will be added here.

Note: Most modules have their own unit tests which can be built
within the src/ folder.

usage:
$ python test_dw.py --dw ../build/dw
"""
import os, sys, argparse, subprocess

def check_help_section(command):
    """Check that there is some consistency in how the help sections
    are written """

    print("")
    print(f"Checking the help section for")
    print(f"$ {' '.join(command)}")

    res = subprocess.run(command, capture_output=True)
    if res.returncode != 0:
        print(f"Non-zero return code: {res.returncode}")
        return

    assert(len(res.stderr) == 0)
    for n, line in enumerate(res.stdout.splitlines()):

        line = line.decode("utf-8")
        line = line.replace('\t', '        ')
        # Not too wide output
        if len(line) >80:
            # tab to 8 whitespaces
            print(f"line {n} is {len(line)} chars, >80")
            print(line)
        # Not to many leading white spaces
        if len(line) > 8:
            if line[0:9] == '         ':
                print(f"Too many initial white spaces on line {n}")
                print(line)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A test script for deconwolf')
    parser.add_argument('--dw', help='path to the dw binary', required=True)

    args = parser.parse_args()

    dw = args.dw

    # Check the consistency of the help sections

    check_help_section([dw, '--help'])
    subcmds = ['maxproj', 'merge', 'dots', 'psf',
               'psf-STED', 'nuclei', 'background',
               'tif2npy', 'npy2tif'];
    for subcmd in subcmds:
        check_help_section([dw, subcmd, '--help'])

    # Check that tif2npy and npy2tif works

    # Check if CPU and GPU implementations give similar results

    # Check that tiling works
    # at least that it does not crash
