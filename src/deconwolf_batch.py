#!/bin/python3
import sys
import glob
import os


def im_to_psf(imfile):
    # Map image file name to psf file name
    imfile = os.path.basename(imfile)
    imfile = imfile[0:-4]
    parts = imfile.split('_')
    if(len(parts) == 2):
        return(f"PSF_{parts[0]}.tif")
    else:
        return f""


def usage():
    print("Usage:")
    print(f"{sys.argv[0]} image_folder psf_folder deonwolf_settings")
    print("Example, using default settings:")
    print(f"{sys.argv[0]} image_folder psf_folder ''")
    print("Example:")
    print(f'{sys.argv[0]} ./ ../psf/ "--iter 40 --tilesize 512"')
    print(f"which you could `| bash -c` or direct to a file ...")


if __name__ == '__main__':
    if(len(sys.argv) != 4):
        usage()
        sys.exit(1)

    imd = sys.argv[1]
    psd = sys.argv[2]
    dcw = sys.argv[3]
    images = glob.glob(imd + '*.tif')

    for im in images:
        psf_file = im_to_psf(im)
        if(len(psf_file)>0):
            psf_file = os.path.join(psd, psf_file)
            print(f"deconwolf {dcw} {im} {psf_file}")

    sys.exit(0)
