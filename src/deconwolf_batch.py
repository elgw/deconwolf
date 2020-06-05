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
    if(len(sys.argv) <= 3):
        usage()
        sys.exit(1)

    imd = sys.argv[1]  # image folder
    psd = sys.argv[2]  # psf folder
    dcw = ""  # default options to deconwolf
    if len(sys.argv) > 3:
        dcw = sys.argv[3]  # common options to deconwolf

    images_tif = glob.glob(imd + '*.tif')
    images_tiff = glob.glob(imd + '*.tiff')
    images = images_tif + images_tiff
    if len(images) == 0:
        print(f"No files ending with .tif, or .tiff in {imd}")
        sys.exit(1)

    psf_images = glob.glob(psd + '*.tif') + glob.glob(psd + '*.tiff')
    if len(psf_images) == 0:
        print(f"No files ending with .tif or .tiff in {psd}")
        sys.exit(1)

    for im in images:
        psf_file = im_to_psf(im)
        if(len(psf_file) > 0):
            psf_file = os.path.join(psd, psf_file)
            print(f"deconwolf {dcw} {im} {psf_file}")

    sys.exit(0)
