import os
import sys
import glob

def getConfig(NA=1.45, lam=780, method="BW", resz=200,
              resxy=130.0, NZ = 181, NXY = 181, nIM=1.51,
              nSample=1.33,
              wDist=150, # um
              ppos=2000): # nm
    # NA: numerical aperture
    # lam: lambda (nm)
    # method: RW or BW (or GL)
    # resxy: lateral pixel size (nm)
    # resz: axial pixel size (nm)
    # NXY: Number of pixels in x and y
    # NZ: number of pixels in z
    # nIM: n for the immersion oil
    # nSample only used for GL
    # wDist only for GL -- working distance of the objective
    # ppos only for GL -- particle position

    config = f"""
    #PSFGenerator
    #Wed May 13 15:44:13 CEST 2020
    psf-Circular-Pupil-defocus=100.0
    psf-RW-NI={nIM}
    psf-Oriented-Gaussian-axial=Linear
    ResLateral={resxy}
    psf-Defocus-DBot=30.0
    psf-Oriented-Gaussian-focus=0.0
    psf-TV-NI=1.5
    psf-Astigmatism-axial=Linear
    psf-RW-accuracy=Best
    psf-Cardinale-Sine-axial=Linear
    psf-Lorentz-axial=Linear
    psf-Gaussian-defocus=100.0
    psf-BW-accuracy=Best
    psf-Cardinale-Sine-focus=0.0
    psf-GL-ZPos={ppos}
    psf-Astigmatism-focus=0.0
    psf-Lorentz-focus=0.0
    psf-Koehler-dMid=3.0
    psf-TV-TI=150.0
    psf-Circular-Pupil-axial=Linear
    psf-Koehler-dTop=1.5
    psf-Koehler-n1=1.0
    psf-Lorentz-defocus=100.0
    psf-Koehler-n0=1.5
    psf-Circular-Pupil-focus=0.0
    psf-Koehler-dBot=6.0
    psf-Defocus-ZI=2000.0
    psf-VRIGL-NS2=1.4
    psf-VRIGL-NS1=1.33
    psf-TV-ZPos=2000.0
    psf-Double-Helix-defocus=100.0
    psf-Cosine-defocus=100.0
    Scale=Linear
    psf-BW-NI={nIM}
    psf-GL-NS={nSample}
    Lambda={lam}
    PSF-shortname={method}
    psf-GL-NI={nIM}
    psf-Astigmatism-defocus=100.0
    ResAxial={resz}
    psf-GL-TI={wDist}
    psf-Double-Helix-axial=Linear
    LUT=Grays
    psf-VRIGL-RIvary=Linear
    psf-Double-Helix-focus=0.0
    NZ={NZ}
    NY={NXY}
    NX={NXY}
    psf-VRIGL-accuracy=Best
    psf-Gaussian-axial=Linear
    psf-VRIGL-ZPos=2000.0
    psf-Gaussian-focus=0.0
    psf-Cardinale-Sine-defocus=100.0
    psf-VRIGL-NI=1.5
    psf-VRIGL-NG=1.5
    psf-Oriented-Gaussian-defocus=100.0
    NA={NA}
    Type=32-bits
    psf-VRIGL-TI=150.0
    psf-Defocus-DMid=1.0
    psf-VRIGL-TG=170.0
    psf-Defocus-K=275.0
    psf-Cosine-axial=Linear
    psf-Cosine-focus=0.0
    psf-TV-NS=1.0
    psf-GL-accuracy=Best
    psf-Defocus-DTop=30.0"""
    return config

# Common settings
ss_BS1_100x = {'run' : 1, # set to 0 for dry run
      'NA' : 1.45,
      'resxy' : 130.0,
      'resz' : 200,
      'nIM' : 1.51,
      'NZ': 181,
      'NXY' : 181,
      'method': "BW",
      'oversample': 0,
      'nSample' : 1.33,
      'wDist': 130, # um
      'ppos' : 2000,
      'mag' : '100x'}

ss_BS1_60x = {'run' : 1, # set to 0 for dry run
      'NA' : 1.4,
      'resxy' : 216.6,
      'resz' : 300,
      'nIM' : 1.515,
      'NZ': 181,
      'NXY' : 181,
      'method': "BW",
      'oversample': 0,
      'nSample' : 1.33,
      'wDist': 130, # um (Su said 140 ...)
      'ppos' : 2000,
      'mag' : '60x'}

ss_BS2 = {'run' : 1, # set to 0 for dry run
      'NA' : 1.4,
      'resxy' : 108.3,
      'resz' : 600,
      'nIM' : 1.515,
      'NZ': 181,
      'NXY' : 181,
      'method': "BW",
      'oversample': 0,
      'nSample' : 1.33,
      'wDist': 130, # um (Su said 140 ...)
      'ppos' : 2000,
      'mag' : 'xx'}

ss = ss_BS1_100x

# Based on the dyes
channels = {'dapi': 461,
            'Cy5': 664, # alexa 647
            'a594': 617,
            'ir800': 794, # or 814 if it is alexa790
            'a700':723,
            'a488': 519,
            'tmr': 562}

ss['outfolder'] = f"PSF{ss['method']}_{ss['mag']}_{ss['resz']}/"

if ss['oversample'] > 1:
    print("Warning: You have to bin the PSF afterwards!")
    x = ss['oversample']
    ss['NZ']=x*ss['NZ']
    ss['NXY']=x*ss['NXY']
    ss['resxy'] = ss['resxy']/x
    ss['resz'] = ss['resz']/x
    ss['outfolder'] = 'PSF3X/'

print("Generating PSFs using the following settings:")
for s, v in ss.items():
    print(f"{s}={v}")

print("Channels and wavelengths")
for cc, lam in channels.items():
    print(f"channel: {cc} lambda: {lam}")


try:
    os.mkdir(ss['outfolder'])
except FileExistsError:
    print(f"{ss['outfolder']} did already exist")
    # ok, if folded did exist

with open(f"{ss['outfolder']}PSF_README.txt", "w") as ofile:
    ofile.write("> PSFs generated from runPSFGenerator.py\n")
    ofile.write("> Generating PSFs using the following settings:\n")
    for s, v in ss.items():
        ofile.write(f"{s}={v}\n")

    print("> Channels and wavelengths\n")
    for cc, lam in channels.items():
        ofile.write(f"channel: {cc} lambda: {lam}\n")


for cc, lam in channels.items():
    print(f"channel: {cc} lambda: {lam}")
    cfname = f"{ss['outfolder']}config_{cc}.txt"
    config = getConfig(NA=ss['NA'], resxy=ss['resxy'], resz=ss['resz'],
                       nIM=ss['nIM'], NZ=ss['NZ'], NXY=ss['NXY'],
                       lam=lam, method=ss['method'],
                       nSample=ss['nSample'], wDist=ss['wDist'],
                       ppos=ss['ppos']);
    print(f"Writing config to {cfname}")
    with open(cfname, "w") as cf:
        cf.write(config)

    cmd = f"java -cp PSFGenerator.jar PSFGenerator {cfname}"
    if ss['run']:
        os.system(cmd)
    else:
        print(f"# {cmd}")

    # PSFGenerator always outputs in the current folder
    outfile = f"PSF {ss['method']}.tif"
    psffile = f"{ss['outfolder']}PSF_{cc}.tif"

    if ss['run']:
        os.rename(outfile, psffile)
    else:
        print(f"# mv {outfile} {psffile}")

