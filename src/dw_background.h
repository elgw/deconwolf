#pragma once
#define dw_module_background

/* Average a set of images to estimate vignetting. Finish with a
 * Gaussian low pass filter.
 *
 * It might be better to use a parametric model, than to use low pass
 * filtering. See for example:
 *
 * Bal A, Palus H. Image Vignetting Correction Using a Deformable
 * Radial Polynomial Model. Sensors (Basel). 2023 Jan
 * 19;23(3):1157. doi: 10.3390/s23031157. PMID: 36772200; PMCID:
 * PMC9921563.
 **/

int
dw_background(int argc, char ** argv);
