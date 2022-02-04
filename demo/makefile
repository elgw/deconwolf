threads=8

dw_dapi_001.tif: PSF_dapi.tif
	dw --threads $(threads) --iter 20 dapi_001.tif PSF_dapi.tif

PSF_dapi.tif:
	dw_bw --threads $(threads) --resxy 130 --resz 300 --lambda 461 --NA 1.45 --ni 1.512 PSF_dapi.tif

clean:
	rm dw_dapi_001.tif
	rm dw_dapi_001.tif.log.txt
	rm PSF_dapi.tif
	rm PSF_dapi.tif.log.txt
