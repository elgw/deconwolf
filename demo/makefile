dw_dapi_001.tif: PSF_dapi.tif
	dw maxproj dapi_001.tif
	dw --iter 40 dapi_001.tif PSF_dapi.tif --overwrite --noplan
	dw maxproj dw_dapi_001.tif

# No need to specify --nslice
PSF_dapi.tif:
	dw bw --resxy 130 --resz 300 --lambda 461 --NA 1.45 --ni 1.512 PSF_dapi.tif --overwrite --nslice 79

clean:
	rm -f dw_dapi_001.tif
	rm -f dw_dapi_001.tif.log.txt
	rm -f PSF_dapi.tif
	rm -f PSF_dapi.tif.log.txt
	rm -f max_dapi_001.tif
	rm -f max_dw_dapi_001.tif
