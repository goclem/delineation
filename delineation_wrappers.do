/* 
# !/usr/local/stata16
# Description: Stata wrappers for the delineation scripts
# Author: Clement Gorin
# Contact: gorinclem@gmail.com
# Version: 2021.10.26
*/

* Wrapper for delineation_raster.R
cap: program drop delineation_raster
program define delineation_raster
	syntax, [density(string) unlivable(string) outdir(string) tmpdir(string) nboots(integer 100) niter(integer 1) bandwidth(integer 15) quantile(integer 95) replace(integer 1) joinall(integer 0) joinunl(integer 0) filter(integer 0) workers(integer -1) memory(integer -1) seed(integer 1) help pause]
	if "pause" != "" & c(os) == "Windows" local pause "& `pause'"
	if "`help'" == "help" {
		!"${Rscript}" delineation_raster.R --help `pause'
	}
	else {
		!"${Rscript}" delineation_raster.R --density="`density'" --unlivable="`unlivable'" --outdir="`outdir'" --tmpdir="`tmpdir'" --nboots=`nboots' --niter=`niter' --bandwidth=`bandwidth' --quantile=`quantile' --replace=`replace' --joinall=`joinall' --joinunl=`joinunl' --filter=`filter'  --workers=`workers' --seed=`seed' `pause'
	}
end

* Wrapper for delineation_building.R
cap: program drop delineation_building
program define delineation_building
	syntax, [density(string) reference(string) unlivable(string) outdir(string) tmpdir(string) nboots(integer 100) niter(integer 1) bandwidth(integer 15) quantile(integer 95) replace(integer 1) joinall(integer 0) joinunl(integer 0) filter(integer 0) workers(integer -1) memory(integer -1) seed(integer 1) help pause]
	if "pause" != "" & c(os) == "Windows" local pause "& `pause'"
	if "`help'" == "help" {
		!"${Rscript}" delineation_building.R --help `pause'
	}
	else {
		!"${Rscript}" delineation_building.R --density="`density'" --reference="`reference'" --unlivable="`unlivable'" --outdir="`outdir'" --tmpdir="`tmpdir'" --nboots=`nboots' --niter=`niter' --bandwidth=`bandwidth' --quantile=`quantile' --replace=`replace' --joinall=`joinall' --joinunl=`joinunl' --filter=`filter'  --workers=`workers' --seed=`seed' `pause'
	}
end
