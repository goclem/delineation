/* 
# !/usr/local/stata16
# Description: Wrapper for delineation.R
# Author: Clement Gorin
# Contact: gorinclem@gmail.com
# Updated: 28 June 2021
*/

* Program
cap: program drop delineation
program define  delineation
	syntax, [density(string) unlivable(string) outdir(string) tmpdir(string) nboots(integer 100) niter(integer 1) bandwidth(integer 15) quantile(integer 95) replace(integer 1) joinall(integer 0) joinunl(integer 0) filter(integer 0) workers(integer -1) memory(integer -1) seed(integer 1) help pause]
	if "pause"!="" local pause "& `pause'"
	if "`help'"=="help" {
		!"${Rscript}" delineation.R --help & pause
	}
	else {
		!"${Rscript}" delineation.R --density="`density'" --unlivable="`unlivable'" --outdir="`outdir'" --tmpdir="`tmpdir'" --nboots=`nboots' --niter=`niter' --bandwidth=`bandwidth' --quantile=`quantile' --replace=`replace' --joinall=`joinall' --joinunl=`joinunl' --filter=`filter'  --workers=`workers' --seed=`seed' `pause'
	}
end
