/* 
# !/usr/local/stata16
# Description: R utilities for Stata
# Author: Clement Gorin
# Contact:gorinclem@gmail.com
# Updated: 2021/01/25
*/

* Globals
global Rscript //[path to Rscript]
global script  //[path to delineation script]

* Wrapper fordelineation.R
cap: program drop delineation
program define  delineation
	syntax, [density(string) unlivable(string) outdir(string) tmpdir(string) nboots(integer 100) bandwidth(integer 15) quantile(integer 95) replace(integer 1) joinall(integer 1) joinunl(integer 2) filter(integer 1) workers(integer -1) memory(integer -1) seed(integer 1) help pause]
	if "pause"!="" local pause "& `pause'"
	if "`help'"=="help" {
		!"${Rscript}" "${script}" --help & pause
	}
	else {
		!"${Rscript}" "${script}" --density="`density'" --unlivable="`unlivable'" --outdir="`outdir'" --tmpdir="`tmpdir'" --nboots=`nboots' --bandwidth=`bandwidth' --quantile=`quantile' --replace=`replace' --joinall=`joinall' --joinunl=`joinunl' --filter=`filter'  --workers=`workers' --seed=`seed' `pause'
	}
end
