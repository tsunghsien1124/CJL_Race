global outputdir "/home/sylee/now/ananth/natnurt/final_real/results"
global statadir "/home/sylee/now/ananth/natnurt/stata_final_real/dtafiles"

*-----------------------------------------------------------------------------
* load data
*-----------------------------------------------------------------------------

import delimited "$outputdir/hhidgen1.out", delimiter(" ",collapse) clear
keep abi col hh_* ss_*
foreach var of varlist abi-ss_11{
	rename `var' `var'_p
}

g hhid = _n
save "$statadir/benchmark.dta", replace

import delimited "$outputdir/hhidgen2.out", delimiter(" ",collapse) clear
keep abi col ll_0 ll_1 ll_2 hh_* ss_*
foreach var of varlist abi-ss_11{
	rename `var' `var'_k
}

g hhid = _n
merge 1:1 hhid using "$statadir/benchmark.dta", nogen
save "$statadir/benchmark.dta", replace

import delimited "$outputdir/hhidgen3.out", delimiter(" ",collapse) clear
keep abi col hh_1 hh_2 hh_3 hh_4 ss_4
foreach var of varlist abi-ss_4{
	rename `var' `var'_g
}

g hhid = _n
merge 1:1 hhid using "$statadir/benchmark.dta", nogen
save "$statadir/benchmark.dta", replace

import delimited "$outputdir/lifetime.out", delimiter(" ",collapse) clear
drop v1

g hhid = _n
merge 1:1 hhid using "$statadir/benchmark.dta", nogen
save "$statadir/benchmark.dta", replace

*-----------------------------------------------------------------------------
* create new benchmark with ranks
*-----------------------------------------------------------------------------

order hhid abi*
foreach var of varlist lfkavge-ss_11_p{
	sort `var'
	g `var'_r = _n/_N*100
}

*-----------------------------------------------------------------------------
* replace abi wt actual values
*-----------------------------------------------------------------------------
scalar sig_a = 0.30259820598356931
scalar mu_a= 0.82887605607062897

scalar agrid1=0.5163
scalar agrid2=0.7921
scalar agrid3=1.2151
foreach gen in _k _p _g {
	g abi`gen'i = abi`gen'
	replace abi`gen' = log(agrid1) if abi`gen'i==1
	replace abi`gen' = log(agrid2) if abi`gen'i==2
	replace abi`gen' = log(agrid3) if abi`gen'i==3
}

*-----------------------------------------------------------------------------
save "$statadir/benchmark.dta", replace
