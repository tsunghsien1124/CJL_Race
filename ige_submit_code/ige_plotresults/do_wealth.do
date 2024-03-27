global outputdir "/home/sylee/now/ananth/natnurt/final_real/results"
global statadir "/home/sylee/now/ananth/natnurt/stata_final_real/dtafiles"

*-----------------------------------------------------------------------------

import delimited "$outputdir/hhidgen1.out", delimiter(" ",collapse) clear
keep ss*
foreach var of varlist ss*{
	rename `var' `var'_p
}
egen pavgwlth = rmean(ss*)
save "$statadir/wealth.dta", replace

import delimited "$outputdir/hhidgen2.out", delimiter(" ",collapse) clear
keep ss*
foreach var of varlist ss*{
	rename `var' `var'_k
}
egen kavgwlth = rmean(ss*)
merge 1:1 _n using "$statadir/wealth.dta", nogen
save "$statadir/wealth.dta", replace

import delimited "$outputdir/hhidgen2.out", delimiter(" ",collapse) clear
keep ss*
foreach var of varlist ss*{
	rename `var' `var'_k
}
egen kavgwlth = rmean(ss*)
merge 1:1 _n using "$statadir/wealth.dta", nogen
save "$statadir/wealth.dta", replace

import delimited "$outputdir/lifetime.out", delimiter(" ",collapse) clear
keep *kwlth *pwlth
merge 1:1 _n using "$statadir/wealth", nogen
save "$statadir/wealth.dta", replace
*-----------------------------------------------------------------------------
*-----------------------------------------------------------------------------
g logp = .
g logk = .

cap prog drop welast
prog welast

	replace logk = log(`1')
	replace logp = log(`2')

	reg logk logp

end
*------------------

log using "$logdir/wealth.log", replace
	welast kavg lfpwlth
	welast ss_4_k lfpwlth
	welast ss_6_k lfpwlth

	welast kavg pavg
	welast ss_4_k pavg
	welast ss_6_k pavg

	g temp = ss_9_p+ss_10_p
	welast ss_6_k temp
log close
*-----------------------------------------------------------------------------
*-----------------------------------------------------------------------------
drop log*

g prank = .
g krank = .

cap prog drop wrrcorr
prog wrrcorr

	sort `1'
	replace krank = _n/_N*100

	sort `2'
	replace prank = _n/_N*100

	correlate krank prank

end
*------------------

log using "$logdir/wealth.log", append
	wrrcorr kavg lfpwlth
	wrrcorr ss_4_k lfpwlth
	wrrcorr ss_6_k lfpwlth

	wrrcorr kavg pavg
	wrrcorr ss_4_k pavg
	wrrcorr ss_6_k pavg

	replace temp = ss_9_p+ss_10_p
	wrrcorr ss_6_k temp
log close
