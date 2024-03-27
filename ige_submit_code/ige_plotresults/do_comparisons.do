capture log close
set matsize 800
set more off
set scheme s1mono

*-----------------------------------------------------------------------------
* hello
*-----------------------------------------------------------------------------
scalar NI = 120000

global outputdir "/home/sylee/now/ananth/natnurt/final_real/results"
global statadir "/home/sylee/now/ananth/natnurt/stata_final_real/dtafiles"
global logdir "/home/sylee/now/ananth/natnurt/stata_final_real/logfiles"
global graphdir "/home/sylee/now/ananth/natnurt/stata_final_real/graphs"
*-----------------------------------------
import delimited "$outputdir/bt_ss.out", delimiter(" ",collapse) clear
drop v1
foreach var of varlist *child *parent{
	rename `var' `var'_ss
}
replace abp_ss = log(abp_ss)
save "$statadir/bthm_dist.dta", replace
*-----------------------------------------
import delimited "$outputdir/bt_dist.out", delimiter(" ",collapse) clear
drop v1
foreach var of varlist *child *parent{
	rename `var' `var'_bt
}
merge 1:1 _n using "$statadir/bthm_dist.dta", nogen
save "$statadir/bthm_dist.dta", replace
*-----------------------------------------
import delimited "$outputdir/netincp.out", delimiter(" ",collapse) clear
drop v1

g netincp = inc12
forval j = 11(-1)4{
	replace netinc = inc`j' + netinc/(1.04)^6
}
keep netinc inc*

merge 1:1 _n using "$statadir/bthm_dist.dta", nogen
save "$statadir/bthm_dist.dta", replace
*-----------------------------------------
import delimited "$outputdir/netinck.out", delimiter(" ",collapse) clear
drop v1

g netinck = inc12
forval j = 11(-1)4{
	replace netinc = inc`j' + netinc/(1.04)^6
}
keep netinc

merge 1:1 _n using "$statadir/bthm_dist.dta", nogen
save "$statadir/bthm_dist.dta", replace
*-----------------------------------------
use hhid abi_p* lfp* ss*p* lfk* hh*k* hh*p* ss_4_k* using $statadir/benchmark.dta, clear
sort hhid
merge 1:1 _n using "$statadir/bthm_dist.dta", nogen
*---------------

*---------------
foreach var of varlist schild* hchild* netinc* *parent* resource*{
	sort `var'
	g `var'_r = _n/_N*100
}
*---------------
save "$statadir/bthm_dist.dta", replace


*-----------------------------------------
replace lfpavge = log(lfpavge)
replace lfkavge = log(lfkavge)
*-----------------------------------------------------------------------------
log using "$logdir/comparisons.log", replace
*-----------------------------------------------------------------------------
* total igc's
	correlate lfpavge_r lfkavge_r
	correlate hparent_bt_r hchild_bt_r
	correlate hparent_ss_r hchild_ss_r

	g learnpbt = log(hparent_bt)
	g learnkbt = log(hchild_bt)
	reg learnkbt learnpbt

	g learnpss = log(hparent_ss)
	g learnkss = log(hchild_ss)
	reg learnkss learnpss

	correlate schild_bt_r ss_4_p_r
	correlate schild_ss_r sparent_ss_r

	g lwlthpbt = log(ss_4_p)
	g lwlthkbt = log(schild_bt)
	reg lwlthkbt lwlthpbt

	g lwlthpss = log(sparent_ss)
	g lwlthkss = log(schild_ss)
	reg lwlthkss lwlthpss

log close


*-----------------------------------------------------------------------------
* igc's by intergenerational transfers
*-----------------------------------------------------------------------------
local temp = 0
forval i = 0/6{
	local temp = `temp' + 1/(1.04)^`i'
}
summ schild_bt
local bt = r(mean)/`temp'
summ schild_ss
local ss = r(mean)/`temp'

g byte upper_iv = ss_4_k>1
g byte upper_bt = schild_bt>`bt'
g byte upper_ss = schild_ss>`ss'

preserve
	collapse (sum) upper*
	local bm = 1-upper_iv/120000
	local bt = 1-upper_bt/120000
	local ss = 1-upper_ss/120000
restore

//sort ss_4_k
//g byte upper_iv = _n>_N/2

//sort schild_bt
//g byte upper_bt = _n>_N/2

//sort schild_ss
//g byte upper_ss = _n>_N/2

log using "$logdir/comparisons.log", append
* parents bequests and ige
	disp `bm'
	reg lfkavge lfpavge if !upper_iv
	reg lfkavge lfpavge if upper_iv

	disp `bt'
	reg learnkbt learnpbt if !upper_bt
	reg learnkbt learnpbt if upper_bt

	disp `ss'
	reg learnkss learnpss if !upper_ss
	reg learnkss learnpss if upper_ss

* bequests and earnings
	correlate ss_4_k_r lfkavge_r
	correlate schild_bt_r hchild_bt_r
	correlate schild_ss_r hchild_ss_r

log close

*-----------------------------------------------------------------------------
* igc's by parent's netwlth
*-----------------------------------------------------------------------------
g netwlth = netincp + ss_4_p

sort netwlth
g netwlth_r = _n/_N*100
g byte rich_par = _n>_N/2

sort resource
g resource_bt_r = _n/_N*100
g byte bt_par = _n>_N/2

g resource_ss = (1-0.118)*hparent_ss+sparent_ss
sort resource_ss
g resource_ss_r = _n/_N*100
g byte ss_par = _n>_N/2

log using "$logdir/comparisons.log", append
* parent net wealth and ige
	correlate netwlth_r lfkavge_r
	correlate resource_bt_r hchild_bt_r
	correlate resource_ss_r hchild_ss_r

	reg lfkavge lfpavge if !rich
	reg lfkavge lfpavge if rich

	reg learnkbt learnpbt if !bt_par
	reg learnkbt learnpbt if bt_par

	reg learnkss learnpss if !ss_par
	reg learnkss learnpss if ss_par

//	reg ss_4_k_r netwlth_r lfkavge_r
//!	reg schild_bt_r resource_bt_r hchild_bt_r
//	reg schild_ss_r resource_ss_r hchild_ss_r
log close


*-----------------------------------------------------------------------------
* graph: igc's by parent's netwlth
*-----------------------------------------------------------------------------
* fraction of high ability / a-h correlation by decile
*-----------------------------------------------------------------------------
preserve
	replace netwlth_r = ceil(netwlth_r/10)
	compress netwlth_r
	matrix coeff = J(10,1,.)
	quietly forval i=1/10{
		correlate ss_4_k_r lfkavge_r if netwlth_r==`i'
		matrix coeff[`i',1] = r(rho)
	}
	egen ahigh = count(abi_pi), by(netwlth_r abi_pi)
	replace ahigh = . if abi_pi<3

	collapse ahigh lfpavge_r abi_p lfpavge ss_4_k, by(netwlth_r)
	replace ahigh = ahigh/NI*10
	replace lfpavge_r  = lfpavge_r/100

	svmat coeff
	rename coeff1 benchmark

	keep netwlth_r ahigh lfpavge_r abi_p lfpavge benchmark ss_4_k
	save dtafiles/bthms4.dta, replace
restore
preserve
	replace netwlth_r = ceil(resource_bt_r/10)
	compress netwlth_r

	matrix coeff = J(10,1,.)
	quietly forval i=1/10{
		correlate schild_bt_r hchild_bt_r if netwlth_r==`i'
		matrix coeff[`i',1] = r(rho)
	}

	use dtafiles/bthms4.dta, clear
	svmat coeff
	rename coeff1 BTmodel
	save, replace
restore
preserve
	replace netwlth_r = ceil(resource_ss_r/10)
	compress netwlth_r

	matrix coeff = J(10,1,.)
	quietly forval i=1/10{
		correlate schild_ss_r hchild_ss_r if netwlth_r==`i'
		matrix coeff[`i',1] = r(rho)
	}
	egen ahss = count(abp_ss), by(netwlth_r abp_ss)
	replace ahss = . if abp_ss<0

	collapse ahss hparent_ss_r abp_ss learnpss schild_ss, by(netwlth_r )
	replace ahss = ahss/NI*10
	replace hparent_ss_r  = hparent_ss_r/100

	svmat coeff
	rename coeff1 SSmodel

	keep netwlth_r ahss hparent_ss_r abp_ss learnpss SSmodel schild_ss
	merge 1:1 netwlth_r using dtafiles/bthms4.dta, nogen
	save dtafiles/bthms4.dta, replace

	tw connected benchmark BTmodel SSmodel netwlth_r, ///
		, l1title("Rank Correlation") xtitle("Parents' Net Wealth Decile") ///
		 yline(0) ytitle("") ylabel(,grid) xlabel(1(1)10,grid) legend(cols(3) order (1 "Benchmark" 2 "BT Short" 3 "BT Long"))
	graph export "$graphdir/bthscorr.eps", replace

	summ abi_p
	local bmean = r(mean)
	summ abp_ss
	local smean = r(mean)

	replace abi_p = abi_p-`bmean'
	replace abp_ss = abp_ss-`smean'

	tw connected abi_p abp_ss netwlth_r, ///
		, l1title("Parent's Log Ability") xtitle("Parents' Net Wealth Decile") ///
		 yline(0) ytitle("") ylabel(,grid) xlabel(1(1)10,grid) legend(cols(2) order (1 "Benchmark" 2 "BT Long"))
	graph export "$graphdir/btavgab.eps", replace

	summ lfpavge
	local bmean = r(mean)
	summ learnpss
	local smean = r(mean)

	replace lfpavge = lfpavge-`bmean'
	replace learnpss = learnpss-`smean'

	tw connected lfpavge learnpss netwlth_r, ///
		, l1title("Parent's Lifetime Avg Earnings") xtitle("Parents' Net Wealth Decile") ///
		 yline(0) ytitle("") ylabel(,grid) xlabel(1(1)10,grid) legend(cols(2) order (1 "Benchmark" 2 "BT Long"))
	graph export "$graphdir/btavgearn.eps", replace

	summ ss_4_k
//	local bmean = r(mean)
	summ schild_ss
//	local smean = r(mean)

//	replace ss_4_k = ss_4_k -`bmean'
//	replace schild_ss = schild_ss -`smean'

	tw connected ss_4_k schild_ss netwlth_r, ///
		, l1title("Child's Initial Wealth") xtitle("Parents' Net Wealth Decile") ///
		 yline(0) ytitle("") ylabel(,grid) xlabel(1(1)10,grid) legend(cols(2) order (1 "Benchmark" 2 "BT Long"))
	graph export "$graphdir/btavgs4.eps", replace

restore
*-----------------------------------------------------------------------------*-----------------------------------------------------------------------------

*-----------------------------------------------------------------------------
* wealth when young
*-----------------------------------------------------------------------------
drop netwlth* rich_par

* parent net wealth and child's earnings
g netwlth = inc4+ss_4_p
sort netwlth
g netwlth_r = _n/_N*100
log using "$logdir/comparisons.log", append
		correlate netwlth_r lfkavge_r
log close

*-----------------------------------------------------------------------------
matrix pwkescorr = J(8,4,.)
forval i = 4/8{
	g yvar = inc`i'+ss_`i'_p
//	forval j=4/`i'{
//		replace yvar = yvar - inc`j'/(1.04)^(6*(`j'-4))
//	}
	sort yvar
	replace yvar = _n/_N*100

	correlate yvar hh_4_k_r
	matrix pwkescorr[`i',1] = r(rho)
	correlate yvar ss_4_k_r
	matrix pwkescorr[`i',2] = r(rho)

* parent wlth and ige
	replace yvar = hh_`i'_p
	sort yvar
	replace yvar = _n/_N*100

	correlate yvar hh_4_k_r
	matrix pwkescorr[`i',3] = r(rho)

	correlate yvar ss_4_k_r
	matrix pwkescorr[`i',4] = r(rho)

	drop yvar
}

svmat pwkescorr
rename pwkescorr3 pwkecorr
rename pwkescorr4 pwkscorr
g age = _n*6

preserve
	drop if missing(pwkecorr)

	tw connected pwkecorr pwkscorr age, ///
	, l1title("Rank Correlation") xtitle("Parents' Age") ///
	 yline(0) ytitle("") ylabel(,grid) xlabel(24(6)48,grid) legend(cols(2) order (2 "wt Child's Initial Wealth" 1 "wt Child's Initial Human Capital"))
	graph export "$graphdir/btcomp.eps", replace

restore
*-----------------------------------------------------------------------------
/*-----------------------------------------------------------------------------

preserve
forval j = 4/10{
	egen hkcorr`j' = corr(lfpearn_r hh_`j'_k)
	egen khkcorr`j' = corr(lfkearn_r hh_`j'_k)
}
keep if _n<=2
keep *corr*

forval j = 4/10{
	replace hkcorr`j' = khkcorr`j' if _n==2
}
keep hkcorr*
g byte gen = _n-1

reshape long hkcorr, i(gen) j(age)
reshape wide hkcorr, i(age) j(gen)

replace age = age*6+3
tw (connected hkcorr1 age) (connected hkcorr0 age), xlabel(24(6)66) xtitle("Age") ytitle("Rank Correlation") legend(cols(2) order(1 "with Own" 2 "with Parents'"))
graph export "$graphdir/rrcorr_hh.eps", replace
restore
