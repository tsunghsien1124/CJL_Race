capture log close
set matsize 800
set more off
set scheme s1mono
*-----------------------------------------------------------------------------
* hello
*-----------------------------------------------------------------------------

global outputdir "/home/sylee/now/ananth/natnurt/final_real/results"
global statadir "/home/sylee/now/ananth/natnurt/stata_final_real/dtafiles"
global logdir "/home/sylee/now/ananth/natnurt/stata_final_real/logfiles"

//do do_load_data
use "$statadir/benchmark.dta", clear

*-----------------------------------------------------------------------------
* parametric variances
*-----------------------------------------------------------------------------

preserve
	g ss_4_k1 = log(ss_4_k)
	replace hh_4_k = log(hh_4_k)

	egen shifter = min(ss_4_k)
	replace shifter = -shifter + 1e-16
	replace ss_4_k = log(ss_4_k + shifter)

	egen hmean = mean(hh_4_k)
	egen smean = mean(ss_4_k)
	egen smean1 = mean(ss_4_k1) if !missing(ss_4_k1)

	egen hvar = sd(hh_4_k)
	egen svar = sd(ss_4_k)
	egen svar1 = sd(ss_4_k1) if !missing(ss_4_k1)

	replace hvar = hvar^2
	replace svar = svar^2
	replace svar1 = svar1^2

	log using "$logdir/param_dist.log", replace
		disp(log(mu_a)-0.5*sig_a^2)
		tab hmean
		tab smean
		tab smean1

		disp(sig_a^2)
		tab hvar
		tab svar
		tab svar1

		correlate abi_k hh_4_k
		correlate abi_k ss_4_k
		correlate hh_4_k ss_4_k

		correlate abi_k ss_4_k1
		correlate hh_4_k ss_4_k1

	log close
restore


*-------------------------------
* variance
*-------------------------------

local one = 33.33333334
local two = 66.66666667
local thr = 100.00000001

egen hgcat = cut(hh_4_g_r), at(0 `one' `two' `thr')
egen sgcat = cut(ss_4_g_r), at(0 `one' `two' `thr')
replace hgcat = int(hgcat)/33+1
replace sgcat = int(sgcat)/33+1

egen hkcat = cut(hh_4_k_r), at(0 `one' `two' `thr')
egen skcat = cut(ss_4_k_r), at(0 `one' `two' `thr')
replace hkcat = int(hkcat)/33+1
replace skcat = int(skcat)/33+1

egen h5kcat = cut(hh_5_k_r), at(0 `one' `two' `thr')
egen s5kcat = cut(ss_5_k_r), at(0 `one' `two' `thr')
replace h5kcat = int(h5kcat)/33+1
replace s5kcat = int(s5kcat)/33+1

egen hpcat = cut(hh_4_p_r), at(0 `one' `two' `thr')
egen spcat = cut(ss_4_p_r), at(0 `one' `two' `thr')
replace hpcat = int(hpcat)/33+1
replace spcat = int(spcat)/33+1

egen h5pcat = cut(hh_5_p_r), at(0 `one' `two' `thr')
egen s5pcat = cut(ss_5_p_r), at(0 `one' `two' `thr')
replace h5pcat = int(h5pcat)/33+1
replace s5pcat = int(s5pcat)/33+1


egen h1gcat = cut(hh_1_g_r), at(0 `one' `two' `thr')
egen h2gcat = cut(hh_2_g_r), at(0 `one' `two' `thr')
egen h3gcat = cut(hh_3_g_r), at(0 `one' `two' `thr')
replace h1gcat = int(h1gcat)/33+1
replace h2gcat = int(h2gcat)/33+1
replace h3gcat = int(h3gcat)/33+1

egen lfkecat = cut(lfkearn_r), at(0 `one' `two' `thr')
egen lfkwcat = cut(lfkwlth_r), at(0 `one' `two' `thr')
replace lfkecat = int(lfkecat)/33+1
replace lfkwcat = int(lfkwcat)/33+1

egen lfpecat = cut(lfpearn_r), at(0 `one' `two' `thr')
egen lfpwcat = cut(lfpwlth_r), at(0 `one' `two' `thr')
replace lfpecat = int(lfpecat)/33+1
replace lfpwcat = int(lfpwcat)/33+1

*-------------------------------
* PROGRAM THAT COMPUTES CONDVAR
*-------------------------------
cap prog drop condvar
prog def condvar

	if "`1'"=="all" {
		local ownvars i.col_g i.abi_gi i.hgcat i.sgcat
	}
	else if "`1'"=="nocg" {
		local ownvars i.abi_gi i.hgcat i.sgcat
	}
	else if "`1'"=="noag" {
		local ownvars i.col_g i.hgcat i.sgcat
	}
	else if "`1'"=="nohg" {
		local ownvars i.col_g i.abi_gi i.sgcat
	}
	else if "`1'"=="nosg" {
		local ownvars i.col_g i.abi_gi i.hgcat
	}
	if "`2'"=="k5" {
		local parvars i.col_k i.abi_ki i.h5kcat i.s5kcat i.abi_gi
	}
	else if "`2'"=="k4" {
		local parvars i.col_k i.abi_ki i.hkcat i.skcat
	}
	else if "`2'"=="p4" {
		local parvars i.col_p i.abi_pi i.hpcat i.spcat
	}

	log using "$logdir/lifetime_var.log", append
		reg lfgavge `ownvars' `parvars'
		reg lfgearn `ownvars' `parvars'
		reg lfgwlth `ownvars' `parvars'
*	list
	log close

end
*-------------------------------
log using "$logdir/lifetime_var.log", replace
log close

condvar
condvar all
condvar nocg
condvar noag
condvar nohg
condvar nosg

condvar none k5
condvar none k4
condvar none p4

*-----------------------------------------------------------------------------
* decompose h4,s4
*-----------------------------------------------------------------------------
cap prog drop condvarhs
prog def condvarhs

	if "`2'"=="k5" {
		local parvars i.col_k i.abi_ki i.h5kcat i.s5kcat i.abi_gi
	}
	else if "`2'"=="k4" {
		local parvars i.col_k i.abi_ki i.hkcat i.skcat
	}
	else if "`2'"=="p4" {
		local parvars i.col_p i.abi_pi i.hpcat i.spcat
	}

	log using "$logdir/hs4_var.log", append
		reg `1' `parvars'
	log close
end

log using "$logdir/hs4_var.log", replace
log close

condvarhs hh_4_g
condvarhs hh_4_g k5
condvarhs hh_4_g k4
condvarhs hh_4_g p4

condvarhs ss_4_g
condvarhs ss_4_g k5
condvarhs ss_4_g k4
condvarhs ss_4_g p4



*-----------------------------------------------------------------------------
* regressions
*-----------------------------------------------------------------------------
* PROGRAM THAT COMPUTES LIFETIME REGRESSIONS
*-------------------------------

cap prog drop lifetime
prog def lifetime

	local depvar `1'
	log using "$logdir/`1'_reg.log", replace

		reg `depvar' i.col_g i.abi_gi hh_4_g_r ss_4_g_r
		reg `depvar' i.col_k i.abi_ki hh_5_k_r ss_5_k_r i.abi_gi

*--------------
* categorical fixed effects
*--------------

		reg `depvar' i.col_g i.abi_gi i.hgcat i.sgcat
		reg `depvar' i.col_k i.abi_ki i.h5kcat i.s5kcat i.abi_gi

	log close

end
*-------------------------------
preserve
replace lfgearn = log(lfgearn)
replace lfgwlth = log(lfgwlth)
lifetime lfgearn
lifetime lfgwlth
restore

lifetime lfgearn_r
lifetime lfgwlth_r

*-----------------------------------------------------------------------------
* regressions with hh_4_g
*-----------------------------------------------------------------------------

cap prog drop hh4g
prog def hh4g

	local depvar `1'
	log using "$logdir/`1'_reg.log", replace

		reg `depvar' hh_3_g_r
		reg `depvar' hh_2_g_r
		reg `depvar' hh_1_g_r

		reg `depvar' hh_3_g_r hh_2_g_r
		reg `depvar' hh_3_g_r hh_2_g_r hh_1_g_r

		reg `depvar' i.col_k i.abi_ki hh_5_k_r ss_5_k_r i.abi_gi
		reg `depvar' i.col_k i.abi_ki i.h5kcat i.s5kcat i.abi_gi

* now just resources and abilities

		reg `depvar' lfkearn_r
		reg `depvar' lfkwlth_r

		reg `depvar' i.lfkecat
		reg `depvar' i.lfkwcat

* now parent states and grandparents

		reg `depvar' lfkearn_r lfpearn_r
		reg `depvar' lfkwlth_r lfpwlth_r

		reg `depvar' i.lfkecat i.lfpecat
		reg `depvar' i.lfkwcat i.lfpwcat

	log close
end
*-------------------------------
hh4g hh_4_g_r
hh4g ss_4_g_r

*-----------------------------------------------------------------------------
* fraction of inter-vivos by ability differences
*-----------------------------------------------------------------------------
scalar RR = 1/(1.04)^(6*5)
preserve
	g ivfrac = ss_4_g*RR/lfkwlth
	g lffrac = lfgearn/lfkearn
	g time = (ll_0+ll_1+ll_2)/3
	g inith = hh_1_g

	collapse ss_4_g lfkwlth lfgearn lfkearn ivfrac time inith, by(abi_gi abi_ki)
	replace ss_4_g = ss_4_g*RR/lfkwlth
	replace lfgearn = lfgearn/lfkearn
	drop lfkearn lfkwlth

	log using "$logdir/abi_rel.log", replace
	list
	log close
restore


*--------------
* r-r ige
*--------------
log using "$logdir/rrige.log", replace
	reg lfgavge_r lfkavge_r
	predict lfgkeres

	reg lfgavge_r lfpavge_r
	reg lfgavge_r lfkavge_r lfpavge_r

	reg lfgkeres lfpavge_r

	reg lfgavge_r lfkwlth_r
	predict lfgkwres

	reg lfgavge_r lfpwlth_r
	reg lfgavge_r lfkwlth_r lfpwlth_r

	reg lfgkwres lfpwlth_r
log close
