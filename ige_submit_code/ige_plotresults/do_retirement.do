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
global graphdir "/home/sylee/now/ananth/natnurt/stata_final_real/graphs"

//do do_load_data
scalar RR = 1/(1.04)^6
scalar pdv = 1+RR+RR^2+RR^3+RR^4+RR^5+RR^6

log using "$logdir/retirement.log", replace
log close

*-----------------------------------------------------------------------------
* program: with or without zeros
*-----------------------------------------------------------------------------
cap prog drop retirement
prog def retirement

	use lfkearn ss_4_g ss_`2'_k using $statadir/benchmark.dta, clear
	replace lfkearn = lfkearn/pdv/6
	summ(lfkearn)
	local lfkmean = r(mean)
	g retwlth = ss_`2'_k + ss_4_g/RR^(`2'-9)
	replace retwlth = retwlth/`lfkmean'
	if `1'==1 {
		sort lfkearn
		drop if retwlth<=0
	}

	log using "$logdir/retirement.log", append
	disp(`1')
	correlate lfkearn retwlth
	log close

*-----------
* lifetime income deciles
*-----------

	sort lfkearn
	g lfkearn_r = _n/_N*100

	g reverse_r = 100-lfkearn_r
	egen earndec = cut(reverse_r), at(0 10 20 30 40 50 60 70 80 90 100)
	drop reverse_r

	replace earndec = int(earndec/10)
	replace earndec = 10-earndec

* indiv saving rate
	g saving = retwlth*RR^(`2'-5)/lfkearn

	summ saving
	local scv = r(sd)/r(mean)
	log using "$logdir/retirement.log", append
	disp(`1')
	disp(`scv')
	log close


*-----------
* total saving rate
*-----------
	preserve
		collapse ///
			(mean) lfmean=lfkearn ssmean=retwlth ratemean=saving ///
			(p50) lfmed=lfkearn ssmed=retwlth ratemed=saving

		g tmean = ssmean*RR^(`2'-5)/lfmean
		g tmed = ssmed*RR^(`2'-5)/lfmed

		save aggstats.dta, replace

	restore

*-----------
* dec saving rate and gini's
*-----------
	preserve
		ineqdec0 retwlth, bygroup(earndec)
		matrix gini = J(10,1,.)
		forval dec = 1/10{
			matrix gini[`dec',1] = r(gini_`dec')
		}
		matrix gini = gini \ r(gini)

		collapse ///
			(mean) lfmean=lfkearn ssmean=retwlth ratemean=saving ///
			(p50) lfmed=lfkearn ssmed=retwlth ratemed=saving ///
			, by(earndec)

		g tmean = ssmean*RR^(`2'-5)/lfmean
		g tmed = ssmed*RR^(`2'-5)/lfmed

		append using aggstats
		rm aggstats.dta

		svmat gini
		rename gini1 gini

		save "$statadir/retirement`1'.dta", replace

		log using "$logdir/retirement.log", append
		disp(`1')
		list if missing(earndec)
		summ gini if _n<11
		log close

		tw connected gini earndec, ///
			l1title("Gini") xtitle("Lifetime Earnings Decile") ///
			ytitle("") yscale(range(0 1)) ylabel(0(0.1)1,grid) xlabel(1(1)10,grid)
		graph export "$graphdir/gini`1'.eps", replace

		tw connected ssmean ssmed earndec, ///
			l1title("Mean/Median Wealth") xtitle("Lifetime Earnings Decile") ///
			ytitle("") ylabel(,grid) xlabel(1(1)10,grid) legend(order (1 "Means" 2 "Medians"))
		graph export "$graphdir/meanmed`1'.eps", replace

	restore

end
*-----------------------------------------------------------------------------
*-----------------------------------------------------------------------------

retirement 0 11
retirement 1 11
