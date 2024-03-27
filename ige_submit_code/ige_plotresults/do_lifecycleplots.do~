global outputdir "/home/sylee/now/ananth/natnurt/final_real/results"
global graphdir "/home/sylee/now/ananth/natnurt/stata_final_real/graphs"

*-----------------------------------------------------------------------------

scalar colshare = 0.4798591
scalar eprem = 1.519012

scalar delk = 0.067
scalar alph = 0.33
scalar ces = 1-1/1.441

scalar r0 = 0.04
scalar ww = 0.71711392481767999
scalar tfp = 1.6737677842781544

scalar Windx = (1+r0+delk)^6-1
scalar Windx = (1-alph)*(alph/Windx)^(alph/(1-alph))
scalar Windx = Windx*tfp

scalar ups = colshare/(1-colshare) *eprem
scalar ups = ww*ces *ups^(1-ces)
scalar ups = 1/(1+ups)

scalar wskill1 = (1-ups)^(1/(1-ces)) *ww^(ces/(ces-1))
scalar wskill1 = ups^(1/(1-ces)) + wskill1
scalar wskill1 = Windx *wskill1^((1-ces)/ces)
scalar wskill2 = wskill1 *ww

*-----------------------------------------------------------------------------

import delimited "$outputdir/hhidgen2.out", delimiter(" ",collapse) clear
keep col ee* hh* ll*

forval j = 4/9{
	g nn_`j' = ee_`j'/hh_`j'
	replace nn_`j' = nn_`j'/wskill1 if col==1
	replace nn_`j' = nn_`j'/wskill2 if col==2

	replace nn_`j' = 1-nn_`j'
}

forval j = 5/7{
	local k = `j'-5
	rename ll_`k' ll_`j'
}

forval j = 5/7{
	replace nn_`j' = nn_`j'-ll_`j'
}
g nn_10 = 0


collapse ee* hh* nn* ll*, by(col)
reshape long ee_ hh_ nn_ ll_, i(col) j(age)
keep if age>3

scalar norm = hh[1]
replace hh = hh/norm

replace age = age*6+3

tw (connected hh age if col==1) (connected hh age if col==2), xlabel(24(6)66) xtitle("Age") ytitle("Human Capital") legend(cols(2) order (1 "High School" 2 "College"))
graph export "$graphdir/hc_by_col.eps", replace

tw (connected nn ll age if col==1, mstyle(p1 p1) lstyle(p1 p1)) (connected nn ll age if col==2, mstyle(p2 p2) lstyle(p2 p2)), xlabel(24(6)66) xtitle("Age") ytitle("Own Time Investments") legend(cols(2) order (1 "High School" 3 "College"))
graph export "$graphdir/nnll_by_col.eps", replace
