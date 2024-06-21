version 18 // needed for the xthdidregress function

ssc install did_imputation, replace
ssc install binscatter, replace
ssc install estout, replace
net install simulate2, from("https://janditzen.github.io/simulate2/") replace
net install parallel, from("https://raw.github.com/gvegayon/parallel/master/") replace
net install did2s, replace from("https://raw.githubusercontent.com/kylebutts/did2s_stata/main/ado/")

cd ""

gl num_replications = 1000

* Creating data
{
clear
set seed 123
set obs 20000
g id = mod(_n, 1000)+1
sort id
g time = mod(_n, 20)+1
xtset id time 

preserve
keep if time == 1
sort id
g cohort = runiformint(1, 15)
replace cohort = 0 if cohort<10
tempfile to_merge
save `to_merge'
restore

merge m:1 id using `to_merge', nogen

fvset base 0 cohort
sort id time

g fe_unit = id/1000
g fe_time = time/20

g exposure_time = time - cohort if cohort > 0
su exposure_time
di r(min)
g exposure_time_adj = exposure_time - r(min)
replace exposure_time     = 999 if cohort == 0
replace exposure_time_adj = 999 if cohort == 0
local base = -r(min)-1
di `base'
fvset base `base' exposure_time_adj

g treated = cohort>0
g post_   = cond(!cohort, time >= 15, time >= cohort)
g treat   = cond(!cohort, 0, time >= cohort)

gl effect = 5

g treat_effect_1 = $effect * treat // constant, homogeneous

su cohort if treat
g treat_effect_2 =  ($effect + cohort - r(mean)) * treat // constant, heterogenous

su exposure_time if treat
g treat_effect_3 = ($effect + exposure_time - r(mean)) * treat // dynamic, homogeneous
g treat_effect_5 = ($effect - exposure_time + r(mean)) * treat // dynamic, homogeneous

g cohorttime = (cohort) + ((cohort/5-2) * exposure_time)
su cohorttime if treat
g treat_effect_4 = ($effect + cohorttime - r(mean)) * treat // dynamic, heterogenous
g treat_effect_6 = ($effect - cohorttime + r(mean)) * treat // dynamic, heterogenous

su treat_effect_* if treat

g epsilon = rnormal()

g x = rnormal()
egen x_avg = mean(x), by(cohort)
g x_demean = x - x_avg

forval i = 1/6 {
g y_`i' = x + fe_unit + fe_time + treat_effect_`i' + epsilon
}

save reg_data, replace
}

* True ATT Plots
{
use reg_data, replace
collapse treat_effect_* if treated, by(cohort time)

g label = "True Treatment Effect"
g scatter = 5 if time == 5

* graphs
{
twoway ///
(line treat_effect_1 time if cohort == 10) ///
(line treat_effect_1 time if cohort == 11) ///
(line treat_effect_1 time if cohort == 12) ///
(line treat_effect_1 time if cohort == 13) ///
(line treat_effect_1 time if cohort == 14) ///
(line treat_effect_1 time if cohort == 15) ///
(scatter scatter time, mcol(none) mlabel(label) mlabpos(12) mlabc(black) mlabs(small)) ///
, subti("1: Homogeneous, Static") xti("") yti("Treatment Effect") xlab(none) ///
yline(5, lp(dash) lcol(maroon)) ysc(r(0 15)) ylab(0(5)15, nogrid) ///
legend(off) plotregion(lcol(black)) name(y1, replace)

twoway ///
(line treat_effect_2 time if cohort == 10) ///
(line treat_effect_2 time if cohort == 11) ///
(line treat_effect_2 time if cohort == 12) ///
(line treat_effect_2 time if cohort == 13) ///
(line treat_effect_2 time if cohort == 14) ///
(line treat_effect_2 time if cohort == 15) ///
(scatter scatter time, mcol(none) mlabel(label) mlabpos(12) mlabc(black) mlabs(small)) ///
, subti("2: Heterogeneous, Static") xti("") yti("") ysca(alt) xlab(none) ///
yline(5, lp(dash) lcol(maroon)) ysc(r(0 15)) ylab(0(5)15, nogrid) ///
legend(off) plotregion(lcol(black)) name(y2, replace)

twoway ///
(line treat_effect_3 time if cohort == 10) ///
(line treat_effect_3 time if cohort == 11) ///
(line treat_effect_3 time if cohort == 12) ///
(line treat_effect_3 time if cohort == 13) ///
(line treat_effect_3 time if cohort == 14) ///
(line treat_effect_3 time if cohort == 15) ///
(scatter scatter time, mcol(none) mlabel(label) mlabpos(12) mlabc(black) mlabs(small)) ///
, subti("3: Homogeneous, Dynamic") xti("") yti("Treatment Effect") xlab(none) ///
yline(5, lp(dash) lcol(maroon)) ysc(r(0 15)) ylab(0(5)15, nogrid) ///
legend(off) plotregion(lcol(black)) name(y3, replace)

twoway ///
(line treat_effect_4 time if cohort == 10) ///
(line treat_effect_4 time if cohort == 11) ///
(line treat_effect_4 time if cohort == 12) ///
(line treat_effect_4 time if cohort == 13) ///
(line treat_effect_4 time if cohort == 14) ///
(line treat_effect_4 time if cohort == 15) ///
(scatter scatter time, mcol(none) mlabel(label) mlabpos(12) mlabc(black) mlabs(small)) ///
,subti("4: Heterogeneous, Dynamic") xti("") yti("") ysca(alt) xlab(none) ///
yline(5, lp(dash) lcol(maroon)) ysc(r(0 15)) ylab(0(5)15, nogrid) ///
legend(off) plotregion(lcol(black)) name(y4, replace)

twoway ///
(line treat_effect_5 time if cohort == 10) ///
(line treat_effect_5 time if cohort == 11) ///
(line treat_effect_5 time if cohort == 12) ///
(line treat_effect_5 time if cohort == 13) ///
(line treat_effect_5 time if cohort == 14) ///
(line treat_effect_5 time if cohort == 15) ///
(scatter scatter time, mcol(none) mlabel(label) mlabpos(12) mlabc(black) mlabs(small)) ///
,subti("5: Homogeneous, Dynamic (Decreasing)", size(medsmall)) xti("Time") yti("Treatment Effect") ///
yline(5, lp(dash) lcol(maroon)) ysc(r(-2 15)) ylab(0(5)15, nogrid) xlab(, nogrid) ///
legend(off) plotregion(lcol(black)) name(y5, replace)

twoway ///
(line treat_effect_6 time if cohort == 10) ///
(line treat_effect_6 time if cohort == 11) ///
(line treat_effect_6 time if cohort == 12) ///
(line treat_effect_6 time if cohort == 13) ///
(line treat_effect_6 time if cohort == 14) ///
(line treat_effect_6 time if cohort == 15) ///
(scatter scatter time, mcol(none) mlabel(label) mlabpos(12) mlabc(black) mlabs(small)) ///
,subti("6: Heterogeneous, Dynamic (Decreasing)", size(medsmall)) xti("Time") yti("") ysca(alt) ///
yline(5, lp(dash) lcol(maroon)) ysc(r(-2 15)) ylab(0(5)15, nogrid) xlab(, nogrid) ///
legend(off) plotregion(lcol(black)) name(y6, replace)
}

graph combine y1 y2 y3 y4 y5 y6, xcommon imargin(zero) rows(3) graphregion(margin(zero)) ysize(8)

graph export "Table 1 Treatment Effects Chart.png", replace width(6000)
graph export "Table 1 Treatment Effects Chart.jpg", replace width(2000)
}

* Comparing performance of various estimators
{
cap program drop run_regs_fast
program define run_regs_fast, rclass

* Creating data
{
clear
set obs 20000
g id = mod(_n, 1000)+1
sort id
g time = mod(_n, 20)+1
xtset id time 

preserve
keep if time == 1
sort id
g cohort = runiformint(1, 15)
replace cohort = 0 if cohort<10
tempfile to_merge
save `to_merge'
restore

merge m:1 id using `to_merge', nogen

fvset base 0 cohort
sort id time

g fe_unit = id/1000
g fe_time = time/20

g exposure_time = time - cohort if cohort > 0
su exposure_time
di r(min)
g exposure_time_adj = exposure_time - r(min)
replace exposure_time     = 999 if cohort == 0
replace exposure_time_adj = 999 if cohort == 0
local base = -r(min)-1
di `base'
fvset base `base' exposure_time_adj

g treated = cohort>0
g post_   = cond(!cohort, time >= 15, time >= cohort)
g treat   = cond(!cohort, 0, time >= cohort)

gl effect = 5

g treat_effect_1 = $effect * treat // constant, homogeneous

su cohort if treat
g treat_effect_2 =  ($effect + cohort - r(mean)) * treat // constant, heterogenous

su exposure_time if treat
g treat_effect_3 = ($effect + exposure_time - r(mean)) * treat // dynamic, homogeneous
g treat_effect_5 = ($effect - exposure_time + r(mean)) * treat // dynamic, homogeneous

g cohorttime = (cohort) + ((cohort/5-2) * exposure_time)
su cohorttime if treat
g treat_effect_4 = ($effect + cohorttime - r(mean)) * treat // dynamic, heterogenous
g treat_effect_6 = ($effect - cohorttime + r(mean)) * treat // dynamic, heterogenous

g epsilon = rnormal()

g x = rnormal()
egen x_avg = mean(x), by(cohort)
g x_demean = x - x_avg

forval i = 1/6 {
g y_`i' = x + fe_unit + fe_time + treat_effect_`i' + epsilon
}

keep id time cohort exposure_time exposure_time_adj treat treated x x_demean y_1 y_2 y_3 y_4 y_5 y_6

cap drop sa_*
forval c = 10/15 {
forval e = 0/24 {
if `e' != 13 {
g sa_`c'_`e' = cohort == `c' & exposure_time_adj == `e'
}
}
}
}

forval i = 1/6 {
* TWFE
reghdfe y_`i' treat, a(id time)
return scalar  b_twfe_non_`i' = _b[treat]

* TWFE with control
reghdfe y_`i' treat x, a(id time)
return scalar  b_twfe_ctrl_`i' = _b[treat]

* Sun/Abr
{
cap drop sa_effect
reghdfe y_`i' sa_*, a(id time) 

g sa_effect=.
forval c = 10/15 {
forval e = 0/24 {
if `e' != 13 {
replace sa_effect = _b[sa_`c'_`e'] if cohort == `c' & exposure_time_adj == `e' & treat
}
}
}
su sa_effect

return scalar b_sa_non_`i' = r(mean)
}

* Wooldr with control
{
reghdfe y_`i' 1.treat#i.cohort#i.time 1.treat#i.cohort#i.time#c.x_demean, a(i.time##c.x i.cohort##c.x) technique(lsqr)

margins, dydx(treat) subpop(if treat) nose

return scalar  b_wooldr_ctrl_`i' = r(b)[1,1]
}

* Gardner
{
did2s y_`i', first_stage(i.time) second_stage(treat) treatment(treat ) unit(id) cluster(id)

return scalar  b_grdnr_non_`i' = e(b)[1,1]
}

* Gardner with control
{
did2s y_`i', first_stage(i.time x) second_stage(treat) treatment(treat ) unit(id) cluster(id)
return scalar  b_grdnr_ctrl_`i' = e(b)[1,1]
}

preserve
replace cohort = . if cohort == 0
* Borusyak et al.
{
did_imputation y_`i' id time cohort, autosample nose
return scalar  b_borys_non_`i' = e(b)[1, 1]
}

* Borusyak et al. control
{
did_imputation y_`i' id time cohort, controls(x) autosample nose
return scalar  b_borys_ctrl_`i' = e(b)[1, 1]
}
restore

* Callaway/Sant'Anna
{
xthdidregress ra (y_`i') (treat), group(id) vce(cluster id)
estat aggregation
return scalar  b_callwy_non_`i' = r(b)[1, 1]
}

* Callaway/Sant'Anna control
{
xthdidregress ra (y_`i' x) (treat), group(id) vce(cluster id)
estat aggregation
return scalar  b_callwy_ctrl_`i' = r(b)[1, 1]
}
}

eret clear
end
	
psimulate2, seed(123) r(${num_replications}) p(4) saving(monte_carlo_table_fast): run_regs_fast

* Prepare table
{
use monte_carlo_table_fast, replace

g i = _n
reshape long b_borys_ctrl_ b_borys_non_ b_grdnr_ctrl_ b_grdnr_non_ b_wooldr_ctrl_ b_sa_non_ b_twfe_ctrl_ b_twfe_non_ b_callwy_non_ b_callwy_ctrl_, i(i) j(data)

reshape long b_, i(i data) j(name) s

collapse (mean) b_* (count) count_obs = i (sd) sd_=b_, by(name data)
g se_ = sd_/sqrt(count_obs)
drop count_obs

reshape wide b_ sd_ se_, i(name) j(data)

replace name = "Borusyak et al., Control" if name == "borys_ctrl_"
replace name = "Borusyak et al., No Control" if name == "borys_non_"
replace name = "Gardner, Control" if name == "grdnr_ctrl_"
replace name = "Gardner, No Control" if name == "grdnr_non_"
replace name = "Wooldridge, Control" if name == "wooldr_ctrl_"
replace name = "Sun/Abraham, No Control" if name == "sa_non_"
replace name = "TWFE, Control" if name == "twfe_ctrl_"
replace name = "TWFE, No Control" if name == "twfe_non_"
replace name = "Callaway/Sant'Anna, No Control" if name == "callwy_non_"
replace name = "Callaway/Sant'Anna, Control" if name == "callwy_ctrl_"

g       sort = 1  if name == "TWFE, No Control"
replace sort = 2  if name == "Sun/Abraham, No Control"
replace sort = 3  if name == "Gardner, No Control"
replace sort = 4  if name == "Borusyak et al., No Control"
replace sort = 5  if name == "Callaway/Sant'Anna, No Control"
replace sort = 6  if name == "TWFE, Control"
replace sort = 7  if name == "Wooldridge, Control"
replace sort = 8  if name == "Gardner, Control"
replace sort = 9  if name == "Borusyak et al., Control"
replace sort = 10 if name == "Callaway/Sant'Anna, Control"

forval i = 1/6 {
	tostring  b_`i', g( b_format_`i') force format(%05.3f)
	tostring sd_`i', g(sd_format_`i') force format(%05.3f)
	replace sd_format_`i' = "(" + sd_format_`i' + ")"
	
	* statistical significance
	g t_`i'=(b_`i' - ${effect})/sd_`i'
	g stars_`i' = cond(abs(t_`i')>2.326, "***", cond(abs(t_`i')>1.64, "**", cond(abs(t_`i')>1.282, "*", "")))
	replace b_format_`i' = b_format_`i' + stars_`i'
}

order name b_format_* b_1 b_2 b_3 b_4 b_5 b_6 sd_format_* sd_* se_* t_* stars_*

* SDs below estimates
expand 2
sort sort
drop sort
forval i = 1/6 {
	replace b_format_`i' = sd_format_`i' if !mod(_n, 2)
}

export excel "Table 1.xlsx", sh("Stata Output", replace) firstrow(var)
}

* Histogram chart
{
use monte_carlo_table_fast, replace

gl graph_list ""

foreach i in twfe_non sa_non grdnr_non borys_non callwy_non twfe_ctrl wooldr_ctrl grdnr_ctrl borys_ctrl callwy_ctrl {
forval j = 1/6 {
	* labels
	{
	if "`i'" == "twfe_non" {
	local title "TWFE, No Control"
	}
	if "`i'" == "twfe_ctrl" {
	local title "TWFE, Control"
	}
	if "`i'" == "sa_non" {
	local title "Sun/Abraham, No Control"
	}
	if "`i'" == "wooldr_ctrl" {
	local title "Wooldridge, Control"
	}
	if "`i'" == "borys_ctrl" {
	local title "Borusyak et al., Control"
	}
	if "`i'" == "borys_non" {
	local title "Borusyak et al., No Control"
	}
	if "`i'" == "grdnr_ctrl" {
	local title "Gardner, Control"
	}
	if "`i'" == "grdnr_non" {
	local title "Gardner, No Control"
	}
	if "`i'" == "callwy_non" {
	local title "Callaway/Sant'Anna, No Control"
	}
	if "`i'" == "callwy_ctrl" {
	local title "Callaway/Sant'Anna, Control"
	}
	}
	
	twoway hist b_`i'_`j' if b_`i'_`j'>4.25, w(0.01) xsc(r(4 6)) xlab(4(0.5)6, nogrid) xline(5, lcol(maroon) lp(dash)) ylab(, nogrid) percent xti("") plotregion(lcol(black)) name(`i'_`j', replace) ti("`title', `j'", size(medsmall)) yti("")	ylab(none)
	
	gl graph_list "${graph_list} `i'_`j'"
} 
}

graph combine ${graph_list}, xcommon rows(10) name(combinedhist, replace) imargin(0 0 0 0)  graphregion(margin(zero)) ysize(7.5) note("{bf: Notes}: n = ${num_replications} replications. Datasets have 20,000 observations. True ATT is 5, represented by the dashed red line." "TWFE, Sun/Abraham, and Wooldridge estimates use the reghdfe package. Gardner estimates use the did2s package." "Borusyak et al. estimates use the did_imputation package. Monte Carlo simulations run with the simulate2 package.", size(tiny))
graph export "Histograms of Monte Carlo Estimates 1.png", replace width(6000)
graph export "Histograms of Monte Carlo Estimates 1.jpg", replace width(2000)
} 
}

* Point equivalence of estimators
{
use monte_carlo_table_fast, replace

* Sun/Abraham and Callaway/Sant'Anna (No Control)
forval i = 1/6 {
	g d_callwy_sa_`i' = abs(b_callwy_non_`i' - b_sa_non_`i')
	qui su d_callwy_sa_`i'
	di r(max)
}	

* Gardner and Borusyak et al.
forval i = 1/6 {
	g d_borys_gardner_non_`i'  = abs(b_borys_non_`i'  - b_grdnr_non_`i')
	qui su d_borys_gardner_non_`i'
	di r(max)
	g d_borys_gardner_ctrl_`i' = abs(b_borys_ctrl_`i' - b_grdnr_ctrl_`i')
	qui su d_borys_gardner_ctrl_`i'
	di r(max)
}

* Wooldridge (ctrl)
forval i = 1/6 {
	g d_borys_wooldr_ctrl_`i' = abs(b_borys_ctrl_`i' - b_wooldr_ctrl_`i')
	qui su d_borys_wooldr_ctrl_`i'
	di r(max)
}
}

* Comparing performance of various estimators, not-yet-treated group 
{
cap program drop run_regs2_fast
program define run_regs2_fast, rclass

* Creating data
{
clear
set obs 71500
g id = mod(_n, 3575)+1
sort id
g time = mod(_n, 20)+1
xtset id time 

preserve
keep if time == 1
sort id
g cohort = runiformint(1, 15)
replace cohort = 0 if cohort<10
tempfile to_merge
save `to_merge'
restore

merge m:1 id using `to_merge', nogen

fvset base 0 cohort
sort id time

g fe_unit = id/1000
g fe_time = time/20

g exposure_time = time - cohort if cohort > 0
su exposure_time
di r(min)
g exposure_time_adj = exposure_time - r(min)
replace exposure_time     = 999 if cohort == 0
replace exposure_time_adj = 999 if cohort == 0
local base = -r(min)-1
di `base'
fvset base `base' exposure_time_adj

g treated = cohort>0
g post_   = cond(!cohort, time >= 15, time >= cohort)
g treat   = cond(!cohort, 0, time >= cohort)

keep if treated & time < 15

gl effect = 5

g treat_effect_1 = $effect * treat // constant, homogeneous

su cohort if treat
g treat_effect_2 =  ($effect + cohort - r(mean)) * treat // constant, heterogenous

su exposure_time if treat
g treat_effect_3 = ($effect + exposure_time - r(mean)) * treat // dynamic, homogeneous
g treat_effect_5 = ($effect - exposure_time + r(mean)) * treat // dynamic, homogeneous

g cohorttime = (cohort) + ((cohort/5-2) * exposure_time)
su cohorttime if treat
g treat_effect_4 = ($effect + cohorttime - r(mean)) * treat // dynamic, heterogenous
g treat_effect_6 = ($effect - cohorttime + r(mean)) * treat // dynamic, heterogenous

g epsilon = rnormal()

g x = rnormal()
egen x_avg = mean(x), by(cohort)
g x_demean = x - x_avg

forval i = 1/6 {
g y_`i' = x + fe_unit + fe_time + treat_effect_`i' + epsilon
}

keep id time cohort exposure_time exposure_time_adj treat treated x x_demean y_1 y_2 y_3 y_4 y_5 y_6 treat_effect_*

cap drop sa_*
forval c = 10/15 {
forval e = 0/24 {
if `e' != 13 {
g sa_`c'_`e' = cohort == `c' & exposure_time_adj == `e'
}
}
}
}

return scalar n = _N

forval i = 1/6 {

* TWFE
reghdfe y_`i' treat, a(id time)
return scalar  b_twfe_non_`i' = _b[treat]

* TWFE with control
reghdfe y_`i' treat x, a(id time)
return scalar  b_twfe_ctrl_`i' = _b[treat]

* Sun/Abr
{
cap drop sa_effect
reghdfe y_`i' sa_*, a(id time) 

g sa_effect=.
forval c = 10/15 {
forval e = 0/24 {
if `e' != 13 {
replace sa_effect = _b[sa_`c'_`e'] if cohort == `c' & exposure_time_adj == `e' & treat
}
}
}
su sa_effect

return scalar b_sa_non_`i' = r(mean)
}

* Wooldr with control
{
reghdfe y_`i' 1.treat#i.cohort#i.time 1.treat#i.cohort#i.time#c.x_demean, a(i.time##c.x i.cohort##c.x) technique(lsqr)

margins, dydx(treat) subpop(if treat) nose

return scalar  b_wooldr_ctrl_`i' = r(b)[1,1]
}

* Gardner
{
did2s y_`i', first_stage(i.time) second_stage(treat) treatment(treat ) unit(id) cluster(id)

return scalar  b_grdnr_non_`i' = e(b)[1,1]
}

* Gardner with control
{
did2s y_`i', first_stage(i.time x) second_stage(treat) treatment(treat ) unit(id) cluster(id)
return scalar  b_grdnr_ctrl_`i' = e(b)[1,1]
}

preserve
replace cohort = . if cohort == 0
* Borusyak et al.
{
did_imputation y_`i' id time cohort, autosample nose
return scalar  b_borys_non_`i' = e(b)[1, 1]
}

* Borusyak et al. control
{
did_imputation y_`i' id time cohort, controls(x) autosample nose
return scalar  b_borys_ctrl_`i' = e(b)[1, 1]
}
restore

* Callaway/Sant'Anna
{
xthdidregress ra (y_`i') (treat), group(id) vce(cluster id) controlgroup(notyet)
estat aggregation
return scalar  b_callwy_non_`i' = r(b)[1, 1]
}

* Callaway/Sant'Anna control
{
xthdidregress ra (y_`i' x) (treat), group(id) vce(cluster id) controlgroup(notyet)
estat aggregation
return scalar  b_callwy_ctrl_`i' = r(b)[1, 1]
}
}
eret clear
end

psimulate2, seed(123) r(${num_replications}) p(4) saving(monte_carlo_table2_fast): run_regs2_fast

* Prepare table
{
use monte_carlo_table2_fast, replace

su n

g i = _n
reshape long b_borys_ctrl_ b_borys_non_ b_grdnr_ctrl_ b_grdnr_non_ b_wooldr_ctrl_ b_sa_non_ b_twfe_ctrl_ b_twfe_non_ b_callwy_non_ b_callwy_ctrl_, i(i) j(data)

reshape long b_, i(i data) j(name) s

collapse (mean) b_* (count) count_obs = i (sd) sd_=b_, by(name data)
g se_ = sd_/sqrt(count_obs)
drop count_obs

reshape wide b_ sd_ se_, i(name) j(data)

replace name = "Borusyak et al., Control" if name == "borys_ctrl_"
replace name = "Borusyak et al., No Control" if name == "borys_non_"
replace name = "Gardner, Control" if name == "grdnr_ctrl_"
replace name = "Gardner, No Control" if name == "grdnr_non_"
replace name = "Wooldridge, Control" if name == "wooldr_ctrl_"
replace name = "Sun/Abraham, No Control" if name == "sa_non_"
replace name = "TWFE, Control" if name == "twfe_ctrl_"
replace name = "TWFE, No Control" if name == "twfe_non_"
replace name = "Callaway/Sant'Anna, No Control" if name == "callwy_non_"
replace name = "Callaway/Sant'Anna, Control" if name == "callwy_ctrl_"

g       sort = 1  if name == "TWFE, No Control"
replace sort = 2  if name == "Sun/Abraham, No Control"
replace sort = 3  if name == "Gardner, No Control"
replace sort = 4  if name == "Borusyak et al., No Control"
replace sort = 5  if name == "Callaway/Sant'Anna, No Control"
replace sort = 6  if name == "TWFE, Control"
replace sort = 7  if name == "Wooldridge, Control"
replace sort = 8  if name == "Gardner, Control"
replace sort = 9  if name == "Borusyak et al., Control"
replace sort = 10 if name == "Callaway/Sant'Anna, Control"

forval i = 1/6 {
	tostring  b_`i', g( b_format_`i') force format(%05.3f)
	tostring sd_`i', g(sd_format_`i') force format(%05.3f)
	replace sd_format_`i' = "(" + sd_format_`i' + ")"
	
	* statistical significance
	g t_`i'=(b_`i' - ${effect})/sd_`i'
	g stars_`i' = cond(abs(t_`i')>2.326, "***", cond(abs(t_`i')>1.64, "**", cond(abs(t_`i')>1.282, "*", "")))
	replace b_format_`i' = b_format_`i' + stars_`i'
}

order name b_format_* b_1 b_2 b_3 b_4 b_5 b_6 sd_format_* sd_* se_* t_* stars_*

* SDs below estimates
expand 2
sort sort
drop sort
forval i = 1/6 {
	replace b_format_`i' = sd_format_`i' if !mod(_n, 2)
}


* true effect row
local n = _N+1
set obs `n'
replace name = "True Average Effect" if mi(name)
g sort = cond(name != "True Average Effect", _n, 0)
sort sort
drop sort
forval i = 1/6 {
	replace b_format_`i' = "$effect" if name == "True Average Effect"
}

export excel "Table 2.xlsx", sh("Stata Output", replace) firstrow(var)
}

* Histogram chart
{
use monte_carlo_table2_fast, replace

gl graph_list ""

foreach i in twfe_non sa_non grdnr_non borys_non callwy_non twfe_ctrl wooldr_ctrl grdnr_ctrl borys_ctrl callwy_ctrl {
forval j = 1/6 {
	* labels
	{
	if "`i'" == "twfe_non" {
	local title "TWFE, No Control"
	}
	if "`i'" == "twfe_ctrl" {
	local title "TWFE, Control"
	}
	if "`i'" == "sa_non" {
	local title "Sun/Abraham, No Control"
	}
	if "`i'" == "wooldr_ctrl" {
	local title "Wooldridge, Control"
	}
	if "`i'" == "borys_ctrl" {
	local title "Borusyak et al., Control"
	}
	if "`i'" == "borys_non" {
	local title "Borusyak et al., No Control"
	}
	if "`i'" == "grdnr_ctrl" {
	local title "Gardner, Control"
	}
	if "`i'" == "grdnr_non" {
	local title "Gardner, No Control"
	}
	if "`i'" == "callwy_non" {
	local title "Callaway/Sant'Anna, No Control"
	}
	if "`i'" == "callwy_ctrl" {
	local title "Callaway/Sant'Anna, Control"
	}
	}
	
	twoway hist b_`i'_`j' , w(0.01) xsc(r(4 6)) xlab(4(0.5)6, nogrid) xline(5, lcol(maroon) lp(dash)) ylab(, nogrid) percent xti("") plotregion(lcol(black)) name(`i'_`j', replace) ti("`title', `j'", size(medsmall)) yti("")	ylab(none)
	
	gl graph_list "${graph_list} `i'_`j'"
} 
}

graph combine ${graph_list}, xcommon rows(10) name(combinedhist, replace) imargin(zero)  graphregion(margin(zero)) ysize(7.5) note("{bf: Notes}: n = ${num_replications} replications. TWFE, Sun/Abraham, and Wooldridge estimates use the reghdfe package. Gardner estimates use the did2s package." "Borusyak et al. estimates use the did_imputation package. Monte Carlo simulations run with the simulate2 package.", size(tiny))
graph export "Histograms of Monte Carlo Estimates 2.png", replace width(6000)
graph export "Histograms of Monte Carlo Estimates 2.jpg", replace width(2000)
} 
}

* TWFE bias by % of never-treated
{
cap program drop sim_unit
program sim_unit, rclass
clear
set obs 20000
g id = mod(_n, 1000)+1
sort id
g time = mod(_n, 20)+1

gl pc_never_treated = runiform(0, 0.9)

preserve
keep if time == 1
sort id
g prob = runiform()
g cohort = 0 if prob < $pc_never_treated
replace cohort = runiformint(10, 15) if mi(cohort)
tempfile to_merge
save `to_merge'
restore

merge m:1 id using `to_merge', nogen

g treat   = cond(!cohort, 0, time >= cohort)

count if cohort == 0
local num = r(N)
count if !treat
di 100*`num'/r(N)
return scalar pc_never_treated = 100*`num'/r(N)

g fe_unit = id/1000
g fe_time = time/20

g exposure_time = time - cohort if cohort > 0
replace exposure_time     = 999 if cohort == 0


gl effect = 5
g treat_effect_1 = $effect * treat // constant, homogeneous
su cohort if treat
g treat_effect_2 =  ($effect + cohort - r(mean)) * treat // constant, heterogenous
su exposure_time if treat
g treat_effect_3 = ($effect + exposure_time - r(mean)) * treat // dynamic, homogeneous
g treat_effect_5 = ($effect - exposure_time + r(mean)) * treat // dynamic, homogeneous

g cohorttime = (cohort) + ((cohort/5-2) * exposure_time)
su cohorttime if treat
g treat_effect_4 = ($effect + cohorttime - r(mean)) * treat // dynamic, heterogenous
g treat_effect_6 = ($effect - cohorttime + r(mean)) * treat // dynamic, heterogenous
g epsilon = rnormal()

forval i = 1/6 {
g y_`i' = fe_unit + fe_time + treat_effect_`i' + epsilon
reghdfe y_`i' treat, a(id time)
return scalar coef_`i' = _b[treat]
return scalar bias_`i' = _b[treat]-$effect
}
eret clear
end

gl num_replications = 10000
psimulate2, seed(123) r(${num_replications}) p(4) saving(twfe_bias_by_pc_never_treated): sim_unit

use twfe_bias_by_pc_never_treated, clear

binscatter coef_1 coef_2 coef_3 coef_4 coef_5 coef_6 pc_never_treated, line(connect) yti("Estimated Treatment Effect") xti("% Never Treated") legend(order(1 "1: Homogeneous, Static" 2 "2: Heterogeneous, Static" 3 "3: Homogeneous, Dynamic" 4 "4: Heterogeneous, Dynamic" 5 "5: Homogeneous, Dynamic (Decreasing)" 6 "6: Heterogeneous, Dynamic (Decreasing)") pos(6) rows(3)) plotregion(lcol(black)) yline(5, lcol(maroon) lp(dash)) note("{bf: Notes}: n = ${num_replications} replications. TWFE estimates use the reghdfe package. Monte Carlo simulations run with the simulate2 package." "% Never Treated = (Count Never-Treated) / (Count All Untreated)." "This chart made with the binscatter package.") colors(stc1 stc2 stc3 stc4 stc5 stc6)

graph export "TWFE Bias by Percent Never Treated.png", replace width(6000)
graph export "TWFE Bias by Percent Never Treated.jpg", replace width(2000)

twoway ///
(scatter coef_1 pc_never_treated, msize(tiny) msymbol(smcircle_hollow) mcol(%50)) ///
(scatter coef_2 pc_never_treated, msize(tiny) msymbol(smcircle_hollow) mcol(%50)) ///
(scatter coef_3 pc_never_treated, msize(tiny) msymbol(smcircle_hollow) mcol(%50)) ///
(scatter coef_4 pc_never_treated, msize(tiny) msymbol(smcircle_hollow) mcol(%50)) ///
(scatter coef_5 pc_never_treated, msize(tiny) msymbol(smcircle_hollow) mcol(%50)) ///
(scatter coef_6 pc_never_treated, msize(tiny) msymbol(smcircle_hollow) mcol(%50)) ///
(scatteri . ., mcol(stc1)) ///
(scatteri . ., mcol(stc2)) ///
(scatteri . ., mcol(stc3)) ///
(scatteri . ., mcol(stc4)) ///
(scatteri . ., mcol(stc5)) ///
(scatteri . ., mcol(stc6)) ///
, yti("Estimated Treatment Effect") xti("% Never Treated") legend(order(7 "1: Homogeneous, Static" 8 "2: Heterogeneous, Static" 9 "3: Homogeneous, Dynamic" 10 "4: Heterogeneous, Dynamic" 11 "5: Homogeneous, Dynamic (Decreasing)" 12 "6: Heterogeneous, Dynamic (Decreasing)") pos(6) rows(3)) plotregion(lcol(black)) yline(5, lcol(maroon) lp(dash)) note("{bf: Notes}: n = ${num_replications} replications. TWFE estimates use the reghdfe package. Monte Carlo simulations run with the simulate2 package." "% Never Treated = (Count Never-Treated) / (Count All Untreated).")

graph export "TWFE Bias by Percent Never Treated Scatter.png", replace width(6000)
graph export "TWFE Bias by Percent Never Treated Scatter.jpg", replace width(2000)
}

* TWFE Bias by % balanced
{
cap program drop sim_unit
program sim_unit, rclass
clear

gl obs_count 50000
set obs $obs_count
local ids = ${obs_count}/20
g id = mod(_n, `ids')+1
sort id
g time = mod(_n, 20)+1

preserve
keep if time == 1
sort id
g cohort = runiformint(1, 15)
replace cohort = 0 if cohort<10
tempfile to_merge
save `to_merge'
restore
merge m:1 id using `to_merge', nogen

g fe_unit = id/1000
g fe_time = time/20
g treat = cond(!cohort, 0, time >= cohort)

gl cutoff = runiform(0, 0.49)
return scalar cutoff = ${cutoff}

gen rd =uniform()
preserve
keep id time rd
gsort id rd
gen x=(rd<${cutoff})
egen xx  = min(time) if x== 1, by(id)
egen xx2 = max(time) if x== 1, by(id)
egen drop_t  = mean(xx), by(id)
egen drop_t2 = mean(xx2), by(id)

* prevent singletons
replace drop_t  = 19 if drop_t == 20
replace drop_t2 = drop_t + 1 if drop_t == drop_t2

replace drop_t = 0 if mi(drop_t)
replace drop_t2 = 0 if mi(drop_t2)
drop if !inrange(time, drop_t, drop_t2)  // create late entrant in the data
keep id time
tempfile a
save `a', replace
restore
merge 1:1 id time using `a', keep(3)

return scalar pc_balanced = 100*(_N/${obs_count})

gl effect = 5
g treat_effect_1 = $effect * treat // constant, homogeneous
su cohort if treat
g treat_effect_2 =  ($effect + cohort - r(mean)) * treat // constant, heterogenous
g exposure_time = time - cohort if cohort > 0
replace exposure_time     = 999 if cohort == 0
su exposure_time if treat
g treat_effect_3 = ($effect + exposure_time - r(mean)) * treat // dynamic, homogeneous
g treat_effect_5 = ($effect - exposure_time + r(mean)) * treat // dynamic, homogeneous

g cohorttime = (cohort) + ((cohort/5-2) * exposure_time)
su cohorttime if treat
g treat_effect_4 = ($effect + cohorttime - r(mean)) * treat // dynamic, heterogenous
g treat_effect_6 = ($effect - cohorttime + r(mean)) * treat // dynamic, heterogenous
g epsilon = rnormal()
forval i = 1/6 {
g y_`i' = fe_unit + fe_time + treat_effect_`i' + epsilon
reghdfe y_`i' treat, a(id time)
return scalar coef_`i' = _b[treat]
return scalar bias_`i' = _b[treat]-$effect
}
eret clear
end

gl num_replications = 10000
psimulate2, seed(123) r(${num_replications}) p(4) saving(twfe_bias_by_pc_balanced) seedsave(twfe_bias_by_pc_balanced_seed): sim_unit

use twfe_bias_by_pc_balanced, clear

format bias_* coef_* %03.1f

binscatter coef_1 coef_2 coef_3 coef_4 coef_5 coef_6 pc_balanced if pc_balanced>=10, line(connect) yti("Estimated Treatment Effect") xti("% Balanced") legend(order(1 "1: Homogeneous, Static" 2 "2: Heterogeneous, Static" 3 "3: Homogeneous, Dynamic" 4 "4: Heterogeneous, Dynamic" 5 "5: Homogeneous, Dynamic (Decreasing)" 6 "6: Heterogeneous, Dynamic (Decreasing)") pos(6) rows(3)) plotregion(lcol(black)) yline(5, lcol(maroon) lp(dash)) note("{bf: Notes}: n = 10000 replications. TWFE estimates use the reghdfe package." "Monte Carlo simulations run with the simulate2 package. % Balanced defined as (# observations)/(#units * #times)." "This chart made with the binscatter package.") xlab(10(10)90) colors(stc1 stc2 stc3 stc4 stc5 stc6)

graph export "TWFE Bias by Percent Balanced Crop.png", replace width(6000)
graph export "TWFE Bias by Percent Balanced Crop.jpg", replace width(2000)

twoway ///
(scatter coef_1 pc_balanced, msize(tiny) msymbol(smcircle_hollow) mcol(%50)) ///
(scatter coef_2 pc_balanced, msize(tiny) msymbol(smcircle_hollow) mcol(%50)) ///
(scatter coef_3 pc_balanced, msize(tiny) msymbol(smcircle_hollow) mcol(%50)) ///
(scatter coef_4 pc_balanced, msize(tiny) msymbol(smcircle_hollow) mcol(%50)) ///
(scatter coef_5 pc_balanced, msize(tiny) msymbol(smcircle_hollow) mcol(%50)) ///
(scatter coef_6 pc_balanced, msize(tiny) msymbol(smcircle_hollow) mcol(%50)) ///
(scatteri . ., mcol(stc1)) ///
(scatteri . ., mcol(stc2)) ///
(scatteri . ., mcol(stc3)) ///
(scatteri . ., mcol(stc4)) ///
(scatteri . ., mcol(stc5)) ///
(scatteri . ., mcol(stc6)) ///
 if pc_balanced>=10, yti("Estimated Treatment Effect") xti("% Balanced") ///
legend(order(7 "1: Homogeneous, Static" 8 "2: Heterogeneous, Static" 9 "3: Homogeneous, Dynamic" 10 "4: Heterogeneous, Dynamic" 11 "5: Homogeneous, Dynamic (Decreasing)" 12 "6: Heterogeneous, Dynamic (Decreasing)") pos(6) rows(3)) xlab(10(10)90) ///
plotregion(lcol(black)) yline(5, lcol(maroon) lp(dash)) note("{bf: Notes}: n = 10000 replications. TWFE estimates use the reghdfe package." "Monte Carlo simulations run with the simulate2 package. % Balanced defined as (# observations)/(#units * #times).")

graph export "TWFE Bias by Percent Balanced Scatter Crop.png", replace width(6000)
graph export "TWFE Bias by Percent Balanced Scatter Crop.jpg", replace width(2000)
}

* Comparing performance of various estimators by % balanced
{
cap program drop run_regs_fast
program define run_regs_fast, rclass
* create data
{
clear
gl obs_count 50000
set obs $obs_count
local ids = ${obs_count}/20
g id = mod(_n, `ids')+1
sort id
g time = mod(_n, 20)+1
xtset id time

preserve
keep if time == 1
sort id
g cohort = runiformint(1, 15)
replace cohort = 0 if cohort<10
tempfile to_merge
save `to_merge'
restore
merge m:1 id using `to_merge', nogen

g fe_unit = id/1000
g fe_time = time/20
g treated = cohort > 0
g treat = cond(!cohort, 0, time >= cohort)

gl cutoff = runiform(0, 0.49)
return scalar cutoff = ${cutoff}

gen rd =uniform()
preserve
keep id time rd
gsort id rd
gen x=(rd<${cutoff})
egen xx  = min(time) if x== 1, by(id)
egen xx2 = max(time) if x== 1, by(id)
egen drop_t  = mean(xx), by(id)
egen drop_t2 = mean(xx2), by(id)

* prevent singletons
replace drop_t  = 19 if drop_t == 20
replace drop_t2 = drop_t + 1 if drop_t == drop_t2

replace drop_t = 0 if mi(drop_t)
replace drop_t2 = 0 if mi(drop_t2)
drop if !inrange(time, drop_t, drop_t2)  // create late entrant in the data
keep id time
tempfile a
save `a', replace
restore
merge 1:1 id time using `a', keep(3)

return scalar n = _N
return scalar pc_balanced = 100*(_N/${obs_count})

gl effect = 5
g treat_effect_1 = $effect * treat // constant, homogeneous
su cohort if treat
g treat_effect_2 =  ($effect + cohort - r(mean)) * treat // constant, heterogenous

g exposure_time = time - cohort if cohort > 0
su exposure_time
di r(min)
g exposure_time_adj = exposure_time - r(min)
replace exposure_time     = 999 if cohort == 0
replace exposure_time_adj = 999 if cohort == 0

su exposure_time if treat
g treat_effect_3 = ($effect + exposure_time - r(mean)) * treat // dynamic, homogeneous
g treat_effect_5 = ($effect - exposure_time + r(mean)) * treat // dynamic, homogeneous

g cohorttime = (cohort) + ((cohort/5-2) * exposure_time)
su cohorttime if treat
g treat_effect_4 = ($effect + cohorttime - r(mean)) * treat // dynamic, heterogenous
g treat_effect_6 = ($effect - cohorttime + r(mean)) * treat // dynamic, heterogenous
g epsilon = rnormal()

g x = rnormal()
egen x_avg = mean(x), by(cohort)
g x_demean = x - x_avg

forval i = 1/6 {
g y_`i' = x + fe_unit + fe_time + treat_effect_`i' + epsilon
}

keep id time cohort exposure_time exposure_time_adj treat treated x x_demean y_1 y_2 y_3 y_4 y_5 y_6

cap drop sa_*
forval c = 10/15 {
forval e = 0/24 {
if `e' != 13 {
g sa_`c'_`e' = cohort == `c' & exposure_time_adj == `e'
}
}
}
}

forval i = 1/6 {
* TWFE
reghdfe y_`i' treat, a(id time)
return scalar  b_twfe_non_`i' = _b[treat]

* TWFE with control
reghdfe y_`i' treat x, a(id time)
return scalar  b_twfe_ctrl_`i' = _b[treat]

* Sun/Abr
{
cap drop sa_effect
reghdfe y_`i' sa_*, a(id time) 

g sa_effect=.
forval c = 10/15 {
forval e = 0/24 {
if `e' != 13 {
replace sa_effect = _b[sa_`c'_`e'] if cohort == `c' & exposure_time_adj == `e' & treat
}
}
}
su sa_effect

return scalar b_sa_non_`i' = r(mean)
}

* Wooldr with control
{
reghdfe y_`i' 1.treat#i.cohort#i.time 1.treat#i.cohort#i.time#c.x_demean, a(i.time##c.x i.cohort##c.x) technique(lsqr)

margins, dydx(treat) subpop(if treat) nose

return scalar  b_wooldr_ctrl_`i' = r(b)[1,1]
}

* Gardner
{
did2s y_`i', first_stage(i.time) second_stage(treat) treatment(treat ) unit(id) cluster(id)

return scalar  b_grdnr_non_`i' = e(b)[1,1]
}

* Gardner with control
{
did2s y_`i', first_stage(i.time x) second_stage(treat) treatment(treat ) unit(id) cluster(id)
return scalar  b_grdnr_ctrl_`i' = e(b)[1,1]
}

preserve
replace cohort = . if cohort == 0
* Borusyak et al.
{
did_imputation y_`i' id time cohort, autosample nose
return scalar  b_borys_non_`i' = e(b)[1, 1]
}

* Borusyak et al. control
{
did_imputation y_`i' id time cohort, controls(x) autosample nose
return scalar  b_borys_ctrl_`i' = e(b)[1, 1]
}
restore

* Callaway/Sant'Anna
{
xthdidregress ra (y_`i') (treat), group(id) vce(cluster id)
estat aggregation
return scalar  b_callwy_non_`i' = r(b)[1, 1]
}

* Callaway/Sant'Anna control
{
xthdidregress ra (y_`i' x) (treat), group(id) vce(cluster id)
estat aggregation
return scalar  b_callwy_ctrl_`i' = r(b)[1, 1]
}
}

eret clear
end

gl num_replications 1000
psimulate2, seed(123) r(${num_replications}) p(4) saving(monte_carlo_table_unbalanced): run_regs_fast

* Combined scatter chart
{
use monte_carlo_table_unbalanced, replace

gl graph_list ""

foreach i in twfe_non sa_non grdnr_non borys_non callwy_non twfe_ctrl wooldr_ctrl grdnr_ctrl borys_ctrl callwy_ctrl {
forval j = 1/6 {
	* labels
	{
	if "`i'" == "twfe_non" {
	local title "TWFE, No Control"
	}
	if "`i'" == "twfe_ctrl" {
	local title "TWFE, Control"
	}
	if "`i'" == "sa_non" {
	local title "Sun/Abraham, No Control"
	}
	if "`i'" == "wooldr_ctrl" {
	local title "Wooldridge, Control"
	}
	if "`i'" == "borys_ctrl" {
	local title "Borusyak et al., Control"
	}
	if "`i'" == "borys_non" {
	local title "Borusyak et al., No Control"
	}
	if "`i'" == "grdnr_ctrl" {
	local title "Gardner, Control"
	}
	if "`i'" == "grdnr_non" {
	local title "Gardner, No Control"
	}
	if "`i'" == "callwy_non" {
	local title "Callaway/Sant'Anna, No Control"
	}
	if "`i'" == "callwy_ctrl" {
	local title "Callaway/Sant'Anna, Control"
	}
	}
	
	twoway scatter b_`i'_`j' pc_balanced if pc_balanced>10, msize(tiny) yline(5, lcol(maroon) lp(dash)) ysc(r(3 7)) ylab(3(1)7, nogrid) xti("% Balanced") xlab(10(10)90) plotregion(lcol(black)) name(`i'_`j', replace) ti("`title', `j'", size(medsmall)) yti("")
	
	gl graph_list "${graph_list} `i'_`j'"
} 
}

graph combine ${graph_list}, xcommon ycommon rows(10) name(combinedscatter, replace) imargin(zero)  graphregion(margin(zero)) ysize(7.5) note("{bf: Notes}: n = ${num_replications} replications. True ATT is 5, represented by the dashed red line." "TWFE, Sun/Abraham, and Wooldridge estimates use the reghdfe package. Gardner estimates use the did2s package." "Borusyak et al. estimates use the did_imputation package. Monte Carlo simulations run with the simulate2 package." "% Balanced defined as (# observations)/(#units * #times).", size(tiny))
graph export "Scatters by Percent Balanced Common Scale.png", replace width(6000)
graph export "Scatters by Percent Balanced Common Scale.jpg", replace width(2000)

* relative scale
{
gl graph_list ""

foreach i in twfe_non sa_non grdnr_non borys_non callwy_non twfe_ctrl wooldr_ctrl grdnr_ctrl borys_ctrl callwy_ctrl {
forval j = 1/6 {
	* labels
	{
	if "`i'" == "twfe_non" {
	local title "TWFE, No Control"
	}
	if "`i'" == "twfe_ctrl" {
	local title "TWFE, Control"
	}
	if "`i'" == "sa_non" {
	local title "Sun/Abraham, No Control"
	}
	if "`i'" == "wooldr_ctrl" {
	local title "Wooldridge, Control"
	}
	if "`i'" == "borys_ctrl" {
	local title "Borusyak et al., Control"
	}
	if "`i'" == "borys_non" {
	local title "Borusyak et al., No Control"
	}
	if "`i'" == "grdnr_ctrl" {
	local title "Gardner, Control"
	}
	if "`i'" == "grdnr_non" {
	local title "Gardner, No Control"
	}
	if "`i'" == "callwy_non" {
	local title "Callaway/Sant'Anna, No Control"
	}
	if "`i'" == "callwy_ctrl" {
	local title "Callaway/Sant'Anna, Control"
	}
	}
	
	qui su b_`i'_`j'  if pc_balanced>10
	local min_dev = abs(5-r(min))
	local max_dev = abs(5-r(max))
	local dev = cond(`max_dev' > `min_dev', `max_dev', `min_dev')
	local chart_min = 5-`dev'
	local chart_max = 5+`dev'
	
	twoway scatter b_`i'_`j' pc_balanced if pc_balanced>10, msize(tiny) ysc(r(`chart_min' `chart_max')) yline(5, lcol(maroon) lp(dash)) ylab(#5, nogrid) xti("% Balanced") xlab(10(10)90) plotregion(lcol(black)) name(`i'_`j', replace) ti("`title', `j'", size(medsmall)) yti("")
	
	gl graph_list "${graph_list} `i'_`j'"
} 
}

graph combine ${graph_list}, xcommon rows(10) name(combinedscatter, replace) imargin(zero)  graphregion(margin(zero)) ysize(7.5) note("{bf: Notes}: n = ${num_replications} replications. True ATT is 5, represented by the dashed red line." "TWFE, Sun/Abraham, and Wooldridge estimates use the reghdfe package. Gardner estimates use the did2s package." "Borusyak et al. estimates use the did_imputation package. Monte Carlo simulations run with the simulate2 package." "% Balanced defined as (# observations)/(#units * #times).", size(tiny))
graph export "Scatters by Percent Balanced Relative Scale.png", replace width(6000)
graph export "Scatters by Percent Balanced Relative Scale.jpg", replace width(2000)
}
} 

* Combined binscatter chart
{
use monte_carlo_table_unbalanced, replace

gl graph_list ""

foreach i in twfe_non sa_non grdnr_non borys_non callwy_non twfe_ctrl wooldr_ctrl grdnr_ctrl borys_ctrl callwy_ctrl {
forval j = 1/6 {
	* labels
	{
	if "`i'" == "twfe_non" {
	local title "TWFE, No Control"
	}
	if "`i'" == "twfe_ctrl" {
	local title "TWFE, Control"
	}
	if "`i'" == "sa_non" {
	local title "Sun/Abraham, No Control"
	}
	if "`i'" == "wooldr_ctrl" {
	local title "Wooldridge, Control"
	}
	if "`i'" == "borys_ctrl" {
	local title "Borusyak et al., Control"
	}
	if "`i'" == "borys_non" {
	local title "Borusyak et al., No Control"
	}
	if "`i'" == "grdnr_ctrl" {
	local title "Gardner, Control"
	}
	if "`i'" == "grdnr_non" {
	local title "Gardner, No Control"
	}
	if "`i'" == "callwy_non" {
	local title "Callaway/Sant'Anna, No Control"
	}
	if "`i'" == "callwy_ctrl" {
	local title "Callaway/Sant'Anna, Control"
	}
	}
	
	 binscatter b_`i'_`j' pc_balanced if pc_balanced>10, line(connect) yline(5, lcol(maroon) lp(dash)) ysc(r(3 7)) ylab(3(1)7, nogrid) xti("% Balanced") xlab(10(10)90) plotregion(lcol(black)) name(`i'_`j', replace) ti("`title', `j'", size(medsmall)) yti("")
	
	gl graph_list "${graph_list} `i'_`j'"
} 
}

graph combine ${graph_list}, xcommon ycommon rows(10) name(combinedscatter, replace) imargin(zero)  graphregion(margin(zero)) ysize(7.5) note("{bf: Notes}: n = ${num_replications} replications. True ATT is 5, represented by the dashed red line." "TWFE, Sun/Abraham, and Wooldridge estimates use the reghdfe package. Gardner estimates use the did2s package." "Borusyak et al. estimates use the did_imputation package. Monte Carlo simulations run with the simulate2 package." "This chart made with the binscatter package." "% Balanced defined as (# observations)/(#units * #times).", size(tiny))
graph export "Binscatters by Percent Balanced Common Scale.png", replace width(6000)
graph export "Binscatters by Percent Balanced Common Scale.jpg", replace width(2000)

* relative scale
{
gl graph_list ""

foreach i in twfe_non sa_non grdnr_non borys_non callwy_non twfe_ctrl wooldr_ctrl grdnr_ctrl borys_ctrl callwy_ctrl {
forval j = 1/6 {
	* labels
	{
	if "`i'" == "twfe_non" {
	local title "TWFE, No Control"
	}
	if "`i'" == "twfe_ctrl" {
	local title "TWFE, Control"
	}
	if "`i'" == "sa_non" {
	local title "Sun/Abraham, No Control"
	}
	if "`i'" == "wooldr_ctrl" {
	local title "Wooldridge, Control"
	}
	if "`i'" == "borys_ctrl" {
	local title "Borusyak et al., Control"
	}
	if "`i'" == "borys_non" {
	local title "Borusyak et al., No Control"
	}
	if "`i'" == "grdnr_ctrl" {
	local title "Gardner, Control"
	}
	if "`i'" == "grdnr_non" {
	local title "Gardner, No Control"
	}
	if "`i'" == "callwy_non" {
	local title "Callaway/Sant'Anna, No Control"
	}
	if "`i'" == "callwy_ctrl" {
	local title "Callaway/Sant'Anna, Control"
	}
	}
	
	qui su b_`i'_`j'  if pc_balanced>10
	local min_dev = abs(5-r(min))
	local max_dev = abs(5-r(max))
	local dev = cond(`max_dev' > `min_dev', `max_dev', `min_dev')
	local chart_min = 5-`dev'
	local chart_max = 5+`dev'
	
	binscatter b_`i'_`j' pc_balanced if pc_balanced>10, line(connect) ysc(r(`chart_min' `chart_max')) yline(5, lcol(maroon) lp(dash)) ylab(#5, nogrid) xti("% Balanced") xlab(10(10)90) plotregion(lcol(black)) name(`i'_`j', replace) ti("`title', `j'", size(medsmall)) yti("")
	
	gl graph_list "${graph_list} `i'_`j'"
} 
}

graph combine ${graph_list}, xcommon rows(10) name(combinedscatter, replace) imargin(zero)  graphregion(margin(zero)) ysize(7.5) note("{bf: Notes}: n = ${num_replications} replications. True ATT is 5, represented by the dashed red line." "TWFE, Sun/Abraham, and Wooldridge estimates use the reghdfe package. Gardner estimates use the did2s package." "Borusyak et al. estimates use the did_imputation package. Monte Carlo simulations run with the simulate2 package." "This chart made with the binscatter package." "% Balanced defined as (# observations)/(#units * #times).", size(tiny))
graph export "Binscatters by Percent Balanced Relative Scale.png", replace width(6000)
graph export "Binscatters by Percent Balanced Relative Scale.jpg", replace width(2000)
}
} 


}

* Checking for correlation in generating unbalanced data
{
cap program drop create_unbalanced
program define create_unbalanced, rclass
* create data
{
clear
gl obs_count 50000
set obs $obs_count
local ids = ${obs_count}/20
g id = mod(_n, `ids')+1
sort id
g time = mod(_n, 20)+1
xtset id time

preserve
keep if time == 1
sort id
g cohort = runiformint(1, 15)
replace cohort = 0 if cohort<10
tempfile to_merge
save `to_merge'
restore
merge m:1 id using `to_merge', nogen

g fe_unit = id/1000
g fe_time = time/20
g treated = cohort > 0
g treat = cond(!cohort, 0, time >= cohort)

gl cutoff = runiform(0, 0.49)
return scalar cutoff = ${cutoff}

gen rd =uniform()
preserve
keep id time rd
gsort id rd
gen x=(rd<${cutoff})
egen xx  = min(time) if x== 1, by(id)
egen xx2 = max(time) if x== 1, by(id)
egen drop_t  = mean(xx), by(id)
egen drop_t2 = mean(xx2), by(id)

* prevent singletons
replace drop_t  = 19 if drop_t == 20
replace drop_t2 = drop_t + 1 if drop_t == drop_t2

replace drop_t = 0 if mi(drop_t)
replace drop_t2 = 0 if mi(drop_t2)
drop if !inrange(time, drop_t, drop_t2)  // create late entrant in the data
keep id time
tempfile a
save `a', replace
restore
merge 1:1 id time using `a', keep(3)

return scalar n = _N
return scalar pc_balanced = 100*(_N/${obs_count})

gl effect = 5
g treat_effect_1 = $effect * treat // constant, homogeneous
su cohort if treat
g treat_effect_2 =  ($effect + cohort - r(mean)) * treat // constant, heterogenous

g exposure_time = time - cohort if cohort > 0
su exposure_time
di r(min)
g exposure_time_adj = exposure_time - r(min)
replace exposure_time     = 999 if cohort == 0
replace exposure_time_adj = 999 if cohort == 0

su exposure_time if treat
g treat_effect_3 = ($effect + exposure_time - r(mean)) * treat // dynamic, homogeneous
g treat_effect_5 = ($effect - exposure_time + r(mean)) * treat // dynamic, homogeneous

g cohorttime = (cohort) + ((cohort/5-2) * exposure_time)
su cohorttime if treat
g treat_effect_4 = ($effect + cohorttime - r(mean)) * treat // dynamic, heterogenous
g treat_effect_6 = ($effect - cohorttime + r(mean)) * treat // dynamic, heterogenous
g epsilon = rnormal()

g x = rnormal()
egen x_avg = mean(x), by(cohort)
g x_demean = x - x_avg

forval i = 1/6 {
g y_`i' = x + fe_unit + fe_time + treat_effect_`i' + epsilon
}

keep id time cohort exposure_time exposure_time_adj treat treated x x_demean y_1 y_2 y_3 y_4 y_5 y_6

cap drop sa_*
forval c = 10/15 {
forval e = 0/24 {
if `e' != 13 {
g sa_`c'_`e' = cohort == `c' & exposure_time_adj == `e'
}
}
}
}
eret clear
end

set seed 123

forval i = 1/100 {
	create_unbalanced
	g iteration = `i'
	save combined_data`i', replace
}

clear
forval i = 1/100 {
	append using combined_data`i'
	erase combined_data`i'.dta
}
save combined_data, replace

use combined_data, replace
contract time treat
twoway (hist time if !treat [fw=_freq], percent discrete col(blue%50)) (hist time if treat [fw=_freq], percent discrete col(red%50)), legend(order(1 "Untreated Observations" 2 "Treated Observations") pos(6) rows(1)) plotregion(lcol(black)) xti("Time Period in Data") note("{bf: Note}: 100 randomly unbalanced datasets were generated and stacked.")

graph export "Average Distribution of Data in Unbalanced Simulations by Treated.png", replace width(6000)
graph export "Average Distribution of Data in Unbalanced Simulations by Treated.jpg", replace width(2000)

use combined_data, replace
contract time
twoway (hist time[fw=_freq], percent discrete col(blue%50)), plotregion(lcol(black)) xti("") note("{bf: Note}: 100 randomly unbalanced datasets were generated and stacked.")

graph export "Average Distribution of Data in Unbalanced Simulations.png", replace width(6000)
graph export "Average Distribution of Data in Unbalanced Simulations.jpg", replace width(2000)

use combined_data, replace
pwcorr time treat, sig
}

* Sensitivities on TWFE bias
{
cap program drop sim_unit_sensitivity
program sim_unit_sensitivity, rclass
	clear
local times = runiformint(1, 10)*5
local units = runiformint(5, 10)*500
local n = `times'*`units'

return scalar n = `n'
return scalar count_units = `units'
return scalar count_times = `times'

set obs `n'
g id = mod(_n, `units')+1
sort id
g time = mod(_n, `times')+1

gl pc_never_treated = runiform()

local cohort_start = runiformint(`times'/5+3, `times'-1)
local cohort_end   = runiformint(`cohort_start'+1, `times')
return scalar cohort_start = `cohort_start'
return scalar cohort_end   = `cohort_end'

preserve
keep if time == 1
sort id
g prob = runiform()
g cohort = 0 if prob < $pc_never_treated
replace cohort = runiformint(`cohort_start', `cohort_end') if mi(cohort)
tempfile to_merge
save `to_merge'
restore

merge m:1 id using `to_merge', nogen

count if cohort == 0
return scalar pc_never_treated = 100*r(N)/_N

g fe_unit = id/1000
g fe_time = time/20

g exposure_time = time - cohort if cohort > 0
replace exposure_time = 999 if cohort == 0

g treat = cond(!cohort, 0, time >= cohort)

* unbalanced
{
gl cutoff = runiform(0.05, 0.2)
return scalar cutoff = ${cutoff}

gen rd =uniform()
preserve
keep id time rd
gsort id rd
gen x=(rd<${cutoff})
egen xx  = min(time) if x== 1, by(id)
egen xx2 = max(time) if x== 1, by(id)
egen drop_t  = mean(xx), by(id)
egen drop_t2 = mean(xx2), by(id)

* prevent singletons
replace drop_t  = 19 if drop_t == 20
replace drop_t2 = drop_t + 1 if drop_t == drop_t2

replace drop_t = 0 if mi(drop_t)
replace drop_t2 = 0 if mi(drop_t2)

drop if !inrange(time, drop_t, drop_t2)  // create late entrant in the data
keep id time
tempfile a
save `a', replace
restore
merge 1:1 id time using `a', keep(3)

return scalar pc_balanced = 100*(_N/`n')
}

su treat
return scalar pc_treated = r(mean) * 100

gl effect = runiformint(-10, 10)
return scalar true_effect = $effect

* individual-specific treatment effect modifier
if runiform() < 0.25 {
	return scalar effect_noise_flag = 1
	if $effect == 0 {
	gl effect_noise_parameter = abs(runiform(1/10, abs(1/2))) // required or the below rule will give undefined value
	}
	else {
	gl effect_noise_parameter = abs(runiform(${effect}/10, abs(${effect}/2)))
	}
	
	return scalar effect_noise_sd = ${effect_noise_parameter}
	g treat_modifier = rnormal(0, ${effect_noise_parameter})
}
else {
	return scalar effect_noise_flag = 0
	g treat_modifier = 0
	return scalar effect_noise_sd = 0
}

g treat_effect_1 = ($effect + treat_modifier) * treat // constant, homogeneous

su cohort if treat
g treat_effect_2 =  (($effect + treat_modifier) + cohort - r(mean)) * treat // constant, heterogenous

su exposure_time if treat
g treat_effect_3 = (($effect + treat_modifier) + exposure_time - r(mean)) * treat // dynamic, homogeneous
g treat_effect_5 = (($effect + treat_modifier) - exposure_time + r(mean)) * treat // dynamic, homogeneous

g cohorttime = (cohort) + ((cohort/5-2) * exposure_time)
su cohorttime if treat
g treat_effect_4 = (($effect + treat_modifier) + cohorttime - r(mean)) * treat // dynamic, heterogenous
g treat_effect_6 = (($effect + treat_modifier) - cohorttime + r(mean)) * treat // dynamic, heterogenous

g epsilon = rnormal()

* covariate modifiers
if runiform() < 0.25 {
	return scalar correlation_effect_flag = 1
	gl x_correlation_parameter = rnormal()
	return scalar x_correlation_parameter = ${x_correlation_parameter}
	g correlation_effect = ${x_correlation_parameter} * (treat == 1)
}
else {
	return scalar correlation_effect_flag = 0
	return scalar x_correlation_parameter = 0
	g correlation_effect = 0
}
gl x_size = runiformint(-5, 5)
return scalar x_size = ${x_size}
g x = rnormal(${x_size}) + correlation_effect

forval i = 1/6 {
g y_`i' = x + fe_unit + fe_time + treat_effect_`i' + epsilon

cap reghdfe y_`i' treat, a(id time)
if _rc {
return scalar regression_error = 1
}
else {
return scalar regression_error = 0
su treat_effect_`i' if treat
return scalar b_`i' = _b[treat]-r(mean)
}
}

eret clear
end

psimulate2, seed(123) r(50000) p(4) saving(twfe_bias_sensitivities_new): sim_unit_sensitivity

use twfe_bias_sensitivities_new, replace

* absolute bias
forval i = 1/6 {
	g abs_b_`i' = abs(b_`i')
	label var abs_b_`i' "Data `i'"
}

replace count_units = count_units / 1000

g pc_never_treated2 = pc_never_treated^2
g pc_balanced2 = pc_balanced^2

g longest_exposure = count_times - cohort_start

* labels
{
label var count_units  "Count of Units (1000s)"
label var count_times  "Count of Times"
label var cohort_start  "First Treatment Time"
label var cohort_end  "Last Treatment Time"
label var pc_never_treated  "\% Never-Treated"
label var pc_never_treated2 "\% Never-Treated Squared"
label var pc_balanced  "\% Balanced"
label var pc_balanced2 "\% Balanced Squared"
label var true_effect  "True Effect Size"
label var effect_noise_sd  "Individual-Specific Treatment Noise"
label var x_size "Covariate Size"
label var x_correlation_parameter  "Covariate Treatment Factor"
label var longest_exposure  "Longest Treatment Exposure Time"
label var pc_treated  "\% Actively Treated"
}

preserve
forval i = 1/6 {
	label var abs_b_`i' "Absolute Bias, Treatment Effect `i'"
}

estpost su abs_b_* count_units count_times pc_balanced pc_treated pc_never_treated cohort_start cohort_end longest_exposure true_effect effect_noise_sd x_correlation_parameter x_size if !regression_error, d
est store a

esttab a using "sensitivity summary stats.tex", replace ///
collabels(n Mean Std.Dev. Median Min Max) ///
cells("count(fmt(%7.0fc)) mean(fmt(3)) sd(fmt(3)) p50(fmt(3)) min(fmt(3)) max(fmt(3))") label noobs nomtitles nonumbers
restore

eststo clear

forval i = 1/6 {
reg abs_b_`i' pc_never_treated pc_never_treated2 count_units cohort_start longest_exposure pc_balanced pc_balanced2 true_effect effect_noise_sd x_correlation_parameter x_size	

if `i' == 1 {
	local title "1: Homogeneous, Static"
}
if `i' == 2 {
	local title "2: Heterogeneous, Static"
}
if `i' == 3 {
	local title "3: Homogeneous, Dynamic"
}
if `i' == 4 {
	local title "4: Heterogeneous, Dynamic"
}
if `i' == 5 {
	local title "5: Homogeneous, Dynamic (Decreasing)"
}
if `i' == 6 {
	local title "6: Heterogeneous, Dynamic (Decreasing)"
}
eststo, title("`title'")
}

esttab using "sensitivity regression 1.tex", ///
cells(b(star fmt(3) label(Coefficient)) se(par fmt(3) label(SE))) ///
stats(N r2_a F p, fmt() ///
labels("Row Count" "Adjusted R Squared" "F Statistic" "Chi Squared p-value"))  ///
label legend replace ti("Regression Against Absolute Bias") varlabels(_cons Constant) nonumbers

eststo clear

forval i = 1/6 {
reg abs_b_`i' (c.pc_never_treated c.pc_never_treated2)##(c.count_units c.cohort_start c.longest_exposure c.pc_balanced c.pc_balanced2 c.true_effect c.effect_noise_sd c.true_effect#c.effect_noise_sd c.x_correlation_parameter c.x_size c.x_correlation_parameter#c.x_size) 

if `i' == 1 {
	local title "1: Homogeneous, Static"
}
if `i' == 2 {
	local title "2: Heterogeneous, Static"
}
if `i' == 3 {
	local title "3: Homogeneous, Dynamic"
}
if `i' == 4 {
	local title "4: Heterogeneous, Dynamic"
}
if `i' == 5 {
	local title "5: Homogeneous, Dynamic (Decreasing)"
}
if `i' == 6 {
	local title "6: Heterogeneous, Dynamic (Decreasing)"
}
eststo, title("`title'")
}

esttab using "sensitivity regression 2.tex", ///
cells(b(star fmt(3) label(Coefficient)) se(par fmt(3) label(SE))) ///
stats(N r2_a F p, fmt() ///
labels("Row Count" "Adjusted R Squared" "F Statistic" "Chi Squared p-value"))  ///
label legend replace ti("Regression Against Absolute Bias") varlabels(_cons Constant) nonumbers
}
