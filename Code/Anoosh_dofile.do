********************************************
*      EIO PS1: Static Demand for diff'd products
********************************************

*Loading the data
cd "C:\Studies_TSE\Emp IO\PS1" //path to folder w dta
global output	"C:\Studies_TSE\Emp IO\PS1\output" // folder where I output tables to

 foreach name in USCPI UShouseholds gasprice{
 	import delimited using `name'.csv
	sort year
	save `name'.dta, replace
	clear
        }
		
import delimited carpanel.csv
sort year
	merge m:1 year using USCPI.dta, keepusing(cpi) nogen keep(master match)
	merge m:1 year using UShouseholds.dta, keepusing(nb_hh) nogen keep(master match)
	merge m:1 year using gasprice.dta, keepusing(gasprice) nogen keep(master match)
save carpanel.dta, replace

*Label the variables
label variable name model.name
label variable id "model id"
label variable yr "model year abbreviated"
label variable cy cylinders
label variable dr "number of doors"
label variable at "Automatic transmission"
label variable ps "Power steering"
label variable air "Air conditioning"
label variable drv "Front-wheel drive"
label variable p "nominal price"
label variable wt Weight
label variable dom "domestic brand"
label variable disp "engine displacement"
label variable hp "Horse power"
label variable lng length
label variable wdt width
label variable wb Wheelbase
label variable mpg "miles per gallon"
label variable q quantity
label variable firmids "car brand"
label variable euro "European"
label variable reli "reliability index"
label variable dfi "direct fuel injection"
label variable hp2wt "Horsepower/weight"
label variable size "Size"
label variable japan "Japanese"
label variable year "model year"
label variable cat "Size category"

*CPI adjusted car price and gas price:
gen padj = p*(cpi/100)
label variable padj "Price"
gen padj_th = padj/1000
label variable padj_th "Price (in 000s)"

gen gaspadj = gasprice*(cpi/100)
*Real dollars per mile:
gen dpm = gaspadj/mpg
label variable dpm "Dollars per mile"

*Compute potential market
gen msize = nb_hh*1000 //nb_hh in 000s
label variable msize "Market size"


********************************************
*               EXERCISE 1: LOGIT
********************************************

tab firmids, gen(firmids) // for brand FEs
tab year, gen(year) // year FEs
tab dr, gen(door) // dummies for num of doors 
	forvalues i = 4(-1)2{
		local j = `i'+1
		rename door`i' door`j' //rename acc to PS instructions
	}

// Market share variables of products + the outside good
gen sj = q/msize // where msize = Potential market
egen totq = sum(q), by(year) // gives Qm = Total Qty of cars purchased
gen lsj0 = log(q/(msize-totq)) // gives the dependent variable ~ log(sjm/s0m)

*** Generating instruments
//we should take into account possible endogeneity issues so let's construct BLP-type instruments

	gen con=1 // sums over constant is the same as taking counts in BLP as a measure of market competitiveness
foreach var of varlist con dpm door1 door3 door4 door5 at ps air drv wt hp2wt hp euro japan size wb{
	bysort year: egen sum1_`var' = sum(`var') 
	bysort year firmids: egen sum2_`var' = sum(`var')
	bysort year cat: egen sum3_`var' = sum(`var')
	bysort year cat firmids: egen sum4_`var' = sum(`var')
					
	gen i1_`var' = sum2_`var' - `var' 	// this is BLP's first instrument
	gen i2_`var' = sum1_`var' - sum2_`var' // this is BLP's second instrument
	gen i3_`var' = sum4_`var' - `var' 	// BLP's third
	gen i4_`var' = sum3_`var' - sum4_`var' // BLP's fourth 
		
	label var i1_`var' "BLP1 sum of `var' from products of the same firm (Logit)"
	label var i2_`var' "BLP2 sum of `var' from products of the other firms (Logit)"
	label var i3_`var' "BLP3 sum of `var' from products of the same firm, in the same size category (NL)"
	label var i4_`var' "BLP4 sum of `var' from products of other firms, in the same size category (NL)"
	}
drop sum*


* OLS w/o brand FEs
eststo: reg lsj0 padj_th dpm door3 door4 door5 at ps air drv wt hp2wt hp euro japan size wb , r
	estadd loc fe1 "" 

* OLS w/ brand FEs
eststo: reg lsj0 padj_th dpm door3 door4 door5 at ps air drv wt hp2wt hp euro japan size wb firmids2-firmids19, r
	estadd loc fe1 "$\checkmark$" 

* IV w/o brand FEs
eststo: ivregress 2sls lsj0  dpm door3 door4 door5 at ps air drv wt hp2wt hp euro japan size wb (padj_th =  i1* i2*), robust
	estadd loc fe1 "" 

* IV w/ brand FEs
eststo: ivregress 2sls lsj0 euro japan size dpm door3 door4 door5 at ps air drv wt hp2wt hp  wb firmids2-firmids19 (padj_th =  i1* i2*),  robust
	estadd loc fe1 "$\checkmark$" 

estat endogenous
estat overid
* !! Seems we need better instruments...

esttab * using "$output\logit.tex", tex title("Logit") ///
	replace se star(* 0.10 ** 0.05 *** 0.01) label scalars(r2) mtitle ///	
	drop(_cons door* firmids*) stats(N r2 fe1, label("Number of observations" "R-squared" "Brand FE")) ///
	noobs nonotes addnotes("Notes:" ///
	"Standard errors are heteroskedasticity-robust" ///
    "* $ p < 0.05$, ** $ p < 0.01$, *** $ p < 0.001$" ///
    ) ///
	mlabels("OLS" "" "IV" "" ) 

eststo clear

********************************************
*               EXERCISE 2: NESTED LOGIT
********************************************
ssc install weakiv // for the Weak IV test

// Generating qty variables by size category g in the nested structure
bysort year cat: egen qg = sum(q) // level of group: the car size category

// Generating share variables for product j by group
gen sjg = q/qg // market share in size category
label variable sjg "market share in size category"
// and taking logs
gen lsjg = log(sjg)
label variable lsjg "Log(market share in size category)"


* OLS w/o brand FEs
eststo: reg lsj0 padj_th lsjg dpm door3 door4 door5 at ps air drv wt hp2wt hp euro japan size wb , r
	estadd loc fe1 "" 

* OLS w/ brand FEs
eststo: reg lsj0 padj_th lsjg dpm door3 door4 door5 at ps air drv wt hp2wt hp euro japan size wb firmids2-firmids19, r
	estadd loc fe1 "$\checkmark$" 

* IV w/o brand FEs
eststo: ivreg2 lsj0 dpm door3 door4 door5 at ps air drv wt hp2wt hp euro japan size wb (padj_th lsjg =  i3* i4*), robust
	estadd loc fe1 "" 		
	estadd scalar hansenj = e(jp)
qui weakiv
	estadd scalar weakclr = `r(clr_p)' :est3

* IV w/ brand FEs
eststo: ivreg2 lsj0 dpm door3 door4 door5 at ps air drv wt hp2wt hp euro japan size wb firmids2-firmids19  (padj_th lsjg =  i3* i4*),  robust
	estadd loc fe1 "$\checkmark$" 
	estadd scalar hansenj = e(jp)
qui weakiv
	estadd scalar weakclr = `r(clr_p)' :est4

	
esttab * using "$output\nl.tex", tex title("Nested Logit") ///
	replace b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) label scalars(r2) mtitle ///	
	drop(_cons door* firmids*) stats(N r2 fe1 weakclr hansenj, label("Number of observations" "R-squared" "Brand FE" "CLR Weak IV" "Hansen J statistic")) ///
	noobs nonotes addnotes("Notes:" ///
	"Standard errors are heteroskedasticity-robust" ///
    "* $ p < 0.05$, ** $ p < 0.01$, *** $ p < 0.001$" ///
    ) ///
	mlabels("OLS" "" "IV" "" ) 

	eststo clear

********************************************
*               EXERCISE 3: RANDOM COEFF LOGIT
********************************************
*Install nec packages
ssc install ranktest
ssc install avar
ssc install rcl
ssc install blp

* Generating a Gandhi-Houde type "differentiation IV"
bysort year: egen count = sum(con) 
replace count = count-1 // Total num of products less own

	* using size
	gen size_sqr = size^2
	bysort year: egen sum1_size = sum(size) 
	bysort year: egen sum_size_sqr = sum(size_sqr) 

	gen oth_size = sum1_size - size 	// this is sum of size from all other products
	gen oth_size_sqr = sum_size_sqr - size_sqr 	// this is sum of size squared from all other products
	gen diff_iv_size = count*size_sqr -2*size*oth_size + oth_size_sqr
	label var oth_size "Sum of size from all other products"
	label var oth_size_sqr "Sum of size squared from all other products"
	label variable diff_iv_size "Size - differentiation IV"

	* using displacement
	gen disp_sqr = disp^2
	bysort year: egen sum1_disp = sum(disp) 
	bysort year: egen sum_disp_sqr = sum(disp_sqr) 

	gen oth_disp = sum1_disp - disp 	// this is sum of disp from all other products
	gen oth_disp_sqr = sum_disp_sqr - disp_sqr 	// this is sum of disp squared from all other products
	gen diff_iv_disp = count*disp_sqr -2*disp*oth_disp + oth_disp_sqr
	label var oth_disp "Sum of disp from all other products"
	label var oth_disp_sqr "Sum of disp squared from all other products"
	label variable diff_iv_disp "Disp - differentiation IV"


rcl sj dpm door3 door4 door5 at ps air drv wt hp2wt hp euro japan wb  (padj_th = diff_iv_size diff_iv_disp  i3* i4*) , market(year) rc(size) msize(msize) gmm2s
// gives sigma = 1.885 & significant
// q: whether or not to also include size in the list of exo variables?

blp sj dpm door3 door4 door5 at ps air drv wt hp2wt hp euro japan wb, stochastic(size) endog(padj_th size = diff_iv_size diff_iv_disp i3_con i3_dpm i3_door1 i3_door3 i3_door4 i3_door5 i3_at i3_ps i3_air i3_drv i3_wt i3_hp2wt i3_hp i3_euro i3_japan i3_size i3_wb i4_con i4_dpm i4_door1 i4_door3 i4_door4 i4_door5 i4_at i4_ps i4_air i4_drv i4_wt i4_hp2wt i4_hp i4_euro i4_japan i4_size i4_wb) markets(year) 
//^this comman requires that the coeff on size also have a mean utility component i.e. !=0
// which does not match the thing Ana has prescribed in the PS



********************************************
*               EXERCISE 4: MARKUPS, ELASTICITIES
********************************************
rcl sj dpm door3 door4 door5 at ps air drv wt hp2wt hp euro japan wb  (padj_th = diff_iv_size diff_iv_disp  i3* i4*) , market(year) rc(size) msize(msize) gmm2s msimulation(firmids) onlymc
// gives pre-merger markups, marg cost

********************************************
*               EXERCISE 5: MERGER SIMULATION
********************************************
gen firmids_post = firmids 
replace firmids_post=7 if firmids==17 // Fiat becomes VW 
rcl sj dpm door3 door4 door5 at ps air drv wt hp2wt hp euro japan wb  (padj_th = diff_iv_size diff_iv_disp  i3* i4*) , market(year) rc(size) msize(msize) gmm2s msimulation(firmids firmids_post)

gen dp = __p_post - padj_th

***    END	***