	program define belieffit, eclass sortpreserve
		version 17.0
		// ylog Xvars
		syntax varlist(min=1) [if] [in] [aw fw iw pw], ///
			[ ABSORB(varlist fv) VCE(varlist) REPs(integer 200) ///
			  TOLerance(real 1e-6) MAXITer(integer 100) ///
			  FEPREFIX(name) VARPREFIX(name) ///
			  DF(string) AUX(string) replace prednames(string) calc(string)]
		
		quietly {
		marksample touse
		gettoken Ylog Xvars : varlist

		if `reps'==1 {
			noi di as error "reps() must be 0 or >1."
			exit 301
		}
		if "`prednames'"!="" {
			loc numnames: word count `prednames'
			if `numnames'!=2 {
				noi di as error "prednames() requires two variable names - one for the corrected mean, one for the predicted variance."
			}
			foreach nm in `prednames' {
				cap confirm variable `nm'
				if _rc!=0 &"`replace'"=="" {
					noi di as error "Variable `nm' already exist."
					exit 301
				}
				else if _rc==0&"`replace'"!="" {
					drop `nm'
				}
			}
		}
		if ("`feprefix'"  == "") local feprefix  "fe"
		if ("`varprefix'" == "") local varprefix "S"

		if ("`df'"=="") local df "adj"
		local df = lower("`df'")
		if !inlist("`df'","adj","none") {
			di as err "df() must be adj or none"
			exit 198
		}

		// user weights
		tempvar userwt
		if ("`weight'" != "") gen double `userwt' `exp' if `touse'
		else                   gen double `userwt' = 1     if `touse'
		
		//Check exsting variable name conflict
		if "`absorb'"!="" {
			loc i=0
			foreach str in `absorb' `absorb' {
				loc ++i
				cap confirm variable `feprefix'`i'
				if _rc==0&"`replace'"=="" {
					noi di as error "Variable `feprefix'`i' already exist in data, use option replace to overwrite."
					exit 301
				}
				else if _rc==0&"`replace'"!="" {
					drop `feprefix'`i'
				}
			}
			if `reps'>0 {
				forvalues j=1/`i' {
					forvalues k=1/`j' {
						cap confirm variable `varprefix'`j'`k'
						if _rc==0&"`replace'"=="" {
							noi di as error "Variable `varprefix'`j'`k' already exist in data, use option replace to overwrite."
							exit 301
						}
						else if _rc==0&"`replace'"!="" {
							drop  `varprefix'`j'`k'
						}
					}
				}
			}
		}
		

		// ---- Baseline (warm start): point estimates + FE cols (for K) + aux blocks ----
		qui belieffit_engine `Ylog' `Xvars' if `touse' [aw=`userwt'], ///
			absorb(`absorb') tol(`tolerance') maxiter(`maxiter') ///
			feprefix(`feprefix') `ppmlfromopt' aux(`aux') prednames(`prednames') calc(`calc')

		
		tempname b_base
		matrix `b_base' = e(b)
		

		// If reps<=0: post baseline and exit (no bootstrap)
		if (`reps' <= 0) {
			ereturn clear
			ereturn post `b_base'
			ereturn scalar reps = 0
			ereturn local cmd "belieffit"
			ereturn local absorb "`absorb'"
			ereturn local feprefix "`feprefix'"
			ereturn local varprefix "`varprefix'"
			ereturn local regressors "`Xvars'"
			ereturn local depvar "`Ylog'"
			ereturn local vcetype "None"
			noi ereturn display
			exit
		}

		// ---- Count FE columns produced by engine (fe1..feKtotal) ----
		// (Ktotal = Kmean + Kvar = 2 * #absorbs when savefe(final) and absorb() non-empty)
		quietly ds `feprefix'*, has(type numeric)
		local K : word count `r(varlist)'

		// Initialize S*_* (FE covariance across reps) if FE present
		if (`K'>0) {		
			// Running means Mj and accumulators SSj_k, plus output Sj_k
			forvalues j=1/`K' {
				tempvar M`j'
				qui gen double `M`j'' = 0 if `touse'
			}
			forvalues j=1/`K' {
				forvalues k=1/`j' {
					tempvar SS`j'`k'
					qui gen double `SS`j'`k''     = 0 if `touse'
					qui gen double `varprefix'`j'`k' = 0 if `touse'
				}
			}
		}
		else {
			di as txt "(note: no FE columns found; will still compute coefficient bootstrap SEs)"
		}

		
		tempname b_bs V_bs
		
		// Map clusters to ids once (product exponential bootstrap)
		local clids
		if ("`vce'" != "") {
			foreach cv of varlist `vce' {
				tempvar cl`cv'
				qui egen long `cl`cv'' = group(`cv') if `touse'
				local clids `clids' `cl`cv''
			}
		}

		// variance of FE setup
		
		tempvar R_row
		gen long `R_row' = 0 if `touse'  // per-row count of successful updates
		
		// make one gid per absorb term; interactions -> group(var1 var2 ...)
		local gidlist
		local gidx = 0
		foreach term of local absorb {
			// strip fv prefixes and parentheses
			local clean : subinstr local term "i." "", all
			local clean : subinstr local clean "c." "", all
			local clean : subinstr local clean "b." "", all
			local clean : subinstr local clean "o." "", all
			local clean : subinstr local clean "("  "", all
			local clean : subinstr local clean ")"  "", all

			// split by '#'
			local vars ""
			local piece `"`clean'"'
			while (strpos("`piece'", "#")) {
				local p = strpos("`piece'", "#")
				local left  = substr("`piece'", 1, `p'-1)
				local right = substr("`piece'", `p'+1, .)
				local vars `vars' `left'
				local piece "`right'"
			}
			// last piece
			local vars `vars' `piece'
			// collapse doubled spaces
			local vars : list retokenize vars

			// create the group id for this absorb term (with missing kept as its own cell)
			tempvar gid_`++gidx'
			egen long `gid_`gidx'' = group(`vars'), missing
			local gidlist `gidlist' `gid_`gidx''
		}
		tempvar ok
		gen byte `ok'=1
		
		tempname fe


		// ---- Bootstrap loop ----
		loc r=1
		loc iter=0
		loc maxiter=1000
		while  (`r' <= `reps')&(`iter'<`maxiter') {
			loc ++iter
			if `iter'==1 noi _dots 0, title("Performing bootstrap repetitions") reps(`reps')
			
			// product-Exp weights
			tempvar bootwt	
			qui gen double `bootwt' = 1 if `touse'
			if ("`vce'" != "") {
				local idx = 0
				foreach cv of varlist `vce' {
					local ++idx
					tempvar gw`idx'
					bysort `: word `idx' of `clids'' : gen double `gw`idx'' = cond(_n==1, max(1e-5,-ln(runiform())), .) if `touse'
					by `: word `idx' of `clids'' : replace `gw`idx'' = `gw`idx''[1] if `touse'
					qui replace `bootwt' = `bootwt' * `gw`idx'' if `touse'
					drop `gw`idx''
				}
			}
			else {
				qui replace `bootwt' = -ln(runiform()) if `touse'
			}
			tempvar wrep
			qui gen double `wrep' = `bootwt' * `userwt' if `touse'
			
			
			//CHECK FOR COLLINEARITY
			// --- compute Kish effective N per absorb group and drop ~singletons ---
			// Kish effective N = (Σw)^2 / Σw^2 by joint FE cell
			/*
			if (`K'>0) {
				tempvar w2 Neff sumw sumw2
				gen double `w2' = `wrep'^2 if `touse'

				bys `gidlist': egen double `sumw'  = total(`wrep')   if `touse'
				bys `gidlist': egen double `sumw2' = total(`w2')     if `touse'
				gen double `Neff' = (`sumw'^2) / `sumw2'            if `touse'

				// drop cells with Neff <= 1 + tol
				local tol_neff = 1e-6
				replace `ok' = (`Neff' > 1 + `tol_neff') if `touse'
			} 
			else replace `ok'=1
			*/
			
			
			// Refit engine with replication weights; pass through noppmlfrom + aux()
			capture belieffit_engine `Ylog' `Xvars' if `touse' /*&`ok'*/ [aw=`wrep'], ///
				absorb(`absorb') tol(`tolerance') maxiter(`maxiter') ///
				feprefix(`fe') `ppmlfromopt' aux(`aux')
			
			if (_rc) {
				continue
			}
			
			loc ++r
			mat `b_bs'=nullmat(`b_bs')\e(b)
			   
			   
			   // == per-row Welford updates of covariance across reps (lower triangle) ==

			if (`K' > 0) {

				// STEP 0: increment per-row sample size for rows we're updating this rep
				// right now ok is always ==1 because you commented Kish. If you reintroduce Kish,
				// keep the &`ok' condition below.
				replace `R_row' = `R_row' + 1 if `touse' /* & `ok' */

				// STEP 1: compute deltas using the OLD means (M# before updating them this rep)
				// dj = x_j - M_j(old)
				forvalues j = 1/`K' {
					tempvar dj`j'
					gen double `dj`j'' = (`fe'`j' - `M`j'') if `touse' /* & `ok' */
				}

				// STEP 2: update means to NEW means using the NEW count R_row
				// M_j(new) = M_j(old) + dj / R_row
				forvalues j = 1/`K' {
					replace `M`j'' = `M`j'' + `dj`j'' / `R_row' if `touse' /* & `ok' */
				}

				// STEP 3: update cross-product accumulators
				// SS_{j,k} += dj_j * (x_k - M_k(new))
				forvalues j = 1/`K' {
					forvalues k = 1/`j' {
						tempvar devk
						gen double `devk' = (`fe'`k' - `M`k'') if `touse' /* & `ok' */
						replace `SS`j'`k'' = `SS`j'`k'' + `dj`j'' * `devk' ///
							if `touse' /* & `ok' */
						drop `devk'
					}
				}

				// STEP 4: produce covariance estimates per row
				forvalues j = 1/`K' {
					forvalues k = 1/`j' {
						if ("`df'"=="adj") {
							replace `varprefix'`j'`k' = ///
								cond(`R_row'>1, `SS`j'`k'' / (`R_row'-1), 0) ///
								if `touse'
						}
						else {
							replace `varprefix'`j'`k' = ///
								cond(`R_row'>0, `SS`j'`k'' / `R_row', 0) ///
								if `touse'
						}
					}
				}

				// STEP 5: clean up temp deltas
				forvalues j = 1/`K' {
					drop `dj`j''
				}
			}


			noi _dots `r' 0
		}
		

		// ---- Bootstrap coef covariance & mean ----
		mata: st_matrix("`V_bs'",variance(st_matrix("`b_bs'")))
		local names: colfullnames `b_bs'
		mat colnames `V_bs'=`names'
		mat rownames `V_bs'=`names'
		
		// ---- Post results ----	
		ereturn clear
		ereturn post `b_base' `V_bs'
		ereturn scalar reps = `reps'
		ereturn local cmd "belieffit"
		ereturn local absorb "`absorb'"
		ereturn local feprefix "`feprefix'"
		ereturn local varprefix "`varprefix'"
		ereturn local regressors "`Xvars'"
		ereturn local depvar "`Ylog'"
		ereturn local vcetype "Bootstrap"
		}
		ereturn display
	end

	program define belieffit_engine, eclass sortpreserve
		version 17.0
		// Usage:
		// belieffit_engine ylog Xvars [if] [in] [weights], absorb(fevars)
		//     [ tol(#) maxiter(#) feprefix(name) savefe(final|none) noppmlfrom aux(string asis) ]
		syntax varlist(min=1) [if] [in] [aw fw iw pw], ///
			[ ABSORB(varlist fv) TOLerance(real 1e-6) MAXITer(integer 100) ///
			  FEPREFIX(name)  AUX(string asis) prednames(string) calc(string)]

		marksample touse
		gettoken Ylog Xvars : varlist
		
		if ("`weight'" != "") local W "[`weight'`exp']"   // e.g. [aw=_wrep] or [pw=_wrep]

		// defaults
		if ("`feprefix'" == "") local feprefix "fe"

		// weights passthrough
		if ("`weight'" != "") local wopt "[`weight'`exp']"

		// deps
		 capture which reghdfe
			if (_rc) {
				di as err "reghdfe required: ssc install reghdfe"
				exit 499
			}
		if ("`absorb'" == "") loc absorbstr noabsorb
		else loc absorbstr absorb(`absorb')

		capture which ppmlhdfe
		if (_rc) {
			di as err "ppmlhdfe required: ssc install ppmlhdfe"
			exit 499
		}

		// FE dimensions (K)
		local Kfe = 0
		if ("`absorb'" != "") {
			capture unab _ablist : `absorb'
			local Kfe : word count `_ablist'
		}

		// working vars
		tempvar tau2 mhat varhat res res2
		quietly gen double `tau2'   = 0 if `touse'
		quietly gen double `mhat'   = . if `touse'
		quietly gen double `varhat' = . if `touse'

		// clean possible leftovers
		forvalues k = 1/`=2*`Kfe'' {
			capture drop `feprefix'`k'
		}

		// loop controls
		local tol   = `tolerance'
		local maxit = `maxiter'
		local iter  = 0
		local maxchg = 1

		// keep last ppmlhdfe coefficients for from()
		tempname b_mu b_var b_ppml
		local have_ppml 0


		// ---------- outer fixed-point loop ----------
		while (`iter' < `maxit' & `maxchg' > `tol') {
			local ++iter

			// (a) adjust log posterior means
			quietly replace `mhat' = `Ylog' - 0.5*`tau2' if `touse'
				
			// (b) mean equation

			quietly reghdfe `mhat' `Xvars' [aw`exp'] if `touse', `absorbstr' residuals(`res')

			// (c) variance equation (PPML on squared residuals)
			quietly gen double `res2' = max(`res'^2, 1e-8) if `touse'

			// build ppml command (with from() unless noppmlfrom or first call)
		   quietly  ppmlhdfe `res2' `Xvars' [pw`exp'] if `touse', `absorbstr' tol(1e-6) d

			tempvar muhat
			quietly predict double `muhat' if e(sample), mu
			quietly replace `varhat' = `muhat' if `touse'

			// (d) convergence
			tempvar diff
			quietly gen double `diff' = abs(`tau2' - `varhat') if `touse'
			quietly summarize `diff' if `touse', meanonly
			local maxchg = r(max)
			quietly replace `tau2' = `varhat' if `touse'
			drop `diff' `res' `res2'
			
		}

		// (e) save FE from final iteration if requested + predictions
		
		if "`absorb'"!="" loc absorbstr absorb(`absorb', savefe)
		
		// Mean: pred_mean`
		   quietly reghdfe `mhat' `Xvars' [aw`exp'] if `touse', `absorbstr' residuals(`res')
			matrix `b_mu' = e(b)
			// (c) variance equation (PPML on squared residuals)
			quietly gen double `res2' = max(`res'^2, 1e-8) if `touse'

		   if "`absorb'"!="" {
			local idx = 0
			foreach v of varlist __hdfe*  {
				local ++idx
				capture drop `feprefix'`idx'
				rename `v' `feprefix'`idx'
			}
			}

		// Variance final: pred_var (+ FE if requested)
		quietly ppmlhdfe  `res2' `Xvars' [pw`exp'] if `touse',  `absorbstr' tol(1e-6) d 
		matrix `b_var' = e(b)

		if "`absorb'"!="" {
			foreach v of varlist __hdfe*  {
				local ++idx
				capture drop `feprefix'`idx'
				quietly rename `v' `feprefix'`idx'
			}
		}
		
		if "`prednames'"!="" {
			tokenize `prednames'
			predict double `2' if `touse', mu
			gen double `1'=`Ylog'-0.5*`2' if `touse'
		}
		
		drop `res' `res2'

		// (f) assemble main b
		tempname b
		matrix coleq `b_mu'  = mean
		matrix coleq `b_var' = variance
		matrix `b' = `b_mu', `b_var'

		// g) Calc - calculate treatment intensities
		if "`calc'"!="" {
			local calc = subinstr("`calc'", "{W}", "`W'", .)
			local calc = subinstr("`calc'", "{mean}", "`mhat'", .)
			local calc = subinstr("`calc'", "{var}", "`varhat'", .)
			`calc'
			noi su M I S
			pause
			}
		
		// (h) AUX: run extra commands and append their e(b)
		// aux() may contain multiple blocks separated by ';'
		// Each block may optionally start with: name(eqname) <space> <command>
		local auxspec `"`aux'"'
		local auxcount = 0
		if ("`auxspec'" != "") {
			// consume blocks until auxspec is empty
			while (length("`auxspec'") > 0) {
				// split on first ';'
				local pos = strpos(`"`auxspec'"', ";")
				if (`pos'==0) {
					local this `"`auxspec'"'
					local auxspec ""
				}
				else {
					local this    = substr(`"`auxspec'"', 1, `pos'-1)
					local auxspec = substr(`"`auxspec'"', `pos'+1, .)
				}
				local this = trim(`"`this'"')
				if ("`this'"=="") continue

				// optional leading name()
				local eqname ""
				if (substr(`"`this'"',1,5)=="name(") {
					local endn = strpos(`"`this'"', ")")
					if (`endn'>5) {
						local eqname = substr(`"`this'"', 6, `endn'-6)
						// drop leading name(...) and any following space
						local this = trim(substr(`"`this'"', `endn'+1, .))
					}
				}

				// run block verbatim; must produce e(b)
				local ++auxcount
				local this = subinstr("`this'", "{W}", "`W'", .)
				local this = subinstr("`this'", "{mean}", "`mhat'", .)
				local this = subinstr("`this'", "{var}", "`varhat'", .)
				`this'
				tempname b_aux
				matrix `b_aux' = e(b)
				if ("`eqname'"=="") local eqname = "aux`auxcount'"
				matrix coleq `b_aux' = `eqname'
				matrix `b' = `b', `b_aux'
			}
		}
		

		
		// (h) post
		ereturn clear
		ereturn post `b'
		ereturn scalar iterations = `iter'
		ereturn scalar maxdiff    = `maxchg'
		ereturn local cmd "belieffit_engine"
		ereturn local depvar "`Ylog'"
		ereturn local absorb "`absorb'"
		ereturn local regressors "`Xvars'"
		ereturn local feprefix "`feprefix'"
		// leaves: pred_mean, pred_var, and optionally fe1..fe(2K) in data
		ereturn display
	end


