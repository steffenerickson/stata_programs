

// TO DO -> Need to estimate standard errors for just the effects 



version 18.0

// ---------------------------------------------------------------------------//
// Stata Code 
// ---------------------------------------------------------------------------//

program flowcheck, rclass 
	syntax, ESTimates(string) z(varlist max = 1)[level1(integer 1) level2(integer 2)]
	
	python: from __main__ import SymbolicEquations
	
	estimates restore `estimates'
	if  (`:word count `e(lyvars)'' == 0) {
		local option = 0 
	}
	else {
		local option = 1
	}
	collectsem, option(`option') z(`z')
	// Main Mata routine
	mata:  results = main(lxvars,lyvars,oyvars,params,params_real,glvls,cmdline,(`level1',`level2'))
	if  (`option' == 0) {	// Model with exogenous latent variables only 
		setupexogenousonly
		python: SymbolicEquations.wrapper(["Lambdax1","Lambdax2","kappa1","kappa2"],option=0) 	
		testnlroutine , eqxonxi(`eq0')
		tempname diff
		mat `diff' = xonxi - upsilon_x
		mat tablexonxi = upsilon_x , xonxi , `diff' , r(resultxonxi)
		
		mata return_rownames(results.lambdax_bag[1...,1])
		local varnames `disp_rownames'
		local rownames
		foreach x of local varnames {
			local rowname   "`x'"
			local rownames `" `rownames' "`rowname'" "'
		}
		
		mat colnames tablexonxi = "Obs" "MI" "MI - obs" "Pval" "Chi^2" "df"
		mat rownames tablexonxi = `rownames'
		matlist tablexonxi ,title("xonxi") lines(oneline)
		return matrix tablexonxi = tablexonxi
	}
	else { // Full SEM model 
		setupfullsem
		python: SymbolicEquations.wrapper(["Lambdax1","Lambdax2","Lambday1","Lambday2","kappa1","kappa2","Gamma1","Gamma2","Beta1","Beta2","alpha1","alpha2"],option=1) 
		testnlroutine , eqxonxi(`eq0') eqyonxi(`eq1') eqyoneta(`eq2') eqyonxieta(`eq3') 
		local names xonxi yonxi yoneta yonxieta
		forvalues i = 1/4 {
			tempname diff 
			local name `:word `i' of `names''
			if `i' == 1 local observed upsilon_x
			else local observed upsilon_y
			mat `diff' = `name' - `observed'
			mat table`name' = `observed' , `name' ,`diff' , r(result`name')
			
			if `i' == 1 mata return_rownames(results.lambdax_bag[1...,1])
			else mata return_rownames(results.lambday_bag[1...,1])
			local varnames `disp_rownames'
			local rownames
			foreach x of local varnames {
				local rowname   "`x'"
				local rownames `" `rownames' "`rowname'" "'
			}
			
			mat colnames table`name' = "Obs" "MI" "MI - obs" "Pval" "Chi^2" "df"
			mat rownames table`name' = `rownames'
			matlist table`name' , title("`name'") lines(oneline)
			return matrix table`name' = table`name'	
		}
	}
end 

program collectsem, nclass 
	syntax, option(real) z(varlist max=1)
	// Command line and parameter matrices (string and real)
	local cmdline = e(cmdline)
	mata  cmdline = st_local("cmdline")
	local colnames : colfullnames e(b)
	mata  params = tokens(st_local("colnames"))'
	matrix params_real = e(b)'
	mata  params_real = st_matrix("params_real")
	
	// List of observed and latent variables
	local lxvars = e(lxvars)
	mata  lxvars = tokens(st_local("lxvars")) 
	local oyvars = e(oyvars)
	mata  oyvars = tokens(st_local("oyvars"))
	// Set Beta to missing if odel contains only exogenous latent variables. This option controls processing in the mata code
	if  (`option') == 0 {
		mata lyvars = J(1,0,"")
	}
	else {
		local lyvars = e(lyvars)
		mata  lyvars = tokens(st_local("lyvars"))
	}
	// Group levels
	levelsforsem `z', local(glvls)
	local glvls = r(glvls)
	mata  glvls = tokens(st_local("glvls")) 
end 

program levelsforsem, rclass  
	syntax varlist(max=1), Local(string)
	qui levelsof `varlist', local(temp)
	forvalues i = 1/`:word count `temp'' {
		local lvl = `:word `i' of `temp''
		if (`i' == 1) local new `lvl'bn.`varlist'
		else local new `lvl'.`varlist'     
		local `local' : list `local' | new 
	}
	return local `local'  ``local''
end

// These two programs set symbolic mata matrices so they can be called by the python class "SymbolicEquations"
// They also call mata matrices into stata that are used for the "testnlroutine" routine 

program setupexogenousonly, nclass 
	mata  upsilon_x_string = asarray(results.intercepts_string,"upsilon_x")
	mata: upsilon_x_real = asarray(results.observed_effects,"upsilon_x")
	mata: xonxi_real = asarray(results.mi_effects ,"mi_effects_xonxi")
	mata: Lambdax1 = asarray(results.parammatrices_string,"Lambdax1")
	mata: Lambdax2 = asarray(results.parammatrices_string,"Lambdax2")
	mata: kappa1 = asarray(results.parammatrices_string,"kappa")[1...,1]
	mata: kappa2 = asarray(results.parammatrices_string,"kappa")[1...,2]
	mata: st_matrix("upsilon_x",upsilon_x_real)
	mata: st_matrix("xonxi", xonxi_real)
end 

program setupfullsem, nclass
	mata  upsilon_x_string = asarray(results.intercepts_string,"upsilon_x")
	mata  upsilon_y_string = asarray(results.intercepts_string,"upsilon_y")
	mata: upsilon_x_real = asarray(results.observed_effects ,"upsilon_x")
	mata: upsilon_y_real = asarray(results.observed_effects ,"upsilon_y")
	mata: xonxi_real = asarray(results.mi_effects ,"mi_effects_xonxi")
	mata: yonxi_real = asarray(results.mi_effects ,"mi_effects_yonxi")
	mata: yoneta_real = asarray(results.mi_effects ,"mi_effects_yoneta")
	mata: yonxieta_real = asarray(results.mi_effects ,"mi_effects_yonxieta")
	mata: Beta1 = asarray(results.parammatrices_string,"Beta1")
	mata: Beta2 = asarray(results.parammatrices_string,"Beta2")
	mata: Gamma1 = asarray(results.parammatrices_string,"Gamma1")
	mata: Gamma2 = asarray(results.parammatrices_string,"Gamma2")
	mata: Lambdax1 = asarray(results.parammatrices_string,"Lambdax1")
	mata: Lambdax2 = asarray(results.parammatrices_string,"Lambdax2")
	mata: Lambday1 = asarray(results.parammatrices_string,"Lambday1")
	mata: Lambday2 = asarray(results.parammatrices_string,"Lambday2")
	mata: kappa1 = asarray(results.parammatrices_string,"kappa")[1...,1]
	mata: kappa2 = asarray(results.parammatrices_string,"kappa")[1...,2]
	mata: alpha1 = asarray(results.parammatrices_string,"alpha")[1...,1]
	mata: alpha2 = asarray(results.parammatrices_string,"alpha")[1...,2]
	mata {
		if (Beta1 == J(0,0,.)) Beta1 = Beta2 = J(cols(Lambday1),cols(Lambday1),0)
	}
	mata: st_matrix("upsilon_x",upsilon_x_real)
	mata: st_matrix("upsilon_y",upsilon_y_real)
	mata: st_matrix("xonxi", xonxi_real)
	mata: st_matrix("yonxi", yonxi_real)
	mata: st_matrix("yoneta", yoneta_real)
	mata: st_matrix("yonxieta", yonxieta_real)
		
end 

program testnlroutine, rclass 
	syntax, eqxonxi(string) [eqyonxi(string) eqyoneta(string) eqyonxieta(string)]
	tempname tempmat resultxonxi resultyonxi resultyoneta resultyonxieta
	
	mata: obseqs = upsilon_x_string[1...,2] :+ "-" :+ upsilon_x_string[1...,1]
	mata: mieqs = ustrsplit(st_local("eqxonxi"),",")'
	mata: st_local("numrows",strofreal(rows(mieqs)))
	forvalues i = 1/`numrows' {
		mata st_local("mieq`i'",mieqs[`i'])
		mata st_local("obseq`i'",obseqs[`i'])
		qui testnl `obseq`i'' = `mieq`i''
		mat `tempmat' = r(p),r(chi2),r(df)
		mat `resultxonxi ' = (nullmat(`resultxonxi') \ `tempmat')
	}
	return matrix resultxonxi  = `resultxonxi' 
	if "`eqyonxi'" != "" {		
		// eqyonxi
		mata: obseqs = upsilon_y_string[1...,2] :+ "-" :+ upsilon_y_string[1...,1]
		mata: mieqs1 = ustrsplit(st_local("eqyonxi"),",")'
		mata: st_local("numrows",strofreal(rows(mieqs1)))
		forvalues i = 1/`numrows' {
			mata st_local("mieq`i'",mieqs1[`i'])
			mata st_local("obseq`i'",obseqs[`i'])
			qui testnl `obseq`i'' = `mieq`i''
			mat `tempmat' = r(p),r(chi2),r(df)
			mat `resultyonxi' = (nullmat(`resultyonxi') \ `tempmat')
		}
		return matrix resultyonxi  = `resultyonxi' 
		
		// eqyoneta
		mata: obseqs = upsilon_y_string[1...,2] :+ "-" :+ upsilon_y_string[1...,1]
		mata: mieqs2 = ustrsplit(st_local("eqyoneta"),",")'
		mata: st_local("numrows",strofreal(rows(mieqs2)))
		forvalues i = 1/`numrows' {
			mata st_local("mieq`i'",mieqs2[`i'])
			mata st_local("obseq`i'",obseqs[`i'])
			qui testnl `obseq`i'' = `mieq`i''
			mat `tempmat' = r(p),r(chi2),r(df)
			mat `resultyoneta' = (nullmat(`resultyoneta') \ `tempmat')
		}
		return matrix resultyoneta  = `resultyoneta' 

		// eqyonxieta
		mata: obseqs = upsilon_y_string[1...,2] :+ "-" :+ upsilon_y_string[1...,1]
		mata: mieqs3 = ustrsplit(st_local("eqyonxieta"),",")'
		mata: st_local("numrows",strofreal(rows(mieqs3)))
		forvalues i = 1/`numrows' {
			mata st_local("mieq`i'",mieqs3[`i'])
			mata st_local("obseq`i'",obseqs[`i'])
			qui testnl `obseq`i'' = `mieq`i''
			mat `tempmat' = r(p),r(chi2),r(df)
			mat `resultyonxieta' = (nullmat(`resultyonxieta') \ `tempmat')
		}
		return matrix resultyonxieta  = `resultyonxieta' 
	}
end 

// ---------------------------------------------------------------------------//
// Mata Code 
// ---------------------------------------------------------------------------//

// ---------------------------------------------------------------------------//
// Structures 
// ---------------------------------------------------------------------------//
mata
// Define structures for the SEM problem and derived results
struct myproblem  {
	string vector   lxvars           		// Latent exogenous variables
	string vector   lyvars           		// Latent endogenous variables
	string vector   oxvars           		// Observed exogenous variables
	string vector   oyvars           		// Observed endogenous variables
	string vector   params           		// Model parameters (strings)
	real   vector   params_real				// Model parameters (real)
	string vector   glvls            		// Group levels
	string scalar   cmdline          		// Command line (model specification)
	real   vector   comparegroups    		// Indices of groups to compare
	real   scalar   higherorder   	 		// Higher order factor model method option 
	string scalar   higherorderitem  		// Higher order factor model achor item 
	struct derived  scalar d         		// Nested struct for derived results
}

struct derived { 
	string matrix   pathcoefs        		// Selected path coefficients
	string matrix   intercepts       		// Selected intercept parameters
	string matrix   beta_bag         		// Endogenous-to-endogenous path coefs
	string matrix   gamma_bag        		// Exogenous-to-endogenous path coefs
	string matrix   lambdax_bag				// Exogenous loadings
	string matrix   lambday_bag      		// Endogenous loadings
	transmorphic    anchorkey       		// Anchors for latent variables
	transmorphic    parammatrices_string 	// Stores symbolic coefficients
	transmorphic    parammatrices_real 		// Stores real coefficients
	transmorphic    intercepts_string 		// Stores symbolic intercepts
	transmorphic    intercepts_real 		// Stores real intercepts
	transmorphic    observed_effects		// Observed Effects 
	transmorphic    mi_effects 				// Model implied effects
}

// Initialize the derived structure with empty matrices and a new anchorkey asarray
void initialize_objects(struct derived scalar d)
{
	d.pathcoefs = J(0,0,"") 
	d.intercepts = J(0,0,"") 
	d.anchorkey = asarray_create() 
	d.beta_bag = J(0,0,"") 
	d.gamma_bag = J(0,0,"") 
	d.lambdax_bag = J(0,0,"") 
	d.lambday_bag = J(0,0,"") 
	d.parammatrices_string = asarray_create()
	d.parammatrices_real = asarray_create()	
	d.intercepts_string = asarray_create()	
	d.intercepts_real = asarray_create()	
	d.observed_effects = asarray_create()	
	d.mi_effects = asarray_create()	
}

// ---------------------------------------------------------------------------//
// Main entry point: initialize and run the routines to process SEM results
// ---------------------------------------------------------------------------//

struct derived main(lxvars,lyvars,oyvars,params,params_real,glvls,cmdline,comparegroups)
{
    struct myproblem scalar pr
    initialize_objects(pr.d)
	
	// Store Inputs 
    // Store input arguments into the myproblem structure
    pr.lxvars  = lxvars
    pr.lyvars  = lyvars
    pr.oyvars  = oyvars
    pr.params  = params
	pr.params_real  = params_real
    pr.glvls   = glvls 
    pr.cmdline = cmdline 
    pr.comparegroups = comparegroups 
	
	// Sub-routines 
	get_coefficients(pr) 
	get_anchorkey(pr)
	group_pathcoefs(pr)
	group_intercepts(pr)
	get_latent_means(pr)
	get_Lambdax(pr)
	if (pr.lyvars != J(1,0,"")) {
		get_Lambday(pr)
		get_Beta(pr)
		get_Gamma(pr)
	}
	fill_pathcoef_matrices(pr)
	fill_intercepts(pr)
	observed_effects(pr)
	mi_effects_xonxi(pr)
	if (pr.lyvars != J(1,0,"")) {
		mi_effects_yonxi(pr)
		mi_effects_yoneta(pr)
		mi_effects_yonxieta(pr)
	}
	
	// return the derived structure
	return(pr.d)
}

// ---------------------------------------------------------------------------//
// Code to parse SEM output and organize parameters for later steps 
// ---------------------------------------------------------------------------//

// get_coefficients: Extracts path coefficients and intercepts from the raw params
void get_coefficients(struct myproblem scalar pr)
{
    real   matrix  s1, s2a, s2b, s3a, s3b
	string matrix  params2, params3
    
    // Remove means, variances, and covariances from params (we only want regression-like coefs)
    s1 = rowsum(strpos(pr.params,("var","mean","cov")))
    params2 = select(pr.params, opposite_mask(s1))
    
    // Separate coefficients by the specified group comparisons
    s2a  = strpos(params2,pr.glvls[pr.comparegroups[1]]) 
    s2b  = strpos(params2,pr.glvls[pr.comparegroups[2]]) 
    params3 = select(params2, s2a), select(params2, s2b)
    
    // Separate intercepts and path coefficients using presence of '#'
    s3a = rowsum(strpos(params3,"#"))
    s3b = opposite_mask(s3a)
    pr.d.pathcoefs  = select(params3, s3a)
    pr.d.intercepts = select(params3, s3b)
}
void get_anchorkey(struct myproblem scalar pr)
{
    real   matrix  s1
    real   scalar pos,i,j
    string scalar cmdline2,matched, anchor, latentvar 
    string matrix equations,equations2,check1,latentvars 
    
    // Extract portion of cmdline from "(" to "," to isolate model equations
    cmdline2 = substr(pr.cmdline,strpos(pr.cmdline,"("),strpos(pr.cmdline,",")-strpos(pr.cmdline,"("))
    equations = ustrsplit(cmdline2, "\)\s*\(")'
    equations = subinstr(subinstr(equations,"(",""),")","")
    
    // Check which equations relate to oyvars
    check1 = J(rows(equations),1,pr.oyvars)
    s1 = rowsum(strpos(equations, check1)) // Creating a mask
    equations2 = select(equations, s1)

    // For each equation, look for patterns like x1@1 indicating anchors
	latentvars = pr.lyvars,pr.lxvars
	
    for(i=1;i<=rows(equations2);i++) {
        if (regexm(equations2[i], "([A-Za-z0-9]+@1)")) {
            matched = regexs(1)
            pos = strpos(matched,"@") - 1 
            anchor = substr(matched,1,pos)
        }
        // Determine which latent variable this anchor belongs to
        for (j=1;j<=length(latentvars);j++) {
            if (strpos(equations2[i], latentvars[j]) > 0) {
                latentvar  = latentvars[j]
                asarray(pr.d.anchorkey,latentvar,anchor)
                break
            }	
        }		
    }
}
// group_pathcoefs: Separate structural coefficients (beta, gamma) from loadings (lambda)
void group_pathcoefs(struct myproblem scalar pr)
{
    real   matrix res,s1,s1a,s1b,s2,s2a,s2b
    real   scalar i 
    string matrix check1,structcoefs,loadings

	if (pr.lyvars == J(1,0,"")) { // used when there are only exogenous latent variables 
		pr.d.lambdax_bag = pr.d.pathcoefs
	}
	else {
		// Structural Coefficients 
		check1 = J(rows(pr.d.pathcoefs),1,pr.oyvars)
		s1 = opposite_mask(rowmissing(editvalue(indexnot(check1, pr.d.pathcoefs[1...,1]),0,.))) // need to get rid of this and replace with the method in get_intercepts , i think the code will break when you use nocapslatent method
		structcoefs = select(pr.d.pathcoefs,s1)
		res = J(rows(structcoefs),length(pr.lxvars),.)
		for (i=1;i<=length(pr.lxvars);i++) res[1...,i] = strpos(structcoefs[1...,1],pr.lxvars[i])
		s1a = rowsum(res)
		s1b = opposite_mask(s1a)
		pr.d.gamma_bag = select(structcoefs,s1a)
		pr.d.beta_bag = select(structcoefs,s1b)
		
		// Loadings
		s2 = opposite_mask(s1)
		loadings = select(pr.d.pathcoefs,s2)
		res = J(rows(loadings),length(pr.lyvars),.)
		for (i=1;i<=length(pr.lyvars);i++) res[1...,i] = strpos(loadings[1...,1],pr.lyvars[i])
		s2a = rowsum(res)
		s2b = opposite_mask(s2a)
		pr.d.lambday_bag = select(loadings,s2a)
		pr.d.lambdax_bag = select(loadings,s2b)

	}
}
// group_intercepts: Organize intercept parameters into y and x vectors 
void group_intercepts(struct myproblem scalar pr)
{
    real   matrix pos1,pos2,s1,s2,res
    string matrix ys,xs,upsilon_y,upsilon_x
	
	pos1 = strpos(pr.d.lambdax_bag[1...,1],":")
	xs = substr(pr.d.lambdax_bag[1...,1],J(rows(pr.d.lambdax_bag),1,1),pos1:-1) 
	res = J(rows(pr.d.intercepts),length(xs),.)
	for (i=1;i<=length(xs);i++) res[1...,i] = strpos(pr.d.intercepts[1...,1],xs[i])
	s1 = rowsum(res)
	upsilon_x = select(pr.d.intercepts,s1)
	upsilon_x = "_b[" :+ upsilon_x  :+ "]"
	asarray(pr.d.intercepts_string,"upsilon_x",upsilon_x)

	if (pr.lyvars != J(1,0,"")) {
		s2 = opposite_mask(s1)
		upsilon_y = select(pr.d.intercepts,s2)
		upsilon_y = "_b[" :+ upsilon_y  :+ "]"
		asarray(pr.d.intercepts_string,"upsilon_y",upsilon_y)
	}

}

// ---------------------------------------------------------------------------//
// Code to create symbolic parameter matrices to estimate total effects 
// ---------------------------------------------------------------------------//
void get_latent_means(struct myproblem scalar pr) 
{
	string matrix x_anchors,y_anchors,check,intercepts,kappa,alpha 
	real matrix s 
	string scalar key, anchor 
	real scalar i 
	
	x_anchors = J(1,length(pr.lxvars),"")
	for (i=1;i<=length(pr.lxvars);i++) {
		for (loc=asarray_first(pr.d.anchorkey); loc!=NULL; loc=asarray_next(pr.d.anchorkey, loc)) {
			key = asarray_key(pr.d.anchorkey, loc)
			anchor = asarray_contents(pr.d.anchorkey, loc)
			if (pr.lxvars[i] == key) x_anchors[i] = anchor 
		}
	}
	
	intercepts = asarray(pr.d.intercepts_string, "upsilon_x")
	check = J(rows(intercepts),1,x_anchors)
    s = rowsum(strpos(intercepts[1...,1], check)) 
    kappa = select(intercepts, s)
	asarray(pr.d.parammatrices_string,"kappa",kappa)
	
	if (pr.lyvars != J(1,0,"")) {
		y_anchors = J(1,length(pr.lyvars),"")
		for (i=1;i<=length(pr.lyvars);i++) {
			for (loc=asarray_first(pr.d.anchorkey); loc!=NULL; loc=asarray_next(pr.d.anchorkey, loc)) {
				key = asarray_key(pr.d.anchorkey, loc)
				anchor = asarray_contents(pr.d.anchorkey, loc)
				if (pr.lyvars[i] == key) y_anchors[i] = anchor 
			}
		}
		intercepts = asarray(pr.d.intercepts_string, "upsilon_y")
		check = J(rows(intercepts),1,y_anchors)
		s = rowsum(strpos(intercepts[1...,1], check)) 
		alpha = select(intercepts, s)
		asarray(pr.d.parammatrices_string,"alpha",alpha)
	}
}
void get_Lambdax(struct myproblem scalar pr) 
{
	real matrix   B, pos,s
	real scalar   i,j,k
	string matrix xobserved, Lambdax1,Lambdax2
	
	Lambdax1 = J(rows(pr.d.lambdax_bag),length(pr.lxvars),"")
	Lambdax2 = J(rows(pr.d.lambdax_bag),length(pr.lxvars),"")
	pos = strpos(pr.d.lambdax_bag[1...,1],":"):-1
	xobserved = substr(pr.d.lambdax_bag[1...,1],J(rows(pos),1,1),pos)
	for(i=1;i<=rows(Lambdax1);i++) {
		for(j=1;j<=cols(Lambdax1);j++) {
			for(k=1;k<=rows(pr.d.lambdax_bag);k++) {
				B = strpos(pr.d.lambdax_bag[k,1],(xobserved[i],pr.lxvars[j]))
				if (rowmissing(editvalue(B,0,.)) == 0 & B[1] < B[2]) {
					Lambdax1[i,j] = "_b[" + pr.d.lambdax_bag[k,1] + "]"
					Lambdax2[i,j] = "_b[" + pr.d.lambdax_bag[k,2] + "]"
				}
			}
		}	
	}
	// Remove duplicate rows caused by indicators having factor complexity > 1 
	s = J(rows(Lambdax1),1,1)
	for (i=2;i<=rows(Lambdax1);i++) {
		j = i - 1 
		if (invtokens(Lambdax1[i,1...])  == invtokens(Lambdax1[j,1...])) s[i] = 0
	}
	asarray(pr.d.parammatrices_string,"Lambdax1",select(Lambdax1,s))
	asarray(pr.d.parammatrices_string,"Lambdax2",select(Lambdax2,s))
}
void get_Lambday(struct myproblem scalar pr) 
{
	real matrix   B, pos,s
	real scalar   i,j,k
	string matrix yobserved,Lambday1,Lambday2
	
	Lambday1 = J(rows(pr.d.lambday_bag),length(pr.lyvars),"")
	Lambday2 = J(rows(pr.d.lambday_bag),length(pr.lyvars),"")
	pos = strpos(pr.d.lambday_bag[1...,1],":"):-1
	yobserved = substr(pr.d.lambday_bag[1...,1],J(rows(pos),1,1),pos)
	for(i=1;i<=rows(Lambday1);i++) {
		for(j=1;j<=cols(Lambday1);j++) {
			for(k=1;k<=rows(pr.d.lambday_bag);k++) {
				B = strpos(pr.d.lambday_bag[k,1],(yobserved[i],pr.lyvars[j]))
				if (rowmissing(editvalue(B,0,.)) == 0 & B[1] < B[2]) {
					Lambday1[i,j] = "_b[" + pr.d.lambday_bag[k,1] + "]"
					Lambday2[i,j] = "_b[" + pr.d.lambday_bag[k,2] + "]"
				}
			}
		}	
	}
	// Remove duplicate rows caused by indicators having factor complexity > 1 
	s = J(rows(Lambday1),1,1)
	for (i=2;i<=rows(Lambday1);i++) {
		j = i - 1 
		if (invtokens(Lambday1[i,1...])  == invtokens(Lambday1[j,1...])) s[i] = 0
	}
	asarray(pr.d.parammatrices_string,"Lambday1",select(Lambday1,s))
	asarray(pr.d.parammatrices_string,"Lambday2",select(Lambday2,s))
}
void get_Beta(struct myproblem scalar pr) 
{
	real matrix B 
	real scalar i,j,k
	string matrix Beta1,Beta2
	
	Beta1 = J(length(pr.lyvars),length(pr.lyvars),"")
	Beta2 = J(length(pr.lyvars),length(pr.lyvars),"")
	for(i=1;i<=rows(Beta1);i++) {
		for(j=1;j<=cols(Beta1);j++) {
			if (pr.lyvars[i] != pr.lyvars[j]) {
				for(k=1;k<=rows(pr.d.beta_bag);k++) {
					B = strpos(pr.d.beta_bag[k,1],(pr.lyvars[i],pr.lyvars[j]))
					if (rowmissing(editvalue(B,0,.)) == 0 & B[1] < B[2]) {
						Beta1[i,j] = "_b[" + pr.d.beta_bag[k,1] + "]"
						Beta2[i,j] = "_b[" + pr.d.beta_bag[k,2] + "]"
					}
				}
			}
		}	
	}
	asarray(pr.d.parammatrices_string,"Beta1",Beta1)
	asarray(pr.d.parammatrices_string,"Beta2",Beta2)
}
void get_Gamma(struct myproblem scalar pr) 
{
	real matrix B 
	real scalar i,j,k
	string matrix Gamma1,Gamma2
	
	Gamma1 = J(length(pr.lyvars),length(pr.lxvars),"")
	Gamma2 = J(length(pr.lyvars),length(pr.lxvars),"")
	for(i=1;i<=rows(Gamma1);i++) {
		for(j=1;j<=cols(Gamma1);j++) {
			for(k=1;k<=rows(pr.d.gamma_bag);k++) {
				B = strpos(pr.d.gamma_bag[k,1],(pr.lyvars[i],pr.lxvars[j]))
				if (rowmissing(editvalue(B,0,.)) == 0 & B[1] < B[2]) {
					Gamma1[i,j] = "_b[" + pr.d.gamma_bag[k,1] + "]"
					Gamma2[i,j] = "_b[" + pr.d.gamma_bag[k,2] + "]"
				}
			}
		}	
	}
	asarray(pr.d.parammatrices_string,"Gamma1",Gamma1)
	asarray(pr.d.parammatrices_string,"Gamma2",Gamma2)
}

// ---------------------------------------------------------------------------//
// Code to fill symbolic parameter matrices to estimate total effects 
// ---------------------------------------------------------------------------//
void fill_pathcoef_matrices(struct myproblem scalar pr)
{
	string matrix tempmat
	real matrix result
	string scalar key 
	real scalar i,j,k,tempval,loc

	for (loc=asarray_first(pr.d.parammatrices_string); loc!=NULL; loc=asarray_next(pr.d.parammatrices_string, loc)) {
		key = asarray_key(pr.d.parammatrices_string, loc)
		tempmat = asarray_contents(pr.d.parammatrices_string, loc)
		result = J(rows(tempmat),cols(tempmat),.)
		for(i=1;i<=rows(tempmat);i++) {
			for(j=1;j<=cols(tempmat);j++) {
				if (tempmat[i,j] == "") tempval = 0 
				else {
					for (k=1;k<=length(pr.params);k++) {
						if (tempmat[i,j] == "_b[" + pr.params[k] + "]") {
							tempval = pr.params_real[k]
							break 
						}
					}
				}
				result[i,j] = tempval
			}
		}
		asarray(pr.d.parammatrices_real,key,result)
	}
}
void fill_intercepts(struct myproblem scalar pr)
{
	string matrix tempmat
	real matrix result
	string scalar key 
	real scalar i,j,k,tempval,loc

	for (loc=asarray_first(pr.d.intercepts_string); loc!=NULL; loc=asarray_next(pr.d.intercepts_string, loc)) {
		key = asarray_key(pr.d.intercepts_string, loc)
		tempmat = asarray_contents(pr.d.intercepts_string, loc)
		result = J(rows(tempmat),cols(tempmat),.)
		for(i=1;i<=rows(tempmat);i++) {
			for(j=1;j<=cols(tempmat);j++) {
				if (tempmat[i,j] == "") tempval = 0 
				else {
					for (k=1;k<=length(pr.params);k++) {
						if (tempmat[i,j] == "_b[" + pr.params[k] + "]") {
							tempval = pr.params_real[k]
							break 
						}
					}
				}
				result[i,j] = tempval
			}
		}
		asarray(pr.d.intercepts_real,key,result)
	}
}

//----------------------------------------------------------------------------//
// Calculate Effects
//----------------------------------------------------------------------------//
void observed_effects(struct myproblem scalar pr)
{
	real matrix tempmat,result,effects 	
	string scalar key
	real scalar loc

	for (loc=asarray_first(pr.d.intercepts_real); loc!=NULL; loc=asarray_next(pr.d.intercepts_real, loc)) {
		key = asarray_key(pr.d.intercepts_string, loc)
		tempmat = asarray_contents(pr.d.intercepts_real, loc)
		result = tempmat[1...,2] - tempmat[1...,1]
		asarray(pr.d.observed_effects,key,result)
	}
}
void mi_effects_xonxi(struct myproblem scalar pr)
{
	real matrix kappa1,kappa2,Lambdax1,Lambdax2,result,effects 	
	
	 Lambdax1 = asarray(pr.d.parammatrices_real,"Lambdax1")
     Lambdax2 = asarray(pr.d.parammatrices_real,"Lambdax2")
	 kappa1 = asarray(pr.d.parammatrices_real,"kappa")[1...,1]
	 kappa2 = asarray(pr.d.parammatrices_real,"kappa")[1...,2]
	 
	 result = J(rows(Lambdax1),2,.)
	 result[1...,1] = Lambdax1*kappa1
	 result[1...,2] = Lambdax2*kappa2
	 effects = result[1...,2] - result[1...,1]
	 
	 asarray(pr.d.mi_effects,"mi_effects_xonxi",effects)
}
void mi_effects_yonxi(struct myproblem scalar pr)
{
	real matrix kappa1,kappa2,Lambday1,Lambday2,
	            Beta1,Beta2,Gamma1,Gamma2,result, effects 	
	Lambday1 = asarray(pr.d.parammatrices_real,"Lambday1")
	Lambday2 = asarray(pr.d.parammatrices_real,"Lambday2")	
	Beta1 = asarray(pr.d.parammatrices_real,"Beta1")
	Beta2 = asarray(pr.d.parammatrices_real,"Beta2")
	if (Beta1 == J(0,0,.)) Beta1 = Beta2 = J(cols(Lambday1),cols(Lambday1),0)
	Gamma1 = asarray(pr.d.parammatrices_real,"Gamma1")
	Gamma2 = asarray(pr.d.parammatrices_real,"Gamma2")
	kappa1 = asarray(pr.d.parammatrices_real,"kappa")[1...,1]
	kappa2 = asarray(pr.d.parammatrices_real,"kappa")[1...,2]
	
	result = J(rows(Lambday1),2,.)
	result[1...,1] = Lambday1 * luinv(I(rows(Beta1)) - Beta1) * Gamma1 * kappa1
	result[1...,2] = Lambday2 * luinv(I(rows(Beta2)) - Beta2) * Gamma2 * kappa2 
	effects = result[1...,2] - result[1...,1]
	
	asarray(pr.d.mi_effects,"mi_effects_yonxi",effects)
}
void mi_effects_yoneta(struct myproblem scalar pr)
{
	real matrix kappa1,kappa2, alpha1,alpha2,Lambday1,Lambday2,
	            Beta1,Beta2,Gamma1,Gamma2,result,effects 	 	
	Lambday1 = asarray(pr.d.parammatrices_real,"Lambday1")
	Lambday2 = asarray(pr.d.parammatrices_real,"Lambday2")	
	Beta1 = asarray(pr.d.parammatrices_real,"Beta1")
	Beta2 = asarray(pr.d.parammatrices_real,"Beta2") 
	if (Beta1 == J(0,0,.)) Beta1 = Beta2 = J(cols(Lambday1),cols(Lambday1),0)
	Gamma1 = asarray(pr.d.parammatrices_real,"Gamma1")
	Gamma2 = asarray(pr.d.parammatrices_real,"Gamma2")
	alpha1 = asarray(pr.d.parammatrices_real,"alpha")[1...,1]
	alpha2 = asarray(pr.d.parammatrices_real,"alpha")[1...,2]
	kappa1 = asarray(pr.d.parammatrices_real,"kappa")[1...,1]
	kappa2 = asarray(pr.d.parammatrices_real,"kappa")[1...,2]
	
	result = J(rows(Lambday1),2,.)
	result[1...,1] = (Lambday1 * alpha1) - (Lambday1 * luinv(I(rows(Beta1)) - Beta1) * Gamma1 * kappa1)
	result[1...,2] = (Lambday2 * alpha2) - (Lambday2 * luinv(I(rows(Beta2)) - Beta2) * Gamma2 * kappa2) 
	effects = result[1...,2] - result[1...,1]
	
	asarray(pr.d.mi_effects,"mi_effects_yoneta",effects)
}
void mi_effects_yonxieta(struct myproblem scalar pr)
{
	real matrix alpha1,alpha2,Lambday1,Lambday2,result, effects 	 	
	Lambday1 = asarray(pr.d.parammatrices_real,"Lambday1")
	Lambday2 = asarray(pr.d.parammatrices_real,"Lambday2")	
	alpha1 = asarray(pr.d.parammatrices_real,"alpha")[1...,1]
	alpha2 = asarray(pr.d.parammatrices_real,"alpha")[1...,2]

	result = J(rows(Lambday1),2,.)
	result[1...,1] = (Lambday1 * alpha1) 
	result[1...,2] = (Lambday2 * alpha2) 
	
	//result[1...,1] = Lambday1 * luinv(I(rows(Beta1)) - Beta1) * (alpha1 + Gamma1 * kappa1)
	//result[1...,2] = Lambday2 * luinv(I(rows(Beta2)) - Beta2) * (alpha2 + Gamma2 * kappa2)
	
	effects = result[1...,2] - result[1...,1]
	
	asarray(pr.d.mi_effects,"mi_effects_yonxieta",effects)
}

// Helper functions 
real vector opposite_mask(real vector mask) return(mm_cond(mask:>0,mask:-mask,mask:+1))

void return_rownames(string matrix strmat)
{
	real   matrix pos1, pos2,s
	string matrix sub1, sub2, names
	
	pos1 = strpos(strmat,":")
	pos2 = strrpos(strmat,".")
	sub1 = substr(strmat,1,pos1:-1)
	sub2 = substr(strmat,pos2:+1,.)
	names = sub2 :+ ":" :+ sub1
	s = J(rows(sub1),1,1)
	for (i=2;i<=rows(sub1);i++) {
		j = i - 1 
		if (invtokens(sub1[i,1...])  == invtokens(sub1[j,1...])) s[i] = 0
	}
	st_local("disp_rownames",invtokens(select(names,s)'))
}

end 


//----------------------------------------------------------------------------//
// Python code 
//----------------------------------------------------------------------------//
python 

class SymbolicEquations:
    import numpy as np
    from sympy import Symbol, Matrix, eye
    from sfi import Mata, Macro

    @classmethod
    def wrapper(cls, matnames, option=0):
        dict_matrices = cls.create_symbolic_matrices(matnames)

        if option == 0:
            eq0 = cls.equations_xonxi(dict_matrices)
            eqlist = [eq0]
        else:
            eq0 = cls.equations_xonxi(dict_matrices)
            eq1 = cls.equations_yonxi(dict_matrices)
            eq2 = cls.equations_yoneta(dict_matrices)
            eq3 = cls.equations_yonxieta(dict_matrices)
            eqlist = [eq0, eq1, eq2, eq3]

        cls.return_as_macro(eqlist)

    @classmethod
    def create_symbolic_matrices(cls, matnames):
        matrices = {}
        for x in range(len(matnames)):
            matrix = cls.np.matrix(cls.Mata.get(f"{matnames[x]}"))
            variables = {}
            for index, element in cls.np.ndenumerate(matrix):
                var_name = f"{matnames[x]}{index[0]}_{index[1]}"
                if element == "":
                    variables[var_name] = 0
                else:
                    variables[var_name] = cls.Symbol(f"{element}")

            rows, cols = matrix.shape
            matrices[f"{matnames[x]}"] = cls.Matrix(
                rows, cols, lambda i, j: variables[f"{matnames[x]}{i}_{j}"]
            )
        return matrices

    @classmethod
    def equations_xonxi(cls, matrices):
        eq = (
            matrices["Lambdax2"] * matrices["kappa2"]
            - matrices["Lambdax1"] * matrices["kappa1"]
        )
        return eq

    @classmethod
    def equations_yonxi(cls, matrices):
        I = cls.eye(matrices["Beta1"].shape[0])
        eq = (
            matrices["Lambday2"]
            * (I - matrices["Beta2"]).inv()
            * matrices["Gamma2"]
            * matrices["kappa2"]
            - matrices["Lambday1"]
            * (I - matrices["Beta1"]).inv()
            * matrices["Gamma1"]
            * matrices["kappa1"]
        )
        return eq

    @classmethod
    def equations_yoneta(cls, matrices):
        I = cls.eye(matrices["Beta1"].shape[0])
        eq = (
            (matrices["Lambday2"] * matrices["alpha2"])
            - (
                matrices["Lambday2"]
                * (I - matrices["Beta2"]).inv()
                * matrices["Gamma2"]
                * matrices["kappa2"]
            )
            - (
                (matrices["Lambday1"] * matrices["alpha1"])
                - (
                    matrices["Lambday1"]
                    * (I - matrices["Beta1"]).inv()
                    * matrices["Gamma1"]
                    * matrices["kappa1"]
                )
            )
        )
        return eq

    @classmethod
    def equations_yonxieta(cls, matrices):
        eq = (
            matrices["Lambday2"] * matrices["alpha2"]
            - matrices["Lambday1"] * matrices["alpha1"]
        )
        return eq

    @classmethod
    def return_as_macro(cls, eqlist):
        for i, eq in enumerate(eqlist):
            nested_list = eq.tolist()
            string = ", ".join(
                map(str, [item for sublist in nested_list for item in sublist])
            )
            cls.Macro.setLocal(f"eq{i}", string)
end 

