
* ---------------------------------------------------------------------------- * 
*! Multivariate Generalizability Study Decompostion Command
*! With some post estimation-mata functions that allow you to 
*! calculate confidence intervals and some G-theory coefficients (like relative profile variance)


//----------------------------------------------------------------------------//
// Main Stata Routine
//----------------------------------------------------------------------------//

clear all 

program mvgstudy, rclass
	syntax anything [if] [in] 
	
	marksample touse 
	tempname tempframe 
	
	// 1: Set up 
	
	*Initialize mvgstudy() class in mata 
	mata c = mvgstudy()
	* Parse syntax expression
	parse_equation `anything'
	local tempvarlist `r(vars)'
	local effects `r(effects)'
	foreach var of varlist `tempvarlist' {
		local varlist : list varlist | var
	}
	* Lengths for Loops
	local lengthdepvars	   :  list sizeof local(varlist)
	local lengtheffectlist :  list sizeof local(effects)
	* Facets
	mata facets = uniqrows(tokens(subinstr(subinstr(invtokens(st_local("effects")),"|"," "), "#"," "))')
	mata facets = invtokens(facets')
	mata st_local("facets",facets)
	* Facet Levels
	capture matrix drop facetlevels 
	foreach facet of local facets {
		qui levelsof `facet' if `touse' == 1
		matrix facetlevels = (nullmat(facetlevels) \ `r(r)')
	}
	* Find Residual 
	mata residual = c.sort_desc_length(tokens(st_local("effects")))[1]
	mata st_local("residual",residual)
	local effects2 : list effects - residual
	
	// 2: Manova and Collect Results
	* Run manova
	frame put `varlist' `facets' if `touse' == 1, into(`tempframe')
	frame `tempframe' {
		qui manova `varlist' = `effects2'
	}
	estimates store manova_estimates
	recovermanova, residual(`residual')
	mata sscplist = asarray_create() 
	forvalues i = 1/`r(j)' {
		mata asarray(sscplist ,"`i'", st_matrix("r(sscp`i')"))
	}
	mata df = st_matrix("r(df)")
	mata effects = tokens(st_local("effects"))
	mata facets = tokens(st_local("facets"))
	mata facetlevels = st_matrix("facetlevels")
	mata c.init_inputs(sscplist,effects,facets,facetlevels,df)
	mata for (loc=asarray_first(c.sscplist); loc!=NULL; loc=asarray_next(c.sscplist, loc)) asarray_contents(c.sscplist, loc)
	mata mata drop sscplist effects facets facetlevels df
	

	// 3: Mata MvGstudy Main Routine 
	mata c.main_routine()
	mata for (loc=asarray_first(c.covcomps); loc!=NULL; loc=asarray_next(c.covcomps, loc)) st_matrix(asarray_key(c.covcomps, loc),asarray_contents(c.covcomps, loc))
	mata st_matrix("df",c.df)
	mata st_matrix("P",c.P)
	* Display results 
	foreach x of local varlist {
		local name   "`x'"
		local names `" `names' "`name'" "'
	}
	forvalues i = 1/`lengtheffectlist' {
		matrix rownames emcp`i' = `names'
		matrix colnames emcp`i' = `names'
		matlist emcp`i' , twidth(20) title("`:word `i' of `effects'' Component Covariance Matrix")
	}
	*Return results 
	forvalues i = 1/`lengtheffectlist' { 
		return matrix emcp`i' = emcp`i'
	}
	return matrix P = P 
	return matrix df = df 
	
end 

//----------------------------------------------------------------------------//
// Stata Sub Routines
//----------------------------------------------------------------------------//

program parse_equation, rclass
	syntax anything 
	
	tempname str pos i 
	mata `str' = st_local("anything")
	mata `str' = tokens(subinstr(subinstr(`str',"(",""),")",""))
	mata `pos' = strpos(`str',"=")
	mata `i' = 1
	mata while (`pos'[`i'] < 1) ++`i'
	mata st_local("vars", invtokens(`str'[1,1..`i' - 1]))
	mata st_local("effects",invtokens(`str'[1,`i' + 1..cols(`str')]))
	return local vars `vars' 
	return local effects `effects'
end

program recovermanova, rclass 
	syntax,  RESIDUAL(string)
	
	* Count the number of effects
	local a = e(cmdline) 
	local check = strrpos("`a'","=")
	local newlist = substr("`a'",`check'+1,.)
	local length : list sizeof local(newlist)
	local j = `length' + 1
	forvalues i = 1/`j' {
		local name sscp`i'
		local sscp_tempmats : list sscp_tempmats | name
	}
	tempname tempmat df `sscp_tempmats'
	
	* Degrees of freedom 
	forvalues i = 1/`length' {
		matrix `tempmat' = `e(df_`i')'
		matrix `df' = (nullmat(`df') \ `tempmat')
	}
	matrix `df' = `df' \ e(df_r)
	
	* Sums of squared cross product matrices 
	forvalues i = 1/`length' {
		matrix `sscp`i'' = e(H_`i')
		local++i
	}
	matrix `sscp`j'' = e(E)
	
	* Return 
	return matrix df = `df'
	return scalar j = `j'
	forvalues i = 1/`j' {
		return matrix sscp`i' = `sscp`i''
	}
end 

//----------------------------------------------------------------------------//
// Mata Code
//----------------------------------------------------------------------------//

	
mata

// Class Definition 

class mvgstudy 
{
	//----------------------------------------------------//
	public:
	
	void init_inputs()
	void main_routine()
	string matrix sort_desc_length()
	string matrix sort_asc_length()
	real matrix satterthwaitecitable()
	
	string vector   effects,facets
	real vector     facetlevels,flproducts, df
	transmorphic    sscplist,mscplist,covcomps,citables
	real matrix     P, compsstacked

	//----------------------------------------------------//
	private:
	
	// covariance comp estimation 
	void get_facetlevelproducts()
	void get_p()
	void get_mscplist()
	void get_compsstacked() 
	void emcp()
	void emcpmatrixprocedure()
	void emcpscalarprocedure()
	
	// post estimation
	real matrix satterthwaite()
	real scalar relativevariance()
	real scalar profilevar()
	real scalar errorprofilevar()
	
}

// ---------------------------------------------------------------------------//
//Public Routines
// ---------------------------------------------------------------------------//

void mvgstudy::init_inputs(transmorphic sscplist,
                           string vector effects,
						   string vector facets,
						   real vector facetlevels,
						   real df)
{
	this.sscplist = sscplist
	this.effects = effects
	this.facets = facets
	this.facetlevels = facetlevels
	this.df = df
}

void mvgstudy::main_routine()
{
	get_facetlevelproducts()
	get_p()
	get_mscplist()
	get_compsstacked() 
	emcp()
}
string matrix mvgstudy:: sort_desc_length(string vector M)
{
	real   scalar 		i,j
	string scalar		temp 
	string matrix 		res
	
	res = M 
	for (i=2;i<=length(M);i++){
		temp = res[i]
		j = i - 1
		while (j >= 1 & strlen(temp) > strlen(res[j])) {
			res[j + 1] = res[j]
			j--
			if (j < 1) break 
		}
		res[j+1] = temp 
	}
	return(res)
}

string matrix mvgstudy:: sort_asc_length(string vector M)
{
	real   scalar 		i,j
	string scalar		temp 
	string matrix 		res
	
	res = M 
	for (i=2;i<=length(M);i++){
		temp = res[i]
		j = i - 1
		while (j >= 1 & strlen(temp) < strlen(res[j])) {
			res[j + 1] = res[j]
			j--
			if (j < 1) break 
		}
		res[j+1] = temp 
	}
	return(res)
}

// ---------------------------------------------------------------------------//
//Private Routines
// ---------------------------------------------------------------------------//


// (1) covariance component calculation routine 

void mvgstudy::get_facetlevelproducts()
{
	real scalar 		flp,x,i,j

	flproducts = J(length(effects),1,.)
	for(i=1;i<=length(effects);i++) {
		flp = 1 
		for(j=1;j<=length(facets);j++) {
			if (strpos(effects[i],facets[j]) == 0) x = facetlevels[j]
			else x = 1
			flp = flp * x 
		}
		flproducts[i] = flp 
	}
}
void mvgstudy::get_p()
{
	string rowvector 	effects_stripped
	real scalar 		row, col 

	
	effects_stripped = tokens(subinstr(subinstr(invtokens(effects),"|",""), "#",""))
	P = J(length(effects_stripped),length(effects_stripped),0)
	for(col=1;col<=length(effects_stripped);col++){
		row = col
		while (row >= 1) {
			P[row, col] = flproducts[col]
			--row
		}
	}
	for(row=1;row<=length(effects_stripped);row++) {
		for(col=1;col<=length(effects_stripped)-1;col++) {
			if (strpos(effects_stripped[col],effects_stripped[row])== 0) P[row,col] = 0 
		}
	}
}
transmorphic matrix emcpmatrixprocedure(real matrix stacked, real matrix p, string matrix names)
{
	transmorphic matrix 	emcp
	real scalar 			skip 
	real matrix 			rangemat, starting, ending, inv_p,temp,res
	real scalar 			i, j
	
	emcp = asarray_create()   
	skip = cols(stacked)
	rangemat = range(0,rows(stacked),skip)
	starting = rangemat[1..(rows(rangemat)-1)]:+1
	ending = rangemat[2..rows(rangemat)]
	inv_p = luinv(p)
	for(i=1;i<=cols(p);i++){
		res = J(cols(stacked),cols(stacked),0)
		for(j=1;j<=cols(p);j++) {
			temp = stacked[starting[j]..ending[j],1..cols(stacked)]
			res = res + (temp * inv_p[i,j])
		}
		asarray(emcp,names[i],res)
	}
	return(emcp)
}
void mvgstudy::get_mscplist()
{
	
	real matrix temp

	mscplist = asarray_create() 
	for(i=1;i<=length(effects);i++) {
		temp = asarray(sscplist, strofreal(i)):/df[i]
		asarray(mscplist ,strofreal(i), temp)
	}
}
void mvgstudy::get_compsstacked() 
{
	
	real scalar cols 
	
	cols = cols(asarray(mscplist,"1"))
	compsstacked = J(0,cols,.)
	for(i=1;i<=length(effects);i++) {
		compsstacked = compsstacked \ asarray(mscplist, strofreal(i))
		
	}
}
void mvgstudy::emcp()
{
	if (cols(compsstacked) == 1) emcpscalarprocedure()
	else emcpmatrixprocedure()

}
void mvgstudy:: emcpscalarprocedure()
{
	real vector ems
	covcomps = asarray_create()   
	ems = luinv(P) * compsstacked
	for (i=1;i<=rows(ems);i++) {
		name = "emcp" + strofreal(i) 
		asarray(covcomps,name,ems[i])
	}
}
void mvgstudy:: emcpmatrixprocedure()
{
	real scalar 			skip 
	real matrix 			rangemat, starting, ending, inv_p,temp,res
	real scalar 			i, j
	
	covcomps = asarray_create()   
	skip = cols(compsstacked)
	rangemat = range(0,rows(compsstacked),skip)
	starting = rangemat[1..(rows(rangemat)-1)]:+1
	ending = rangemat[2..rows(rangemat)]
	inv_p = luinv(P)
	for(i=1;i<=cols(P);i++){
		res = J(cols(compsstacked),cols(compsstacked),0)
		for(j=1;j<=cols(P);j++) {
			temp = compsstacked[starting[j]..ending[j],1..cols(compsstacked)]
			res = res + (temp * inv_p[i,j])
		}
		name = "emcp" + strofreal(i) 
		asarray(covcomps,name,res)
	}
}


// (2) Reliability Projections 


// (3) Profile Variances 


// (4) Disattenuated Correlations 


// (5) Partial Correlations 


end 
	
	
	
	
	
asd 
/*
	
	
	
	
	
	

capture program drop emcpmatrixprocedure
program emcpmatrixprocedure, rclass 
	syntax,  DF(string) P(string) 
	
	capture matrix drop stacked
	forvalues i = 1/`:colsof(`p')' {
		scalar a = df[`i',1]
		mata df = st_numscalar("a")
		mata temp = st_matrix("sscp`i'")
		mata mscp`i' = temp :/df
		mata st_matrix("mscp`i'",mscp`i')
		matrix stacked = (nullmat(stacked) \ mscp`i')
		local emcpname emcp`i'
		local emcpnames : list emcpnames | emcpname
	}
	if `:colsof(stacked)' == 1 {
		mata ems = luinv(st_matrix("`p'")) * st_matrix("stacked") 
		mata st_matrix("ems",ems)
		forvalues i = 1/`:colsof(`p')'{
			matrix comp`i' = ems[`i',1]
			return matrix comp`i' = comp`i'
		}
	}
	else {
		mata A = emcpmatrixprocedure(st_matrix("stacked"),st_matrix("`p'"),tokens(st_local("emcpnames")))	
		mata for (loc=asarray_first(A); loc!=NULL; loc=asarray_next(A, loc)) st_matrix(asarray_key(A, loc),asarray_contents(A, loc))
		forvalues i = 1/`:colsof(`p')'{
			return matrix comp`i' = emcp`i'
		}
	}
	
	
	
	
	
	
	
	
	mata 
	
		for (loc=asarray_first(sscplist); loc!=NULL; loc=asarray_next(sscplist, loc)) {
		temp = asarray_contents(sscplist, loc)
		temp
	}
	
	
	end 
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	matrix df = r(df)
	matrix flproducts = r(flproducts)
	forvalues i = 1/`lengtheffectlist' {
		matrix sscp`i' = r(sscp`i')
		local mat sscp`i'
		local matrixlist : list matrixlist | mat 
	}	
	* (3) Create P matrix (EMCP equation matrix)
	createP, effects(`effects') flproducts(flproducts)
	matrix P = r(P)
	* (4) Expected Mean Cross Products (covariance components)
	emcpmatrixprocedure, df(df) p(P) 		   
	forvalues i = 1/`lengtheffectlist' { 
		matrix covcomps`i' = r(comp`i')
	}
	// Display Results ----------
	foreach x of local varlist {
		local name   "`x'"
		local names `" `names' "`name'" "'
	}
	forvalues i = 1/`lengtheffectlist' {
		matrix rownames covcomps`i' = `names'
		matrix colnames covcomps`i' = `names'
		matlist covcomps`i' , twidth(30) title("`:word `i' of `effects'' Component Covariance Matrix")
	}
	// Return results ----------
	forvalues i = 1/`lengtheffectlist' { 
		return matrix covcomps`i' = covcomps`i'
	}
	return matrix P = P 
	return matrix df = df 
	return matrix flproducts = flproducts 		
end 

cap prog drop parse_equation
program parse_equation, rclass
	syntax anything 
	
	tempname str pos i 
	mata `str' = st_local("anything")
	mata `str' = tokens(subinstr(subinstr(`str',"(",""),")",""))
	mata `pos' = strpos(`str',"=")
	mata `i' = 1
	mata while (`pos'[`i'] < 1) ++`i'
	mata st_local("vars", invtokens(`str'[1,1..`i' - 1]))
	mata st_local("effects",invtokens(`str'[1,`i' + 1..cols(`str')]))
	return local vars `vars' 
	return local effects `effects'
end

capture program drop recovermanova
program recovermanova, nclass 
	syntax,  FACETS(string) FACETLEVELS(string) EFFECTS(string) RESIDUAL(string)
	
	* Count the number of effects
	local a = e(cmdline) 
	local check = strrpos("`a'","=")
	local newlist = substr("`a'",`check'+1,.)
	local length : list sizeof local(newlist)
	local j = `length' + 1
	forvalues i = 1/`j' {
		local name sscp`i'
		local sscp_tempmats : list sscp_tempmats | name
	}
	
	tempname tempmat df flproducts `sscp_tempmats'
	
	* Degrees of freedom 
	forvalues i = 1/`length' {
		matrix `tempmat' = `e(df_`i')'
		matrix `df' = (nullmat(`df') \ `tempmat')
	}
	matrix `df' = `df' \ e(df_r)
	
	* Sums of squared cross product matrices 
	forvalues i = 1/`length' {
		matrix `sscp`i'' = e(H_`i')
		local++i
	}
	matrix `sscp`j'' = e(E)
	
	* Send to Mata 
	mata sscplist = asarray_create() 
	forvalues i = 1/`j' {
		mata asarray(sscplist ,"sscp`i'", st_matrix("`sscp`i''"))
	}
	mata effects = tokens(st_local("effects"))
	mata facets = tokens(st_local("facets"))
	mata facetlevels = st_matrix("`facetlevels'")
	mata flproducts = get_facetlevelproducts(effects,facets,facetlevels)
	mata df = st_matrix("`df'")

end 

capture program drop createP
program createP, rclass 
	syntax,   EFFECTS(string) FLPRODUCTS(string)
	
	mata st_matrix("P",get_p(tokens(st_local("effects")),st_matrix("`flproducts'")))
	foreach x of local effects {
		local name   "`x'"
		local names `" `names' "`name'" "'
	}
	matrix rownames P = `names'
	matrix colnames P = `names'
	return matrix P = P
end 

capture program drop emcpmatrixprocedure
program emcpmatrixprocedure, rclass 
	syntax,  DF(string) P(string) 
	
	capture matrix drop stacked
	forvalues i = 1/`:colsof(`p')' {
		scalar a = df[`i',1]
		mata df = st_numscalar("a")
		mata temp = st_matrix("sscp`i'")
		mata mscp`i' = temp :/df
		mata st_matrix("mscp`i'",mscp`i')
		matrix stacked = (nullmat(stacked) \ mscp`i')
		local emcpname emcp`i'
		local emcpnames : list emcpnames | emcpname
	}
	if `:colsof(stacked)' == 1 {
		mata ems = luinv(st_matrix("`p'")) * st_matrix("stacked") 
		mata st_matrix("ems",ems)
		forvalues i = 1/`:colsof(`p')'{
			matrix comp`i' = ems[`i',1]
			return matrix comp`i' = comp`i'
		}
	}
	else {
		mata A = emcpmatrixprocedure(st_matrix("stacked"),st_matrix("`p'"),tokens(st_local("emcpnames")))	
		mata for (loc=asarray_first(A); loc!=NULL; loc=asarray_next(A, loc)) st_matrix(asarray_key(A, loc),asarray_contents(A, loc))
		forvalues i = 1/`:colsof(`p')'{
			return matrix comp`i' = emcp`i'
		}
	}
end 


mata
string matrix sort_desc_length(string vector M)
{
	real   scalar 		i,j
	string scalar		temp 
	string matrix 		res
	
	res = M 
	for (i=2;i<=length(M);i++){
		temp = res[i]
		j = i - 1
		while (j >= 1 & strlen(temp) > strlen(res[j])) {
			res[j + 1] = res[j]
			j--
			if (j < 1) break 
		}
		res[j+1] = temp 
	}
	return(res)
}

string matrix sort_asc_length(string vector M)
{
	real   scalar 		i,j
	string scalar		temp 
	string matrix 		res
	
	res = M 
	for (i=2;i<=length(M);i++){
		temp = res[i]
		j = i - 1
		while (j >= 1 & strlen(temp) < strlen(res[j])) {
			res[j + 1] = res[j]
			j--
			if (j < 1) break 
		}
		res[j+1] = temp 
	}
	return(res)
}

// Rule π(α̇) = {1 if α = ω; and, otherwise, the product of the sample sizes for all indices in α̇}
real vector get_facetlevelproducts(string vector effects, string vector facets, real vector fls)
{
	real scalar 		flp,x,i,j
	real vector         flps
	
	flps = J(1,length(effects),.)
	for(i=1;i<=length(effects);i++) {
		flp = 1 
		for(j=1;j<=length(facets);j++) {
			if (strpos(effects[i],facets[j]) == 0) x = fls[j]
			else x = 1
			flp = flp * x 
		}
		flps[i] = flp 
	}
	flps = flps'
	return(flps)
	flps
}

// Private - covariance component calculation routine 
void mvgstudy::get_facetlevelproducts()
{
	real scalar 		flp,x,i,j

	flproducts = J(1,length(effects),.)
	for(i=1;i<=length(effects);i++) {
		flp = 1 
		for(j=1;j<=length(facets);j++) {
			if (strpos(effects[i],facets[j]) == 0) x = fls[j]
			else x = 1
			flp = flp * x 
		}
		flproducts[i] = flp 
	}
	flproducts = flproducts'
}










real matrix get_p(string vector effects, real vector flproducts)
{
	string rowvector 	effects_stripped
	real scalar 		row, col 
	real matrix 		P
	
	effects_stripped = tokens(subinstr(subinstr(invtokens(effects),"|",""), "#",""))
	P = J(length(effects_stripped),length(effects_stripped),0)
	for(col=1;col<=length(effects_stripped);col++){
		row = col
		while (row >= 1) {
			P[row, col] = flproducts[col]
			--row
		}
	}
	for(row=1;row<=length(effects_stripped);row++) {
		for(col=1;col<=length(effects_stripped)-1;col++) {
			if (strpos(effects_stripped[col],effects_stripped[row])== 0) P[row,col] = 0 
		}
	}
	return(P)
}

transmorphic matrix emcpmatrixprocedure(real matrix stacked, real matrix p, string matrix names)
{
	transmorphic matrix 	emcp
	real scalar 			skip 
	real matrix 			rangemat, starting, ending, inv_p,temp,res
	real scalar 			i, j
	
	emcp = asarray_create()   
	skip = cols(stacked)
	rangemat = range(0,rows(stacked),skip)
	starting = rangemat[1..(rows(rangemat)-1)]:+1
	ending = rangemat[2..rows(rangemat)]
	inv_p = luinv(p)
	for(i=1;i<=cols(p);i++){
		res = J(cols(stacked),cols(stacked),0)
		for(j=1;j<=cols(p);j++) {
			temp = stacked[starting[j]..ending[j],1..cols(stacked)]
			res = res + (temp * inv_p[i,j])
		}
		asarray(emcp,names[i],res)
	}
	return(emcp)
}

real scalar is_colvec(z) return(anyof(("colvector","scalar"),orgtype(z)))

void row_to_col(real vector v) 
{
    if (is_colvec(v) == 0) v = v'
	else v = v 
}



real matrix satterthwaite(real vector variances,
						  real vector df,
						  real scalar ci_level)
{
	real colvector r, lower, upper 
	real matrix res 

	row_to_col(variances)
	row_to_col(df)
	r = sqrt(df/2)
	upper = variances:* (2*r:^2:/invchi2(df,(1-ci_level)/2))
	lower = variances:* (2*r:^2:/invchi2(df,(1+ci_level)/2))
	res = lower , upper 
	return(res) 
}

real scalar relativevariance(real matrix true, real matrix obs, real colvector obsmeans, real scalar df) 
{
	return(profilevar(true,obsmeans,df)/profilevar(obs,obsmeans,df))
}

real scalar profilevar(real matrix covmat, real colvector obsmeans, real scalar df)
{
	real scalar avgdiag, avgelements, profilevar
	
	avgdiag = sum(diagonal(covmat)) / rows(covmat)
	avgelements = sum(covmat) / (rows(covmat) * rows(covmat))
	profilevar = ((df - 1) / df) * (avgdiag - avgelements) + variance(obsmeans) / rows(obsmeans)
	return(profilevar)
}

real scalar errorprofilevar(real matrix errorcovmat)
{
	real scalar avgdiag, avgelements, profilevar
	
	avgdiag = sum(diagonal(errorcovmat)) / rows(errorcovmat)
	avgelements = sum(errorcovmat) / (rows(errorcovmat) * rows(errorcovmat))
	profilevar = avgdiag - avgelements
	return(profilevar)

}

real matrix off_diag(real matrix cov_mat)
{
	real scalar v, a, b, c, r, i, x
	real matrix res
	
	v = length(cov_mat[1...,1])
	a = v - 1
	b = 1
	for (i=a;i>=1;i--) {
		for (x=1;x<=i;x++) {
			c = v - i
			r = x + b
			if (b == 1 & x == 1) res = cov_mat[r,c]
			else res = (res \ cov_mat[r,c])
		}
		b++
	}
	return(res)	
}

real matrix rebuild_matrix(real matrix cov_mat, 
						    real matrix edited_covs,
					        real matrix variances)
{
	real matrix res 
	real scalar k, g, b, c, r1, r2 , i , x 
	
	res = J(rows(cov_mat),rows(cov_mat), . )
	k = 0
	g = 1
	for (i=(rows(cov_mat));i>=1;i--) {
		b = i - 1
		c = rows(cov_mat) - b 
		
		for (x=1; x <= b; x++) {
			r1 = x + k 
			r2 = x + g
			res[r2,c] = edited_covs[r1,1]
		}
		k = k + b
		g++ 
	}
	_diag(res, variances)
	_makesymmetric(res)
	return(res)
} 

string matrix off_diag_key(real matrix cov_mat)
{
	real   scalar v, a, b, c, r, i, x
	string matrix res
	string scalar str1, str2, str 
	
	v = length(cov_mat[1...,1])
	a = v - 1
	b = 1
	for (i=a;i>=1;i--) {
		for (x=1;x<=i;x++) {
			c = v - i
			r = x + b
			str1 = strofreal(r)
			str2 = strofreal(c)
			str = str2 + " " + str1
			if (b == 1 & x == 1) res = str
			else res = (res \ str)
		}
		b++
	}
	return(res)	
}
end 


