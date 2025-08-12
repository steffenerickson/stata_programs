
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
*! Multivariate Generalizability and Decision Study Decompostion Commands
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
// Main Stata Routine
//----------------------------------------------------------------------------//
program mvgstudy, rclass
	version 19
	syntax anything [if] [in] 
	
	marksample touse 
	tempname tempframe 
	
	* Check if `anything` matches the pattern "( * = * )"
	//if !regexm("`anything'", "^\(\s*[^=]+?\s*=\s*.+?\s*\)$") {
	//    di "Invalid Syntax"
	//	exit
	//} 
	
	// 1: Set up 
	mata c = mvgstudy()
	parse_equation `anything'
	local tempvarlist `r(vars)'
	local effects `r(effects)'
	mata st_local("effects",invtokens(c.sortbylength(tokens(st_local("effects")),0))) // making sure effects are sorted so that the resdiual is on the end of the string list
	foreach var of varlist `tempvarlist' {
		local varlist : list varlist | var
	}
	
	mata facets = uniqrows(tokens(subinstr(subinstr(invtokens(st_local("effects")),"|"," "), "#"," "))')
	mata st_local("facets",invtokens(facets'))
	foreach facet of local facets {
		qui levelsof `facet' if `touse' == 1
		matrix facetlevels = (nullmat(facetlevels) \ `r(r)')
	}
	mata residual = c.sortbylength(tokens(st_local("effects")),1)[1]
	mata st_local("residual",residual)
	local effects2 : list effects - residual
	
	// 2: Manova and collect results
	frame put `varlist' `facets' if `touse' == 1, into(`tempframe')
	frame `tempframe' : qui manova `varlist' = `effects2'
	recovermanova, residual(`residual')
	mata sscplist = asarray_create() 
	forvalues i = 1/`r(j)' {
		mata asarray(sscplist ,"`i'", st_matrix("r(sscp`i')"))
	}
	mata effects = tokens(st_local("effects"))
	mata facets = tokens(st_local("facets"))
	mata varlist = tokens(st_local("varlist"))
	mata facetlevels = st_matrix("facetlevels")
	mata df = st_matrix("r(df)")
	mata c.init_inputs(sscplist,effects,facets,varlist,facetlevels,df)
	mata mata drop sscplist effects facets facetlevels df
	
	// 3: Mata mvgstudy main routine 
	mata c.mvgstudy_main_routine()
	mata for (loc=asarray_first(c.covcomps); loc!=NULL; loc=asarray_next(c.covcomps, loc)) st_matrix(asarray_key(c.covcomps, loc),asarray_contents(c.covcomps, loc))
	mata st_matrix("df",c.df)
	mata st_matrix("P",c.P)
	
	// 4: Results
	local lengtheffectlist : list sizeof local(effects)
	foreach x of local varlist {
		local name   "`x'"
		local names `" `names' "`name'" "'
	}
	forvalues i = 1/`lengtheffectlist' {
		matrix rownames emcp`i' = `names'
		matrix colnames emcp`i' = `names'
		matlist emcp`i' , twidth(20) title("`:word `i' of `effects'' Component")
	}
	forvalues i = 1/`lengtheffectlist' { 
		return matrix emcp`i' = emcp`i'
	}
	return matrix P = P 
	return matrix df = df 
	return local varlist `varlist'
	return local effects `effects'
	
end 

program mvdstudy, rclass 
	version 19
	syntax, Object(string)  Errortype(string) FACETnum(integer) [COMPositeweights(string)]
	
	if ("`errortype'" == "relative") local etype = 0 
	else if ("`errortype'" == "absolute") local etype = 1
	else {
		di "invalid error type"
		exit
	}
	
	if ("`compositeweights'" == "") mata c.init_dstudyinputs("`object'",`etype',`facetnum')
	else mata c.init_dstudyinputs("`object'",`etype',`facetnum', st_matrix("`compositeweights'"))
	
	mata c.mvdstudy_main_routine()
	
	mata {
		stringlist = ""
		for (loc=asarray_first(c.projections); loc!=NULL; loc=asarray_next(c.projections, loc)) {
			stringlist = stringlist + " " + asarray_key(c.projections, loc)
			st_matrix(asarray_key(c.projections, loc),asarray_contents(c.projections, loc))
		}
		colnames = invtokens((tokens(c.projectionmat[rows(c.projectionmat)]) , "Error" , "True" , "Rel."))
		st_local("vars",stringlist)
		st_local("colnames",colnames)
	}
	
	foreach v of local vars {
		local names 
		foreach x of local colnames {
			local name   "`x'"
			local names `" `names' "`name'" "'
		}
		mat colnames `v' = `names'
		matlist `v' , names(c) title("`v'")
		return matrix `v' = `v'
	}
	
end 

//----------------------------------------------------------------------------//
// Stata Sub Routines
//----------------------------------------------------------------------------//
program parse_equation, rclass
	version 19
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
	version 19
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

	//----------//
	public:
	//---------//
	
	// Gstudy routines 
	void init_inputs()
	void mvgstudy_main_routine()
	
	// D Study routines
	void init_dstudyinputs()
	void mvdstudy_main_routine()
	
	// Mata functions used in stata wrapper program 
	string matrix sortbylength()
	real vector flip_select()
	
	// Gstudy objects 
	string vector   effects
	string vector 	facets
	string vector 	varlist
	real vector		facetlevels
	real vector 	flproducts
	real vector 	df
	real matrix     P
	real matrix     compsstacked
	transmorphic  	sscplist
	transmorphic  	mscplist
	transmorphic  	covcomps
	
	// Dstudy objects 
	string scalar objmeasurement
	real scalar errortype
	real scalar projectionnum 
	real colvector weights
	string matrix projectionmat
	transmorphic  	projections

	//----------//
	private:
	//----------//
	
	// Gstudy internal routines
	void get_facetlevelproducts()
	void get_p()
	void get_mscplist()
	void get_compsstacked() 
	void emcp()
	void emcpmatrixprocedure()
	void emcpscalarprocedure()
	
	real scalar relativevariance()
	real scalar profilevar()
	real scalar errorprofilevar()

	// Dstudy internal routines
	void absolute_errors_effects()
	void relative_errors_effects()
	real matrix errorprojections()
	void combine_matrices()
	real matrix rowwise_nonzero_product()
	real matrix compositevariance()
	real vector extract_bracket_numbers()
	
	real vector errorvariances
	real scalar  objmeasurevariance
	
}

//Public Routines
void mvgstudy::init_inputs(transmorphic sscplist,
                           string vector effects,
						   string vector facets,
						   string vector varlist,
						   real vector facetlevels,
						   real df)
{
	this.sscplist = sscplist
	this.effects = effects
	this.facets = facets
	this.varlist = varlist
	this.facetlevels = facetlevels
	this.df = df
}

void mvgstudy::mvgstudy_main_routine()
{
	get_facetlevelproducts()
	get_p()
	get_mscplist()
	get_compsstacked() 
	emcp()
}

void mvgstudy::init_dstudyinputs(string scalar objmeasurement, 
					  real scalar errortype, 
					  real scalar projectionnum, 
					| real vector weights)
{
	this.objmeasurement = objmeasurement
	this.errortype = errortype
	this.projectionnum = projectionnum
	if (args()>3) this.weights = weights 
	else this.weights = J(0,0,.)
}


void mvgstudy::mvdstudy_main_routine()
{
	real vector variances 
	real matrix result

	projections = asarray_create()  
	if (weights == J(0,0,.)) {
		for (i=1;i<=length(varlist);i++) {
			variances = J(length(effects),1,.)
			for (j=1;j<=length(effects);j++) {
				name = "emcp" + strofreal(j) 
				variances[j] = asarray(covcomps, name)[i,i]
			}
			if (errortype == 0) relative_errors_effects(variances)	
			else if (errortype == 1) absolute_errors_effects(variances)
			result = errorprojections()
			asarray(projections,varlist[i],result)
		}
	}
	else {
		variances =  compositevariance()
		if (errortype == 0) relative_errors_effects(variances)	
		else if (errortype == 1) absolute_errors_effects(variances)
		result = errorprojections()
		asarray(projections,"composite",result)
	}
}

string matrix mvgstudy:: sortbylength(string vector M, real scalar descending)
{
	real   scalar 		i,j
	string scalar		temp 
	string matrix 		res
	
	res = M 
	for (i=2;i<=length(M);i++){
		temp = res[i]
		j = i - 1
		if (descending == 1) {
			while (j >= 1 & strlen(temp) > strlen(res[j])) {
				res[j + 1] = res[j]
				j--
				if (j < 1) break 
			}
		}
		else {
			while (j >= 1 & strlen(temp) < strlen(res[j])) {
				res[j + 1] = res[j]
				j--
				if (j < 1) break 
			}
		}
		res[j+1] = temp 
	}
	return(res)
}

real vector mvgstudy::flip_select(real vector mask) return(mm_cond(mask:>0,mask:-mask,mask:+1))


//Private Routines
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
real matrix mvgstudy:: compositevariance()
{
	real matrix M 
	real scalar i 
	string scalar name
	real vector compositevariances
	
	compositevariances = J(length(effects),1,.)
	
	for (i=1;i<=length(effects);i++) {
		name = "emcp" + strofreal(i) 
		M = asarray(covcomps, name)
		compositevariances[i] = weights'*M*weights 
	}
	return(compositevariances)
}

void mvgstudy:: absolute_errors_effects(variances)
{
	real matrix s1,s2
	string matrix erroreffects
	
	s1 = effects':== objmeasurement
	s2 = flip_select(s1)
	objmeasurevariance = select(variances,s1)
	erroreffects   = select(effects',s2)
	errorvariances = select(variances,s2)
	projectionmat = strtrim(subinstr(subinstr(subinstr(erroreffects,"|"," "), "#"," "),objmeasurement," "))
}
void mvgstudy:: relative_errors_effects(variances)
{
	
	real matrix s1,s2,s3
	string matrix erroreffects
	
	s1 = effects':==objmeasurement
	s2 = flip_select(s1)
	objmeasurevariance = select(variances,s1)
	erroreffects   = select(effects',s2)
	errorvariances = select(variances,s2)
	s3 = strpos(erroreffects,objmeasurement)
	erroreffects   = select(erroreffects,s3)
	errorvariances = select(errorvariances,s3)
	projectionmat = strtrim(subinstr(subinstr(subinstr(erroreffects,"|"," "), "#"," "),objmeasurement," "))
}
real matrix mvgstudy::errorprojections()
{
	real vector stringlengths,selectvec,nums,div,temp2
	string vector uniqueffects,group_keys
	string matrix res,temp 
	real matrix first,zero,result,M
	real scalar colspan, difference, i
	transmorphic inner,mats, outerarray
	
	//stringlengths = strlen(projectionmat)
	//selectvec = (stringlengths :== max(stringlengths)[1])
	//uniqueffects = tokens(select(projectionmat,selectvec))
	
	uniqueffects = tokens(projectionmat[rows(projectionmat)])
	colspan = cols(uniqueffects)
	res = J(0,colspan,"")
	for (i=1;i<=rows(projectionmat);i++) {
		difference = colspan - cols(tokens(projectionmat[i]))
		if (difference !=0) temp = tokens(projectionmat[i]), J(1,difference,"")
		else temp = tokens(projectionmat[i])
		res = res \ temp 
	}
	outerarray = asarray_create() 
	for (i=1;i<=cols(uniqueffects);i++) {
		inner = asarray_create() 
		asarray(outerarray,uniqueffects[i],inner)
	}
	for (i=1;i<=cols(uniqueffects);i++) {
		numres = (res:==uniqueffects[i])
		//printf(strofreal(i) + " Matrices")
		for (j=1;j<=projectionnum;j++) {
			M = numres * j 
			asarray(asarray(outerarray,uniqueffects[i]),strofreal(j),M)
		}
	}
	group_keys = asarray_keys(outerarray)
	mats = asarray_create() 
	first = asarray(asarray(outerarray, group_keys[1]), "1")
	zero = first * 0
	combine_matrices(outerarray,"", zero, 1,group_keys,projectionnum,mats)
	result = J(0,cols(first)+3,.)
	
	for (loc=asarray_first(mats); loc!=NULL; loc=asarray_next(mats, loc)) {
		nums = extract_bracket_numbers(asarray_key(mats, loc))
		div = rowwise_nonzero_product(asarray_contents(mats, loc))
		temp2 = nums, colsum(errorvariances:/div),objmeasurevariance, (objmeasurevariance / (objmeasurevariance + colsum(errorvariances:/div)))
		result = result \ temp2 
	}
	
	return(sort(result,(1..(cols(result)-3))))
}
void mvgstudy::combine_matrices(transmorphic outerarray, 
                      string scalar path, 
					  real matrix acc, 
					  real scalar depth,
					  string vector group_keys,
					  real scalar projectionnum,
					  transmorphic mats) 
{
	
	num_groups = length(group_keys)
    if (depth > num_groups) {
        // Base case: all groups processed
        //printf("Combination: %s\n", path)
        //acc
		asarray(mats,path,acc)
		return
    }
    group = group_keys[depth]
    for (i = 1; i <= projectionnum; i++) {
        M = asarray(asarray(outerarray, group), strofreal(i))
        new_path = path + group + "[" + strofreal(i) + "] "
        combine_matrices(outerarray,new_path, acc + M, depth + 1,group_keys, projectionnum, mats)
    }

}
real matrix mvgstudy:: rowwise_nonzero_product(real matrix X) {
    real matrix result 
	real scalar i 
	
	result = J(rows(X), 1, 1) 
    for (i = 1; i <= rows(X); i++) {
        for (j = 1; j <= cols(X); j++) {
            if (X[i, j] != 0) {
                result[i] = result[i] * X[i, j]
            }
        }
    }

    return(result)
}

real vector mvgstudy:: extract_bracket_numbers(string scalar s) 
{
    string scalar pattern 
    real vector results
    real scalar pos 
    real scalar match
	
	pattern = "\[([0-9]+)\]"
	results = J(1, 0, .) 
	pos = 1

    while (regexm(substr(s, pos, .), pattern)) {
        match = strtoreal(regexs(1))
        results = results , match
        // Move position forward to search for next match
        pos = pos + strpos(substr(s, pos, .), "]")
    }
    
    return(results)
}

// (3) Profile Variances 

end 

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Conditional Standard Errors of measurement program 
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
	
program conditionalsem , rclass

	syntax, Variables(varlist) Facets(varlist) OBJ(varlist max=1 ) [COMPositeweights(string)]

	tempname tempframe 
	frame put `facets' `variables' , into(`tempframe')
	quietly {
		frame `tempframe' {
			local pos = strpos("`facets'", "`obj'")
			tuples `facets' , nopython conditionals(!`pos')
			forvalues i = 1/`ntuples' {
				local current `tuple`i''
				local name : subinstr local current " " "", all
				egen group_`name' = group(`current')
				
				preserve
					collapse `variables' , by(`obj' group_`name') 
					mata X = st_data(.,("`obj'","`variables'"))
				restore 
				
				if "`compositeweights'" != "" {
					mata st_matrix("`name'",composite_sem(X,st_matrix("`compositeweights'")))
					preserve
						clear 
						svmat `name'
						local variances
						rename (*) (`obj' comp_`name'_variance)
						tempfile data 
						save `data'
					restore 
				}
				else {
					if (`:word count `variables'' == 1) mata st_matrix("`name'",univariate_sem(X))
					else mata st_matrix("`name'",multivariate_sem(X))
					preserve
						clear 
						svmat `name'
						local variances
						foreach var of local variables {
							local temp `var'_`name'_variance
							local variances : list variances | temp
						}
						rename (*) (`obj' `variances')
						tempfile data 
						save `data'
					restore 
				}
				
				merge m:1 `obj' using `data' 
				drop _merge
				
			}
			if "`compositeweights'" != "" {
				foreach var of local variables  {
					egen `var'_mean = mean(`var') , by(`obj')
				}
				
				local rhs 
				forvalues i = 1/`:rowsof(`compositeweights')/' {
					local var `:word `i' of `variables''
					local w = `compositeweights'[`i',1]
					local term "+`w'*`var'_mean" 
					local rhs : list rhs |term 
				}
				local rhs : subinstr local  rhs "+" ""		
				gen comp_mean = `rhs'
		
				tempvar temp
				egen `temp' = rowmean(comp_*_variance) 
				gen comp_abs = sqrt(`temp')
				collapse comp_mean comp_abs , by (`obj')
			}
			else {
				foreach var of local variables  {
					tempvar temp
					egen `var'_mean = mean(`var') , by(`obj')
					egen `temp' = rowmean(`var'_*_variance) 
					gen `var'_abs = sqrt(`temp')
				}
				collapse *mean *abs , by(`obj')
			}
			mkmat *mean *abs , matrix(result) rownames(`obj')
		}
	}
	return matrix result = result
	
end 

mata 

real matrix composite_sem(real matrix X, real colvector w)
{
	real matrix info,res,subX
	real scalar i, compvariance

	info = panelsetup(X, 1)
	res = J(0,2,.)
	for (i=1; i<=rows(info); i++) {
		subX = panelsubmatrix(X, i, info)
		compvariance = (w' * variance(subX[.,2..cols(subX)]) * w)/ rows(subX)
		res = res \ (subX[1,1], compvariance)
	}
	return(res)
}
real matrix multivariate_sem(real matrix X)
{
	real matrix info,res,subX,variances
	real scalar i 
	
	info = panelsetup(X, 1)
	res = J(0,cols(X),.)
	for (i=1; i<=rows(info); i++) {
		subX = panelsubmatrix(X, i, info)
		variances = (diagonal(variance(subX[.,2..cols(subX)])):/rows(subX))'
		res = res \ (subX[1,1], variances)
	}
	return(res)
}
real matrix univariate_sem(real matrix X)
{
	real matrix info,res,subX
	real scalar i, variance
	
	info = panelsetup(X, 1)
	res = J(0,2,.)
	for (i=1; i<=rows(info); i++) {
		subX = panelsubmatrix(X, i, info)
		variance = variance(subX[.,2]) / rows(subX)
		res = res \ (subX[1,1], variance)
	}
	return(res)
}

end 
	
	

