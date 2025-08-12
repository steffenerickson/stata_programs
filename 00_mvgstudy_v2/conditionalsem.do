
clear all 

cd /Users/steffenerickson/Documents/GitHub/stata_programs/00_mvgstudy_v2
include mvgstudy.ado 





use "/Users/steffenerickson/Box Sync/NSF_DR_K12/measurement/${data}/manova_data.dta", clear 
rename (task rater coaching) (t r treat)
recode time (2 =1)
drop if (r == 3) | (treat == 2) | (x1 == . & x2 == . & x3 == . & x4 == . & x5 == .)
encode participantid , gen(id)
egen p = group (id site semester treat)
drop x6
keep if time == 1 
collapse x1 x2 x3 x4 x5, by(p t r)
keep p t  r x*
order p t  r x*
sort p t r 
//egen counts = count(x1) , by(p)
//keep if counts == 4
//drop counts 
egen p2 = group(p)
replace p = p2 
drop p2

	
	
conditionalsem, variables(x1 x2 x3 x4 x5) facets(p t r) obj(p)
matrix semmat = r(result) 
local colnames : colfullnames semmat 
preserve
clear 
svmat semmat 
rename (*) (`colnames')
forvalues i = 1/5{
		twoway  (scatter x`i'_abs x`i'_mean) (lowess x`i'_abs x`i'_mean) , name(g`i', replace) legend(off)
}
graph combine g1 g2 g3 g4 g5 , ycommon xcommon 
restore
	

	
matrix weights = J(5,1,1/5)
conditionalsem, variables(x1 x2 x3 x4 x5) facets(p t r) obj(p) comp(weights)
matrix semmat = r(result) 
local colnames : colfullnames semmat 
preserve
clear 
svmat semmat 
rename (*) (`colnames')
twoway (scatter comp_abs comp_mean) (lowess comp_abs comp_mean) , legend(off) name(g1, replace)
restore


matrix weights = J(5,1,1/5)

matrix weights = (.05,.35,.35,.2,.05)'
conditionalsem, variables(x1 x2 x3 x4 x5) facets(p t r) obj(p) comp(weights)
matrix semmat = r(result) 
local colnames : colfullnames semmat 
preserve
clear 
svmat semmat 
rename (*) (`colnames')
twoway (scatter comp_abs comp_mean) (lowess comp_abs comp_mean) , legend(off) name(g2, replace)
restore
	
graph combine g1 g2  , ycommon xcommon 

	
	
	
	

	twoway (scatter x1_abs x1_mean) (qfit x1_abs x1_mean) ///
	       (scatter x2_abs x2_mean) (qfit x2_abs x2_mean) ///
		   (scatter x3_abs x3_mean) (qfit x3_abs x3_mean) ///
		   (scatter x4_abs x4_mean) (qfit x4_abs x4_mean) ///
		   (scatter x5_abs x5_mean) (qfit x5_abs x5_mean) 


		   
	



forvalues i = 1/5{
		twoway  (scatter x`i'_abs x`i'_mean) (qfit x`i'_abs x`i'_mean) , name(g`i', replace) legend(off)
}
graph combine g1 g2 g3 g4 g5 , ycommon xcommon 


/*

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


conditionalsem, variables(x1 x2 x3 x4 x5) facets(p t r) obj(p)

	
matrix weights = J(5,1,1/5)
conditionalsem, variables(x1 x2 x3 x4 x5) facets(p t r) obj(p) comp(weights)

set trace off
conditionalsem, variables(x1 ) facets(p t r) obj(p) 

	



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



	local facets p t r
	local variables x1 x2 x3 x4 x5
	local obj p




















use mvgstudyexampledata.dta, clear
mata data = st_data(.,"*")
mata info = panelsetup(data, 1)
mata 
res = J(0,11,.)
for (i=1; i<=rows(info); i++) {
    X = panelsubmatrix(data, i, info)
	vars = (diagonal(variance(X[.,3..7])):/rows(X))'
	means = mean(X[.,3..7])
	id = X[1,1]
	rows = id, means, vars
	res = res \ rows
}

st_matrix("res",res)
end 
clear 
svmat res
sort res3
twoway (scatter res8 res3) (qfit res8 res3)



use mvgstudyexampledata.dta, clear
mata data = st_data(.,"*")
mata info = panelsetup(data, 1)
mata 
res = J(0,3,.)
weights = J(5,1,1/5)
for (i=1; i<=rows(info); i++) {
    X = panelsubmatrix(data, i, info)
	compvar = (weights' * variance(X[.,3..7]) * weights)/ rows(X)
	compvar = sqrt(compvar)
	means = mean(mean(X[.,3..7])')
	id = X[1,1]
	rows = id, means, compvar
	res = res \ rows
}

st_matrix("res",res)
end 
clear 
svmat res
sort res2
twoway (scatter res3 res2) (qfit res3 res2)



use "/Users/steffenerickson/Box Sync/NSF_DR_K12/measurement/${data}/manova_data.dta", clear 
rename (task rater coaching) (t r treat)
recode time (2 =1)
drop if (r == 3) | (treat == 2) | (x1 == . & x2 == . & x3 == . & x4 == . & x5 == .)
encode participantid , gen(id)
egen p = group (id site semester treat)
drop x6
keep if time == 1 
collapse x1 x2 x3 x4 x5, by(p t r)
keep p t  r x*
order p t  r x*
sort p t r 
//egen counts = count(x1) , by(p)
//keep if counts == 4
//drop counts 
egen p2 = group(p)
replace p = p2 
drop p2


mata data = st_data(.,"*")
mata info = panelsetup(data, 1)
mata 
res = J(0,3,.)
weights = J(5,1,1/5)
for (i=1; i<=rows(info); i++) {
    X = panelsubmatrix(data, i, info)
	compvar = (weights' * variance(X[.,4..8]) * weights)/ rows(X)
	compvar = sqrt(compvar)
	means = mean(mean(X[.,4..8])')
	id = X[1,1]
	rows = id, means, compvar
	res = res \ rows
}

st_matrix("res",res)
end 
clear 
svmat res
sort res2
twoway (scatter res3 res2) (qfit res3 res2)


use "/Users/steffenerickson/Box Sync/NSF_DR_K12/measurement/${data}/manova_data.dta", clear 
rename (task rater coaching) (t r treat)
recode time (2 =1)
drop if (r == 3) | (treat == 2) | (x1 == . & x2 == . & x3 == . & x4 == . & x5 == .)
encode participantid , gen(id)
egen p = group (id site semester treat)
drop x6
keep if time == 1 
collapse x1 x2 x3 x4 x5, by(p t)
keep p t  x*
order p t  x*
sort p t 
//egen counts = count(x1) , by(p)
//keep if counts == 4
//drop counts 
egen p2 = group(p)
replace p = p2 
drop p2


mata data = st_data(.,"*")
mata info = panelsetup(data, 1)
mata 
res = J(0,3,.)
weights = J(5,1,1/5)
for (i=1; i<=rows(info); i++) {
    X = panelsubmatrix(data, i, info)
	compvar = (weights' * variance(X[.,3..7]) * weights)/ rows(X)
	compvar = sqrt(compvar)
	means = mean(mean(X[.,3..7])')
	id = X[1,1]
	rows = id, means, compvar
	res = res \ rows
}

st_matrix("res",res)
end 
clear 
svmat res
sort res2
twoway (scatter res3 res2) (qfit res3 res2)



use "/Users/steffenerickson/Box Sync/NSF_DR_K12/measurement/${data}/manova_data.dta", clear 
rename (task rater coaching) (t r treat)
recode time (2 =1)
drop if (r == 3) | (treat == 2) | (x1 == . & x2 == . & x3 == . & x4 == . & x5 == .)
encode participantid , gen(id)
egen p = group (id site semester treat)
drop x6
keep if time == 1 
collapse x1 x2 x3 x4 x5, by(p t r)
keep p t  r x*
order p t  r x*
sort p t r 
//egen counts = count(x1) , by(p)
//keep if counts == 4
//drop counts 
egen p2 = group(p)
replace p = p2 
drop p2



	local facets p t r
	local variables x1 x2 x3 x4 x5
	local obj p

	tempname tempframe 
	frame put `facets' `variables' , into(`tempframe')
	frame `tempframe' {
		local pos = strpos("`facets'", "`obs'")
		tuples `facets' , nopython conditionals(!`pos')
		forvalues i = 1/`ntuples' {
			local current `tuple`i''
			local name : subinstr local current " " "", all
			egen group_`name' = group(`current')
			
			preserve
				collapse `variables' , by(`obj' `groupvar') 
				mata X = st_data(.,("`obj'","`variables'"))
			restore 
			
			if "`compositeweights'" != "" {
				mata st_matrix("`name'",composite_sem(X,st_matrix("`compositeweights'")))
			}
			else if `:word count `variables'' > 1 {
				mata st_matrix("`name'",multivariate_sem(X))
			}
			else {
				mata st_matrix("`name'",univariate_sem(X))
			}
			
			preserve
				clear 
				svmat `name'
				tempfile data 
				save `data'
			restore 
			
			merge m:1 `obj' using `data' 
			drop _merge
			

		}
		
			
	

	}
	
	mata 
	
	real matrix composite_sem(real matrix X, real colvector w)
	{
		real matrix info,res,subX
		real scalar i, compvariance, mu 
		
		info = panelsetup(X, 1)
		res = J(0,3,.)
		for (i=1; i<=rows(info); i++) {
			subX = panelsubmatrix(X, i, info)
			compvariance = (w' * variance(X[.,2...]) * w)/ rows(subX)
			mu = mean(mean(subX[.,2...])')
			res = res \ subX[1,1], means, compvariance
		}
	}
	real matrix multivariate_sem(real matrix X)
	{
		real matrix info,res,subX, mu, variances
		real scalar i 
		
		info = panelsetup(X, 1)
		colnum = (cols(X)-1) * 2 + 1
		res = J(0,colnum,.)
		for (i=1; i<=rows(info); i++) {
			subX = panelsubmatrix(X, i, info)
			variances = (diagonal(variance(X[.,2...])):/rows(subX))'
			mu = mean(subX[.,2...])
			res = res \ subX[1,1], mu, variances
		}
	}
	real matrix univariate_sem(real matrix X)
	{
		real matrix info,res,subX
		real scalar i, variance, mu
		
		info = panelsetup(X, 1)
		res = J(0,3,.)
		for (i=1; i<=rows(info); i++) {
			subX = panelsubmatrix(X, i, info)
			variance = variance(X[.,2]) / rows(subX)
			mu = mean(subX[.,2])
			res = res \ subX[1,1], mu, variance
		}
	}
	
	end 






