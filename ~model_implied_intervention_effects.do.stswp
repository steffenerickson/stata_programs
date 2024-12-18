*! Version 2 : Added Higher order factor model routine

mata mata clear


cap prog drop levelsforsem 
prog levelsforsem, rclass  
	syntax varlist(max=1), Local(string)
	qui levelsof `varlist', local(temp)
	forvalues i = 1/`:word count `temp'' {
		local lvl = `:word `i' of `temp''
		if (`i' == 1) local new `lvl'bn.`varlist'
		else          local new `lvl'.`varlist'
		local `local' : list `local' | new 
	}
	return local `local'  ``local''
end

/*
cap prog drop structual_analysis
program structual_analysis, nclass
	syntax, Model(string) Var(varlist(max=1)) Levels(string) ///
	[Higherorder(real 0) HIgherorderitem(string)]
	
	levelsforsem t, local(glvls)
	local glvls = r(glvls)
	
	estimates restore `model'
	local cmdline = e(cmdline)
	local colnames : colfullnames e(b)
	local lxvars = e(lxvars)
	local lyvars = e(lyvars)
	local oxvars = e(oxvars)
	local oyvars = e(oyvars)
	
	mata  glvls = tokens(st_local("glvls")) 
	mata  cmdline = st_local("cmdline")
	mata  lxvars = tokens(st_local("lxvars")) 
	mata  lyvars = tokens(st_local("lyvars"))
	mata  oxvars = tokens(st_local("oxvars")) 
	mata  oyvars = tokens(st_local("oyvars")) 
	mata  params = tokens(st_local("colnames"))'
	mata  levels = st_matrix("`levels'")
	
	if (`higherorder' != 0) mata  results = main(lxvars,lyvars,oyvars,params,glvls,cmdline,levels,`higherorder',`higherorderitem')
	else mata  results = main(lxvars,lyvars,oyvars,params,glvls,cmdline,levels)

end 
*/

mata
// Define structures for the SEM problem and derived results
	struct myproblem  {
	string vector   lxvars           	// Latent exogenous variables
	string vector   lyvars           	// Latent endogenous variables
	string vector   oxvars           	// Observed exogenous variables
	string vector   oyvars           	// Observed endogenous variables
	string vector   params           	// Model parameters (strings)
	real   vector   params_real         // Model parameters (strings)
	string vector   glvls            	// Group levels
	string scalar   cmdline          	// Command line (model specification)
	real   vector   comparegroups    	// Indices of groups to compare
	real   scalar   higherorder   	 	// Higher order factor model method option 
	string scalar   higherorderitem  	// Higher order factor model achor item 
	struct derived  scalar d         	// Nested struct for derived results
}

struct derived { 
	string 			matrix  pathcoefs        		// Selected path coefficients
	string 			matrix  intercepts       		// Selected intercept parameters
	transmorphic    		anchorkey        		// Anchors for latent variables
	string 			matrix  beta_bag         		// Endogenous-to-endogenous path coefs
	string 			matrix  gamma_bag        		// Exogenous-to-endogenous path coefs
	string 			matrix  lambdax_bag		
	string 			matrix  lambdax_a_bag    		// Anchored exogenous loadings
	string 			matrix  lambdax_na_bag   		// Non-anchored exogenous loadings
	string 			matrix  lambday_bag      		// Endogenous loadings
	string 			matrix  upsilon_y_bag    		// Intercepts for endogenous indicators
	string 			matrix  upsilon_xa_bag   		// Intercepts for anchored exogenous indicators
	string 			matrix  upsilon_xna_bag  		// Intercepts for non-anchored exogenous indicators
	transmorphic    		parammatrices_string 	// stores the following 
	string 			matrix  Beta1
	string 			matrix  Beta2
	string 			matrix  Gamma1
	string 			matrix  Gamma2
	string 			matrix  Lambday1
	string 			matrix  Lambday2
	string 			matrix  Lambdaxna1
	string 			matrix  Lambdaxna2
	transmorphic    		parammatrices_real 		// stores the following 
	real   			matrix  Beta1_real
	real   			matrix  Beta2_real
	real   			matrix  Gamma1_real
	real   			matrix  Gamma2_real
	real   			matrix  Lambday1_real
	real   			matrix  Lambday2_real
	real   			matrix  Lambdaxna1_real
	real   			matrix  Lambdaxna2_real	
	
	// Maybe should write an additonal routine to collect variances and covariances
}

// Initialize the derived structure with empty matrices and a new anchorkey asarray
void initialize_objects(struct derived scalar d)
{
	d.pathcoefs = J(0,0,"") 
	d.intercepts = J(0,0,"") 
	d.anchorkey  = asarray_create() 
	d.beta_bag  = J(0,0,"") 
	d.gamma_bag = J(0,0,"") 
	d.lambdax_bag = J(0,0,"") 
	d.lambdax_a_bag = J(0,0,"") 
	d.lambdax_na_bag = J(0,0,"") 
	d.lambday_bag  = J(0,0,"") 
	d.upsilon_y_bag = J(0,0,"") 
	d.upsilon_xa_bag = J(0,0,"") 
	d.upsilon_xna_bag = J(0,0,"") 
	d.parammatrices_string = asarray_create()
	d.Beta1 = J(0,0,"") 
	d.Beta2 = J(0,0,"") 
	d.Gamma1 = J(0,0,"") 
	d.Gamma2 = J(0,0,"") 
	d.Lambday1 = J(0,0,"") 
	d.Lambday2 = J(0,0,"") 
	d.Lambdaxna1 = J(0,0,"") 
	d.Lambdaxna2 = J(0,0,"") 
	d.parammatrices_real = asarray_create()	
	d.Beta1_real = J(0,0,"")
	d.Beta2_real = J(0,0,"")
	d.Gamma1_real = J(0,0,"")
	d.Gamma2_real = J(0,0,"")
	d.Lambday1_real = J(0,0,"")
	d.Lambday2_real = J(0,0,"")
	d.Lambdaxna1_real = J(0,0,"")
	d.Lambdaxna2_real = J(0,0,"")	 
}

// Main entry point: initialize and run the routines to process SEM results
struct derived main(lxvars,lyvars,oyvars,params,params_real,glvls,cmdline,comparegroups, | higherorder, higherorderitem)
{
    struct myproblem scalar pr
    initialize_objects(pr.d)
	
	if (args() <= 8) higherorder = 0 
    // Store input arguments into the myproblem structure
    pr.lxvars  = lxvars
    pr.lyvars  = lyvars
    pr.oyvars  = oyvars
    pr.params  = params
	pr.params_real  = params_real
    pr.glvls   = glvls 
    pr.cmdline = cmdline 
    pr.comparegroups = comparegroups 
	pr.higherorder = higherorder
	pr.higherorderitem = higherorderitem 

	get_coefficients(pr) 
	if (args() <= 8) get_anchorkey(pr)
	group_pathcoefs(pr)
	group_intercepts(pr)
	get_Beta(pr)
	get_Gamma(pr)
	get_Lambday(pr)
	if (length(pr.d.lambdax_a_bag) < length(pr.d.lambdax_bag)) get_Lambdaxna(pr)
	fill_matrices(pr)
	//Additional steps (e.g., build_matrices/equations(), standard errors) could follow
	return(pr.d)
}

// get_coefficients: Extracts path coefficients and intercepts from the raw params
void get_coefficients(struct myproblem scalar pr)
{
    real   matrix  s1, s2a, s2b, s3, s4a, s4b, hold1, hold2
    real scalar i 
    string matrix  params2, params3, params4
    
    // Remove means, variances, and covariances from params (we only want regression-like coefs)
    s1 = rowsum(strpos(pr.params,("var","mean","cov")))
    s1 = mm_cond(s1:>0,s1:-s1,s1:+1)
    params2 = select(pr.params, s1)
    
    // Separate coefficients by the specified group comparisons
    s2a  = strpos(params2,pr.glvls[pr.comparegroups[1]]) 
    s2b  = strpos(params2,pr.glvls[pr.comparegroups[2]]) 
    params3 = select(params2, s2a), select(params2, s2b)
    
    // Identify and remove endogenous intercept parameters (not needed for structural model)
    hold1 = J(rows(params3),length(pr.lyvars),.)
    for (i=1;i<=length(pr.lyvars);i++) hold1[1...,i] = rowsum(strpos(params3,pr.lyvars[i]))
    hold1 = rowsum(hold1)
    hold2 = rowsum(strpos(params3,"#"))
    s3 = J(rows(params3),1,0)
    for (i=1;i<=length(s3);i++) if (hold1[i] > 0 & hold2[i] == 0) s3[i] = 1
    s3 = mm_cond(s3:>0,s3:-s3,s3:+1)
    params4 = select(params3, s3)
    
    // Separate intercepts and path coefficients using presence of '#'
    s4a = rowsum(strpos(params4,"#"))
    s4b = mm_cond(s4a:>0,s4a:-s4a,s4a:+1)
    pr.d.pathcoefs  = select(params4, s4a)
    pr.d.intercepts = select(params4, s4b)
}

// get_anchorkey: Identify anchors (fixed factor loadings set at 1) from the command line
void get_anchorkey(struct myproblem scalar pr)
{
    real  matrix  s1,s2
    real scalar pos,i,j
    string scalar cmdline2,matched, anchor, latentvar 
    string matrix equations,equations2,equations3,check1,check2 
    
    // Extract portion of cmdline from "(" to "," to isolate model equations
    cmdline2 = substr(pr.cmdline,strpos(pr.cmdline,"("),strpos(pr.cmdline,",")-strpos(pr.cmdline,"("))
    equations = ustrsplit(cmdline2, "\)\s*\(")'
    equations = subinstr(subinstr(equations,"(",""),")","")
    
    // Check which equations relate to oyvars
    check1 = J(rows(equations),1,pr.oyvars)
    s1 = rowsum(strpos(equations, check1)) // Creating a mask
    equations2 = select(equations, s1)

    // Now filter down to those involving lxvars
    check2 = J(rows(equations2),1,pr.lxvars)
    s2 = rowsum(strpos(equations2, check2))
    equations3 = select(equations2, s2)
    
    // For each equation, look for patterns like LV@1 indicating anchors
    for(i=1;i<=rows(equations3);i++) {
        if (regexm(equations3[i], "([A-Za-z0-9]+@1)")) {
            matched = regexs(1)
            pos = strpos(matched,"@") - 1 
            anchor = substr(matched,1,pos)
        }
        // Determine which latent variable this anchor belongs to
        for (j=1;j<=length(pr.lxvars);j++) {
            if (strpos(equations3[i], pr.lxvars[j]) > 0) {
                latentvar  = pr.lxvars[j]
                asarray(pr.d.anchorkey,latentvar,anchor)
                break
            }	
        }		
    }
}

// group_pathcoefs: Separate structural coefficients (beta, gamma) from loadings (lambda)
void group_pathcoefs(struct myproblem scalar pr)
{
    real matrix res,a,s1,s1a,s1b,s2,s2a,s2b,s2c,s2d
    real scalar i 
    string matrix check1,structcoefs,loadings

    // Identify structural coefficients: these relate to oyvars
    check1 = J(rows(pr.d.pathcoefs),1,pr.oyvars)
    s1 = rowmissing(editvalue(indexnot(check1, pr.d.pathcoefs[1...,1]),0,.))
    s1 = mm_cond(s1:>0,s1:-s1,s1:+1)
    structcoefs = select(pr.d.pathcoefs,s1)
    
    // Determine which structural coefs pertain to endogenous variables (lyvars)
    res = J(rows(structcoefs),length(pr.lyvars),.)
    for (i=1;i<=length(pr.lyvars);i++) {
        a = strpos(structcoefs[1...,1],pr.lyvars[i]) 
        res[1...,i] = mm_cond(a:>1,a:-a:+ 1,a:-a)
    }
    s1a = rowsum(res)
    
    // Beta coefs: among endogenous vars
    pr.d.beta_bag = select(structcoefs,s1a)
    // Gamma coefs: exogenous to endogenous (complement of beta)
    s1b = mm_cond(s1a:>0,s1a:-s1a,s1a:+1)
    pr.d.gamma_bag = select(structcoefs,s1b)
    
    // Identify loadings: complement of structural coefs
    s2 = mm_cond(s1:>0,s1:-s1,s1:+1)
    loadings = select(pr.d.pathcoefs,s2)
	
	if (pr.higherorder != 0) { // If a higher-order model is specified, use specialized grouping routines. Note that there are are no non-anchored exogenous items 
		// Split into anchor item and loadings are for endogenous indicators (lyvars) 
		s2a = strpos(loadings[1...,1],pr.higherorderitem)
		pr.d.lambdax_a_bag = select(loadings,s2a)
		
		s2b = mm_cond(s2a:>0,s2a:-s2a,s2a:+1)
		pr.d.lambday_bag  = select(loadings,s2b)
	}
	else {
		// Determine which loadings are for endogenous indicators (lyvars)
		res = J(rows(loadings),length(pr.lyvars),.)
		for (i=1;i<=length(pr.lyvars);i++) {
			a = strpos(loadings[1...,1],pr.lyvars[i]) 
			res[1...,i] = mm_cond(a:>1,a:-a:+ 1,a:-a)
		}
		s2a = rowsum(res)
		// Lambday bag: loadings for endogenous indicators
		pr.d.lambday_bag = select(loadings,s2a)
		
		// Exogenous loadings: opposite mask
		s2b = mm_cond(s2a:>0,s2a:-s2a,s2a:+1)
		pr.d.lambdax_bag = select(loadings,s2b)
		// Determine anchored vs non-anchored among exogenous loadings
		res = J(rows(pr.d.lambdax_bag),length(pr.lxvars),.)
		for (i=1;i<=length(pr.lxvars);i++) {
			res[1...,i] = strpos(pr.d.lambdax_bag[1...,1],asarray(pr.d.anchorkey,pr.lxvars[i])) 
		}
		s2c = rowsum(res)
		// Anchored exogenous loadings
		pr.d.lambdax_a_bag = select(pr.d.lambdax_bag,s2c)
	
		// Non-anchored exogenous loadings
		s2d = mm_cond(s2c:>0,s2c:-s2c,s2c:+1)
		if (length(pr.d.lambdax_a_bag) < length(pr.d.lambdax_bag)) pr.d.lambdax_na_bag = select(pr.d.lambdax_bag,s2d)
	}
}

// group_intercepts: Organize intercept parameters into y, x anchored, and x non-anchored sets
void group_intercepts(struct myproblem scalar pr)
{
    real matrix pos1,pos2,pos3,s1,s2,s3
    string matrix check1,check2,check3,ys,xas,xnas
    
    // Extract variable names from lambday_bag for endogenous intercepts
    pos1 = strpos(pr.d.lambday_bag[1...,1],":")
    ys = substr(pr.d.lambday_bag[1...,1],J(rows(pr.d.lambday_bag),1,1),pos1:-1) 
    check1 = J(rows(pr.d.intercepts),1,ys')
    s1 = rowmissing(editvalue(indexnot(check1, pr.d.intercepts[1...,1]),0,.)) 
    pr.d.upsilon_y_bag = select(pr.d.intercepts,s1)
    
    // Anchored exogenous intercepts
    pos2 = strpos(pr.d.lambdax_a_bag[1...,1],":")
    xas = substr(pr.d.lambdax_a_bag[1...,1],J(rows(pr.d.lambdax_a_bag),1,1),pos2:-1) 
    check2 = J(rows(pr.d.intercepts),1,xas')
    s2 = rowmissing(editvalue(indexnot(check2, pr.d.intercepts[1...,1]),0,.)) 
    pr.d.upsilon_xa_bag  = select(pr.d.intercepts,s2)
    
	if (pr.d.lambdax_na_bag != J(0,0,"")) {
		// Non-anchored exogenous intercepts
		pos3 = strpos(pr.d.lambdax_na_bag[1...,1],":")
		xnas = substr(pr.d.lambdax_na_bag[1...,1],J(rows(pr.d.lambdax_na_bag),1,1),pos3:-1) 
		check3 = J(rows(pr.d.intercepts),1,xnas')
		s3 = rowmissing(editvalue(indexnot(check3, pr.d.intercepts[1...,1]),0,.)) 
		pr.d.upsilon_xna_bag  = select(pr.d.intercepts,s3)
	}
}
void get_Beta(struct myproblem scalar pr) 
{
	real matrix B 
	real scalar i,j,k
	
	pr.d.Beta1 = J(length(pr.lyvars),length(pr.lyvars),"")
	pr.d.Beta2 = J(length(pr.lyvars),length(pr.lyvars),"")
	for(i=1;i<=rows(pr.d.Beta1);i++) {
		for(j=1;j<=cols(pr.d.Beta1);j++) {
			if (pr.lyvars[i] != pr.lyvars[j]) {
				for(k=1;k<=rows(pr.d.beta_bag);k++) {
					B = strpos(pr.d.beta_bag[k,1],(pr.lyvars[i],pr.lyvars[j]))
					if (rowmissing(editvalue(B,0,.)) == 0 & B[1] < B[2]) {
						pr.d.Beta1[i,j] = "_b[" + pr.d.beta_bag[k,1] + "]"
						pr.d.Beta2[i,j] = "_b[" + pr.d.beta_bag[k,2] + "]"
						asarray(pr.d.parammatrices_string,"Beta1",pr.d.Beta1)
						asarray(pr.d.parammatrices_string,"Beta2",pr.d.Beta2)
					}
				}
			}
		}	
	}
}
void get_Gamma(struct myproblem scalar pr) 
{
	real matrix B 
	real scalar i,j,k
	
	pr.d.Gamma1 = J(length(pr.lyvars),length(pr.lxvars),"")
	pr.d.Gamma2 = J(length(pr.lyvars),length(pr.lxvars),"")
	for(i=1;i<=rows(pr.d.Gamma1);i++) {
		for(j=1;j<=cols(pr.d.Gamma1);j++) {
			for(k=1;k<=rows(pr.d.gamma_bag);k++) {
				B = strpos(pr.d.gamma_bag[k,1],(pr.lyvars[i],pr.lxvars[j]))
				if (rowmissing(editvalue(B,0,.)) == 0 & B[1] < B[2]) {
					pr.d.Gamma1[i,j] = "_b[" + pr.d.gamma_bag[k,1] + "]"
					pr.d.Gamma2[i,j] = "_b[" + pr.d.gamma_bag[k,2] + "]"
					asarray(pr.d.parammatrices_string,"Gamma1",pr.d.Gamma1)
					asarray(pr.d.parammatrices_string,"Gamma2",pr.d.Gamma2)
				}
			}
		}	
	}
}
void get_Lambday(struct myproblem scalar pr) 
{
	real matrix   B, pos
	real scalar   i,j,k
	string matrix yobserved 
	
	pr.d.Lambday1 = J(rows(pr.d.lambday_bag),length(pr.lyvars),"")
	pr.d.Lambday2 = J(rows(pr.d.lambday_bag),length(pr.lyvars),"")
	pos = strpos(pr.d.lambday_bag[1...,1],":"):-1
	yobserved = substr(pr.d.lambday_bag[1...,1],J(rows(pos),1,1),pos)
	for(i=1;i<=rows(pr.d.Lambday1);i++) {
		for(j=1;j<=cols(pr.d.Lambday1);j++) {
			for(k=1;k<=rows(pr.d.lambday_bag);k++) {
				B = strpos(pr.d.lambday_bag[k,1],(yobserved[i],pr.lyvars[j]))
				if (rowmissing(editvalue(B,0,.)) == 0 & B[1] < B[2]) {
					pr.d.Lambday1[i,j] = "_b[" + pr.d.lambday_bag[k,1] + "]"
					pr.d.Lambday2[i,j] = "_b[" + pr.d.lambday_bag[k,2] + "]"
					asarray(pr.d.parammatrices_string,"Lambday1",pr.d.Lambday1)
					asarray(pr.d.parammatrices_string,"Lambday2",pr.d.Lambday2)
				}
			}
		}	
	}
}

void get_Lambdaxna(struct myproblem scalar pr) 
{
	real matrix   B, pos
	real scalar   i,j,k
	string matrix xobserved 
	
	pr.d.Lambdaxna1 = J(rows(pr.d.lambdax_na_bag),length(pr.lxvars),"")
	pr.d.Lambdaxna2 = J(rows(pr.d.lambdax_na_bag),length(pr.lxvars),"")
	pos = strpos(pr.d.lambdax_na_bag[1...,1],":"):-1
	xobserved = substr(pr.d.lambdax_na_bag[1...,1],J(rows(pos),1,1),pos)
	for(i=1;i<=rows(pr.d.Lambdaxna1);i++) {
		for(j=1;j<=cols(pr.d.Lambdaxna1);j++) {
			for(k=1;k<=rows(pr.d.lambdax_na_bag);k++) {
				B = strpos(pr.d.lambdax_na_bag[k,1],(xobserved[i],pr.lxvars[j]))
				if (rowmissing(editvalue(B,0,.)) == 0 & B[1] < B[2]) {
					pr.d.Lambdaxna1[i,j] = "_b[" + pr.d.lambdax_na_bag[k,1] + "]"
					pr.d.Lambdaxna2[i,j] = "_b[" + pr.d.lambdax_na_bag[k,2] + "]"
					asarray(pr.d.parammatrices_string,"Lambdaxna1",pr.d.Lambdaxna1)
					asarray(pr.d.parammatrices_string,"Lambdaxna2",pr.d.Lambdaxna2)
				}
			}
		}	
	}
}
void fill_matrices(struct myproblem scalar pr)
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
end












/*

Python code to write equations

from sympy import Symbol, Matrix, eye, simplify

# Define symbolic variables for matrix entries using Symbol (not symbols)
x2_F2 = Symbol("_b[x2:0bn.t#c.F2]")
x3_F3 = Symbol("_b[x3:0bn.t#c.F3]")
F2_F1 = Symbol("_b[F2:0bn.t#c.F1]")
F3_F2 = Symbol("_b[F3:0bn.t#c.F2]")

# Define Lambda_y as a symmetric matrix
Lambda_y = Matrix([
    [x2_F2, 0],
    [0, x3_F3]
])

# Define Gamma1 as a column vector
Gamma1 = Matrix([
    [F2_F1],
    [0]
])

# Define Beta1 matrix
Beta1 = Matrix([
    [0, 0],
    [F3_F2, 0]
])

# Create identity matrix (I)
I = eye(2)

try:
    # Perform the operation: Lambda_y * (I - Beta1)^(-1) * Gamma1
    result = Lambda_y * (I - Beta1).inv() * Gamma1
    result = simplify(result)

    print("Result of Lambda_y * (I - Beta1)^-1 * Gamma1:")
    print(result)

except Exception as e:
    print(f"An error occurred: {e}")


	
from sympy import Symbol, Matrix, eye, simplify

# Define symbolic variables for Lambda_y entries (non-zero entries)
s1_r01_l02_s1 = Symbol("_b[s1_r01_l02:0bn.t#c.S1]")
s1_r02_l01_s1 = Symbol("_b[s1_r02_l01:0bn.t#c.S1]")
s1_r02_l02_s1 = Symbol("_b[s1_r02_l02:0bn.t#c.S1]")

s2_r01_l01_s2 = Symbol("_b[s2_r01_l01:0bn.t#c.S2]")
s2_r01_l02_s2 = Symbol("_b[s2_r01_l02:0bn.t#c.S2]")
s2_r02_l01_s2 = Symbol("_b[s2_r02_l01:0bn.t#c.S2]")
s2_r02_l02_s2 = Symbol("_b[s2_r02_l02:0bn.t#c.S2]")

s3_r01_l01_s3 = Symbol("_b[s3_r01_l01:0bn.t#c.S3]")
s3_r01_l02_s3 = Symbol("_b[s3_r01_l02:0bn.t#c.S3]")
s3_r02_l01_s3 = Symbol("_b[s3_r02_l01:0bn.t#c.S3]")
s3_r02_l02_s3 = Symbol("_b[s3_r02_l02:0bn.t#c.S3]")

s4_r01_l01_s4 = Symbol("_b[s4_r01_l01:0bn.t#c.S4]")
s4_r01_l02_s4 = Symbol("_b[s4_r01_l02:0bn.t#c.S4]")
s4_r02_l01_s4 = Symbol("_b[s4_r02_l01:0bn.t#c.S4]")
s4_r02_l02_s4 = Symbol("_b[s4_r02_l02:0bn.t#c.S4]")

s5_r01_l01_s5 = Symbol("_b[s5_r01_l01:0bn.t#c.S5]")
s5_r01_l02_s5 = Symbol("_b[s5_r01_l02:0bn.t#c.S5]")
s5_r02_l01_s5 = Symbol("_b[s5_r02_l01:0bn.t#c.S5]")
s5_r02_l02_s5 = Symbol("_b[s5_r02_l02:0bn.t#c.S5]")

s6_r01_l01_s6 = Symbol("_b[s6_r01_l01:0bn.t#c.S6]")
s6_r01_l02_s6 = Symbol("_b[s6_r01_l02:0bn.t#c.S6]")
s6_r02_l01_s6 = Symbol("_b[s6_r02_l01:0bn.t#c.S6]")
s6_r02_l02_s6 = Symbol("_b[s6_r02_l02:0bn.t#c.S6]")

# Define Gamma1 entries
S1_F1 = Symbol("_b[S1:0bn.t#c.F1]")
S2_F1 = Symbol("_b[S2:0bn.t#c.F1]")
S3_F1 = Symbol("_b[S3:0bn.t#c.F1]")
F2_F1 = Symbol("_b[F2:0bn.t#c.F1]")

# Define Beta1 entries
S4_F2 = Symbol("_b[S4:0bn.t#c.F2]")
S5_F2 = Symbol("_b[S5:0bn.t#c.F2]")
S6_F2 = Symbol("_b[S6:0bn.t#c.F2]")

# Construct Lambda_y (23 x 7)
# Rows are indexed from 1 to 23 and columns from 1 to 7 as shown in the given matrix
Lambda_y_data = [
 [s1_r01_l02_s1,              0,              0,              0,              0,              0, 0],
 [s1_r02_l01_s1,              0,              0,              0,              0,              0, 0],
 [s1_r02_l02_s1,              0,              0,              0,              0,              0, 0],
 [0,              s2_r01_l01_s2,              0,              0,              0,              0, 0],
 [0,              s2_r01_l02_s2,              0,              0,              0,              0, 0],
 [0,              s2_r02_l01_s2,              0,              0,              0,              0, 0],
 [0,              s2_r02_l02_s2,              0,              0,              0,              0, 0],
 [0,                           0, s3_r01_l01_s3,              0,              0,              0, 0],
 [0,                           0, s3_r01_l02_s3,              0,              0,              0, 0],
 [0,                           0, s3_r02_l01_s3,              0,              0,              0, 0],
 [0,                           0, s3_r02_l02_s3,              0,              0,              0, 0],
 [0,                           0,              0, s4_r01_l01_s4,              0,              0, 0],
 [0,                           0,              0, s4_r01_l02_s4,              0,              0, 0],
 [0,                           0,              0, s4_r02_l01_s4,              0,              0, 0],
 [0,                           0,              0, s4_r02_l02_s4,              0,              0, 0],
 [0,                           0,              0,              0, s5_r01_l01_s5,              0, 0],
 [0,                           0,              0,              0, s5_r01_l02_s5,              0, 0],
 [0,                           0,              0,              0, s5_r02_l01_s5,              0, 0],
 [0,                           0,              0,              0, s5_r02_l02_s5,              0, 0],
 [0,                           0,              0,              0,              0, s6_r01_l01_s6, 0],
 [0,                           0,              0,              0,              0, s6_r01_l02_s6, 0],
 [0,                           0,              0,              0,              0, s6_r02_l01_s6, 0],
 [0,                           0,              0,              0,              0, s6_r02_l02_s6, 0]
]

Lambda_y = Matrix(Lambda_y_data)

# Construct Gamma1 (7 x 1)
Gamma1 = Matrix([
 [S1_F1],
 [S2_F1],
 [S3_F1],
 [0],
 [0],
 [0],
 [F2_F1]
])

# Construct Beta1 (7 x 7)
Beta1 = Matrix([
 [0, 0, 0, 0, 0, 0, 0],
 [0, 0, 0, 0, 0, 0, 0],
 [0, 0, 0, 0, 0, 0, 0],
 [0, 0, 0, 0, 0, 0, S4_F2],
 [0, 0, 0, 0, 0, 0, S5_F2],
 [0, 0, 0, 0, 0, 0, S6_F2],
 [0, 0, 0, 0, 0, 0, 0]
])

# Create identity matrix of the same dimension as Beta1 (7 x 7)
I = eye(7)

try:
    # Compute Lambda_y * (I - Beta1)^-1 * Gamma1
    result = Lambda_y * (I - Beta1).inv() * Gamma1
    result = simplify(result)

    print("Result of Lambda_y * (I - Beta1)^-1 * Gamma1:")
    print(result)
except Exception as e:
    print(f"An error occurred: {e}")
	
	
Result of Lambda_y * (I - Beta1)^-1 * Gamma1:
Matrix([[_b[S1:0bn.t#c.F1]*_b[s1_r01_l02:0bn.t#c.S1]], 
[_b[S1:0bn.t#c.F1]*_b[s1_r02_l01:0bn.t#c.S1]], 
[_b[S1:0bn.t#c.F1]*_b[s1_r02_l02:0bn.t#c.S1]], 
[_b[S2:0bn.t#c.F1]*_b[s2_r01_l01:0bn.t#c.S2]], 
[_b[S2:0bn.t#c.F1]*_b[s2_r01_l02:0bn.t#c.S2]], 
[_b[S2:0bn.t#c.F1]*_b[s2_r02_l01:0bn.t#c.S2]], 
[_b[S2:0bn.t#c.F1]*_b[s2_r02_l02:0bn.t#c.S2]], 
[_b[S3:0bn.t#c.F1]*_b[s3_r01_l01:0bn.t#c.S3]], 
[_b[S3:0bn.t#c.F1]*_b[s3_r01_l02:0bn.t#c.S3]], 
[_b[S3:0bn.t#c.F1]*_b[s3_r02_l01:0bn.t#c.S3]], 
[_b[S3:0bn.t#c.F1]*_b[s3_r02_l02:0bn.t#c.S3]], 
[_b[F2:0bn.t#c.F1]*_b[S4:0bn.t#c.F2]*_b[s4_r01_l01:0bn.t#c.S4]], 
[_b[F2:0bn.t#c.F1]*_b[S4:0bn.t#c.F2]*_b[s4_r01_l02:0bn.t#c.S4]], 
[_b[F2:0bn.t#c.F1]*_b[S4:0bn.t#c.F2]*_b[s4_r02_l01:0bn.t#c.S4]], 
[_b[F2:0bn.t#c.F1]*_b[S4:0bn.t#c.F2]*_b[s4_r02_l02:0bn.t#c.S4]], 
[_b[F2:0bn.t#c.F1]*_b[S5:0bn.t#c.F2]*_b[s5_r01_l01:0bn.t#c.S5]], 
[_b[F2:0bn.t#c.F1]*_b[S5:0bn.t#c.F2]*_b[s5_r01_l02:0bn.t#c.S5]], 
[_b[F2:0bn.t#c.F1]*_b[S5:0bn.t#c.F2]*_b[s5_r02_l01:0bn.t#c.S5]], 
[_b[F2:0bn.t#c.F1]*_b[S5:0bn.t#c.F2]*_b[s5_r02_l02:0bn.t#c.S5]], 
[_b[F2:0bn.t#c.F1]*_b[S6:0bn.t#c.F2]*_b[s6_r01_l01:0bn.t#c.S6]], 
[_b[F2:0bn.t#c.F1]*_b[S6:0bn.t#c.F2]*_b[s6_r01_l02:0bn.t#c.S6]], 
[_b[F2:0bn.t#c.F1]*_b[S6:0bn.t#c.F2]*_b[s6_r02_l01:0bn.t#c.S6]], 
[_b[F2:0bn.t#c.F1]*_b[S6:0bn.t#c.F2]*_b[s6_r02_l02:0bn.t#c.S6]]])


























