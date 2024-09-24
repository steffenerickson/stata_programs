//----------------------------------------------------------------------------//
// Title   : Permutations and Combinations
// Author  : Steffen Erickson 
// Date    : 8/10/24
//----------------------------------------------------------------------------//

//-----------------------//
// Permutations
//-----------------------//
capture program drop permutations
program permutations, rclass 
	syntax,  N(integer) K(integer)
	tempname tempframe a v 
	mata st_matrix("`a'", J(1,`k',range(1,`n',1)))
	mkf `tempframe'
	frame `tempframe' {
		svmat `a' 
		fillin *
		drop _fillin
		qui ds
		tokenize `r(varlist)'
		local x = 2
		local k_1 = `k' - 1
		forvalues y = 1/`k_1' {
			forvalues z = 1/`y' {
				capture drop if ``z'' == ``x''
			}
			local++x	
		}
		mkmat *, matrix(`v') rowprefix(perm)
	}
	return matrix permutations = `v'
end 

//-----------------------//
// Combinations 
//-----------------------//

capture program drop combinations 
program combinations , rclass 
	syntax , N(integer) K(integer)
	tempname tempframe a v
	mata st_matrix("`a'", J(1,`k',range(1,`n',1)))
	mkf `tempframe'
	frame `tempframe' {
		svmat `a' //, name(`v')
		fillin *
		drop _fillin
		qui ds 
		tokenize `r(varlist)'
		local i = 1
		local j = 1
		local count = `:word count `r(varlist)''
		while (`j' < `count') {
			local j = `i' + 1
			drop if ``i'' >= ``j''
			local++i
		}
		mkmat * , matrix(`v') rowprefix(comb)
	}
	return matrix combinations  = `v'
end 


//permutations , n(6) k(4)
//mat li r(permutations)
//
//combinations , n(6) k(4)
//mat li r(combinations)



/* 
clear all 
mata a = J(1,3,range(1,4,1))
mata a
mata st_matrix("a",a)
svmat a 
set trace on 
fillin *

