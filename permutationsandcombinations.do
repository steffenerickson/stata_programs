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
	local k_1 = `k' - 1
	mata st_matrix("`a'",range(1,`n',1))
	mkf `tempframe'
	frame `tempframe' {
		svmat `a', name(`v')
		forvalues i = 1/`k_1' {
			local j = `i' + 1
			gen `v'`j' = `v'`i'
		}
		fillin *
		drop _fillin
		qui ds
		tokenize `r(varlist)'
		local x = 2
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
	local k_1 = `k' - 1
	mata st_matrix("`a'",range(1,`n',1))
	mkf `tempframe'
	frame `tempframe' {
		svmat `a', name(`v')
		forvalues i = 1/`k_1' {
			local j = `i' + 1
			gen `v'`j' = `v'`i'
		}
		fillin *
		drop _fillin
		qui ds 
		tokenize `r(varlist)'
		local i = 1
		local j = 1
		while (`j' < `:word count `r(varlist)'') {
			local j = `i' + 1
			drop if ``i'' >= ``j''
			local++i
		}
		mkmat * , matrix(`v') rowprefix(comb)
	}
	return matrix combinations  = `v'
end 


permutations , n(3) k(2)
mat li r(permutations)

combinations , n(3) k(2)
mat li r(combinations)
