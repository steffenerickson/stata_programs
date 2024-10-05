//----------------------------------------------------------------------------//
// MvGstudy Demo 
// Steffen Erickson, 10/5/24
//----------------------------------------------------------------------------//

// Commands 
clear all
include mvgstudy.ado
program stubs_for_reshape, rclass 
	syntax varlist, Start(integer) End(integer)
	tempname a b 
	mata `a' = tokens(st_local("varlist"))
	mata `b' = J(1,cols(`a'),"")
	mata for(i=1;i<=(cols(`a'));i++) `b'[1,i] = substr(`a'[1,i], `start', `end')
	mata `b' = uniqrows(`b'')'
	mata st_local("stubs",invtokens(`b'))
	return local stubs `stubs' 
end 

// Data
input x111 x112 x121 x122 x131 x132 x211 x212 x221 x222 x231 x232
1 3 3 1 3 2 1 2 2 2 2 1 1
2 2 1 1 1 1 1 1 1 1 1 1 1
3 1 2 2 2 1 1 1 1 2 2 1 1
4 1 2 2 2 1 1 1 2 2 2 2 1
end 
gen p = _n 
stubs_for_reshape x* , s(1) e(3)
reshape long `r(stubs)' , i(p) j(r)
stubs_for_reshape x* , s(1) e(2)
reshape long `r(stubs)' , i(p r) j(t)

// Method
mvgstudy (x* = p t r|t p#t p#(r|t))
matrix T = r(covcomps1)
matrix E = r(covcomps2) + r(covcomps3) + r(covcomps4) + r(covcomps5)
matrix w = (.5\.5)
matrix var_t = w'*T*w
matrix var_e = w'*E*w
di var_t[1,1] / (var_t[1,1] + var_e[1,1])














