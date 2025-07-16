

clear all 

clear 
frame reset 
mkf designmatrix 
qui frame designmatrix :crosseddesign a b c, facetlevels(30 4 10)
qui frame designmatrix : gen y = runiform() * a + runiform() * b + runiform() * c + runiform()
qui frame designmatrix : manova y = a##b##c
ereturn list 
qui frame designmatrix : gen randnum = runiform()
qui frame designmatrix : sort randnum
qui frame designmatrix : gen n = _n
qui frame designmatrix : keep if n < 800
qui frame designmatrix : manova y = a c b 

qui frame designmatrix : fvexpand a##b##c
return list

mkmat i.a##i.b##i.c, matrix(X)
mat list X

cap program drop expandvarlist
program expandvarlist

	syntax varlist(fv)
	di "`varlist'"
	
end 
frame designmatrix : expandvarlist a#(b c)


frame designmatrix : expandvarlist a##b##c





frame change designmatrix


tab c b


mata : permn(10,10)



clear 
frame reset 
mkf designmatrix 
frame designmatrix: crosseddesign a b , facetlevels(100 100)
frame designmatrix: randobselecter b, nobs(10)
frame designmatrix: tab a b
