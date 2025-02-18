*! Program to compute standard errors for weights of variables used to form an emergent (composite) latent variable via the Henseler–Ogasawara specification
*! Version 1
*! Steffen Erickson

version 18 

program weightSEs, rclass
	syntax , Matrix(string)
	
	tempname temp result x rn cn weqs pos 
	
	python: from __main__ import compositeweights
	
	* create symbol matrix 
	local rownames: rowfullnames `matrix'
	local colnames: colfullnames `matrix' 
	mata `x' = st_matrix("`matrix'")
	mata `rn' = tokens(st_local("rownames"))
	mata `cn' = tokens(st_local("colnames"))
	mata symbolmat = symbolmat(`x',`rn',`cn')
	
	* symbolically invert symbol matrix in python 
	quietly python: compositeweights("symbolmat")
	
	* use equations to compute SEs via delta method
	mata: `weqs' = ustrsplit(st_local("eq"),",")
	mata: st_local("numrows",strofreal(cols(`weqs')))
	forvalues i = 1/`numrows' {
		mata st_local("weq`i'",`weqs'[`i'])
		quietly nlcom `weq`i''
		mat `temp' = r(table)[1..6,1]'
		mat `result' = (nullmat(`result') \ `temp')
	}
	
	* results 
	mata `pos' = strpos(`rn'[1],":")
	mata st_local("vars",invtokens(substr(`rn' ,`pos'+1)))
	foreach x of local vars {
		local name   "`x'"
		local names `" `names' "`name'" "'
	}
	matrix rownames `result' = `names'
	matlist `result'
	return matrix table = `result'
end

mata 
string matrix symbolmat(real matrix X,
                       string vector rn,
					   string vector cn)
{
	real scalar   	pos1,pos2
	string vector 	rn2,cn2 
	string matrix 	symbolmat
	
	pos1 = strpos(rn[1],":")
	pos2 = strpos(cn[1],":")
	rn2  = substr(rn ,pos1+1)
	cn2  = substr(cn ,pos2+1)
	
	symbolmat = J(rows(X),cols(X),"")
	for(i=1;i<=rows(X);i++){
		for(j=1;j<=cols(X);j++){
			if (X[i,j] != 0) symbolmat[i,j] = "_b[" + rn2[i] + ":" + cn2[j] + "]"
		}
	}
		
	return(symbolmat)
}
end 

python 

import numpy as np
from sympy import Symbol, Matrix
from sfi import Macro, Mata

class compositeweights:	
	def __init__(self, matrixname):
		self._matrix = np.matrix(Mata.get(f"{matrixname}"))
		self._symbolic_matrix = self._createsymbolic()
		self._composite_weights = self._compweights()
		self._return_as_macro()

	def _createsymbolic(self):
		variables = {}
		for index, element in np.ndenumerate(self._matrix):
			var_name = f"{index[0]}_{index[1]}"
			if element == "":
				variables[var_name] = 0
			else:
				variables[var_name] = Symbol(str(element))
		rows, cols = self._matrix.shape
		return Matrix(rows, cols, lambda i, j: variables[f"{i}_{j}"])

	def _compweights(self):
		return self._symbolic_matrix.inv()[0, :]

	def _return_as_macro(self):
		nested_list = self._composite_weights.tolist()
		flattened = [item for sublist in nested_list for item in sublist]
		string = ", ".join(map(str, flattened))
		Macro.setLocal("eq", string)
end 




