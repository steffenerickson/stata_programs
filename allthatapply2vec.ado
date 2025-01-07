//----------------------------------------------------------------------------//
// Title   : allthatapply2vec
// Purpose : command to vectorize variable containing lists of strings 
// Author  : Steffen Erickson 
// Date    : 5/27/24

//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
// Stata command 
//----------------------------------------------------------------------------//
capture program drop allthatapply2vec
program allthatapply2vec, nclass
	syntax varlist(min=1 max=1 string), STUB(string) REGularexpression(string)
	
	mata: object = driver(st_sdata(.,"`varlist'"),"`regularexpression'")
	mata: responses = object.responses 
	mata: matchmatrix = object.matchmatrix
	getmata (`stub'*) = matchmatrix
	mata: st_local("len",strofreal(length(object.responses)))
	forvalues i = 1/`len' {
		mata: st_local("varlabel",object.responses[`i'])
		label variable `stub'`i' "`varlabel'"
	}
end 

//----------------------------------------------------------------------------//
// Mata routine 
//----------------------------------------------------------------------------//

mata: 
mata clear 

//------------ Structures -----------------//
struct myproblem  {
	string colvector 		stringvec
	string scalar 			regularexpression
	struct derived_objects 	scalar dobj
}
struct derived_objects { 
	real scalar 	max
	string matrix 	stringmat
	string vector  	responses
	real matrix  	matchmatrix
}
//------------ Main routine -----------------//

struct derived_objects driver(string vector stringvec, string scalar regularexpression)
{
	struct myproblem scalar pr
	
	initialize_objects(pr.dobj)
	row_to_col(stringvec) 
	pr.stringvec = stringvec
	pr.regularexpression = regularexpression
	mainrountine(pr)
	return(pr.dobj)
}
void initialize_objects(struct derived_objects scalar dobj)
{
	dobj.max 		 = .
	dobj.stringmat 	 = J(0,0,"")
	dobj.responses	 = J(1,0,"")
	dobj.matchmatrix = J(0,0,.) 
}
void mainrountine(struct myproblem scalar pr)
{
	maxnumofstrings(pr)   // Determines the numbers of columns for the stringmat 
	splitstringmatrix(pr) // create a matrix where each row is a vector of string responses for an obseravtion 
	responsematrix(pr)    // Creates a vector of unique responses and a matrix that indicates if the observation gave the response 
}

//------------ Sub routines -----------------//

real scalar is_colvec(z) return(anyof(("colvector","scalar"),orgtype(z)))

void row_to_col(string vector v) 
{
	if (is_colvec(v) == 0) v = v'
}

void maxnumofstrings(struct myproblem scalar pr)
{
	real scalar		i,x
	real matrix		numofstrings 
	
	numofstrings = J(0,1,.)
	for (i=1;i<=rows(pr.stringvec);i++) {
		x = cols(ustrsplit(pr.stringvec[i], pr.regularexpression))
		numofstrings = numofstrings \ x
	}
	pr.dobj.max = max(numofstrings)
}

void splitstringmatrix(struct myproblem scalar pr)
{		
	real scalar		i, empty
	string matrix	temp,temp2
	
	pr.dobj.stringmat = J(0,pr.dobj.max,"")
	for (i=1;i<=rows(pr.stringvec);i++) {
		empty = pr.stringvec[i] == ""
		temp  = (empty == 1  ? J(1,pr.dobj.max,""): ustrsplit(pr.stringvec[i],pr.regularexpression)) 
		temp2 = (cols(temp) == pr.dobj.max ? temp : (temp , J(1,pr.dobj.max - cols(temp),"")))
		pr.dobj.stringmat = pr.dobj.stringmat \ temp2
	}
}

void responsematrix(struct myproblem scalar pr)
{
	real scalar		i,j,n, match
	real matrix		hasmatch,checkmatch
	
	pr.dobj.responses = select(uniqrows(vec(pr.dobj.stringmat))', uniqrows(vec(pr.dobj.stringmat))' :!= "")
	pr.dobj.matchmatrix = J(rows(pr.dobj.stringmat),0,.)
	for (n=1;n<=cols(pr.dobj.responses);n++) {
		hasmatch = J(0,1,.)
		for (i=1;i<=rows(pr.dobj.stringmat);i++) {
			checkmatch = J(0,1,.)
			for (j=1;j<=cols(pr.dobj.stringmat);j++) {
				checkmatch = checkmatch \ (pr.dobj.stringmat[i,j] == pr.dobj.responses[n] ? 1 : 0) 
			}
			match = max(checkmatch)
			hasmatch = hasmatch \ match 
		}
		pr.dobj.matchmatrix = pr.dobj.matchmatrix,hasmatch
	}
}

end



