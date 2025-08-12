


// Composite Reliability 

mata objmeasurement = "p"
mata effects = c.effects
mata varlist = c.varlist
mata covcomps = c.covcomps 
mata w = J(length(varlist),1,1/length(varlist))
mata compositevariances = J(length(effects),1,.)
mata 
for(i=1;i<=length(effects);i++) {
	name = "emcp" + strofreal(i) 
	M = asarray(covcomps, name)
	M
	compositevariances[i] = w'*M*w
}
end 
mata compositevariances

// ---------------------------------------------------------------------------//

mata s1 = effects':==objmeasurement
mata s2 = c.flip_select(s1)
mata select(compositevariances,s1)
mata erroreffects   = select(effects',s2)
mata errorvariances = select(compositevariances,s2)
mata s3   = strpos(erroreffects ,objmeasurement)
mata relerroreffects   = select(erroreffects,s3)
mata relerrorvariances = select(errorvariances,s3)
mata a = subinstr(subinstr(subinstr(erroreffects,"|"," "), "#"," "),objmeasurement," ")
mata a = strtrim(a) 


mata 
stringlengths = strlen(a)
selectvec = (stringlengths :== max(stringlengths))
uniqueffects = tokens(select(a,selectvec))
colspan = cols(uniqueffects)
res = J(0,colspan,"")
for (i=1;i<=rows(a);i++) {
	difference = colspan - cols(tokens(a[i]))
	if (difference !=0) temp = a[i], J(1,difference,"")
	else temp = tokens(a[i])
	res = res \ temp 
}

res
uniqueffects
relerrorvariances
end 

mata 
outerarray = asarray_create() 
for (i=1;i<=cols(uniqueffects);i++) {
	inner = asarray_create() 
	asarray(outerarray,uniqueffects[i],inner)
}
end 

mata asarray(outerarray,"t")
mata asarray(outerarray,"r")
mata
for (i=1;i<=cols(uniqueffects);i++) {
	numres = (res:==uniqueffects[i])
	printf(strofreal(i) + " Matrices")
	for (j=1;j<=6;j++) {
		M = numres * j 
		asarray(asarray(outerarray,uniqueffects[i]),strofreal(j),M)
	}
}
end 

// Step 2: Recursive function to add combinations
mata mata drop combine_matrices()

// Step 1: Get list of group names (e.g., "1", "2", ..., "K")
mata group_keys = asarray_keys(outerarray)
mata mats = asarray_create() 
mata 
void combine_matrices(transmorphic outerarray, 
                      string scalar path, 
					  real matrix acc, 
					  real scalar depth,
					  string vector group_keys,
					  real scalar projectnum,
					  transmorphic mats) 
{
	
	num_groups = length(group_keys)
    if (depth > num_groups) {
        // Base case: all groups processed
        printf("Combination: %s\n", path)
        acc
		asarray(mats,path,acc)
		return
    }
    group = group_keys[depth]
    for (i = 1; i <= projectnum; i++) {
        M = asarray(asarray(outerarray, group), strofreal(i))
        new_path = path + group + "[" + strofreal(i) + "] "
        combine_matrices(outerarray,new_path, acc + M, depth + 1,group_keys, projectnum, mats)
    }

}

real matrix rowwise_nonzero_product(real matrix X) {
    real matrix result 
	
	
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

real vector extract_bracket_numbers(string scalar s) 
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


end 

mata 
// Step 3: Start recursion with zero matrix
// Grab any one matrix to get the correct shape
first = asarray(asarray(outerarray, group_keys[1]), "1")
zero = first * 0
combine_matrices(outerarray,"", zero, 1,group_keys,6,mats)

res = J(0,cols(first)+1,.)
for (loc=asarray_first(mats); loc!=NULL; loc=asarray_next(mats, loc)) {
	nums = extract_bracket_numbers(asarray_key(mats, loc))
	div = rowwise_nonzero_product(asarray_contents(mats, loc))
	temp = nums, colsum(errorvariances:/div) 
	res = res \ temp 
}
end


  r[1] t[1] 
       1   2
    +---------+
  1 |  1   0  |
  2 |  1   1  |
    +---------+
                 1
    +---------------+
  1 |  .0117101651  |
  2 |  .0461092863  |
    +---------------+








