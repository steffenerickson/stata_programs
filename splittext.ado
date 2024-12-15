program splittext, nclass 
	syntax varlist(max=1), Id(varlist) Regex(string)
	tempname hold data res 
	
	mkf `hold' 
	quietly {
		sort `id'
		egen id = group(`id')
		mata `data' = st_sdata(.,"`varlist'")	
		mata `res'  = splittext(`data',"`regex'") 
		frame `hold'{
			getmata (text id)  = `res'
			destring id , replace
			bysort id: gen linenum = _n
			tempfile expandedtext
			save `expandedtext'
		}
		
		drop `varlist'
		merge 1:m id using `expandedtext'
		drop _merge 
	}
end 

mata 
string matrix splittext(string matrix text, string scalar regex)
{
	string matrix temp1,temp2,res
	real scalar i 
	
	res = J(0,2,"")
	for (i=1;i<=rows(text);i++) {
		temp1 = ustrsplit(text[i],regex)'
		temp2 = temp1,J(rows(temp1),1,strofreal(i))
		res = res \ temp2
	}
	return(res) 
}
end 
