




clear all 
use http://www.stata-press.com/data/r13/dfex, clear
include "/Users/steffenerickson/Documents/GitHub/stata_programs/00_mvgstudy/mvgstudy_rewrite.do"
gen year = year(dofm(month))
tostring year, gen(ytemp)
encode ytemp, gen(y)
gen m = month(dofm(month))

global varlist  unemp   hours   inc96   ipman  
foreach v of global varlist {
	tempvar temp 
	egen `temp' = std(`v')
	replace `v' = `temp'
}
mvgstudy ( $varlist = m y m#y)






clear all 
include "/Users/steffenerickson/Documents/GitHub/stata_programs/00_mvgstudy/mvgstudy_rewrite.do" 
use https://www.stata-press.com/data/r19/sem_2fmmby, clear 



rename (phyab* appear* peerrel* parrel*) (x1* x2* x3* x4*)


foreach var of varlist x* {
	tempvar temp 
	egen `temp' = std(`var')
	replace `var' = `temp'
	
}



alpha x1*
gen p = _n
reshape long x1 x2 x3 x4 , i(grade p) j(item)
mvgstudy (x1 = p p#item) 
mvdstudy, obj(p) errortype(relative) facetnum(10)
mvgstudy (x1 x2 x3 x4 = p p#item) 





mvgstudy (x1 = p item p#item) 
mvgstudy (x2 = p item p#item) 
mvgstudy (x3 = p item p#item) 
mvgstudy (x4 = p item p#item) 
mvdstudy, obj(p) errortype(relative) facetnum(4)


preserve
tempvar temp 
egen `temp' = mean(x4) , by(p)
replace x4 = `temp' if item == 4
mvgstudy (x1 x2 x3 x4 = p p#item) 
restore


use https://www.stata-press.com/data/r19/automiss, clear 
alpha price headroom rep78 trunk weight length turn displacement , std
di  r(alpha) 
keep price headroom rep78 trunk weight length turn displacement
gen p = _n 
foreach var in price headroom rep78 trunk weight length turn displacement {
	tempvar temp 
	egen `temp' = std(`var')
	replace `var' = `temp'
	rename `var' y`var'
}
reshape long y, i(p) j(item) string
encode item, gen(i)

tempvar temp 
egen `temp' = mean(y) , by(i)
replace y = `temp' if y == . 

mvgstudy (y = p i#p)
mvdstudy, obj(p) errortype(relative) facetnum(8)

