/*
clear all 
include "/Users/steffenerickson/Documents/GitHub/stata_programs/00_mvgstudy/mvgstudy_rewrite.do"
use "/Users/steffenerickson/Box Sync/NSF_DR_K12/measurement/${data}/manova_data.dta", clear 
rename (task rater coaching) (t r treat)
recode time (2 =1)
drop if (r == 3) | (treat == 2) | (x1 == . & x2 == . & x3 == . & x4 == . & x5 == .)
encode participantid , gen(id)
egen p = group (id site semester treat)
drop x6
keep if time == 1 
keep p t r x*
order p t r x*
sort p t r 
egen counts = count(x1) , by(p)
keep if counts == 6
drop counts 
sort p t r 
egen p2 = group(p)
replace p = p2 
drop p2
collapse x1 x2 x3 x4 x5, by(p t)
label variable p "pst"
label variable t "task"
label variable x1 "Objective"
label variable x2 "Unpacking"
label variable x3 "Self-Instruction"
label variable x4 "Self-Regulation"
label variable x5 "Ending"
save /Users/steffenerickson/Documents/GitHub/stata_programs/00_mvgstudy/mvgstudyexampledata, replace 
*/


clear all 
include "/Users/steffenerickson/Documents/GitHub/stata_programs/00_mvgstudy/mvgstudy_rewrite.do"
use "/Users/steffenerickson/Box Sync/NSF_DR_K12/measurement/${data}/manova_data.dta", clear 
rename (task rater coaching) (t r treat)
recode time (2 =1)
drop if (r == 3) | (treat == 2) | (x1 == . & x2 == . & x3 == . & x4 == . & x5 == .)
encode participantid , gen(id)
egen p = group (id site semester treat)
drop x6
keep if time == 1 
collapse x1 x2 x3 x4 x5, by(p t)
keep p t  x*
order p t  x*
sort p t 
egen counts = count(x1) , by(p)
keep if counts == 4
drop counts 
egen p2 = group(p)
replace p = p2 
drop p2




label variable p "pst"
label variable t "task"
label variable x1 "Objective"
label variable x2 "Unpacking"
label variable x3 "Self-Instruction"
label variable x4 "Self-Regulation"
label variable x5 "Ending"
save /Users/steffenerickson/Documents/GitHub/stata_programs/00_mvgstudy/mvgstudyexampledata, replace 




use "/Users/steffenerickson/Documents/GitHub/stata_programs/00_mvgstudy/mvgstudyexampledata.dta", clear 

mvgstudy (x* = p t r|t p#t r|p#t) 
mat w1 = J(5,1,1/5)
pcamat r(emcp1) , n(64)
mat w2 = e(L)[1...,1]
mvdstudy, obj(p) errortype(relative) facet(30) comp(w1)
mat rel = r(composite)

preserve
	clear 
	svmat rel
	rename (*) (r t error true rel)
	
	twoway (scatter rel t if r == 1 , connect(l)) ///
	       (scatter rel t if r == 2 , connect(l)) ///
		   (scatter rel t if r == 3 , connect(l)) ///
		   (scatter rel t if r == 4 , connect(l)) ///
		   (scatter rel t if r == 5 , connect(l)) ///
		   (scatter rel t if r == 6 , connect(l)) ///
		   (scatter rel t if r == 7 , connect(l)) ///
		   (scatter rel t if r == 8 , connect(l)) , legend(off)  
restore





collapse x* , by(p)
graph box x*


reshape long x,i(p t r) j(i)

mvgstudy (x = p t i r|t p#t p#i i#t r|p#t r|p#i r|p#t#i) 
mvdstudy p#i, errortype(relative) facet(8) 
mat rel = r(x)
preserve

	clear 
	svmat rel
	rename (*) (r t i error true rel)
	scatter (rel i)
	
	//graph3d rel r t  
		   
restore

graph3d x y z



















collapse x*, by(time t p)
//set varabbrev off
mvgstudy (x* = p t p#t) if time == 1
mvdstudy p,  e(absolute) facet(1)
mvdstudy p,   errortype(absolute) facet(1)


mvgstudy (x* = p t p#t) if time == 1
mvdstudy p,  e(absolute) facet(1)
mvdstudy p,   errortype(absolute) facet(1)










collapse x*, by(time t p)
drop x6
//set varabbrev off
mvgstudy (x* = p t p#t) if time == 1
mvdstudy p,  e(absolute) facet(1)
set trace on 
mvdstudy p,   errortype(absolute) facet(3)
