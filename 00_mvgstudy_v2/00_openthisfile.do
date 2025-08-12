*! MVGSTUDY PACKAGE 
*! Steffen Erickson

clear all 

cd /Users/steffenerickson/Documents/GitHub/stata_programs/00_mvgstudy_v2
include mvgstudy.ado 


help mvgstudy
help mvdstudy 



use mvgstudyexampledata.dta, clear



info = panelsetup(V1, 1)
 
 
 
