{smcl}
{* *! version 1.0.0  4aug2025}{...}
{viewerdialog mvgstudy "dialog _mvgstudy"}{...}
{viewerjumpto "Syntax" "mvgstudy##syntax"}{...}
{viewerjumpto "Description" "mvgstudy##description"}{...}
{viewerjumpto "Options" "mvgstudy##options"}{...}
{viewerjumpto "Examples" "mvgstudy##examples"}{...}
{viewerjumpto "Stored results" "mvgstudy##results"}{...}
{viewerjumpto "Reference" "mvgstudy##reference"}{...}
{vieweralsosee "[MV] manova" "help manova"}{...}

{p2col:{bf:mvgstudy}} Decision study for univariate and multivariate multifaceted designs  {p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 18 2}
{cmd:mvdstudy} {cmd:,} {opt o:bject}{cmd:(}{it:string}{cmd:)} {opt e:rrortype}{cmd:(}{it:string}{cmd:)} {opt facet:num}{cmd:(}{it:integer}{cmd:)} [{opt comp:ositeweights}{cmd:(}{it:string}{cmd:)}]


{marker description}{...}
{title:Description}


{pstd}{opt mvgstudy} is a post estimation command for the {helpb mvgstudy} command. {p_end}


{marker options}{...}
{title:Options}

{phang}
{opt o:bject} is required and specifies the object of measurement for the decision study.

{phang}
{opt e:rrortype} is required and must be specified as either relative or absolute. See Brennan (2001) for a discussion of error types. Relative errors are used to calculate generalizability coefficients, while absolute errors are used to calculate phi coefficients.

{phang}
{opt facet:num} is required and specifies the number of levels per facet over which you want to calculate reliability coefficients.

{phang}
{opt comp:ositeweights} is optional and specifies a column vector of weights used to calculate composite scores. This option only applies in multivariate designs. If specified, {cmd:mvdstudy} calculates reliability coefficients for the composite scores. If the design is multivariate and the option is not specified, {cmd:mvdstudy} calculates reliability coefficients by item.


{marker results}{...}
{title: Stored Results}

{pstd}
{cmd:mvgtudy} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:r(namelist)}} D-study projection matrix. The last column contains the reliability coefficient; the second-to-last column contains the true score variance for a single score; the third-to-last column contains the error variance (either relative or absolute); and the remaining columns indicate the number of levels per facet.}



{marker examples}{...}
{title:Examples}

{hline}

{pstd}
Composite Score Analysis.
{p_end}

{p 8 10 2}
{stata use mvgstudyexampledata.dta, clear}
{p_end}

{pstd} Data Structure{p_end}

{p 8 10 2}
{stata table p t if p <= 20}
{p_end}

{pstd} Run mvgstudy {p_end}

{p 8 10 2}
{stata mvgstudy (x1 x2 x3 x4 x5 = p t p#t)}
{p_end}

{pstd} Run mvdstudy {p_end}

{p 8 10 2}
{stata matrix weights = J(5,1,1/5)}
{p_end}
{p 8 10 2}
{stata mvdstudy , o(p) e(relative) facet(10) comp(weights)}
{p_end}
{p 8 10 2}
{stata matrix rel = r(composite)}
{p_end}

{pstd} Plot Results {p_end}

{p 8 10 2}
{stata clear}
{p_end}
{p 8 10 2}
{stata svmat rel}
{p_end}
{p 8 10 2}
{stata rename (*) (t error true rel)}   
{p_end}
{p 8 10 2}
{stata scatter rel t, connect(l) ytitle(G coefficient) xtitle(tasks) title(Composite)}
{p_end}

{hline}

{pstd}Sub Score Analysis{p_end}

{pstd}Setup{p_end}
{p 8 10 2}
{stata use mvgstudyexampledata.dta, clear}
{p_end}

{pstd}Data Structure{p_end}

{p 8 10 2}
{stata table p t if p <= 20}
{p_end}

{pstd} Run mvgstudy {p_end}

{p 8 10 2}
{stata mvgstudy (x1 x2 x3 x4 x5 = p t p#t)}
{p_end}

{pstd} Run mvdstudy {p_end}

{p 8 10 2}
{stata mvdstudy, o(p) e(relative) facet(10)}
{p_end}
{p 8 10 2}
{stata matrix rel = (J(10,1,1), r(x1)) \ (J(10,1,2), r(x2)) \ (J(10,1,3),r(x3)) \ (J(10,1,4), r(x4)) \ (J(10,1,5), r(x5))}
{p_end}

{pstd} Plot Results {p_end}
{p 8 10 2}
{stata clear}
{p_end}
{p 8 10 2}
{stata svmat rel}
{p_end}
{p 8 10 2}
{stata rename (*) (i t error true rel)}   
{p_end}
{p 8 10 2}
{stata scatter rel t if i == 1, connect(l) ytitle(G coefficient) xtitle(tasks) title(x1) name(g1,replace) }
{p_end}
{p 8 10 2}
{stata scatter rel t if i == 2, connect(l) ytitle(G coefficient) xtitle(tasks) title(x2) name(g2,replace)}
{p_end}
{p 8 10 2}
{stata scatter rel t if i == 3, connect(l) ytitle(G coefficient) xtitle(tasks) title(x3) name(g3,replace)}
{p_end}
{p 8 10 2}
{stata scatter rel t if i == 4, connect(l) ytitle(G coefficient) xtitle(tasks) title(x4) name(g4,replace)}
{p_end}
{p 8 10 2}
{stata scatter rel t if i == 5, connect(l) ytitle(G coefficient) xtitle(tasks) title(x5) name(g5,replace)}
{p_end}
{p 8 10 2}
{stata graph combine g1 g2 g3 g4 g5 , altshrink ycommon}
{p_end}


{marker reference}{...}
{title:Reference}

{marker HTF}{...}
{phang} 
Brennan, R. L. (2001). Generalizability theory. Springer. https://doi.org/10.1007/978-1-4757-3456-0
{p_end}


