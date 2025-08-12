{smcl}
{* *! version 1.0.0  4aug2025}{...}
{viewerdialog mvgstudy "dialog _mvgstudy"}{...}
{viewerjumpto "Syntax" "mvgstudy##syntax"}{...}
{viewerjumpto "Description" "mvgstudy##description"}{...}
{viewerjumpto "Options" "mvgstudy##postestimation"}{...}
{viewerjumpto "Examples" "mvgstudy##examples"}{...}
{viewerjumpto "Stored results" "mvgstudy##results"}{...}
{viewerjumpto "Reference" "mvgstudy##reference"}{...}
{vieweralsosee "[MV] manova" "help manova"}{...}

{p2col:{bf:mvgstudy}} Variance and covariance component estimation for multifaceted designs  {p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 18 2}
{cmd:mvgstudy} ({depvarlist} {cmd:=} {it:termlist}) {ifin} 

{p 8 10 2}
where {it:termlist} is a factor variable list of nested and crossed facets with the following features:

{p 8 10 2}
o Variables are assumed to be categorical

{p 8 10 2}
o The {cmd:#} symbol is used to indicate crossing

{p 8 10 2}
o The {cmd:|} symbol (indicating nesting) may be used in place of the
    {cmd:#} symbol and indicates that the term to the left of the {cmd:|} is nested within the term to the right.
{p_end}

{marker description}{...}
{title:Description}


{pstd}{opt mvgstudy} is a wrapper program for the {helpb manova} command. It calculates variance and covariance components using the SSCP matrices generated from a manova decomposition, following the method outlined in Brennan, 2001. The command can calculate variance and covariance components for any nested and crossed multifaceted designs. When {depvarlist} contains multiple variables, the command generates covariance component matrices corresponding to each effect in {it:termlist}. When {depvarlist} contains one variable, the command generates variance components corresponding to each effect in {it:termlist}. The command assumes that the design is balanced. {p_end}


{marker results}{...}
{title: Stored Results}

{pstd}
{cmd:mvgtudy} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:r(df)}} degrees of freedom corresponding to the facets used in {it:termlist} {p_end}
{synopt:{cmd:r(P)}}  upper triangular matrix of EMCP equation coefficients used to calculate emcp matrices (covariance components) using the SSCP matrices {p_end}
{synopt:{cmd:r(emcpi)}} covariance component matrix corresponding to the i_th score effect in termlist


{marker postestimation}{...}
{title:Post Estimation Commands}

{pstd}The variance and covariance components estimated using the {bf:mvgstudy} command can be used in many subsequent analyses. For example, they can inform a decision study (D-study), where a researcher projects the reliability of scores averaged over conditions that vary in facet-level sample sizes. One might ask, for instance, how much more reliable a student's reading achievement score becomes as the number of testing occasions increases. The variance components can also be used to calculate standard errors of measurement, conditional standard errors of measurement, and, in the multivariate case, various indices that characterize the variation of scores within a composite or score profile. {p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Commands}{p_end}
{synopt:{helpb mvdstudy}} Decision study command for univariate and multivariate - multifaceted designs{p_end}


{marker examples}{...}
{title:Examples}

{hline}

{pstd}
Example of a multivariate person-by-task design: each person was assessed on five dimensions of a rubric across three repeated tasks (observations).
{p_end}

{pstd}Setup{p_end}


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


{pstd} Save covariance component matrices {p_end}

{p 8 10 2}
{stata mata p = st_matrix("r(emcp1)")}
{p_end}
{p 8 10 2}
{stata mata t = st_matrix("r(emcp2)")}
{p_end}
{p 8 10 2}
{stata mata pt = st_matrix("r(emcp3)")}
{p_end}

{pstd} Useful statistics using covariance component matrices {p_end}

{pstd} Proportion of true score (co)variance  {p_end}
{p 8 10 2}
{stata "mata p:/(p+t+pt)"}
{p_end}

{pstd} Composite Variances {p_end}
{p 8 10 2}
{stata mata weights = J(5,1,1/5)}
{p_end}
{p 8 10 2}
{stata mata weights'*p*weights}
{p_end}
{p 8 10 2}
{stata mata weights'*t*weights}
{p_end}
{p 8 10 2}
{stata mata weights'*pt*weights}
{p_end}

{hline}


{hline}

{pstd}
Example of a univariate person-by-task design: each person was assessed by two randomly assigned raters using a single average score across three repeated tasks (observations).
{p_end}

{pstd}Setup{p_end}

{p 8 10 2}
{stata use mvgstudyexampledata.dta, clear}
{p_end}

{p 8 10 2}
{stata egen xavg = rowmean(x1 x2 x3 x4 x5)}
{p_end}

{pstd} Run mvgstudy {p_end}

{p 8 10 2}
{stata mvgstudy (xavg = p t p#t)}
{p_end}

{marker reference}{...}
{title:Reference}

{marker HTF}{...}
{phang} 
Brennan, R. L. (2001). Generalizability theory. Springer. https://doi.org/10.1007/978-1-4757-3456-0
{p_end}


