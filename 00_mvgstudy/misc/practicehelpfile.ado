/* START HELP FILE
title[a command to give the parameters of the single stage Gehan design]

desc[
 {cmd:sampsi_gehan} calculates the sample sizes for the first and second stages of the Gehan design
    (1961).
]
opt[beta() specifies the first stage maximum probability of seeing no responses.]
opt[p1() specifies the desired probability of response.]
opt[se()  specifies the desired standard error in the second stage.]
opt[start() specifies the smallest n to start the search from.]
opt[precp() specifies the probability used in the standard error formula.]


opt2[precp() specifies the probability used in the standard error formula.
and I wanted more than one line for the longer  descriptions of the
option precp() later in the help file]

example[
 {stata sampsi_gehan, p1(0.2) beta(0.05) se(0.1) precp(0.4)}
]
author[Prof Adrian Mander]
institute[Cardiff University]
email[mandera@cardiff.ac.uk]

return[n1 The first stage sample size]
return[p1 The interesting p1 ]
return[beta The type 2 error]
return[se Standard error]
return[n2 The second stage sample size]

freetext[]

references[
Gehan, E.A. (1961) The Determination of the Number of Patients Required in a Preliminary and 
Follow-Up Trial of a New Chemotherapeutic Agent. Journal of Chronic Diseases, 13, 346-353.
]

seealso[
{help sampsi_fleming} (if installed)  {stata ssc install sampsi_fleming} (to install this command)

{help simon2stage} (if installed)   {stata ssc install simon2stage} (to install this command)

]

END HELP FILE */

