{smcl}
{* *! version 1.0.3 02jun2026}{...}
{vieweralsosee "reghdfe" "help reghdfe"}{...}
{vieweralsosee "ppmlhdfe" "help ppmlhdfe"}{...}
{viewerjumpto "Syntax" "belieffit##syntax"}{...}
{viewerjumpto "Description" "belieffit##description"}{...}
{viewerjumpto "Options" "belieffit##options"}{...}
{viewerjumpto "Generated variables" "belieffit##generated"}{...}
{viewerjumpto "Stored results" "belieffit##results"}{...}
{viewerjumpto "Examples" "belieffit##examples"}{...}

{title:Title}

{pstd}
{cmd:belieffit} {hline 2} Fit a two-equation belief model with corrected
log means and log-variance predictions
{p_end}

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:belieffit} {it:ylog} [{it:xvars}] [{cmd:if} {it:exp}]
[{cmd:in} {it:range}] [{it:weight}] [{cmd:,} {it:options}]
{p_end}

{synoptset 28 tabbed}{...}
{synopthdr:options}
{synoptline}
{syntab:Model}
{synopt:{cmd:absorb(}{it:varlist}{cmd:)}}absorbed fixed effects{p_end}
{synopt:{cmd:tolerance(}{it:#}{cmd:)}}outer-loop convergence tolerance{p_end}
{synopt:{cmd:maxiter(}{it:#}{cmd:)}}maximum number of outer iterations{p_end}

{syntab:Bootstrap}
{synopt:{cmd:reps(}{it:#}{cmd:)}}multiplier bootstrap repetitions{p_end}
{synopt:{cmd:vce(}{it:varlist}{cmd:)}}cluster variables for bootstrap weights{p_end}
{synopt:{cmd:df(adj)}}use R-1 denominator for FE covariance variables{p_end}
{synopt:{cmd:df(none)}}use R denominator for FE covariance variables{p_end}

{syntab:Generated variables}
{synopt:{cmd:feprefix(}{it:name}{cmd:)}}prefix for saved fixed-effect variables{p_end}
{synopt:{cmd:varprefix(}{it:name}{cmd:)}}prefix for saved FE covariance variables{p_end}
{synopt:{cmd:prednames(}{it:newmean newvar}{cmd:)}}save corrected mean and variance predictions{p_end}
{synopt:{cmd:replace}}overwrite generated variables if they exist{p_end}

{syntab:Advanced}
{synopt:{cmd:aux(}{it:string}{cmd:)}}run additional e-class commands{p_end}
{synopt:{cmd:calc(}{it:string}{cmd:)}}run a custom calculation block{p_end}
{synoptline}

{pstd}
Weights may be {cmd:aweight}, {cmd:fweight}, {cmd:iweight}, or {cmd:pweight}.
{p_end}

{marker description}{...}
{title:Description}

{pstd}
{cmd:belieffit} fits a two-part belief model for a log outcome {it:ylog}.
The first equation is a linear high-dimensional fixed-effect regression for a
corrected posterior mean. The second equation is a PPML high-dimensional
fixed-effect regression for squared residuals, interpreted as the conditional
variance. The two equations are iterated until the predicted variance
stabilizes.
{p_end}

{pstd}
At each iteration, the corrected mean is formed as
{p_end}

{p 12 12 2}
{it:mhat} = {it:ylog} - 0.5*{it:tau2},
{p_end}

{pstd}
where {it:tau2} is updated using the fitted conditional variance from the
PPML residual-squared equation. The final coefficient vector combines the
mean-equation coefficients and the variance-equation coefficients in equations
{cmd:mean:} and {cmd:variance:}.
{p_end}

{pstd}
When {cmd:reps()} is positive, standard errors are obtained from an exponential
multiplier bootstrap. If {cmd:vce()} is supplied, the program draws product
exponential weights at the level of each cluster variable in {cmd:vce()}.
{p_end}

{pstd}
{cmd:belieffit} requires the user-written commands {cmd:reghdfe} and
{cmd:ppmlhdfe}.
{p_end}

{phang2}{cmd:. ssc install reghdfe}{p_end}
{phang2}{cmd:. ssc install ppmlhdfe}{p_end}

{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{cmd:absorb(}{it:varlist}{cmd:)} specifies fixed effects to absorb in both the
corrected-mean equation and the variance equation. Factor-variable notation and
interactions may be used in the form accepted by {cmd:reghdfe} and
{cmd:ppmlhdfe}. When fixed effects are supplied, final fixed-effect
contributions are saved in generated variables.
{p_end}

{phang}
{cmd:tolerance(}{it:#}{cmd:)} sets the convergence tolerance for the outer
fixed-point loop. The default is {cmd:tolerance(1e-6)}.
{p_end}

{phang}
{cmd:maxiter(}{it:#}{cmd:)} sets the maximum number of outer fixed-point
iterations. The default is {cmd:maxiter(100)}.
{p_end}

{dlgtab:Bootstrap}

{phang}
{cmd:reps(}{it:#}{cmd:)} specifies the number of bootstrap repetitions. The
default is {cmd:reps(200)}. Specify {cmd:reps(0)} or a negative number to post
only the baseline coefficient estimates and no bootstrap variance matrix.
{cmd:reps(1)} is not allowed.
{p_end}

{phang}
{cmd:vce(}{it:varlist}{cmd:)} specifies one or more cluster variables for the
multiplier bootstrap. With multiple variables, bootstrap weights are the product
of cluster-level exponential draws across the listed cluster dimensions. If
{cmd:vce()} is omitted, observation-level exponential weights are used.
{p_end}

{phang}
{cmd:df(adj)} or {cmd:df(none)} controls the denominator used for generated
fixed-effect covariance variables. The default is {cmd:df(adj)}. This option
affects the generated covariance variables only, not the coefficient covariance
matrix posted in {cmd:e(V)}.
{p_end}

{dlgtab:Generated variables}

{phang}
{cmd:feprefix(}{it:name}{cmd:)} sets the prefix for generated fixed-effect
variables. The default is {cmd:feprefix(fe)}. With {cmd:absorb()} specified, the
final mean-equation fixed-effect variables are saved first, followed by the
final variance-equation fixed-effect variables.
{p_end}

{phang}
{cmd:varprefix(}{it:name}{cmd:)} sets the prefix for generated fixed-effect
covariance variables from the bootstrap. The default is {cmd:varprefix(S)}.
The generated variables are lower-triangular covariance entries across the
saved fixed-effect columns.
{p_end}

{phang}
{cmd:prednames(}{it:newmean newvar}{cmd:)} requests two generated prediction
variables. The first variable contains the corrected mean. The second variable
contains the PPML-predicted variance. In the current implementation, use
{cmd:replace} with {cmd:prednames()} when creating or overwriting these
variables.
{p_end}

{phang}
{cmd:replace} permits {cmd:belieffit} to drop and overwrite generated
fixed-effect variables, fixed-effect covariance variables, and prediction
variables that would otherwise conflict with existing variable names.
{p_end}

{dlgtab:Advanced}

{phang}
{cmd:aux(}{it:string}{cmd:)} runs one or more additional e-class commands after
the final fit and appends their {cmd:e(b)} vectors to the posted coefficient
vector. Separate multiple auxiliary blocks with semicolons. A block may start
with {cmd:name(}{it:eqname}{cmd:)} to set the equation name in the combined
coefficient vector.
{p_end}

{pmore}
Inside {cmd:aux()}, the placeholder tokens {c -(}W{c )-}, {c -(}mean{c )-},
and {c -(}var{c )-} are replaced by the current weight expression, the corrected
mean variable, and the predicted variance variable, respectively.
{p_end}

{phang}
{cmd:calc(}{it:string}{cmd:)} runs a custom calculation block after the final
fit. As in {cmd:aux()}, {c -(}W{c )-}, {c -(}mean{c )-}, and {c -(}var{c )-}
are substituted before evaluation.
{p_end}

{marker generated}{...}
{title:Generated variables}

{pstd}
If {cmd:absorb()} is specified, {cmd:belieffit} saves fixed-effect contribution
variables using {cmd:feprefix()}. With the default prefix, these are named
{cmd:fe1}, {cmd:fe2}, and so on. The first block corresponds to the final mean
equation and the second block corresponds to the final variance equation.
{p_end}

{pstd}
If bootstrap repetitions are run and fixed-effect columns are present,
{cmd:belieffit} also creates lower-triangular fixed-effect covariance variables.
With the default prefix, these are named {cmd:S11}, {cmd:S21}, {cmd:S22},
{cmd:S31}, and so on, where {cmd:Sjk} is the bootstrap covariance between
fixed-effect columns {it:j} and {it:k} for each observation.
{p_end}

{pstd}
If {cmd:prednames()} is specified, the requested variables contain the corrected
mean and predicted variance.
{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:belieffit} stores the following in {cmd:e()}:
{p_end}

{phang}{cmd:e(reps)}: number of requested bootstrap repetitions; 0 when skipped.{p_end}
{phang}{cmd:e(cmd)}: {cmd:belieffit}.{p_end}
{phang}{cmd:e(depvar)}: dependent log outcome variable.{p_end}
{phang}{cmd:e(regressors)}: right-hand-side variables.{p_end}
{phang}{cmd:e(absorb)}: absorbed fixed effects, if specified.{p_end}
{phang}{cmd:e(feprefix)}: fixed-effect variable prefix.{p_end}
{phang}{cmd:e(varprefix)}: fixed-effect covariance variable prefix.{p_end}
{phang}{cmd:e(vcetype)}: {cmd:Bootstrap} when bootstrap variance is posted;
{cmd:None} when {cmd:reps(0)} is used.{p_end}
{phang}{cmd:e(b)}: combined coefficient row vector. Equations include
{cmd:mean:}, {cmd:variance:}, and any auxiliary equations.{p_end}
{phang}{cmd:e(V)}: bootstrap covariance matrix for {cmd:e(b)}, when
{cmd:reps()} is positive.{p_end}

{marker examples}{...}
{title:Examples}

{phang2}
{cmd:. belieffit logbelief, absorb(person job) reps(500) vce(person job)}
{p_end}

{phang2}
{cmd:. belieffit logbelief x1 x2, absorb(person job) prednames(mhat vhat) }
{cmd:replace reps(0)}
{p_end}

{pstd}
Also run a regression of posterior mean beliefs and posterior variance on gender and include it in
the bootstrap:
{p_end}

{phang2}
{cmd:. belieffit logbelief, absorb(person job) reps(500) vce(person job) }
{cmd:aux(name(mean) regress }{c -(}mean{c )-}{cmd: gender }{c -(}W{c )-}{cmd: ; name(variance) regress }{c -(}var{c )-}{cmd: gender }{c -(}W{c )-}{cmd:)}
{p_end}

{title:Author}

{pstd}
Martin Eckhoff Andresen, Department of Economics, University of Oslo.
Developed for {browse "https://arxiv.org/pdf/2606.02503":Pay Beliefs and the Amenity-Pay Tradeoff}, joint with Manudeep
Bhuller and Alfred Løvgren. 
{p_end}
