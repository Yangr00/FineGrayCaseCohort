# Example data

example.csv is synthetic demonstration data, not patient data or a performance
table. It contains 200 full-cohort subjects and four columns:

* time: observed follow-up time in simulation units.
* status: 0=censored, 1=target event, 2=competing event.
* subcohort: 1=random subcohort member, 0=other cohort member.
* covariate1: simulated numeric covariate.

The file was exported from the archived canonical fixture-generation workflow
in make_fixtures.R. That workflow and the teacher's original scripts are retained
outside the installed package. The supplied CSV has complete measurements; the
vignette shows how to mask measurements not required by the design.

example_ph.csv and example_nonph.csv are synthetic vignette illustrations with
500 subjects each and the same four columns. covariate1 is Uniform(0,1).
The target subdistribution hazard is 4 exp(0.5 z + (-3.5 + theta z)t), with
theta=0 for PH and theta=-12 for non-PH. Both use seed 2026, independent
Uniform(0.05,1.5) censoring and Bernoulli(0.6) random subcohort membership.
Competing-event times conditional on cause 2 have rate exp(-0.5 z).
These are individual illustrations, not estimates of test size or power.
Generation and inverse-CIF validation code are retained outside the package.

example_ph_4000.csv and example_nonph_4000.csv use the same generating settings,
seed and columns as the 500-subject illustrations, with n=4000. They are the
frozen data for the GitHub tutorial and the validated n=4000 plot comparison.
The tutorial adds the same independent two-level factor using seed 2026.
