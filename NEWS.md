# CompetingRiskRPackage 0.3.7

* Use the validated n=4000 examples in the GitHub tutorial and simplify the
  repository README to purpose and downloads. Statistical code is unchanged.

* Reduce covariance memory and repeated risk-set calculations using equivalent
  finite-sum identities; release likelihood caches on return and interruption.
* Make long native calculations responsive to R interrupts.
* Keep invalid residual correlations and their inference unavailable, preserve
  unaffected columns, and warn when selected residual plots cannot be drawn.
* Stabilize covariance, residual and regression calculations under large
  covariate offsets and changes of measurement units.
* Verify the single nlminb result using event-normalized scores and a relative
  parameter correction in standardized coordinates. Preserve the original
  optimizer status and numerical checks in solver_diagnostics.

# CompetingRiskRPackage 0.3.6

* Revise vignette wording to focus on package use, inputs and outputs.
  All example code, statistical implementation and plotting functions are unchanged.

# CompetingRiskRPackage 0.3.5

* Expand formula dots over predictors, excluding named cohort design columns.
* Allow residual plot selection by formula term or coefficient name, preserving
  requested order and expanding factors. Use smaller translucent points.
* Use native plot calls throughout the vignette; retain statistical kernels.

# CompetingRiskRPackage 0.3.4

* Use a multivariable vignette throughout, demonstrate safe selection of all
  covariates and show the complete code for the two-panel residual figure.
* Statistical code and public API are unchanged.

# CompetingRiskRPackage 0.3.3

* Replace the vignette's single residual illustration with paired PH and non-PH
  synthetic examples (n=500 each), using continuous covariates and identical
  plotting settings. Keep the original n=200 worked example and statistical API.

# CompetingRiskRPackage 0.3.2

* Add a concise English PDF vignette with an executed n=200 example, all three
  PH diagnostics, a residual plot and common input and extraction operations.
* Add knitr/rmarkdown as optional vignette-building dependencies. Statistical
  code and the public API are unchanged.

# CompetingRiskRPackage 0.3.1

* Retain beta estimates and descriptive residuals when covariance computation
  fails; expose covariance status and withhold inference instead of aborting.
* Validate complete beta covariance matrices in unit-free coordinates.
* Verify PH regression against the manuscript's estimating equations, numerical
  derivatives, constant/affine Q identities and the full Taylor covariance.
* Document strict variance versus inclusive fitting risk sets; retain teacher
  equations and verify zero/undefined-alpha boundaries independently.

# CompetingRiskRPackage 0.3.0

* Separate input preparation, fitting, PH tests, S3 methods and plotting.
* Let R own the single nlminb invocation; keep cached likelihood evaluation and
  numerical convergence checks in C++.
* Separate beta covariance, residual diagnostics and PH regression kernels;
  preserve their existing risk-set boundaries and equations.
* Reuse one censoring KM distribution throughout each fit.
* Remove numbered development filenames and repeated native matrix conversion.
* Keep fixed-beta testing and legacy field aliases in test helpers only.
* Return large regression influence details only when requested internally.
* Clarify user documentation and retain archived independent reference fixtures.

# CompetingRiskRPackage 0.2.0

* Make cc_fit() the single public fitting entry point, using a covariate formula.
* Add standard coefficient, covariance, residual and summary methods.
* Return full sandwich and information covariance matrices and a unified table
  of cor, rank cor, and joint/component regression proportionality diagnostics.
* Retain one stats::nlminb fit with analytic derivatives; no fallback restart.
* Remove development C++ examples, debug utilities and unused legacy optimizers.
* Preserve teacher-R fixtures and independent regression-score verification.
* Keep fixed-beta hypothesis tests internal, separate from PH diagnostics.
