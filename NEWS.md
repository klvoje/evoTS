# evoTS 1.0.4

# Bug fixes
- Fixed a small bug in sim.multi.OU and sim.multi.URW functions. The latter function now correctly handles unequal time steps by scaling each random increment by its own time interval.
- Fixed a factual error in the @details documentation of opt.decel.single.R.zero.corr, which incorrectly described the model as accelerating rather than decelerating.

# Other changes
- Improved performance of multivariate likelihood functions by pre-computing time-distance matrices and sampling error vectors outside the optimization loop.
- Diagonal elements of the A matrix are now constrained to be positive (> 0) for the "diag", "upper.tri", "lower.tri", ensuring positive definiteness for these model structures. No constraint is applied for the "full" parameterization.
- Improved starting values for multivariate OU model fitting: diagonal A and R elements are initialized from univariate fits and constrained to be positive, off-diagonal R elements are initialized to 0, and trait optima are initialized at the column means of the data.
- When fitting multivariate models with multiple iterations, perturbed starting values are now automatically clamped to respect parameter bounds, reducing the frequency of failed likelihood evaluations.
- Simplified the make.multivar.evoTS function.

