-----------------------------------------
-- Notes on Stan models in this folder --
-----------------------------------------

-------
-- General notes
-------

-- All models here are written to be constrained to two clusters


-------
-- Primary Stan models
-------

(1) cluster_regression_base.stan
	-- simple/base two cluster model
(2) cluster_regression_with_beta.stan
	-- ^^ But with coefficients for group assignment
(3) cluster_regression_with_beta_theta.stan
	-- ^^ Also with coefficients to explain within-group MFI variance


-------
-- Conditional Stan models
	(all of these build upon a base above)
-------

2.1) cluster_regression_with_beta_plate_random.stan
	-- (2) above but with random variation by plate
		TBD



