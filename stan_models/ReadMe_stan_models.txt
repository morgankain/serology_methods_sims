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

(1) cluster_regression_base_1.stan
	-- simple/base two cluster model
(2) cluster_regression_with_beta_1.stan
	-- ^^ But with coefficients for group assignment
(3) cluster_regression_with_beta_theta_1.stan
	-- ^^ Also with coefficients to explain within-group MFI variance


-------
-- Conditional Stan models
	(all of these build upon a base above)
-------

3.1) cluster_regression_with_beta_theta_2.stan
	-- (2) above but with random variation by plate and an additional continuous covariate

1.1) cluster_regression_base_plate_random_1.stan
	-- (1) above but with a random effect by plate


-------
-- Empirical data focused models
-------

1.2) cluster_regression_base_1_test.stan
	-- (1) but playing with generated quantities because of NAs when the data is bad 




