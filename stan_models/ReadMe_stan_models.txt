-----------------------------------------
-- Notes on Stan models in this folder --
-----------------------------------------

-------
-- General notes
-------

-- All models here are written to be constrained to two clusters
-- the "2" at the end of each name (prior to .stan) indicates that the model captures simulations of data complexity "2" (see _targets.R for details) 


-------
-- Stan models used in the publication
-------

publication_model_normal_2.stan
	-- normal normal model

publication_model_lognormal_2.stan
	-- normal lognormal model (normal for seronegative distribution lognormal for seropositive distribution)publication_model_skew_normal_2.stan
	-- normal skew normal model (normal for seronegative distribution skew normal for seropositive distribution)


-------
-- Extra Stan model
-------

cluster_regression_with_beta_theta_2.stan
	-- normal normal example model with random variation by plate and an additional continuous covariate


