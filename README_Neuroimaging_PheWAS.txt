To call the code from within R:
1. Put all .R code in the same folder
2. Edit the paths in Neuroimaging_PheWAS.R to point to where the support functions (PheWAS_models.R and genodata_handling.R) are located
3. Neuroimaging_PheWAS(genofile, imagfile, covfile, out_prefix, min_genetic, max_order, outlier_criteria)
	genofile = input genotype data (.csv format)
	imagfile = input imaging data (.csv format)
	covfile = input covariate data file (.csv format)
	out_prefix = output prefix will generate .csv files and .png with the manhattan plots
	min_genetic = minimum percentage of subjects that have to have each allele to be considered a group. for example in 100 subjects, if 60 have AA, 30 have AG, and 10 have GG & you set this parameter to 10, all the stats will be computed as t-tests between AA and AG because there are not greater than 10% of subjects with GG
	max_order = 0, 1 or 2 refers to how many statistical models (currently linear combinations of covariates) you want to have included. If this parameter is 0, then only the main effects of the genotype will be calculated. If 1, then main effects and first order interactions will be computed, similarly with 2
	outlier_criteria = number of IQR outside of quantile range 0.25-0.75 that is considered an outlier
