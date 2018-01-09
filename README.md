R routines to accompany the manuscript:

H. Michael, L. Tian, and M. Ghebremichael. The ROC Curve for Regularly Measured Longitudinal Biomarkers. (Under revision at Biostatistics.)
Corresponding author: haben.michael@stanford.edu

Contents:

* `ragon.R` has R routines implementing the ROC approximation algorithms discussed in the paper, including:
	1. `roc.1`, `roc.2`, and `roc.3` -- implementing the subject-specific, typical subject, and population ROC curves
	2. `subj.fit` and `pop.fit` -- empirical bayes fits of subject and population parameters
	3. classes `Patient` and `PatientROC` encapsulating patient data and a patient ROC; plot and other routines have been overloaded for these classes
	
* `revision.R` has the master code for the simulations and figures presented in the paper. The simulations in turn are located in the enumerated subdirectories of `simulations/`, as further described in `revision.R`. Also included are *nix and SunGrid scripts to run the simulations in parallel, assuming available installations.

Since the data from the (ongoing) Yale Pediatric HIV Cohort is not publicly available, simulations based on comparable data have been provided. These are described in the Simulations section of the manuscript, and include simulations `g` and `h`. Since the HIV data set is not available, some code, such as simulations `i` and `j`, will not run as provided.
