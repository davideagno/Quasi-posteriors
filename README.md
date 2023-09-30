# Quasi-posteriors

Code and data for the analysis contained in [[1]](#1).
- [example_intro.R](https://github.com/davideagno/HackTheGene/blob/main/BALSAMICO_analysis.R) and [example_intro.stan](https://github.com/davideagno/HackTheGene/blob/main/BALSAMICO_analysis.R) contain the code for replicate Figure 1.





Code and data from the HackTheGene competition. The details are available in [[2]](#2).
- [BALSAMICO_analysis.R](https://github.com/davideagno/HackTheGene/blob/main/BALSAMICO_analysis.R) contains the code for the analysis made with BALSAMICO model of [[1]](#1).
- [GLM_PCA_analysis.R](https://github.com/davideagno/HackTheGene/blob/main/GLM_PCA_analysis.R) contains the code for the analysis made with GLM-PCA model of [[3]](#3).
- [RCF_analysis_tuning.R](https://github.com/davideagno/HackTheGene/blob/main/RCF_analysis_tuning.R) contains the code for the tuning procedure made with RCF.
- [RCF_analysis_best.R](https://github.com/davideagno/HackTheGene/blob/main/RCF_analysis_best.R) contains the code for the imputation of missing values using RCF.
- [counts.csv](https://github.com/davideagno/HackTheGene/blob/main/counts.csv) contains the dataset of the competition - i.e., with the missing values.
- [counts_complete.csv](https://github.com/davideagno/HackTheGene/blob/main/counts_complete.csv) contains the full observed matrix.
- [covariates.csv](https://github.com/davideagno/HackTheGene/blob/main/covariates.csv) contains the covariates for the samples.

### References
<a id="1">[1]</a> 
Abe, K., Hirayama, M., Ohno, K., and Shimamura, T. (2021).
Hierarchical non-negative matrix factorization using clinical information for microbial communities.
BMC genomics, 22(1), 1–17.

<a id="2">[2]</a> 
Agnoletto, D., Collarin, C., Panarotto, A., Stolf, F. (2023+). 
Missing values imputation via Bayesian NMF for genomic data.
Submitted.

<a id="3">[3]</a> 
Townes, F. W., Hicks, S. C., Aryee, M. J., and Irizarry, R. A. (2019).
Feature selection and dimension reduction for single-cell rna-seq based on a multinomial model. 
Genome biology, 20, 1–16.
