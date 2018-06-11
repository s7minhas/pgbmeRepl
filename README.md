This archive includes the files to replicate all analyses in “Modeling Asymmetric Relationships from Symmetric Networks,” Political Analysis.

Replication Instructions
===

The replication archive is organized into four directories:

- main: contains the data files and scripts necessary to reproduce the main results in the paper
- appendix: contains the data files and scripts necessary to reproduce the results in the appendix
- vignette: provides a brief introduction to the P-GBME model we introduce in this paper
- PA_submission: includes .tex and .pdf versions of the paper

The main directory contains four scripts that should be run in the following order: 

- 0_imputeData.R: Imputes missing data among the covariates using `sbgcop`.
    + This should take only a few minutes to run and is run in parallel across two cores.
- 1_createFig1Fig2.R: Run `pgbme` to generate Figures 1 and 2 in the paper
    + Original results from the authors are already included in the main directory, if those files (`results1995.rda` and `results2010.rda`) are deleted, then this script should take approximately one hour and ten minutes to run. 
- 2_inPerfAnalysis.R: Assess performance for in-sample predictions for P-GBME, GBME, and a pooled probit model. Results are saved to `table1_inSamplePerf.txt` and correspond to the two left-most columns in Table 1 of the paper.
    + Original results from the authors are already included in the main directory, if those files (`resultsPGBME_inPerf.rda` and the files in `gbme_inPerf`) are deleted, then this script should take approximately one hour and twenty minutes to run. Also do not delete the gbme_inPerf directory itself, instead just delete the files included therein.
- 3_outPerfAnalysis.R: Assess performance for out-of-sample predictions for P-GBME, GBME, and a pooled probit model. Results are saved to `table1_outSamplePerf.txt` and correspond to the two right-most columns in Table 1 of the paper.
    + Original results from the authors are already included in the main directory, if those files ('resultsPGBME_outPerf_fold[1:5].rda' and the files in `gbme_outPerf/fold[1:5]/`) are deleted, then this script should take approximately one hour and thirty minutes to run. Also do not delete the gbme_outPerf directory or any of the directories labeled fold1, fold2, etc., instead just delete teh files included therein.

The appendix directory contains three scripts that can be run in any order:

- 1_pgbmeSig_figA1.R: Conducts simulation exercise summarized in Figure 1 of the appendix. This script was run on a m5.24xlarge EC2 instance and took approximately 40 minutes to complete using 90 cores in parallel. Results are saved to figureA1.pdf.
- 2_pgbmeCoef_section2.3.R: Estimates PGBME model for every year from 1990 to 2012. This script was run on a m4.10xlarge EC2 instance and took approximately 40 minutes to complete using 23 cores in parallel. Results are saved to figureA2.png and figureA3.pdf
- 3_pgbmeVarK_section2.4.R: Estimates PGBME models varying dimension of multiplicative effect. This script takes approximately one hour to complete and utilizes three cores in parallel. Results are saved to tableA3.txt.

R and package build notes
===

Last, please note the version of each of the libraries that our project relies on (each library was built using R 3.5.0). 

|                     |                   |                          |                      |
|:--------------------|:------------------|:-------------------------|:---------------------|
|abind: (1.4-5)       |foreach: (1.4.4)   |magic: (1.5-8)            |PRROC: (1.3)          |
|caret: (6.0-76)      |ggplot2: (2.2.1)   |magrittr: (1.5)           |RColorBrewer: (1.1-2) |
|coda: (0.19-1)       |gridExtra: (2.2.1) |mnormt: (1.5-5)           |reshape2: (1.4.3)     |
|countrycode: (0.16)  |igraph: (1.2.1)    |msm: (1.6.6)              |sbgcop: (0.975)       |
|data.table: (1.10.4) |latex2exp: (0.4.0) |OptimalCutpoints: (1.1-3) |tidyr: (0.8.0)        |
|doParallel: (1.0.11) |lme4: (1.1-17)     |parallel: (3.5.0)         |verification: (1.42)  |
|dplyr: (0.7.4)       |                   |                          |                      |


Using the model
===

When using this model for your own work, the script that includes the functions you need is in the `main/` directory and is labeled `pgbme.R`. 

If you find any errors or have any further questions, please address them to me via email at minhassh@msu.edu.