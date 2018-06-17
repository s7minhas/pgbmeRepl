This archive includes the files to replicate all analyses in “Modeling Asymmetric Relationships from Symmetric Networks,” Political Analysis.

## Replication Instructions

The replication archive is organized into four directories (note that these files are also available on Github at [https://github/s7minhas/pgbmeRepl](https://github.com/s7minhas/pgbmeRepl)):

- main: contains the data files and scripts necessary to reproduce the main results in the paper
- appendix: contains the data files and scripts necessary to reproduce the results in the appendix
- pgbme: includes the pgbme package which can be installed using the `devtools` package 
- vignette: provides a brief introduction to the P-GBME model we introduce in this paper
- PA_submission: includes .tex and .pdf versions of the paper

#### Setup Information

All of the analysis reported in the manuscript and the appendix is run on a m4.10xlarge EC2 instance. Instructions to set up the ec2 instance are shown below: 

```
sudo adduser minhas # substitute your own username
sudo add-apt-repository 'deb [arch=amd64,i386] https://cran.rstudio.com/bin/linux/ubuntu xenial/'
sudo apt-get update
sudo apt-get install r-base
sudo apt-get install gdebi-core
sudo apt-get install libapparmor1
wget https://download2.rstudio.org/rstudio-server-1.1.442-amd64.deb
sudo gdebi rstudio-server-1.1.442-amd64.deb
sudo rstudio-server verify-installation
sudo apt-get install pkg-config
sudo apt-get install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev libnlopt-dev
sudo apt-get install htop
sudo apt-get install tree
```

Session information for the version of R we ran on this instance is shown below (further information on packages used in the analysis is included at the end of the README): 

```
> sessionInfo()
R version 3.4.4 (2018-03-15)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 16.04.3 LTS

Matrix products: default
BLAS: /usr/lib/libblas/libblas.so.3.6.0
LAPACK: /usr/lib/lapack/liblapack.so.3.6.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

loaded via a namespace (and not attached):
[1] compiler_3.4.4
>
```

#### Running the scripts

The main directory contains four scripts that should be run in the following order (before running any of the scripts please modify the path to match your local environment, the path is specified in line 2 on each of the scripts): 

- 0_imputeData.R: Imputes missing data among the covariates using `sbgcop`.
    + This should take only a few minutes to run and is run in parallel across two cores.
- 1_createFig1Fig2.R: Run `pgbme` to generate Figures 1 and 2 in the paper
    + Original results from the authors are already included in the main directory, if those files (`results1995.rda` and `results2010.rda`) are deleted, then this script should take approximately one hour and ten minutes to run. 
- 2_inPerfAnalysis.R: Assess performance for in-sample predictions for P-GBME, GBME, and a pooled probit model. Results are saved to `table1_inSamplePerf.txt` and correspond to the two left-most columns in Table 1 of the paper.
    + Original results from the authors are already included in the main directory, if those files (`resultsPGBME_inPerf.rda` and the files in `gbme_inPerf`) are deleted, then this script should take approximately one hour and twenty minutes to run. Also do not delete the gbme_inPerf directory itself, instead just delete the files included therein.
- 3_outPerfAnalysis.R: Assess performance for out-of-sample predictions for P-GBME, GBME, and a pooled probit model. Results are saved to `table1_outSamplePerf.txt` and correspond to the two right-most columns in Table 1 of the paper.
    + Original results from the authors are already included in the main directory, if those files ('resultsPGBME_outPerf_folds.rda' and the files in `gbme_outPerf/fold[1:5]/`) are deleted, then this script should take approximately one hour and thirty minutes to run. Also do not delete the gbme_outPerf directory or any of the directories labeled fold1, fold2, etc., instead just delete teh files included therein.

The appendix directory contains three scripts that can be run in any order:

- 1_pgbmeSig_figA1.R: Conducts simulation exercise summarized in Figure 1 of the appendix. This script takes approximately 2 hours to complete using 35 cores in parallel. Results are saved to figureA1.pdf.
- 2_pgbmeCoef_section2.3.R: Estimates PGBME model for every year from 1990 to 2012. This script took approximately 40 minutes to complete using 23 cores in parallel. Results are saved to figureA2.png and figureA3.pdf
- 3_pgbmeVarK_section2.4.R: Estimates PGBME models varying dimension of multiplicative effect. This script takes approximately one hour to complete and utilizes three cores in parallel. Results are saved to tableA3.txt.

#### R package build notes

Last, please note the version of each of the libraries that our project relies on (each library was built using R 3.4.4). 

|                |                 |                   |                 |
|:---------------|:----------------|:------------------|:----------------|
|abind: 1.4-5    |devtools: 1.13.5 |doParallel: 1.0.11 |dplyr: 0.7.5     |
|foreach: 1.4.4  |ggplot2: 2.2.1   |gridExtra: 2.3     |latex2exp: 0.4.0 |
|lme4: 1.1-17    |magic: 1.5-8     |magrittr: 1.5      |mnormt: 1.5-5    |
|msm: 1.6.6      |parallel: 3.4.4  |pgbme: 0.0.0.1     |PRROC: 1.3       |
|reshape2: 1.4.3 |sbgcop: 0.980    |tidyr: 0.8.1       |                 |

#### Using pbgme

When using this model for your own work, see the vignette for getting started. The package can be installed via github as follows: 

```
library(devtools)
devtools::install_github('s7minhas/pgbmeRepl', subdir='pgbme')
```

If you find any errors or have any further questions, please address them to me via email at minhassh@msu.edu.