This repository includes the P-GBME package and replication files for [“Modeling Asymmetric Relationships from Symmetric Networks”](https://arxiv.org/pdf/1711.03838.pdf).

## Abstract

Many bilateral relationships requiring mutual agreement produce observable networks that are symmetric (undirected). However, the unobserved, asymmetric (directed) network is frequently the object of scientific interest. We propose a method that probabilistically reconstructs the latent, asymmetric network from the observed, symmetric graph in a regression-based framework. We apply this model to the bilateral investment treaty network. Our approach successfully recovers the true data generating process in simulation studies, extracts new, politically relevant information about the network structure inaccessible to alternative approaches, and has superior predictive performance.

## Installing P-GBME

P-GBME can be installed via devtools: 

```
library(devtools)
devtools::install_github('s7minhas/pgbmeRepl', subdir='pgbme')
```

A brief vignette is included with the package and can also be viewed [here](https://cdn.rawgit.com/s7minhas/pgbmeRepl/master/pgbme/vignettes/pgbmeVignette.html).

Publication Outlet
---
Forthcoming in Political Analysis
