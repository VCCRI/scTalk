# scTalk

scTalk is an R package for intercellular communication (ligand-receptor) analysis from scRNA-seq data and implements the method described in [Farbehi et al.](https://elifesciences.org/articles/43882). Please read the vignette for an explanation of how to use the software. 

## Installation
To install scTalk, open an R session and use the following commands:

```
install.packages("devtools")
devtools::install_github("VCCRI/scTalk", build = TRUE, build_vignettes = TRUE, build_opts = c("--no-resave-data", "--no-manual"))
```
To access the vignette:

```
library(scTalk)
browseVignettes("scTalk")
```


## Citation
If you find scTalk useful in your research, please cite the following paper:

Farbehi N, Patrick R, Dorison A, Xaymardan M, Janbandhu V, Wystub-Lis K, Ho JWK, Nordon RE, and Harvey RP. Single-cell expression profiling reveals dynamic flux of cardiac stromal, vascular and immune cells in health and injury. eLife, 8:e43882, 2019.

