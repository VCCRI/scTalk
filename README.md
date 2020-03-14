# scTalk

scTalk is an R package for intercellular communication (ligand-receptor) analysis from scRNA-seq data and implements the method described in [Farbehi et al.](https://elifesciences.org/articles/43882).


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
