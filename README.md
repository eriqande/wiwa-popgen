# wiwa-popgen

This is a collection of R scripts and data for the population
genetics analysis of 96 SNPs typed on numerous Wilson's
Warblers.  

We call it *popgen* as distinct from the bioinformatics phase
of the project in which we sifted through a lot of 
next generation sequencing data to design assays for the
96 SNPs.

Note that the order of analyses as presented here is different from
the way we did things, i.e. we first ran **structure** and **geneland**, and 
the results of those analyses informed our choice of reporting units for the 
genetic assignment analyses.  Here we have a script that processes the data
set, does genetic assignment first, and then outputs data files for 
**structure** and also runs geneland.  It doesn't really matter, as this 
still reproduces all of our results.


## Reproducing our results
Make sure that the current working directory is the directory that has the 
file `wiwa-popgen.Rproj` in it.  Then do this in R:
```
source("R-main/wiwa_analysis_main_via_gpiper.R")
```
That will create a directory called `outputs` and will put a lot of output into it. 
Some relevant bits of output are:
* Various plots:
  + Boing
  + Bong
* Various output files:
  + foo
  + bar
* An rda file called `name` that holds some variables that get used in downstream
analyses.