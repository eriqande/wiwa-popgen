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

## TERMS

As a work, at least in part, of the United States Government, the contents
of this repository, except where specifically excepted, is in the
public domain within the United States. The co-authors of this package
who are not Federal employees also place their contributions into the
public domain. Additionally, we waive
copyright and related rights in the work worldwide through the CC0 1.0
Universal public domain dedication.

Note that we include code to download part of a map that was 
ade with Natural Earth. Free vector and raster map data @ naturalearthdata.com. This
is all in the public domain. Many thanks to its creators for putting it in the 
public domain.

### Exceptions to Public Domain
This packages includes code that downloads several files derived from differently-licensed materials.  These materials
are not distributed under the public domain license.  Specifically:

*  The files `wibreed_rast.nc` and `wiwinter_rast.nc` are raster files of the range of Wilson's Warblers in the Americas.
They were generated from shapefiles provided to us by BirdLife International.   (BirdLife International and NatureServe (2012) *Bird species distribution maps of the world*. BirdLife International, Cambridge, UK and NatureServe, Arlington, USA.)  The raster files are provided through download from an external link only for reproduction of the results of the paper we have published on this topic and should not be used for other purposes.  Persons interested in the range map should contact BirdLife International or NatureServe directly.  The full text of the license under which the original shapefiles and meta data were provided to us appears in TERMS-BLI.md.

## Reproducing our results

### Initial Maneuvers
Make sure that the current working directory is the directory that has the 
file `wiwa-popgen.Rproj` in it.  Then source the R code in the file:
```
R-main/01_wiwa_analysis_main_via_gpiper.R
```
That will create a directory called `outputs` and will put a lot of output into it. 
Some relevant bits of output are:

* Various plots.  All saved with a `.pdf` extension.  Some of these were not used in the paper. Some of them were, but
mostly kristen ruegg took elements from them and created plots to her liking in Illustrator.  Not super reproducible, but 
what can you do...
* Various output files. Mostly things that have a `.txt` extension.  Some of these are intermediate files, some are 
outputs of gsi_sim, etc.
* Other files that will be used as input later.  Most notably: 
     + An rda file called `WIWA-main-carryover-variables.Rda` that holds some variables that get used in downstream
analyses.
     + Three files: `wibreed-pipe*.txt` that are used for analysis with https://github.com/eriqande/slg_pipe


### Geneland analysis
Should you choose to do so, you can re-do our Geneland analysis.  We recommend using a fast computer
with at least 10 cores.  **If you don't want to re-run Geneland, that is OK! We provide relevant output so that you can continue with creating the map from our paper if so desired**

We use `mclapply` from the `parallel` package to run Geneland 10 times.  Do this by 
sourcing the R code in:
```
R-main/02_wiwa_popgen_geneland_run.R
```
This creates a boatload of output in a series of files with names beginning with `outputs/GeneLandRun`.  This step may take a
good 12 hours or more.


### Invesigating geneland maximum a posteriori estimates
Once the previous section has finished.  You can source this file:
```
R-main/03_wiwa_popgen_plot_geneland_posterior_modes.R
```
to make plots of the MAP estimates from Geneland to see that the results are pretty similar across most runs. 
If you didn't re-run Geneland, then you won't be able to create all of these plots.

The output plots are put into the directory `./outputs/geneland-posterior-mode-plots`

In our experience there were basically two different modes found in the 10 different runs of the program:
one of the modes had 7 clusters and the other had 6.  We have put a representative plot of each of these in 
the files `PostModeMap-6_clusters.pdf` and   `PostModeMap-7_clusters.pdf` in the directory `example_outputs` 

For our downstream work of making maps and delineating reporting groups we focused on the mode with
6 clusters, as it seemed to correspond well with the *structure* results and the self-assignment capability
of our SNP panel.

We include the spatial posterior probabilities of cluster membership from one of the 6-cluster-mode runs in 
`intermediates/proba.pop.membership.txt`.  This is the output that is used in the next section to make a pretty
map.  Of course, if you want to make your own map from different results that you have generated from
Geneland, you will have to replace `intermediates/proba.pop.membership.txt` with your own output.


### Making a pretty map with color transparency according to geneland posterior probabilities
This is done by sourcing the file:
```
R-main/04_wiwa_popgen_plot_pretty_geneland_map.R
```
into R. By default, this is run with the variables `REGENERATE_BASE_MAP` and `REGENERATE_POLY_RASTS`
set to `FALSE`, which means that, instead, the R-code downloads 3 files of about 200 Mb worth of stuff to
avoid some lengthy computations.  The files it downloads are:  `b3_cropped.nc` (a piece of the Natural Earth
Data base map), `wibreed_rast.nc` (a rasterized version of the breeding range of Wilson's warbler) and `wiwinter_rast.nc`
(a rasterized version of the wintering range).  The latter two files were generated from data provided by BirdLifeInternational.

If you set `REGENERATE_BASE_MAP` and `REGENERATE_POLY_RASTS` to `TRUE` then you will have 
you will have to download a suitable map from http://www.naturalearthdata.com/downloads/10m-raster-data/10m-cross-blend-hypso/
We used the "Cross Blended Hypso with Shaded Relief, Water, and Drainages".  You will also have to get the shapefiles from Wilson's Warblers directly
from BirdLifeInternational (see the **Exceptions to Public Domain** section above). And you will have to put them in the correct location, 
etc.

The main output from this script is:


