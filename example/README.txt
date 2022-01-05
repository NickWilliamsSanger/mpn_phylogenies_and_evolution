This directory contains a reproducible set of analyses for the patient PD6629.

Dependencies:  

rmarkdown
dplyr
data.table
readr
rstan
karyoploteR
kableExtra
metafor

To install the required packages rsimpop and rtreefit:

devtools::install_github("NickWilliamsSanger/rsimpop")
devtools::install_github("NickWilliamsSanger/rtreefit",build_vignettes=TRUE)

To run, start R and cd to the example_code directory (here denoted by "<example_code>"):

* setwd("<example_code>")
* library("rmarkdown")
* render("SingleSampleAnalysisExample.Rmd",output_file = "SingleSampleAnalysisExample_RERUN.html",output_dir="html")

Note that the output_dir should be an existing directory and the output file will be stored there.

A command line alternative is:

R -e "rmarkdown::render('SingleSampleAnalysisExample.Rmd',output_file='SingleSampleAnalysisExample_RERUN.htm',output_dir='html')"

Note that specifying TREEMODEL="nb_tree" will fit a negative binomial model instead of the default "poisson_tree" model.