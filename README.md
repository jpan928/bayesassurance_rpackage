# bayesassurance R package

This R package offers a constructive set of simulation-based functions
and additional features used for determining sample size and assurance
in numerous settings.
These 
suitable for addressing a wide range of clinical trial study design
problems. 

# Setup Instructions

To install the `bayesassurance` package in R, you can compile the package from
source. 

  1. Make sure you have git properly installed.
  2. Run `git clone https://github.com/jpan928/bayesassurance_rpackage.git` from your
  command line and navigate to the cloned project directory via `cd bayesassurance_rpackage`.
  3. Install the package into R using `R CMD INSTALL bayesassurance`. You should now
  be able to run `library(bayesassurance)` and start using the package normally.
  
Alternatively, you can directly build the function in R if you have the
necessary components saved as a compressed file, bayesassurance_0.1.0.tar.gz. 

  1. Extract the files and save to desired location.
  2. Open R and navigate to the corresponding directory where the files are saved.
  3. Install and load in the package `devtools` using `install.packages("devtools")` and 
     `library(devtools)`. 
  4. Load in components of the package using `devtools::load_all()`. 
  5. Install package using `devtools::install()`, which installs the package into
     your R system library. 
  6. Load package using `library("bayesassurance")` and start using package normally. 
  
  
# Vignettes

Vignettes are currently undergoing revisions in preparation of being launched to CRAN. 

