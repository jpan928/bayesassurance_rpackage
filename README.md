# bayesassurance R package

This R package offers a constructive set of simulation-based functions
used for determining sample size and assurance in various settings. 
We hope these functions will be useful for addressing a wide range of 
clinical trial study design problems. 

# Setup Instructions

To install the `bayesassurance` package in R, you can compile the package from
source as the package is not yet available on CRAN. 

## Directly From Github (Mac/Windows)
  1. Open R Studio.
  2. Make sure `devtools` is installed and loaded. If not, run `install.packages("devtools")` and 
  load the package using `library(devtools)` once installation is complete. 
  3. Install the bayesassurance package directly through Github by running
  `devtools::install_github("jpan928/bayesassurance_rpackage")`. 
  You may be asked to install `Rtools`. Please follow the instructions
  on https://cran.rstudio.com/bin/windows/Rtools/. 
  4. Load package using `library(bayesassurance)` and start using package normally. 


Alternatively, you can build the package using the tar.gz file.

## Mac
  1. Download the bayesassurance_0.1.0.tar.gz file. 
  2. In the R prompt, navigate to where this file is stored using `setwd()`. 
  3. Run `install.packages(bayesassurance_0.1.0.tar.gz, repos = NULL, type = "source")`. 
  4. Load package using `library(bayesassurance)` and start using package normally. 

 
## Windows
  1. Download the bayesassurance_0.1.0.tar.gz file.
  2. Open command prompt. 
  3. Identify path of the folder to where R is installed and run `PATH <your.file.path.here>`. 
  An example of this file path is C:\Program Files\R\R-4.1.3\bin\x64. 
  4. On the same command prompt, navigate to the directory containing the `bayesassurance_rpackage` folder.
  (Do not go into this folder.)
  5. Enter `R CMD INSTALL bayesassurance_0.1.0.tar.gz` in the command prompt to install the 
  bayesassurance package. 
  6. Open R Studio and run `library(bayesassurance)` and start using package normally. 
  

# Replication Materials

For JSS reviewers, R scripts containing the necessary code to reproduce figures and examples
in the manuscript can be found under `Replication_Material`. 
  
# Vignettes

Vignettes are currently undergoing revisions and will be available soon. 

