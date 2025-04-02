
## Date Last Modified: 4/2/25
## Author: Adithya Narayanan


## Since covidhm package was authored, many of these packages have undergone updates, 
## and therefore do not really work the way you would intend to if you install from CRAN. 
## Using these specific versions of packages, installed using the code below, 
## the covidhm installation was successful. If you need assistance with this, 
## please reach out to the author.

## Download the covidhm files from https://github.com/biouea/covidhm to install the package. 
## The installation method specifies in the package's repo is a bit outdated and needs specifications of versions of dependencies.


install.packages('devtools')
require(devtools)

directorypath = "/cloud/project/covidhm-master" # point this to your location of the downloaded covidhm package files



# tibble
install_version("dplyr", version = "1.0.4", repos = "http://cran.us.r-project.org")
install_version("purrr", version = "0.3.4", repos = "http://cran.us.r-project.org")
install_version("plyr", version = "1.8.6", repos = "http://cran.us.r-project.org")


# Matrix, MatrixModels, quantreg, magrittr
install_version("sn", version = "2.0.2", repos = "http://cran.us.r-project.org")

# dotCall64, spam
install_version("fields", version = "13.3", repos = "http://cran.us.r-project.org")

install_version("diagram", version = "1.6.5", repos = "http://cran.us.r-project.org")
install_version("shape", version = "1.4.5", repos = "http://cran.us.r-project.org")
install_version("clue", version = "0.3-58", repos = "http://cran.us.r-project.org")

install_version("reshape2", version = "1.4.4", repos = "http://cran.us.r-project.org")

install_version("wfg", version = "0.1", repos = "http://cran.us.r-project.org")


Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE)
devtools::install("covidhm-master",dependencies = TRUE)




### If the following code runs and produces visuals, you should be all set.
## If it yields an error, particularly around mutate, or slice from dplyr, reach out to the author.
setwd(directorypath)
library(covidhm)

#Load association matrices - ensure path is correct for data-raw directory from the covidhm package
load("data-raw/am_list.RData")

#First item in the list is data across all days
m <- am_list[[1]]

#Plot network
plot_network(
  am = m,
  day = 20,
  num.initial.cases = 1,
  prop.asym = 0.4,
  delay_shape =  1,
  delay_scale = 1.4,
  prop.ascertain = 0.9,
  presymrate = 0.2,
  R = 0.8,
  outside = 0.001,
  testing = FALSE,
  s = 333,
  isolation = FALSE,
  secondary = FALSE,
  tracing = FALSE,
  quarantine = FALSE)



library(covidhm)
library(ggplot2)

res <- scenario_sim(net = m, n.sim = 10, num.initial.cases = 1,prop.asym=0.4,
                    prop.ascertain = 0.9, cap_max_days = 70,
                    delay_shape = 1, delay_scale = 1.4, R = 0.8, presymrate = 0.2, scenario = "nothing",
                    testing = FALSE, outside = 0.001, distancing = 0)

# Plot of raw cumulative cases
ggplot(data=res, aes(x=week, y=cumcases,col = sim)) +
  geom_line(show.legend = FALSE, alpha=0.6, aes(group = sim)) +
  scale_y_continuous(name="Weekly number of cases") +
  theme_bw()

