#steps

install.packages(c("devtools", "roxygen2", "testthat", "knitr"))

create_tidy_package(".", copyright_holder = "JNS")
#create_package()

use_mit_license()
#update author name in LICENSE.md

#Rbuild file will ignore the vignette
use_article("GPAbin", "Getting started with GPAbin")

build_vignettes()

document()

install()
check()
build_vignettes()
