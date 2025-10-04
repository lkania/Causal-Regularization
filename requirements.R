options(repos = c(CRAN = "https://cloud.r-project.org/"))

getPackage <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    message(paste("Installing package:", package))
    install.packages(package)
  }
  message(paste("Package", package, "is ready to use."))
}

# computation
getPackage('data.table')
getPackage('limSolve')
getPackage('future.apply')
getPackage('digest')
getPackage('purrr')
getPackage('dplyr')
getPackage('scales')
getPackage('latex2exp')
getPackage('ggplot2')
getPackage('reshape2')
getPackage('ggpubr')
getPackage('ggExtra')
getPackage('grid')
getPackage('gridExtra')
getPackage('ggpubr')
getPackage('gtable')
getPackage('ggnewscale')