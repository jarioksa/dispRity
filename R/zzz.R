## Package version checking
.onAttach <- function(libname = find.package("dispRity"), pkgname = "dispRity") {
    packageStartupMessage(paste0("       --- dispRity package ---\nThis is the CRAN release version (1.0.3) of the package.\nFor more functionalities, news, vignettes and releases,\nvisit https://github.com/TGuillerme/dispRity"))
}