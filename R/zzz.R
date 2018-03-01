## Package version checking
.onLoad <- function(libname = find.package("dispRity"), pkgname = "dispRity") {
    packageStartupMessage(paste0("For the latest news, vignettes and releases,\ncheck on https://github.com/TGuillerme/dispRity"))
}