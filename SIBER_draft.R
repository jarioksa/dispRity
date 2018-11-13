
# remove previously loaded items from the current environment and remove previous graphics.
rm(list=ls())
graphics.off()

# Here, I set the seed each time so that the results are comparable. 
# This is useful as it means that anyone that runs your code, *should*
# get the same results as you, although random number generators change 
# from time to time.
set.seed(1)

library(SIBER)

# Load the viridis package and create a new palette with 3 colours, one for 
# each of the 3 groups we have in this dataset.
library(viridis)
palette(viridis(3))

# load in the included demonstration dataset
data("demo.siber.data")


#
# create the siber object
siber.example <- createSiberObject(demo.siber.data)

par(mfrow=c(1,1))

community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")

# this time we will make the points a bit smaller by 
# cex = 0.5
plotSiberObject(siber.example,
                  ax.pad = 2, 
                  hulls = F, community.hulls.args, 
                  ellipses = F, group.ellipses.args,
                  group.hulls = F, group.hull.args,
                  bty = "L",
                  iso.order = c(1,2),
                  xlab=expression({delta}^13*C~'\u2030'),
                  ylab=expression({delta}^15*N~'\u2030'),
                  cex = 0.5
                  )



# Calculate summary statistics for each group: TA, SEA and SEAc
group.ML <- groupMetricsML(siber.example)
print(group.ML)
#>            1.1       1.2       1.3       2.1       2.2       2.3
#> TA   21.924922 10.917715 17.945127 3.0714363 11.476354 1.4818061
#> SEA   5.783417  3.254484  5.131601 0.8623300  3.458824 0.4430053
#> SEAc  5.989967  3.370715  5.314872 0.8931275  3.582354 0.4588269


# You can add more ellipses by directly calling plot.group.ellipses()
# Add an additional p.interval % prediction ellilpse
plotGroupEllipses(siber.example, n = 100, p.interval = 0.95,
                    lty = 1, lwd = 2)



parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(siber.example, parms, priors)






# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                xlab = c("Community | Group"),
                ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                bty = "L",
                las = 1,
                main = "SIBER ellipses on each group"
                )

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML[3,], col="red", pch = "x", lwd = 2)







# Calculate some credible intervals 
cr.p <- c(0.95, 0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)




true.ellipse <- function(matrix, type) {

    ## Position of the ellipse (the mean of all axis)
    position <- apply(matrix, 2, mean)

    ## Eigen decomposition of the variance/covariance matrix
    decomposition <- eigen(var(matrix))

    ## Getting the semi-major axis
    semi_major <- 1/decomposition$values

    ## Getting the for the first axis
    angle <- asin(decomposition$vectors)
}

ellipse.estimate <- function(matrix, level) {
    ## Measure the ellipse from the variance covariance matrix
    ellipse_coords <- ellipse::ellipse(var(matrix),
                                       centre = apply(matrix, 2, mean),
                                       level = level,
                                       t = sqrt(qchisq(level, 2)),
                                       which = c(1, 2), npoints = 100
                                       )
}


# SIBER::standard.ellipse()
# ???::convexhull()
# SIBER:siber.ellipse()
# siar.density.plot()
# ???::overlap()
