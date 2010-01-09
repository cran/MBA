\name{mba.points}
\alias{mba.points}
\title{Point approximation from bivariate scattered data using multilevel B-splines}

\description{
  The function \code{mba.points} returns points on a surface
  approximated from a bivariate scatter of points using multilevel B-splines.
}

\usage{
mba.points(xyz, xy.est, n = 1, m = 1, h = 8, extend = TRUE,
           verbose = TRUE, ...)
}

\arguments{
  \item{xyz}{a \eqn{n \times 3}{n x 3} matrix or data frame, where \eqn{n} is
    the number of observed points.  The three columns correspond to point x, y, and z
    coordinates.  The z value is the response at the given x, y
    coordinates.  }
  \item{xy.est}{a \eqn{p \times 2}{p x 2} matrix or data frame, where
    \eqn{p} is the number of points for which to estimate a z.  The two
    columns correspond to x, y point coordinates where a z estimate is
    required. }
  \item{n}{initial size of the spline space in the hierarchical
    construction along the x axis. If the rectangular domain is a
    square, n = m = 1 is recommended. If the x axis is k times the length
    of the y axis, n = 1, m = k is recommended. The default is n = 1. }
  \item{m}{initial size of the spline space in the hierarchical
    construction along the y axis. If the y axis is k times the length
    of the x axis, m = 1, n = k is recommended. The default is m = 1.  }
  \item{h}{Number of levels in the hierarchical construction. If, e.g.,
    n = m = 1 and h = 8, the resulting spline surface has a coefficient
    grid of size \eqn{2^h}{2^h} + 3 = 259 in each direction of the
    spline surface. See references for additional information. }
  \item{extend}{if FALSE, points in \code{xy.est} that fall outside of
    the domain defined by \code{xyz} are set to NA with a warning; otherwise, the
    domain is extended to accommodate points in \code{xy.est}
    with a warning. }
  \item{verbose}{if TRUE, warning messages are printed to the screen. }
   \item{...}{currently no additional arguments.  }
}

\value{
  List with 1 component:
  \item{xyz.est}{a \eqn{p \times 3}{p x 3} matrix.  The first two
    columns are \code{xy.est} and the third column is the corresponding z estimates.  }
}

\note{
The function \code{mba.points} relies on the Multilevel B-spline
Approximation (MBA) algorithm.  The underlying code was developed at
SINTEF Applied Mathematics by Dr. Oyvind Hjelle.  Dr. Oyvind Hjelle based the
algorithm on the paper by the originators of Multilevel B-splines:

S. Lee, G. Wolberg, and S. Y. Shin. Scattered data interpolation with
multilevel B-splines. IEEE Transactions on Visualization and Computer
Graphics, 3(3):229-244, 1997.

For additional documentation and references please see:

\url{www.sintef.no/upload/IKT/9011/geometri/MBA/mba_doc/index.html}.

This minor portion of the MBA codebase was ported by Andrew O. Finley \email{finleya@msu.edu}.
}

\examples{
data(LIDAR)

##split the LIDAR dataset into training and validation sets
tr <- sample(1:nrow(LIDAR),trunc(0.5*nrow(LIDAR)))

##look at how smoothing changes z-approximation,
##careful the number of B-spline surface coefficients
##increases at ~2^h in each direction
for(i in 1:10){
  mba.pts <- mba.points(LIDAR[tr,], LIDAR[-tr,c("x","y")], h=i)$xyz.est
  print(sum(abs(LIDAR[-tr,"z"]-mba.pts[,"z"]))/nrow(mba.pts))
}

\dontrun{
##rgl or scatterplot3d libraries can be fun
library(rgl)

##exaggerate z a bit for effect and take a smaller subset of LIDAR
LIDAR[,"z"] <- 10*LIDAR[,"z"]
tr <- sample(1:nrow(LIDAR),trunc(0.99*nrow(LIDAR)))

##get the "true" surface
mba.int <- mba.surf(LIDAR[tr,], 100, 100, extend=TRUE)$xyz.est

open3d()
surface3d(mba.int$x, mba.int$y, mba.int$z)

##now the point estimates
mba.pts <- mba.points(LIDAR[tr,], LIDAR[-tr,c("x","y")])$xyz.est
spheres3d(mba.pts[,"x"], mba.pts[,"y"], mba.pts[,"z"],
          radius=5, color="red")
}
}

\seealso{\code{\link{mba.surf}}}

\keyword{dplot}
\keyword{smooth}
