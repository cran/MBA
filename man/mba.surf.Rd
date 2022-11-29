\name{mba.surf}
\alias{mba.surf}
\title{Surface approximation from bivariate scattered data using multilevel B-splines}

\description{
  The function \code{mba.surf} returns a surface approximated from a
  bivariate scatter of data points using multilevel B-splines.
}

\usage{
mba.surf(xyz, no.X, no.Y, n = 1, m = 1, h = 8, extend=FALSE,
         sp=FALSE, ...)
}

\arguments{
  \item{xyz}{a \eqn{n \times 3}{n x 3} matrix or data frame, where \eqn{n} is
    the number of observed points.  The three columns correspond to point x, y, and z
    coordinates.  The z value is the response at the given x, y
    coordinates.  }
  \item{no.X}{resolution of the approximated surface along the x axis. }
  \item{no.Y}{resolution of the approximated surface along the y axis. }
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
  \item{extend}{if FALSE, a convex hull is computed for the input points
    and all matrix elements in z that have centers outside of this
    polygon are set to \code{NA}; otherwise, all elements in z are given an
    estimated z value. }
  \item{sp}{if TRUE, the resulting surface is returned as a
  \code{SpatialPixelsDataFrame} object; otherwise, the surface is in
  \code{image} format.}
  \item{...}{\code{b.box} is an optional vector to sets the bounding
    box. The vector's elements are minimum x, maximum x, minimum y, and maximum
    y, respectively.  }
}

\value{
  List with 8 component:
  \item{xyz.est}{a list that contains vectors x, y and the \eqn{no.X
      \times no.Y}{no.X x no.Y} matrix z of estimated z-values. }
  \item{no.X}{\code{no.X} from arguments.} 
  \item{no.Y}{\code{no.Y} from arguments.} 
  \item{n}{\code{n} from arguments.} 
  \item{m}{\code{m} from arguments.} 
  \item{h}{\code{h} from arguments.} 
  \item{extend}{\code{extend} from arguments.} 
  \item{sp}{\code{sp} from arguments.}
  \item{b.box}{\code{b.box} defines the bounding box over which z is estimated.}
}

\note{
If \code{no.X != no.Y} then use \code{sp=TRUE} for compatibility with
the \code{image} function.
  
The function \code{mba.surf} relies on the Multilevel B-spline
Approximation (MBA) algorithm.  The underlying code was developed at
SINTEF Applied Mathematics by Dr. Øyvind Hjelle.  Dr. Øyvind Hjelle
based the algorithm on the paper by the originators of Multilevel B-splines:

S. Lee, G. Wolberg, and S. Y. Shin. (1997) Scattered data interpolation with
multilevel B-splines. IEEE Transactions on Visualization and Computer
Graphics, 3(3):229--244.

For additional documentation and references see:

\url{https://www.sintef.no/upload/IKT/9011/geometri/MBA/mba_doc/index.html}.
}

\examples{
\dontrun{
data(LIDAR)

mba.int <- mba.surf(LIDAR, 300, 300, extend=TRUE)$xyz.est

##Image plot
image(mba.int, xaxs="r", yaxs="r")

##Perspective plot
persp(mba.int, theta = 135, phi = 30, col = "green3", scale = FALSE,
      ltheta = -120, shade = 0.75, expand = 10, border = NA, box = FALSE)

##For a good time I recommend using rgl
library(rgl)

ex <- 10
x <- mba.int[[1]]
y <- mba.int[[2]]
z <- ex*mba.int[[3]]
zlim <- range(z)
zlen <- zlim[2] - zlim[1] + 1
colorlut <- heat.colors(as.integer(zlen))
col <- colorlut[ z-zlim[1]+1 ]

open3d()
surface3d(x, y, z, color=col, back="lines")
}
}

\seealso{\code{\link{mba.points}}}

\keyword{dplot}
\keyword{smooth}
