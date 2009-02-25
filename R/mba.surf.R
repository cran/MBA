"mba.surf" <- function(xyz, no.X, no.Y,  n = 1, m = 1, h = 8, extend = FALSE, sp=FALSE, ...){

  ####################################################
  ##Check for unused args
  ####################################################
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
    if(! i %in% formal.args)
      warning("'",i, "' is not a parameter")
  }
  
  if(missing(xyz)){stop("error: xyz matrix or data frame must be specified")}
  if(missing(no.X) || no.X <= 1){stop("error: no.X must be specified as an integer greater than 1")}
  if(missing(no.Y) || no.Y <= 1){stop("error: no.Y must be specified as an integer greater than 1")}
  
  if(!any(is.matrix(xyz), is.data.frame(xyz))){stop("error: xyz must be a matrix or data frame")}
  if(any(ncol(xyz) != 3, nrow(xyz) == 0)){stop("error: xyz must have 3 columns corresponding to x, y, z and at least one row")}

  if(m <= 0){stop("error: m must be a positive integer")}
  if(n <= 0){stop("error: n must be a positive integer")}
  if(h <= 0){stop("error: h must be a positive integer")}

  if(!extend){
    hpts <- chull(xyz[,c(1,2)])
    hpts <- c(hpts, hpts[1])
  }else{
    hpts <- NULL
  }
  
  xyz <- as.matrix(xyz)
  storage.mode(xyz) <- "double"
  
  out <- .Call("MBASurf", xyz, as.integer(no.X), as.integer(no.Y), as.integer(m), as.integer(n), as.integer(h), as.integer(extend), as.integer(hpts))

  if(sp){
    xy <- expand.grid(out[["x"]], out[["y"]])
    grid <- data.frame(z=matrix(out[["z"]], as.integer(length(out[["x"]])*length(out[["y"]])), 1),x=xy[,1], y=xy[,2])
    coordinates(grid) = ~x+y
    gridded(grid) <- TRUE
  }else{
    grid <- out[c("x","y","z")]
  }
  
  out <- list()
  out$xyz.est <- grid
  out$no.X <- no.X
  out$no.Y <- no.Y
  out$n <- n
  out$m <- m
  out$h <- h
  out$extend <- extend
  out$sp <- sp
  out
}
