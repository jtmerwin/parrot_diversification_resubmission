#' Colour volume overlap
#'
#' Compute colour volume overlap in 3D with alphashapes
#'
#' @inheritParams voloverlap
#' @inheritParams tcsshape
#' @param nsamp the number of points to be sampled. Stoddard & Stevens(2011) use
#' around 750,000 points, but more or fewer might be required depending on the
#' degree of overlap.
#'
#' @return a data.frame with the following columns:
#' - `s_in1, s_in2` the number of sampled points that fall within each of the
#' volumes individually.
#' - `s_inboth` the number of sampled points that fall within both volumes.
#' - `s_ineither` the number of points that fall within either of the volumes.
#' - `psmallest` the proportion of points that fall within both volumes divided
#' by the number of points that fall within the smallest volume.
#' - `pboth` the proportion of points that fall within both volumes divided by
#' the total number of points that fall within both volumes.
#'
#' If the Monte Carlo solution is used, a number of points much greater than the default should be
#' considered
#'
#' @examples
#' data(sicalis)
#' tcs.sicalis.C <- subset(colspace(vismodel(sicalis)), "C")
#' tcs.sicalis.T <- subset(colspace(vismodel(sicalis)), "T")
#' tcs.sicalis.B <- subset(colspace(vismodel(sicalis)), "B")
#' overlap3d(tcs.sicalis.T, tcs.sicalis.B, avalue = 1)
#' \donttest{
#' overlap3d(tcs.sicalis.T, tcs.sicalis.C, plot = TRUE, avalue = 0.1)
#' overlap3d(tcs.sicalis.T, tcs.sicalis.C, plot = TRUE, col = seq_len(3), avalue = 0.1)
#' }
#'
#' @importFrom stats runif
#'
#' @export

overlap3d <- function(colsp1, colsp2, plot = FALSE, interactive = TRUE,
                      col = c("blue", "red", "darkgrey"), fill = FALSE,
                      new = TRUE, nsamp = 1000, psize = 0.001, lwd = 1,
                      avalue, ...) {
  
  if (!interactive) {
    warning("interactive = FALSE has not been implemented yet, falling back to",
            " interactive plot.")
  }
  
  dat1 <- colsp1[, c("x", "y", "z")]
  
  dat2 <- colsp2[, c("x", "y", "z")]
  
  shape1 <- alphashape3d::ashape3d(as.matrix(dat1), avalue)
  shape2 <- alphashape3d::ashape3d(as.matrix(dat2), avalue)
  
  vol1 <- alphashape3d::volume_ashape3d(shape1)
  vol2 <- alphashape3d::volume_ashape3d(shape2)
  
  # sample random points
  pmin <- apply(rbind(dat1, dat2), 2, min)
  pmax <- apply(rbind(dat1, dat2), 2, max)
  
  samples <- apply(rbind(pmin, pmax), 2, function(x) runif(nsamp, x[1], x[2]))
  
  invol1 <- alphashape3d::inashape3d(shape1, points = samples)
  invol2 <- alphashape3d::inashape3d(shape2, points = samples)
  
  # how many points are in each category
  
  s_in1 <- sum(invol1)
  s_in2 <- sum(invol2)
  
  s_inboth <- sum(invol1 & invol2)
  
  s_ineither <- sum(invol1 | invol2)
  
  # points in both relative points in smallest
  
  psmallest <- s_inboth / c(s_in1, s_in2)[which.min(c(vol1, vol2))]
  
  # points in both relative to total points in both
  
  pboth <- s_inboth / s_ineither
  
  res <- data.frame(vol1, vol2, s_in1, s_in2, s_inboth, s_ineither, psmallest, pboth)
  
  ############
  # PLOT BEGIN#
  ############
  if (plot) {
    if (length(col) < 3) {
      col <- c(rep(col, 2)[seq_len(2)], "darkgrey")
    }
    # check if rgl is installed and loaded
    if (!requireNamespace("rgl", quietly = TRUE)) {
      stop(dQuote("rgl"), " package needed for interactive plots. Please install it, or use interactive=FALSE.",
           call. = FALSE
      )
    }
    
    if (!isNamespaceLoaded("rgl")) {
      requireNamespace("rgl")
    }
    
    if (new) {
      rgl::open3d(FOV = 1, mouseMode = c("zAxis", "xAxis", "zoom"))
    }
    
    tcsshape(colsp1, col = col[1], fill = FALSE, avalue = avalue)
    tcsshape(colsp2, col = col[2], fill = FALSE, avalue = avalue)
    
    rgl::spheres3d(samples[which(invol1 & !invol2), ],
                   type = "s",
                   lit = FALSE, radius = psize, col = col[1]
    )
    rgl::spheres3d(samples[which(invol2 & !invol1), ],
                   type = "s",
                   lit = FALSE, radius = psize, col = col[2]
    )
    
    if (s_inboth > 0) {
      rgl::spheres3d(samples[which(invol1 & invol2), ],
                     type = "s",
                     lit = FALSE, radius = psize, col = col[3]
      )
    }
  }
  ##########
  # PLOT END#
  ##########
  
  res
}

#' Produces a 3D alphashape in tetrahedral colour space
#'
#' @inheritParams tcsvol
#' @inheritParams tetrashape
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # For plotting
#' data(sicalis)
#' vis.sicalis <- vismodel(sicalis, visual = "avg.uv")
#' tcs.sicalis <- colspace(vis.sicalis, space = "tcs")
#' tcsplot(tcs.sicalis, col = "blue", size = 0.005)
#' tcsshape(tcs.sicalis, avalue = 0.1)
#' }

tcsshape <- function (tcsdata, col = "black", alpha = 0.2, grid.alpha = 1,
                      grid = TRUE, fill = TRUE, lwd = 1, avalue) {
  if (attr(tcsdata, "clrsp") != "tcs") stop("object is not in tetrahedral color space")
  
  # check if rgl is installed and loaded
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop(dQuote("rgl"), " package needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  
  if (!isNamespaceLoaded("rgl")) {
    requireNamespace("rgl")
  }
  
  ashape <- alphashape3d::ashape3d(as.matrix(tcsdata[, c("x", "y", "z")]), 
                                   alpha = avalue)
  tri <- ashape$triang
  vol <- t(tri[tri[, ncol(tri)] %in% c(2,3), c(1, 2, 3)])
  coords <- ashape$x
  
  if (grid) {
    rgl::triangles3d(coords[vol, ],
                     color = col, alpha = grid.alpha, lwd = lwd,
                     front = "lines", back = "lines"
    )
  }
  
  if (fill) {
    rgl::rgl.triangles(coords[vol, ],
                       alpha = alpha, color = col
    )
  }
  
  rgl::material3d(alpha = 1)
}

#' Plot a (possibly concave) tetrahedral colour space
#'
#' Produces a 3D alphashape in tetrahedral colour space when plotting a
#' non-interactive tetrahedral plot.
#'
#' @inheritParams vol
#' @param avalue Value of alpha parameter for the alphashape
#'
#' @examples
#'
#' @examples
#' \dontrun{
#' # For plotting
#' data(sicalis)
#' vis.sicalis <- vismodel(sicalis, visual = "avg.uv")
#' tcs.sicalis <- colspace(vis.sicalis, space = "tcs")
#' plot(tcs.sicalis)
#' tetrashape(tcs.sicalis, avalue = 0.1)
#' }
#'
#' @export

tetrashape <- function(tcsdata, alpha = 0.2, grid = TRUE, fill = TRUE,
                       new = FALSE, avalue, ...) {
  if (!is.null(attr(tcsdata, "clrsp")) && attr(tcsdata, "clrsp") != "tcs")
    stop("object is not in tetrahedral color space")
  
  ashape <- alphashape3d::ashape3d(as.matrix(tcsdata[, c("x", "y", "z")]), 
                                   alpha = avalue)
  
  tri <- ashape$triang
  vol <- t(tri[tri[, ncol(tri)] %in% c(2,3), c(1, 2, 3)])
  coords <- ashape$x
  
  arg <- list(...)
  
  if (new) {
    argempty <- c(list(border = FALSE), arg)
    
    argempty$col <- NULL
    
    if (is.null(argempty$zoom)) {
      argempty$zoom <- 1
    }
    if (is.null(argempty$xlim)) {
      argempty$xlim <- range(tcsdata[, "x"]) / argempty$zoom
    }
    if (is.null(argempty$ylim)) {
      argempty$ylim <- range(tcsdata[, "y"]) / argempty$zoom
    }
    if (is.null(arg$zlim)) {
      argempty$zlim <- range(tcsdata[, "z"]) / argempty$zoom
    }
    
    argempty$zoom <- NULL
    
    if (is.null(argempty$theta)) {
      argempty$theta <- 45
    }
    if (is.null(argempty$phi)) {
      argempty$phi <- 10
    }
    if (is.null(argempty$r)) {
      argempty$r <- 12
    }
    if (is.null(argempty$box)) {
      argempty$box <- FALSE
    }
    
    if (is.null(argempty$margin)) {
      margin <- c(0, 0, 0, 0)
    } else {
      margin <- argempty$margin
      argempty$margin <- NULL
    }
    
    argempty$x <- argempty$xlim
    argempty$y <- argempty$ylim
    argempty$z <- matrix(c(argempty$zlim, argempty$zlim), nrow = 2)
    
    # empty plot
    oPar <- par("mar", "pty")
    on.exit(par(oPar))
    
    par(mar = margin)
    
    P <- do.call(perspbox, argempty)
    
    # Save plot info
    assign("last_plot.tetra", P, envir = .PlotTetraEnv)
  }
  
  
  last_tetraplot <- tryCatch(
    get("last_plot.tetra", envir = .PlotTetraEnv),
    error = function(e) get("last_plot.tetra", envir = pavo:::.PlotTetraEnv)
  )
  
  flatcoords <- data.frame(trans3d(coords[, "x"], coords[, "y"], coords[, "z"], last_tetraplot))
  
  if (is.null(arg$col)) {
    arg$col <- "darkgrey"
  }
  
  darkcolor <- arg$col
  alphacolor <- adjustcolor(arg$col, alpha.f = alpha)
  
  if (fill) {
    arg$border <- NA
    arg$col <- alphacolor
  } else {
    arg$col <- NA
  }
  
  if (grid) {
    arg$border <- darkcolor
  }
  
  # CRAN won't accept triple : arguments and persp.default is not exported,
  # so we need to pass arguments by hand
  perspargs <- c(
    "x", "y", "z", "xlim", "ylim", "zlim", "xlab", "ylab", "zlab",
    "main", "sub", "theta", "phi", "r", "d", "scale", "expand",
    "ltheta", "lphi", "shade", "box", "axes", "nticks", "ticktype", "...", ""
  )
  
  arg[perspargs] <- NULL
  
  
  for (i in seq_len(ncol(vol))) {
    arg$x <- flatcoords[vol[, i], "x"]
    arg$y <- flatcoords[vol[, i], "y"]
    
    do.call(polygon, arg)
  }
}
