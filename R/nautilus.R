#' @title Nautilus Aka Snail Plot
#' @description Draws a modified stars (or radar) plot.
#' @param x numeric matrix of the values to be plotted with observations in rows and variables in columns
#' @param col line colors to be used, e.g. according to a grouping factor, vector of type character and same length as rows in x
#' @param lwd line width, Default: 1
#' @param FUN function to be used for log or log-like transformation, if "lin" no transformation is performed, Default: c("qlog10", "asinh.co", "asinh", "log10", "lin")
#' @param cofactor for asinh.co as e.g. used for (non-negative) cytometry data transformation, Default: 0.2
#' @param seed for randomizing order of plotted lines (to prevent masking), may be turned of by setting to NA, Default: 25
#' @param p quantile of the data which is taken to order variables, e.g. p = 0.5 for medians or p = 1 for maxima, Default: 1
#' @param coef used to remove outliers by IQR rule when set to positive value, e.g. 1.5 as in \code{\link[grDevices]{boxplot.stats}}, Default: 0
#' @param untie this provides an experimental option to reorder fractions of the plot in order to bring correlated variables closer. It uses hierarchical clustering for reordering of at least 3 variables per fraction., Default: 0
#' @param switch.cols logical whether to reverse reordering within each fraction, Default: FALSE
#' @param distfun distance function used in clustering for within fraction reordering, Default: \code{\link[stats]{dist}}
#' @param hclustfun linkage function used in clustering for within fraction reordering, Default: \code{\link[stats]{hclust}}
#' @param NAUT an object as returned from a previous call to the function. This will override calculation of all dimensions of the new plot, including variable ordering. This can come handy if data shall be plotted at same scale but different figures 
#' @param lwd.shell width of the outer line, Default: 0.45
#' @param col.shell color of the outer line, Default: 'black'
#' @param cex character expansion factor for the labels, Default: 0.45
#' @param adj label adjustment, see \code{\link[graphics]{par}}, Default: 0.2
#' @param plot.major logical whether to plot main scale lines, Default: TRUE
#' @param lwd.major.l width of the smallest main scale line, Default: 0.15
#' @param lwd.major.g width of the greatest main scale line. Line widths in-between will be ramped according to these two values, Default: 0.45
#' @param col.major color of main scale lines, Default: 'gray'
#' @param plot.minor logical whether to plot minor scale lines, only effective for log10- or asinh-transformed data, Default: TRUE
#' @param lwd.minor.l width of the smallest minor scale line, Default: 0.05
#' @param lwd.minor.g width of the greatest minor scale line. Line widths in-between will be ramped according to these two values, Default: 0.15
#' @param col.minor color of minor scale lines, Default: 'gray80'
#' @param plot.axis logical whether to plot axis, Default: TRUE
#' @param cex.axis expansion factor for axis labels, Default: 0.35
#' @param padj.axis adjustment for each tick label, see \code{\link[graphics]{axis}}, Default: -4.5
#' @param lwd.major.axis line width for the axis line, Default: 0.75
#' @param lwd.major.ticks line widths the main tick marks, Default: 0.75
#' @param tck.major length of main tick marks, Default: -0.015
#' @param lwd.minor.axis line width for the axis line, Default: 0.05
#' @param lwd.minor.ticks line widths for the minor tick marks, Default: 0.5
#' @param tck.minor length of minor tick marks, Default: -0.0075
#' @param plot.vertical logical whether to plot vertical lines, Default: TRUE
#' @param lwd.vertical width of vertical lines, Default: 0.15
#' @param col.vertical color of vertical lines, Default: 'gray'
#' @return  Invisibly, a list object with following fields:
#' \describe{
#'   \item{\code{FUN}}{transformation used}
#'   \item{\code{SCALES}}{matrix containing plotting coordinates for the shell}
#'   \item{\code{ORD}}{used order of variables}
#'   \item{\code{SHIFT}}{constant to move all values on positive scale}
#'   \item{\code{main.ticks}}{main tick mark coordinates}
#'   \item{\code{minor.ticks}}{minor tick mark coordinates}
#'   \item{\code{log.ticks}}{main tick labels on log scale}
#' }
#' @details Starting with 10 to 15 variables, the stars (or radar) plot can be made to resemble a nautilus (\eqn{Nautilus belauensis}). With more than 50 markers it will however be difficult to grasp the actual expression profiles. Besides the standard log10 transformation, the plotting function offers asinh-based transformations in order to bring the data properly on scale. 
#' @examples 
#' \dontrun{
#' 
#' # mock up data with e.g. 100 samples total, four groups and 30 variables
#' set.seed(25)
#' g <- rep(1:4, each = 25)
#' sc <- 10^sample(seq(0, 3, length.out = 100), 30, TRUE)
#' y <- sapply(sc, function(x) rnorm(100, x, x / 2) + x * rep(sample.int(4), each = 25))
#' colnames(y) <- colnames(y, do.NULL = FALSE, prefix = "VAR")
#' col <- palette()[g + 1] # according to group code
#' 
#' # asinh.co allows display of negative data on "quasi-" logarithmic scale
#' nautilus(y, col = col, FUN = "asinh.co", padj.axis = -5.5)
#' # pair correlated variables
#' nautilus(y, col = col, FUN = "asinh.co", cofactor = 10, padj.axis = -5.5, untie = 4)
#' nautilus(y, col = col, FUN = "asinh.co", cofactor = 10, padj.axis = -5.5, untie = 4, 
#'             switch.cols = TRUE, plot.major = FALSE, plot.minor = FALSE, plot.vertical = FALSE)
#' # default asinh
#' nautilus(y, col = col, FUN = "asinh", padj.axis = -5.5)
#' # quasi-log10 
#' nautilus(y, col = col, FUN = "qlog10", padj.axis = -5.5)
#' 
#' # default log10
#' nautilus(y, col = col, FUN = "log10")
#' 
#' # log-normal data
#' set.seed(25)
#' g <- rep(1:4, each = 25)
#' sc <- 10^sample(seq(-2, 3, length.out = 100), 30, TRUE)
#' y <- sapply(sc, function(x) rnorm(100, x, abs(x))  +  x * rep(sample.int(4), each = 25))
#' colnames(y) <- colnames(y, do.NULL = FALSE, prefix = "VAR")
#' col <- palette()[g + 1]
#' nautilus(y, col = col, FUN = "log10")
#' 
#' # linear scale
#' nautilus(y, col = col, FUN = "lin", adj = 500)
#' 
#' # normal data
#' set.seed(25)
#' g <- rep(1:4, each = 25)
#' sc <- sample(seq(-1, 5, length.out = 100), 30, TRUE)
#' y <- 100 + sapply(sc, function(x) rnorm(100, x, abs(x) / 2)  +  x * rep(sample.int(4), each = 25))
#' colnames(y) <- colnames(y, do.NULL = FALSE, prefix = "VAR")
#' col <- palette()[g + 1]
#' 
#' nautilus(y, col = col, FUN = "lin", p = 0.5, padj.axis = -5.5)
#' 
#' }
#' @seealso 
#'  \code{\link[graphics]{stars}}, \code{\link[graphics]{axis}}
#' @rdname nautilus
#' @export 
#' @import stats
#' @importFrom graphics stars axis text
nautilus <- function(x,
                     col,
                     lwd = 1,
                     FUN = c("qlog10", "asinh.co", "asinh", "log10", "lin"),
                     cofactor = 0.2,
                     seed = 25,
                     p = 1,
                     coef = 0,
                     untie = 0,
                     switch.cols = FALSE,
                     distfun = dist,
                     hclustfun = hclust,
                     # may provide previous plot object
                     NAUT,
                     
                     # tweak outer line
                     lwd.shell = 0.45,
                     col.shell = "black",
                     
                     # label size
                     cex = 0.45,
                     adj = 0.2,
                     
                     plot.major = TRUE,
                     # tweak major scale lines
                     lwd.major.l = 0.15,
                     lwd.major.g = 0.45,
                     col.major = "gray",
                     
                     plot.minor = TRUE,
                     # tweak minor scale lines
                     lwd.minor.l = 0.05,
                     lwd.minor.g = 0.15,
                     col.minor = "gray80",
                     
                     plot.axis = TRUE,
                     # tweak axis and axis ticks
                     cex.axis = 0.35,
                     padj.axis = -4.5,
                     lwd.major.axis = 0.75,
                     lwd.major.ticks = 0.75,
                     tck.major = -0.015,
                     lwd.minor.axis = 0.05,
                     lwd.minor.ticks = 0.5,
                     tck.minor = -0.0075,
                     
                     plot.vertical = TRUE,
                     # tweak vertical scale lines
                     lwd.vertical = 0.15,
                     col.vertical = "gray")

{
  # internal to get correct label locations
  stars.label <- function(x, locations, len, key.loc, cex, adj)
  {
    x <- data.matrix(x)
    key.labels = dimnames(x)[[2L]]
    n.loc <- nrow(x)
    n.seg <- ncol(x)
    locations <-
      cbind(rep.int(locations[1L], n.loc), rep.int(locations[2L], n.loc))
    xloc <- locations[, 1]
    yloc <- locations[, 2]
    angles <-
      seq.int(0, 2 * pi, length.out = n.seg + 1)[-(n.seg + 1)]
    deg <- pi / 180
    
    key.x <- len * cos(angles) + key.loc[1L]
    key.y <- len * sin(angles) + key.loc[2L]
    lab.angl <- angles
    
    label.x <- 1 * (len + adj) * cos(lab.angl) + key.loc[1L]
    label.y <- 1 * (len + adj) * sin(lab.angl) + key.loc[2L]
    for (k in 1L:n.seg) {
      text.adj <- c(
        if (lab.angl[k] < 90 * deg || lab.angl[k] >
            270 * deg)
          0
        else if (lab.angl[k] > 90 * deg &&
                 lab.angl[k] < 250 * deg)
          1
        else
          0.5,
        if (lab.angl[k] <= # 250 270
            90 * deg)
          (1 - lab.angl[k] / (90 * deg)) / 2
        else if (lab.angl[k] <=
                 270 * deg)
          (lab.angl[k] - 90 * deg) / (180 * deg)
        else
          1 -
          (lab.angl[k] - 270 * deg) / (180 * deg)
      )
      graphics::text(
        label.x[k],
        label.y[k],
        labels = key.labels[k],
        cex = cex,
        adj = text.adj
      )
    }
  }
  
  # internal to extract boxplots statistics
  stars.stats <- function (x, coef, p)
  {
    do.out <- TRUE
    if (coef < 0)
      stop("'coef' must not be negative")
    nna <- !is.na(x)
    n <- sum(nna)
    stats <-
      c(stats::fivenum(x, na.rm = TRUE),
        stats::quantile(x, p, na.rm = TRUE))
    iqr <- diff(stats[c(2, 4)])
    if (coef == 0)
      do.out <- FALSE
    else {
      out <- if (!is.na(iqr)) {
        x < (stats[2L] - coef * iqr) | x > (stats[4L] + coef *
                                              iqr)
      }
      else!is.finite(x)
      if (any(out[nna], na.rm = TRUE))
        stats[c(1, 5, 6)] <- c(range(x[!out], na.rm = TRUE),
                               stats::quantile(x[!out], p, na.rm = TRUE))
    }
    list(stats = stats,
         out = if (do.out)
           out & nna
         else
           logical(length(x)))
  }
  
  # quasi-log10
  qlog10 <- function(x)
  {
    asinh(x / 2) / log(10)
  }
  
  # asinh with cofactor, as used for cytometry data
  asinh.co <- function(x)
  {
    asinh(x * cofactor)
  }
  
  # no transform, i.e. linear scale
  lin <- identity
  
  if (!missing(NAUT)) {
    FUN <- NAUT$FUN
    SCALES <- NAUT$SCALES
    ORD <- NAUT$ORD
    SHIFT <- NAUT$SHIFT
    main.ticks <- NAUT$main.ticks
    minor.ticks <- NAUT$minor.ticks
    log.ticks <- NAUT$log.ticks
    
    if (!all(colnames(x)[ORD] == colnames(SCALES)[1:(ncol(SCALES) - 1)])) {
      stop("variables don't match.")
    }
    
  }
  
  re <- ifelse(switch.cols, identity, rev)
  
  # transform to be used
  f <- get(FUN)
  
  mat <- x
  
  # data should be randomized to prevent overplotting
  if (is.numeric(seed)) {
    set.seed(seed)
    ran <- sample.int(nrow(mat))
    col <- col[ran]
    mat <- mat[ran, ]
  }
  
  # permit log
  if (FUN == "log10") {
    mat[mat < 0] <- NA
  }
  
  mat.stats <- apply(mat, 2, stars.stats, coef = coef, p = p)
  
  mat[sapply(mat.stats, '[[', 2)] <- NA # remove outliers
  # remove NA
  if (any(is.na(mat))) {
    mat <- na.omit(mat)
    col <- col[-attr(mat, "na.action")]
  }
  
  if (missing(NAUT)) {
    # get specified quantile, e.g. p = 0.5 for medians or p = 1 for maxima
    colmax <- f(sapply(mat.stats, '[[', 1)[6, ])
    
    # reorder columns according to specified quantile
    ORD <- order(colmax, decreasing = TRUE)
    
    if (untie > 2 & untie < ncol(mat)) {
      n <- ncol(mat)
      ORD <-
        unlist(unname(sapply(split(ORD, c(
          rep(1:floor(n / untie), each = untie), rep(floor(n / untie), n %% untie)
        )),
        function(o) {
          o[re(hclustfun(distfun(t(f(
            mat[, o]
          ))))$order)]
        }, simplify = FALSE)))
    }
    
    tmp.colmax <- colmax[ORD]
    
    # main and minor scale
    if (FUN != "lin") {
      if (FUN == "log10") {
        MIN <- floor(log10(min(mat, na.rm = TRUE)))
        MAX <- ceiling(log10(max(mat, na.rm = TRUE)))
        log.ticks <- MIN:MAX
        main.ticks <- 10 ^ log.ticks
        minor.ticks <-
          unlist(sapply(log.ticks[-1], function(x)
            c(2:9) * 10 ^ (x - 1), simplify = FALSE))
      } else {
        suppressWarnings(MIN <-
                           min(-ceiling(log10(
                             max(-mat[mat < -1], na.rm = TRUE)
                           )), floor(log10(
                             min(mat[mat > 1], na.rm = TRUE)
                           )), na.rm = TRUE))
        suppressWarnings(MAX <-
                           max(-floor(log10(
                             min(-mat[mat < -1], na.rm = TRUE)
                           )), ceiling(log10(
                             max(mat[mat > 1], na.rm = TRUE)
                           )), na.rm = TRUE))
        log.ticks <- MIN:MAX
        log.ticks <- log.ticks[log.ticks != 0]
        near.zero <- mat > -10 & mat < 10
        if (any(near.zero)) {
          lin.ticks <-
            floor(min(mat[near.zero], na.rm = TRUE)):ceiling(max(mat[near.zero], na.rm =
                                                                   TRUE))
          lin.ticks <- lin.ticks[lin.ticks %in% -1:1]
        } else {
          lin.ticks <- NULL
        }
        main.ticks <-
          c(-10 ^ -log.ticks[log.ticks < 0], lin.ticks, 10 ^ log.ticks[log.ticks > 0])
        minor.ticks <- sort(c(unlist(
          sapply((log.ticks[log.ticks < 0]), function(x)
            sign(x) * c(2:9) * 10 ^ (abs(x) - 1), simplify = FALSE)
        ),
        lin.ticks / 2,
        unlist(
          sapply(log.ticks[log.ticks > 0], function(x)
            sign(x) * c(2:9) * 10 ^ (abs(x) - 1), simplify = FALSE)
        )))
        minor.ticks <- 	minor.ticks[minor.ticks > min(main.ticks)]
      }
      SHIFT <- -f(min(main.ticks))
      
      m <-
        matrix(NA,
               nrow = (1 + length(main.ticks) + length(minor.ticks)),
               ncol = length(ORD))
      colnames(m) <- colnames(mat)[ORD]
      SCALES <- data.frame(m, ' ' = 0, check.names = FALSE)
      tmp.outline <-
        seq(f(max(main.ticks)), tmp.colmax[length(tmp.colmax)], length.out = ncol(SCALES) -
              1)
      SCALES[1, ] <-
        c(tmp.outline + max(tmp.colmax - tmp.outline), 0) # outer shell
      SCALES[1 + seq_along(main.ticks), ] <- f(main.ticks) # main scale
      SCALES[1 + length(main.ticks) + seq_along(minor.ticks), ] <-
        f(minor.ticks) # minor scale
      
      
    } else {
      main.ticks <- pretty(range(mat, na.rm = TRUE))
      minor.ticks <- log.ticks <- NA
      MAX <- max(main.ticks)
      MIN <- min(main.ticks)
      SHIFT <- -MIN
      
      m <- matrix(NA, nrow = (1 + length(main.ticks)), ncol = length(ORD))
      colnames(m) <- colnames(mat)[ORD]
      SCALES <- data.frame(m, ' ' = 0, check.names = FALSE)
      tmp.outline <-
        seq(MAX, tmp.colmax[length(tmp.colmax)], length.out = ncol(SCALES) - 1)
      SCALES[1, ] <-
        c(tmp.outline + max(tmp.colmax - tmp.outline), 0) # outer shell
      SCALES[1 + seq_along(main.ticks), ] <-
        f(main.ticks) # main scale only
    }
    
    SCALES <- SCALES + SHIFT
    SCALES[1, ncol(SCALES)] <- 0
    # this creates the actual shape for main and minor scale lines
    for (c in 1:ncol(SCALES)) {
      SCALES[SCALES[, c] > SCALES[1, c], c] <- 0
    }
    
  }
  
  # initialize plot
  graphics::stars(
    SCALES[1, ],
    locations = c(0, 0),
    scale = FALSE,
    radius = FALSE,
    key.loc = NA,
    main = "",
    col.lines = col.shell,
    lwd = 0.005,
    axes = F,
    len = 1.0
  )
  
  if (plot.minor & FUN != "lin") {
    # minor scale
    lwd.tmp <-
      seq(lwd.minor.l, lwd.minor.g, length.out = length(minor.ticks))
    for (i in seq_along(minor.ticks)) {
      graphics::stars(
        SCALES[i + 1 + length(main.ticks), ],
        locations = c(0, 0),
        scale = FALSE,
        radius = FALSE,
        key.loc = NA,
        main = "",
        col.lines = col.minor,
        lwd = lwd.tmp[i],
        axes = F,
        len = 1.0,
        add = T
      )
    }
  }
  if (plot.major) {
    # major scale
    lwd.tmp <-
      seq(lwd.major.l, lwd.major.g, length.out = length(main.ticks))
    for (i in seq_along(main.ticks)) {
      graphics::stars(
        SCALES[i + 1, ],
        locations = c(0, 0),
        scale = FALSE,
        radius = FALSE,
        key.loc = NA,
        main = "",
        col.lines = col.major,
        lwd = lwd.tmp[i],
        axes = F,
        len = 1.0,
        add = T
      )
    }
  }
  if (plot.vertical) {
    # vertical
    graphics::stars(
      SCALES[1, ],
      locations = c(0, 0),
      scale = FALSE,
      radius = TRUE,
      key.loc = NA,
      main = "",
      col.lines = col.vertical,
      lwd = lwd.vertical,
      axes = F,
      len = 1.0,
      add = T
    )
  }
  # final outline
  graphics::stars(
    SCALES[1, ],
    locations = c(0, 0),
    scale = FALSE,
    radius = FALSE,
    key.loc = NA,
    main = "",
    col.lines = col.shell,
    lwd = lwd.shell,
    axes = F,
    len = 1.0,
    add = T
  )
  
  if (plot.axis) {
    # main scale
    if (FUN != "lin") {
      if (FUN == "log10") {
        graphics::axis(
          1,
          pos = c(0, 0),
          at = f(main.ticks[-1]) + SHIFT,
          labels = parse(text = c(
            paste("10", "^", as.character(log.ticks[-1]), sep = "")
          )),
          lwd = lwd.major.axis,
          lwd.ticks = lwd.major.ticks,
          col = col.shell,
          cex.axis = cex.axis,
          padj = padj.axis,
          tck = tck.major
        )
      } else {
        graphics::axis(
          1,
          pos = c(0, 0),
          at = f(main.ticks[-1]) + SHIFT,
          # labels=parse(text=ifelse(log10(abs(main.ticks[-1]))>0,
          # c(paste("10", "^",as.character(log.ticks[-1]), sep="")),
          # main.ticks[-1])),
          labels = main.ticks[-1],
          lwd = lwd.major.axis,
          lwd.ticks = lwd.major.ticks,
          col = col.shell,
          cex.axis = cex.axis,
          padj = padj.axis,
          tck = tck.major
        )
      }
      # minor scale
      graphics::axis(
        1,
        pos = c(0, 0),
        at = f(c(minor.ticks)) + SHIFT,
        labels = FALSE,
        lwd = lwd.minor.axis,
        lwd.ticks = lwd.minor.ticks,
        col = col.shell,
        tck = tck.minor
      )
      
    } else {
      # main linear scale
      graphics::axis(
        1,
        pos = c(0, 0),
        at = f(main.ticks[-1]) + SHIFT,
        labels = main.ticks[-1],
        lwd = lwd.major.axis,
        lwd.ticks = lwd.major.ticks,
        col = col.shell,
        cex.axis = cex.axis,
        padj = padj.axis,
        tck = tck.major
      )
    }
  }
  
  # do transform
  graphics::stars(
    data.frame(f(mat[, ORD]) + SHIFT, ' ' = 0, check.names = FALSE),
    locations = c(0, 0),
    scale = FALSE,
    radius = F,
    key.loc = NA,
    col.lines = col,
    lwd = lwd,
    add = T
  )
  
  # add labels
  stars.label(
    x = SCALES[1, ],
    locations = c(0, 0),
    len = as.numeric(SCALES[1, ]),
    key.loc = c(0, 0),
    cex = cex,
    adj = adj
  )
  
  message("plotted ", nrow(mat), " out of ", nrow(x), " input rows")
  
  invisible(stats::setNames(
    list(FUN, SCALES, ORD, SHIFT, main.ticks, minor.ticks, log.ticks),
    c(
      "FUN",
      "SCALES",
      "ORD",
      "SHIFT",
      "main.ticks",
      "minor.ticks",
      "log.ticks"
    )
  ))
  
}
