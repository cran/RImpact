################################################################################
Metrics <- function(citation.counts, publishing.age = 0, display = TRUE) {
  # Measures scholarly impact using modern citation-based indices.
  #
  # Args:
  #  citation.counts: Number of times each aritcle has been cited. (vector)
  #  publishing.age : Age of the first article author has published. (scalar)
  #  display        : Whether to display metrics (if TRUE, the default) or direct
  #                   output to a file (if FALSE).
  #
  # Returns:
  #   h.index            : h index, the largest number h such that at least h
  #                        articles are cited h times each (Hirsch, 2005).
  #   tapered.h.index    : Tapered h index, credit decreases for citations farther
  #                        from the origin (Anderson, Hankin, & Killworth, 2008).
  #   f.index            : f index, largest value f such that the harmonic mean
  #                        for the f most highly cited articles is at least f
  #                        (Tol, 2009).
  #   t.index            : t index, largest value t such that the geometric mean
  #                        for the t most highly cited articles is at least t
  #                        (Tol, 2009).
  #   g.index            : g index, largest value g such that the mean citations
  #                        for the g most highly cited articles is at least g
  #                        (Egghe, 2006).
  #   hg.index           : hg index, geometric mean of h and g (Alonso,
  #                        Cabrerizo, Herrera-Viedma, & Herrera, 2010).
  #   a.index            : a index, mean citations for papers in Hirsch core
  #                        (Jin, 2006).
  #   m.index            : m index, median citations for papers in Hirsch core
  #                        (Bornmann, Mutz, & Daniel, 2008).
  #   r.index            : r index, square root of citations for papers in Hirsch
  #                        core (Jin, Liang, Rousseau, & Egghe, 2007).
  #   weighted.h.index   : h index weighted by citation impact (Egghe & Rousseau,
  #                        2008).
  #   q2.index           : q2 index, geometric mean of h and m indexes (Cabrerizo,
  #                        Alonso, Herrera-Viedma, & Herrera, 2010).
  #   e.index            : e index, excess citations for papers in Hirsch core
  #                        (Zhang, 2009).
  #   max.product        : Maximum product index, maximum product of article's
  #                        rank and citation count (Kosmulski, 2007).
  #   sqrt.max.product   : Rescales maximum product index from an area to a
  #                        distance measure.
  #   h2.index           : h2 index, analogous to h index with more stringent
  #                        criterion (Kosmulski, 2006).
  #   m.quotient         : m quotient, controlling h index for publishing age
  #                        (Hirsch, 2005).
  #   tapered.m.quotient : Controlling tapered h index for publishing age.
  #
  citation.counts <- sort(citation.counts, decreasing = TRUE)
  n <- length(citation.counts)
  c.total <- sum(citation.counts)
  sqrt.c <- sqrt(c.total)
  m <- c.total / n
  mdn <- median(citation.counts)
  h.index <- h2.index <- g.index <- r0 <- f.index <- t.index <- 1
  tapered.h.index <- 0
  if (c.total > 0) {
    test <- TRUE
    while ((test) & (h.index < n))
      if (citation.counts[h.index + 1] >= (h.index + 1)) {
        h.index <- h.index + 1
      } else {
        test <- FALSE
      }
    for (i in 1:n)
      if (citation.counts[i] > 0)
        for (j in 1:citation.counts[i])
          tapered.h.index <- tapered.h.index + 1 / (2 * max(i, j) - 1)
        test <- TRUE
        while ((test) & (f.index < n))
          if (HarmonicMean(citation.counts[1:(f.index + 1)]) >= (f.index + 1)) {
            f.index <- f.index + 1
          } else {
            test <- FALSE
          }
        test <- TRUE
        while ((test) & (t.index < n))
          if (GeometricMean(citation.counts[1:(t.index + 1)]) >= (t.index + 1)) {
            t.index <- t.index + 1
          } else {
            test <- FALSE
          }
        test <- TRUE
        while ((test) & (g.index < n))
          if (sum(citation.counts[1:(g.index + 1)]) >= (g.index + 1) ^ 2) {
            g.index <- g.index + 1
          } else {
            test <- FALSE
          }
        hg.index <- sqrt(h.index * g.index)
        h.sum <- sum(citation.counts[1:h.index])
        a.index <- h.sum / h.index
        m.index <- median(citation.counts[1:h.index])
        r.index <- sqrt(h.sum)
        test <- TRUE
        while ((test) & (r0 < n))
          if ((sum(citation.counts[1:(r0 + 1)]) / h.index) <=
              citation.counts[r0 + 1]) {
            r0 <- r0 + 1
          } else {
            test <- FALSE
            weighted.h.index <- sqrt(sum(citation.counts[1:r0]))
          }
        q2.index <- sqrt(h.index * m.index)
        e.index <- sqrt(h.sum - h.index ^ 2)
        pr <- citation.counts * 1:n
        max.product <- max(pr)
        sqrt.max.product <- sqrt(max.product)
        test <- TRUE
        while ((test) & (h2.index < n))
          if (citation.counts[h2.index + 1] >= (h2.index + 1) ^ 2) {
            h2.index <- h2.index + 1
          } else {
            test <- FALSE
          }
        if (publishing.age > 0) {
          m.quotient <- h.index / publishing.age
          tapered.m.quotient <- tapered.h.index / publishing.age
        } else {
          m.quotient <- tapered.m.quotient <- 0
        }
  } else {
    h.index <- tapered.h.index <- f.index <- t.index <- g.index <- hg.index <-
      a.index <- m.index <- r.index <- weighted.h.index <- q2.index <-
      e.index <- max.product <- sqrt.max.product <- h2.index <- m.quotient <-
      tapered.m.quotient <- 0
  }
  if (display) {
    cat("\nRank-ordered citation counts:\n")
    cat("Hirsch core: ",sort (citation.counts[1:h.index], decreasing = TRUE),
        "\n")
    cat("Additional: ",sort (citation.counts[h.index + 1:n], decreasing = TRUE),
        "\n")
    cat("\nNumber of articles =", n, "\n")
    cat("Total citations =", c.total, "\n")
    cat("Square root of total citations =", round(sqrt.c, 2), "\n")
    cat("Mean citations =", round(m, 2), "\n")
    cat("Median citations =", mdn, "\n")
    cat("h index =", h.index, "\n")
    cat("tapered h index =", round(tapered.h.index, 2), "\n")
    cat("f index =", round(f.index, 2), "\n")
    cat("t index =", round(t.index, 2), "\n")
    cat("g index =", g.index, "\n")
    cat("hg index =", round(hg.index, 2), "\n")
    cat("a index =", round(a.index, 2), "\n")
    cat("m index =", m.index, "\n")
    cat("r index =", round(r.index, 2), "\n")
    cat("weighted h index =", round(weighted.h.index, 2), "\n")
    cat("q2 index = ", round(q2.index, 2), "\n")
    cat("e index =", round(e.index, 2), "\n")
    cat("maximum product index =", round(max.product), "\n")
    cat("square root of maximum product index =", round(sqrt.max.product, 2), "\n")
    cat("h(2) index =", h2.index, "\n")
    if (publishing.age > 0) {
      cat("m quotient =", round(m.quotient, 2), "\n")
      cat("tapered m quotient =", round(tapered.m.quotient, 2), "\n")
    }
  } else {
    return (c(n, c.total, sqrt.c, m, mdn, h.index, tapered.h.index, f.index,
              t.index, g.index, hg.index, a.index, m.index, r.index,
              weighted.h.index, q2.index, e.index, max.product,
              sqrt.max.product, h2.index, m.quotient, tapered.m.quotient))
  }
}

################################################################################
HarmonicMean <- function(x) {
  # Calculates the harmonic mean (the reciprocal of the arithmetic mean of the
  # reciprocals of n values).
  #
  # Args:
  #   x: Vector of n values whose harmonic mean is to be calculated.
  #
  # Returns:
  #   The harmonic mean of x.
  d <- 0
  for (i in 1:length(x))
    if (x[i] > 0)
      d <- d + 1 / x[i]
    if (d > 0) {
      return(1 / (d / length(x)))
    } else {
      return(0)
    }
}

################################################################################
GeometricMean <- function(x)
  # Calculates the geometric mean (the nth root of the product of n values).
  #
  # x: Vector of n values whose geometric mean is to be calculated.
  #
  # Returns:
  #   The geometric mean of x.
  return(prod(x) ^ (1 / length(x)))



