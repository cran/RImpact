\name{Metrics}
\alias{Metrics}
\title{Metrics
}
\description{Measures scholarly impact using modern citation-based indices.
}
\usage{Metrics(citation.counts, publishing.age = 0, display = TRUE)
}
\arguments{
  \item{citation.counts}{Number of times each aritcle has been cited. (vector)}
  \item{publishing.age}{Age of the first article author has published. (scalar)}
  \item{display}{Whether to display metrics (if TRUE, the default) or direct output to a file (if FALSE).}
}
\value{
  \item{h.index}{h index, the largest number h such that at least h articles are cited h times each (Hirsch, 2005).}
  \item{tapered.h.index}{Tapered h index, credit decreases for citations farther from the origin (Anderson, Hankin, & Killworth, 2008).}
  \item{f.index}{f index, largest value f such that the harmonic mean fo rthe f most highly cited articles is at least f (Tol, 2009).}
  \item{g.index}{g index, larfest value g such that th emean citations for the g most highly cited articles is at least g (Egghe, 2006).}
  \item{hg.index}{hg index, geometric mean of h and g (Alonso, Cabrerizo, Herrera-Viedma, & Herrera, 2010).}
  \item{a.index}{a index, mean citations for the papers in Hirsch core (Jin, 2006).}
  \item{m.index}{m index, median citations for papers in Hirsch core (Bornmann, Mutz, Daniel, 2008).}
  \item{r.index}{r index, square root of citations for papers in Hirsch core (Jin, Liang, Rousseau, Egghe, 2007).}
  \item{weighted.h.index}{h index weighted by citation impact (Egghe & Rousseau, 2008).}
  \item{q2.index}{q2 index, geometric mean of h and m indexes (Cabrerizo, Alonso, Herrera-Videma, & Herrera, 2010).}
  \item{e.index}{e index, excess citations for papers in Hirsch core (Zhang, 2009).}
  \item{max.product}{Maximum product index, maximum product of article's rank and citation count (Kosmulski, 2007).}
  \item{sqrt.max.product}{Rescales maximum product index from an area to a distance measure.}
  \item{h2.index}{h2 index, analogous to h index with more stringent criterion (Kosmulski, 2006).}
  \item{m.quotient}{m quotient, controlling h index for publishing age (Hirsch, 2005).}
  \item{tapered.m.quoteient}{Controlling tapered h index for publishing age.}
}
\references{Ruscio et al. (2012)
}
\author{John Ruscio
}
\examples{
# Running the Metrics program with illustrative data from Ruscio et al. (2012)
x <- c(24, 18, 12, 8, 6, 5, 5, 4, 4, 3, 2, 2, 1, 1, 1, 0, 0, 0)
Metrics(x)
}

