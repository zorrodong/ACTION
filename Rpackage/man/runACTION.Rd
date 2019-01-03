\name{runACTION}
\alias{runACTION}
\title{Main decomposition function in the ACTION framework}
\usage{
runACTION(S_r, min_k, max_k, numThreads)
}
\arguments{
  \item{S_r}{reduced single-cell profile of dimension D x cells (D is the dimension of latent subspace)}

  \item{min_k}{minimum number of archetypes to consider}
  
  \item{max_k}{maximum number of archetypes to consider}
  
  \item{numThreads}{number of threads (for multi-core processing)}
}
\value{
A list of C and H matrices for different values of k:
  C: A list of contribution matrices. Archetypes are computed as S_r \%*\% C.
  H: A list where element H[[i]] is the footprint of archetypes at level i on all cells.
}

}
\description{
It first runs convex NMF to identify a set of promissing archetypal candidates. Then, it runs archetypal analysis starting from the result of convex NMF. Repeats this process for k = k_min, ..., k_max.
}
\examples{
ACTION.results = runACTION(S_r, 2, 20, 8)
C_trace = ACTION.results$C
H_trace = ACTION.results$H
}