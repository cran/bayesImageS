% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getNeighbors.R
\name{getNeighbors}
\alias{getNeighbors}
\title{Get Neighbours of All Vertices of a Graph}
\usage{
getNeighbors(mask, neiStruc)
}
\arguments{
\item{mask}{a vector, matrix, or 3D array specifying vertices within a graph. Vertices of value 1 are within the graph and 0 are not.}

\item{neiStruc}{a scalar, vector of four components, or \eqn{3\times4} matrix
corresponding to 1D, 2D, or 3D graphs. It gives the definition of
neighbours of a graph.
All components of \code{neiStruc} should be positive (\eqn{\ge 0})
even numbers.
For 1D graphs, \code{neiStruc} gives
the number of neighbours of each vertex.
For 2D graphs, \code{neiStruc}[1] specifies
the number of neighbours on vertical direction, \code{neiStruc}[2]
horizontal direction, \code{neiStruc}[3] north-west (NW) to south-east (SE)
diagonal direction, and \code{neiStruc}[4] south-west (SW) to
north-east (NE) diagonal direction. For 3D
graphs, the first row of \code{neiStruc} specifies the number of neighbours on
vertical direction, horizontal direction and two diagonal directions from
the 1-2 perspective, the second row the 1-3 perspective, and the
third row the 2-3 perspective. The index to
perspectives is represented with the leftmost subscript of the
array being the smallest.}
}
\value{
A matrix with each row giving the neighbours of a vertex.
The number of the rows is equal to the number of 
vertices within the graph and the number or columns
is the number of neighbours of each vertex.

For a 1D graph, if each vertex has two neighbours,
The first column are the neighbours on the left-hand side of
corresponding vertices and the second column the right-hand side.
For the vertices on boundaries, missing neighbours are represented by
the number of vertices within a graph plus 1.
When \code{neiStruc} is bigger than 2, The first two columns
are the same as when \code{neiStruc} is equal to 2; the third column
are the neighbours on the left-hand side of the vertices on the first column;
the forth column are the neighbours on the right-hand side of the vertices
on the second column, and so on and so forth. And again for the
vertices on boundaries, their missing neighbours are represented by
the number of vertices within a graph plus 1.

For a 2D graph, the index to vertices is column-wised. For each
vertex, the order of neighbours are as follows. First are those
on the vertical direction, second the horizontal
direction, third the NW to SE diagonal
direction, and forth the SW to NE diagonal
direction. For each direction, the neighbours of every vertex
are arranged in the same way as in a 1D graph.

For a 3D graph, the index to vertices is that
the leftmost subscript of the array moves the fastest.  For each
vertex, the neighbours from the 1-2 perspective
appear first and then the 1-3 perspective and finally
the 2-3 perspective. For each perspective, the neighbours are arranged
in the same way as in a 2D graph.
}
\description{
Obtain neighbours of vertices of a 1D, 2D, or 3D graph.
}
\details{
There could be more than one way to define the same 3D neighbourhood
structure for a graph (see Example 3 for illustration).
}
\examples{
  #Example 1: get all neighbours of a 1D graph.
  mask <- c(0,0,rep(1,4),0,1,1,0,0,1,1,1)
  getNeighbors(mask, neiStruc=2)
  
  #Example 2: get all neighbours of a 2D graph based on neighbourhood structure
  #           corresponding to the second-order Markov random field.
  mask <- matrix(1, nrow=2, ncol=3)
  getNeighbors(mask, neiStruc=c(2,2,2,2))
  
  #Example 3: get all neighbours of a 3D graph based on 6 neighbours structure
  #           where the neighbours of a vertex comprise its available
  #           N,S,E,W, upper and lower adjacencies. To achieve it, there
  #           are several ways, including the two below.
  mask <- array(1, dim=rep(3,3))
  n61 <- matrix(c(2,2,0,0,
                  0,2,0,0,
                  0,0,0,0), nrow=3, byrow=TRUE)
  n62 <- matrix(c(2,0,0,0,
                  0,2,0,0,
                  2,0,0,0), nrow=3, byrow=TRUE)
  n1 <- getNeighbors(mask, neiStruc=n61)
  n2 <- getNeighbors(mask, neiStruc=n62)
  n1 <- apply(n1, 1, sort)
  n2 <- apply(n2, 1, sort)
  all(n1==n2)
  
  #Example 4: get all neighbours of a 3D graph based on 18 neighbours structure
  #           where the neighbours of a vertex comprise its available
  #           adjacencies sharing the same edges or faces.
  #           To achieve it, there are several ways, including the one below.
  
  n18 <- matrix(c(2,2,2,2,
                  0,2,2,2,
                  0,0,2,2), nrow=3, byrow=TRUE)  
  mask <- array(1, dim=rep(3,3))
  getNeighbors(mask, neiStruc=n18)
}
\references{
Winkler, G. (2003)
 "Image Analysis, Random Fields and Markov Chain Monte Carlo Methods: A Mathematical Introduction" (2nd ed.)
 \emph{Springer-Verlag}
 
 Feng, D. (2008)
 "Bayesian Hidden Markov Normal Mixture Models with Application to MRI Tissue Classification"
 \emph{Ph. D. Dissertation, The University of Iowa}
}
\keyword{spatial}
