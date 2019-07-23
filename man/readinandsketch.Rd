\name{readinandsketch}
\alias{readinandsketch}
\title{
Create a sketch from a file containing a (very large) data set
}
\description{
This function calculates a sketch of a file. The sketch can be used to
perform approximate frequentist or Bayesian linear regression. The
sketch is a substitute data set of the same dimension but much smaller
number of observations. The analysis based on the sketch is much faster
and its results are provably close to the results on the original data
set. The file is read in sequentially, making it possible to sketch data
sets that are too large to be loaded into R completely. 
}
\usage{
readinandsketch(file, nrows = 50000, epsilon = NULL, obs_sketch = NULL,
                affine = TRUE, method= c("C", "S", "R"), header = FALSE,
                sep = "", col.names, skip = 0, warn = FALSE, ...)
}
\arguments{
  \item{file}{
  	The name of a file that contains the (large) data set. The data set should consist of both the design matrix X and the vector Y, which contains the values of the dependent variable. The order is arbitrary.
  	}
  \item{nrows}{
  	A positive integer, which controls, how many rows are read into the memory per iteration. Differs from use in \code{\link{read.table}} as the other rows will be read in subsequent iterations. For that reason, \code{\var{nrows}} has to be larger than 0.
	}
  \item{epsilon}{
	Approximation error of the sketch (see Details). Only one of \code{\var{epsilon}} and \code{\var{obs_sketch}} can be used, if both are specified, currently \code{\var{epsilon}} is used and \code{\var{obs_sketch}} is ignored. Possible values for epsilon lie in the interval (0, 0.5].
	}
  \item{obs_sketch}{
  	Desired number of observations of the sketch (see Details). Only one of \code{\var{epsilon}} and \code{\var{obs_sketch}} can be used, if both are specified, currently \code{\var{epsilon}} is used and \code{\var{obs_sketch}} is ignored.
  }
  \item{affine}{
  	Boolean, choose TRUE if your model includes an intercept term and your data set does not contain a corresponding column. The corresponding column will be added as new left-most column of the sketch. If you do not want an added intercept term, choose FALSE.
   } 
  \item{method}{
  	The sketching method to be used. Possible values are "R", "S", and "C". See Details.
   } 
  \item{header}{
  	Boolean, if TRUE, the first line of the file is used as variable names, see \code{\link{read.table}}.
  }
  \item{sep}{
        The field separator character, see \code{\link{read.table}}.
  }
  \item{col.names}{
  	An optional vector containing the variable names, see \code{\link{read.table}}.
  }
  \item{skip}{
  	integer: the number of lines of the data file to skip before beginning to read data, see \code{\link{read.table}}.
  }
  \item{warn}{
  	Boolean, if TRUE show a warning if the sketch will result in a matrix of larger dimension than the original matrix.
   }
  \item{...}{
  	Additional arguments that will be passed on to \code{\link{read.table}}.
  }
}
\details{
This function reads a data set iteratively and calculates/updates a sketch of the read in data set. This sketch can then be used for frequentist or Bayesian linear regression, especially on large data sets. The functionality used here is the same as in \code{\link{sketch}}, but \code{\var{readinandsketch}} can also handle data sets that are too large to be loaded into the working memory.

In principle, \code{\var{nrows}} can be any positive integer value. If using the methods "R" or "C", small integer values will only lead to an increased running time. If using method "S", however, \code{\var{nrows}} has to be at least as large as the number of observations \eqn{k} in the sketch, otherwise there will be an error.

If the number of observations in the data set is a multiple of \code{\var{nrows}}, there will be one additional empty run, where no data is read and a sketch of an empty data set is calculated. This does not influence the resulting sketch.
}
\value{
  Returns a data frame, which contains both the sketched data frame SX and the sketched vector SY. The order of the columns is the same as in the original data set. If affine is TRUE, the corresponding intercept column is added as the new left-most column of the sketch. Please omit the standard intercept term from any models based on sketches in that case.
}
\references{
Geppert, L., Ickstadt, K., Munteanu, A., Quedenfeld, J., & Sohler, C. (2017). Random projections for Bayesian regression. Statistics and Computing, 27(1), 79-101. doi:10.1007/s11222-015-9608-z
}
\author{
LN Geppert, K. Ickstadt, A. Munteanu, J. Quedenfeld, L. Sandig, C. Sohler
}

\seealso{
\code{\link{sketch}}
}
\examples{
\dontrun{
  sketchC = readinandsketch(file.choose(), header = TRUE, sep = '\t',
  nrows = 10000, epsilon = 0.1, method = 'R')
}
}
