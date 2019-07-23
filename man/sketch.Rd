\name{sketch}
\alias{sketch}
\title{
Create a sketch from a matrix, data.frame or file
}
\description{
This function calculates a sketch of a data set, matrix or file. The
sketch can be used to perform approximate frequentist or Bayesian linear
regression. The sketch is a substitute data set of the same dimension
but much smaller number of observations. The analysis based on the
sketch is much faster and its results are provably close to the results
on the original data set. 
}
\usage{
sketch(data, file, epsilon = NULL, obs_sketch = NULL,
       warn = TRUE, affine = TRUE, method= c("R", "S", "C"), ...)}
\arguments{
  \item{data}{
  	Matrix or data.frame which contains both the design matrix X and the vector Y used for the linear regression, in any order. Either \code{\var{data}} or \code{\var{file}} has to be provided.
	}
  \item{file}{
  	The name of a file which contains the matrix or data.frame. Ignored, if \code{\var{data}} is specified. As in \code{\var{data}}, this file should contain both the design matrix X and the vector Y used in the linear regression, in any order. Either \code{\var{file}} or \code{\var{data}} has to be provided. If your data set is very large, consider using \code{\link{readinandsketch}} for more efficient sketching of the data.
  }
  \item{epsilon}{
	Approximation error of the sketch (see Details). Only one of \code{\var{epsilon}} and \code{\var{obs_sketch}} can be used, if both are specified, currently \code{\var{epsilon}} is used and \code{\var{obs_sketch}} is ignored. Possible values for epsilon lie in the interval (0, 0.5].
	}
  \item{obs_sketch}{
  	Desired number of observations of the sketch (see Details). Only one of \code{\var{epsilon}} and \code{\var{obs_sketch}} can be used, if both are specified, currently \code{\var{epsilon}} is used and \code{\var{obs_sketch}} is ignored. If method "C" is chosen, the number is rounded up to the next power of two, see Details.
  }
  \item{warn}{
  	Boolean, if TRUE show a warning if the sketch will result in a matrix of larger dimension than the original matrix. Please note that a sketch with larger dimension than the original matrix will result in an error if method "S" is used.
   }
  \item{affine}{
  	Boolean, choose TRUE if your model includes an intercept term and your data set does not contain a corresponding column. The corresponding column will be added as new left-most column of the sketch. If you do not want an added intercept term, choose FALSE.
   } 
  \item{method}{
  	The sketching method to be used. Possible values are "R", "S", and "C", where "R" is the default. See Details.
   } 
  \item{...}{
    Additional arguments passed on to \code{\link{read.table}} if \code{\var{file}} is specified.
   }
}
\details{
Let X be a \eqn{(n \times d)}{n x d}-matrix and Y a \eqn{(n \times 1)}{n x 1}-vector. This function produces an implicit matrix S and efficiently performs the multiplication, which embeds X and Y into a lower dimension \eqn{k}, with \eqn{k \ll n}{k << n}. The value of k depends on the method used. For "R" and "S", the formula is
\eqn{k = \left \lceil{\frac{d\cdot \ln(d)}{\varepsilon^2}}\right\rceil}{k = ceil(d * ln(d) / eps^2)}, for "C", this changes to \eqn{k = 2^{\lceil\log_2 s \rceil}}{k = 2^ceil(log_2 s)}, where 
\eqn{s = \left \lceil{\frac{d^2}{20\cdot \varepsilon^2}}\right\rceil}{s = ceil(d * d / (20*eps^2))}, that is the smallest power of two larger than s.

The function outputs the sketched data frame and vector SX and SY. These can be used to conduct frequentist or Bayesian linear regression on a smaller data set, saving running time and memory. The results are guaranteed to be close to the results that would have been obtained on the original data set in the sense that the original likelihood is closely approximated by the likelihood on the sketched data set in the case of classical linear regression, or, in the Bayesian case the original posterior is closely approximated by the posterior on the sketched data set.

When using methods "R" and "C", it is possible for the sketch to be of a larger dimension than the original matrix. When using method "S", this will result in an error. Such cases occur when the number of variables is relatively large compared to the number of observations.

For more details, please refer to Geppert et al. (2017) and the references mentioned therein.
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
\code{\link{readinandsketch}}
}
\examples{
# create a small simulated data set
# with 400 observations and
# 4 variables
set.seed(23)
x1 = rnorm(400, 10, 2)
x2 = rnorm(400, 5, 3)
x3 = rnorm(400, -2, 1)
x4 = rnorm(400, 0, 5)
y = 2.4 - 0.6 * x1 + 5.5 * x2 - 7.2 * x3 + 5.7 * x4 + rnorm(400)
# all in one data.frame
data = data.frame(x1, x2, x3, x4, y)

# linear model based on original data set
lm(y ~ ., data = data)

# Calculate a RAD/"R"-sketch with epsilon = 0.2
s1 = sketch(data, epsilon = 0.2, method = 'R', affine = TRUE)
dim(s1)
# very similar results, intercept should be omitted
lm(y ~ . - 1, data = s1)

# use option "obs_sketch" to fix the new number of observations
s2 = sketch(data, obs_sketch = 200, method = 'R', affine = TRUE)
dim(s2)
# some more differences as sketch is smaller
lm(y ~ . - 1, data = s2)

# calculate SRHT/"S"-sketch
s3 = sketch(data, epsilon = 0.2, method = 'S', affine = TRUE)
dim(s3)
lm(y ~ . - 1, data = s3)

# calculate CW/"C"-sketch
s4 = sketch(data, epsilon = 0.2, method = 'C', affine = TRUE)
dim(s4)
# sketch is smaller, because the number of variables is very small
# CW-sketches require a lot more observations compared to RAD/SRHT
# when number of variables increases
lm(y ~ . - 1, data = s4)

# same simulated data set, but with intercept added to data.frame
data2 = data.frame(x0 = 1, x1, x2, x3, x4, y)
lm(y ~ . - 1, data = data2)

# Same as s1, but now option affine = FALSE is adequate
s5 = sketch(data2, epsilon = 0.2, method = 'R', affine = FALSE)
dim(s5)
lm(y ~ . - 1, data = s5)
}
