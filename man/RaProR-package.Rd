\name{RaProR-package} 
\alias{RaProR-package} 
\alias{RaProR} 
\docType{package} 
\title{ 
Create a smaller substitute data set to perform linear regression 
} 
\description{ 
This package can be used to calculate sketches of a data set that can be used to perform approximate classical or Bayesian linear regression. The sketch is a substitute data set of the same dimension but much smaller number of observations. The inference based on the sketch is much faster and is provably close to the exact inference. The calculation is done time- and space-efficiently in C. The two main functions are \code{\link{sketch}} for data sets that fit into the working memory and can be processed at once and \code{\link{readinandsketch}} for data sets that (potentially) do not fit into the working memory and will be read and sketched sequentially blockwise. 
} 
\author{ 
LN Geppert, K. Ickstadt, A. Munteanu, J. Quedenfeld, L. Sandig, C. Sohler 
} 
\references{ 
Geppert, L., Ickstadt, K., Munteanu, A., Quedenfeld, J., Sohler, C. (2017). Random projections for Bayesian regression. Statistics and Computing, 27(1), 79-101. doi:10.1007/s11222-015-9608-z
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

# Calculate an RAD/"R"-sketch with epsilon = 0.2
s1 = sketch(data, epsilon = 0.2, method = 'R', affine = TRUE)
dim(s1)
# very similar results, intercept should be omitted
lm(y ~ . - 1, data = s1)
}
