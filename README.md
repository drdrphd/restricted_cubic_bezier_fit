# restricted_cubic_bezier_fit
R code for performing restricticted cubic Bezier (Hermite spline) regression
Based on Aaron Olsen's [bezier](https://cran.r-project.org/web/packages/bezier/index.html) package (with some code adapted from his)

This is essentially a cubic Hermite spline regression, where the magnitude of the slope vectors is not known. Will likely update code names to reflect this in the future.

## Usage
`restrictedCubicBezierFit <- function(m, ends_from_data = FALSE,
	P0 = NULL, P3 = NULL, slope_start, slope_end,
	na_fill = FALSE, max_iter = 20) `

## Arguments
| Argument | Description |
| ----------- | ----------- |
| __m__ | matrix of datapoints (can be in multiple dimensions -- each column is a dimension). Will use first and last rows as control points P0 and P3 (start and end points). |
| __ends_from_data__ | (Boolean) (if true) Get starting and ending points from data matrix |
| __P0__, __P3__ | (if ends_from_data not specified) Starting and ending points, each as n-dim (x,y,...) vectors. |
| __slope_start__ | Slope vector from starting point |
| __slope_end__ | Slope vector from ending point |
| __na_fill__ | (Boolean) Interpolate missing values in data matrix |
| __max_iter__ | Maximum number of iterations to use when fitting |

## Returns
A fit of a cubic Bezier curve to the data, including:

optim() fit information:
  * fit$par = alpha,beta parameters (magnitude of (P0,P1) and (P2,P3) segments)
  * fit$value = Residual Sum of Squares
  * fit$Points = set of control points for best fit cubic Bezier
  * and information on convergence
	
regression information:
  * fit$residuals
  * fit$rss = Residual Sum of Squares
  * fit$rmse = Root Mean Square Error
  * fit$r2 = R²
  * fit$r2_adj = Adjusted R² -- for two free parameters: magnitudes of both slope vectors

## Details
Splines are used very frequently in computer graphics and for producing interpolations of datasets. Cubic splines are particularly useful in some applications, as they form an interpolation between points. In the case of Cubic Bezier splines where we restrict the slope from (P0,P1) and from (P2,P3) -- equivalent to an Hermite cubic spline from (P0,P1) with slope vectors (m0,m1) where we don't know the magnitude of the slope vector -- this can be used to represent an interpolation between two functions. Cubic Hermite splines are defined on the interval from \[0,1], which can be applied in several dimensions.

For example, when dealing with data representing a proportional switch between categories predicted by another proportion (say: [average perception of a stimulus in noise](https://www.jneurosci.org/content/jneuro/31/17/6339/F8.large.jpg)), you have a switch between y=0 and y=1. The strength of the attraction to those underlying functions can be reasonably represented by the magnitude of the slope vector in the Hermite spline (or the distance from (P0,P1), (P2,P3)). This is particularly relevant for interpolations between functions on different slopes, for example, data that interpolates between y=0 and y=x.

Using this interpolation as a basis function for regression is interpretable as a switching behavior between two wieghted linear functions. Regression fits on test datasets reveal comparable or better RMSE when compared to other classic regression methods (polynomial, logistic, beta), and may be superior in situations where you suspect the underlying switching behavior between known functions and cases where the end points are strictly specified (where polynomial predictions may be off the mark) or beta regressions which are only defined on (0,1) and do not make good predictions without extended methods like zero-one inflated beta regression.

A note on the parameters: the slopes must be given in vector terms in the direction of the spline (that is velocity as t traverses \[0,1]). If your starting slope is m0=\[1,1] and the spline is entirely linear, then the ending slope is m1=\[1,1], not \[-1,-1]. (That is, don't point the end vector in towards the spline.)
