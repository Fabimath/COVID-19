#######################################################################################
##                             B-splines basis functions                             ##          
#######################################################################################

## Arguments
## x = vector of observations
## xl, xr = min(x) and max(x)
## ndx = number of equally spaced knots
## deg = B-spline degree
## max_derivs = maximum order of derivative to compute
## sparse = logical if a sparse representation is required

## Output
## A list with each element refering to derivative of order 0,1,..., max_derivs


basis_array = function(x, xl = min(x), xr = max(x), ndx = 10, deg = 3, max_derivs = 0, sparse = FALSE)
{
  dx = (xr - xl) / ndx
  knots = seq(xl - deg * dx, xr + deg * dx, by = dx)

  B <- list()
  for (i in 1:(max_derivs + 1))
  {
     B[[i]] = splineDesign(knots = knots,
                              x = x,
                              ord = deg + 1,
                              derivs = rep(i - 1, length(x)),
                              outer.ok = TRUE,
                              sparse = sparse)
  }

  B
}
