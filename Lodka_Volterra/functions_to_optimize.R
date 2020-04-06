#######################################################################################################
##                        Objective function lynx-hare dyn                                           ##    
## Compute spline coefficients and optimization criterion for ODE parameters in QL-ODE-P-splines     ##
##                        Lotka-Volterra system of equations                                         ##
#######################################################################################################

## Function arguments:
## theta = current estimates ODE parameters
## tau = current estimate for the precision of measurements
## ck = current estimates spline coefficients
## length_t = length grid B-spline for penalty matrix
## seq_t = grid of time points (B-spline  for penalty)
## B_seq_t = list of B-splines terms for the construction of the penalty
## min_t, max_t = min and max values of observed times
## regul = (optional) extra regularization parameter (default = 0)


#---------------------- function to compute the estimation of the spline coefficients -----------------#

c_hat_function = function(theta, tau, gamma, ck, length_t, seq_t, B_seq_t, regul = 0.0)
{
 pen = pen_fun(theta, ck, length_t, seq_t, B_seq_t, min_t = min_t, max_t = max_t)

 R = pen $ R
 r = pen $ r

 P = as.matrix(gamma * R +
               rbind(cbind(tB_B * tau, tB_B * 0),
                     cbind(tB_B * 0, tB_B * tau) )) +
                     regul * diag(ncol(R))
 p = - gamma * r +
              c(tau * tB_y)

 c_hat = c(solve(P, tol = 1e-50) %*% p)

 return(list(as.numeric(c_hat), pen))
}


#------------------- function to compute the H criterion for estimation of ODE parameters -------------#

H_fun = function(theta, tau, gamma, ck, length_t, seq_t, B_seq_t, regul)
{ 
 ckp1 = c_hat_function(theta, tau, gamma, ck, length_t, seq_t, B_seq_t, regul = 0.0)[[1]]
 
 y_star_1 = (B_0 %*% matrix(ckp1, ncol = 2))[, 1]
 y_star_2 = (B_0 %*% matrix(ckp1, ncol = 2))[, 2]
 y_star = c(y_star_1, y_star_2)

 H =  - 0.5 * tau * sum((c(y_obs) - c(y_star)) ^ 2) 
 return(H)
}
  
