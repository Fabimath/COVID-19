#######################################################################################
##                            Penalty function lynx-hare dyn                         ##                     
##                            Lotka-Volterra ODE equations                           ##
#######################################################################################

## Function arguments:
## theta = current estimates ODE parameters
## ck = current estimates spline coefficients
## length_t = length grid B-spline for penalty matrix
## seq_t = grid of time points (B-spline  for penalty)
## B_seq_t = list of B-splines terms for the construction of the penalty
## min_t, max_t = min and max values of observed times

pen_fun = function(theta, ck, length_t, seq_t, B_seq_t, min_t, max_t)
{
   alpha = theta[1]  
   beta = theta[2]   
   gam = theta[3]    
   delta = theta[4]  

   B_x_0_seq_t = as.matrix(B_seq_t[[1]])
   B_x_1_seq_t = as.matrix(B_seq_t[[2]])

   ck1 = ck[1:M]
   ck2 = ck[-c(1:M)]
  
   B_x_0_seq_t_times_ck1 = B_x_0_seq_t %*% ck1
   B_x_1_seq_t_times_ck1 = B_x_1_seq_t %*% ck1

   B_x_0_seq_t_times_ck2 = B_x_0_seq_t %*% ck2
   B_x_1_seq_t_times_ck2 = B_x_1_seq_t %*% ck2

   w_1_1_1 = rep(1, length_t)
   w_1_1_0 = as.vector(-alpha + beta * B_x_0_seq_t_times_ck2)
   w_1_2_0 = as.vector(beta * B_x_0_seq_t_times_ck1)

   H1 = cbind(w_1_1_1 * B_x_1_seq_t + w_1_1_0 * B_x_0_seq_t,  B_x_0_seq_t * w_1_2_0)

   h1 = - as.vector(B_x_0_seq_t_times_ck1 * B_x_0_seq_t_times_ck2 * beta )

   w_2_2_1 = w_1_1_1
   w_2_2_0 = as.vector(gam - B_x_0_seq_t_times_ck1 * delta )
   w_2_1_0 = as.vector( - delta * B_x_0_seq_t_times_ck2)

   H2 =  cbind(B_x_0_seq_t * w_2_1_0 , w_2_2_1 * B_x_1_seq_t + w_2_2_0 * B_x_0_seq_t)
   h2 =  as.vector(delta * B_x_0_seq_t_times_ck1 * B_x_0_seq_t_times_ck2 )

   R1 = (max_t - min_t) /
          (length_t - 1) *
          (t(H1) %*% H1)

   R2 = (max_t - min_t) /
          (length_t - 1) *
          (t(H2) %*% H2)

   R =  R1 + R2

   r1 = (max_t - min_t) /
          (length_t - 1) *
          (t(H1) %*% h1)
   r2 = (max_t - min_t) /
          (length_t - 1) *
          (t(H2) %*% h2)
   r = c(r1 + r2)

   l1 = (max_t - min_t) /
        (length_t - 1) *
        (t(h1) %*% h1)
   l2 = (max_t - min_t) /
        (length_t - 1) *
        (t(h2) %*% h2)

   return(list(R = R, r = r, l = c(l1 + l2) ))
}

