
###definition of kernel###
f_gk <- function(x, h) {
   gk <- function(t) {
      1 / sqrt(2 * pi) * exp(
         - 1 / 2 * t ^ 2
      )
   }
   kh <- 1 / h * gk(x / h)
   return(kh)
}

###R(K) of gaussian kernel###
rk <- 1 / 2 / sqrt(pi)

###derive LLE###
###x and y: data, t_vec: evaluate points###
f_locpoly_vec <- function(x, y, t_vec) {
   h <- thumbBw(x, y, deg = 1, kernel = gaussK)
   f_hat <- NULL
   for(t in t_vec) {
      x_mat <- x %>% length(.) %>% rep(1, .) %>% cbind(., x - t)
      w_t <- f_gk(x_mat[, 2], h) %>% diag()
      f_hat_i <- solve(
         t(x_mat) %*% w_t %*% x_mat
      ) %*% t(x_mat) %*% w_t %*% y
      f_hat <- c(f_hat, f_hat_i[1, 1])
   }
   return(f_hat)
}


####derive LLE with KernSmooth####
####returned data: phi1, phi2, phi3, t####
f_llr <- function(d, t = NULL, grid = 201L) {
   if (is.null(t)) {
      t <- d[, ncol(d)]
      data <- d[, - ncol(d)]
   } else {
      data <- d
   }
   rgm <- NULL
   h_vec <- NULL
   for(i in 1 : ncol(data)) {
      y <- data[, i]
      h <- thumbBw(t, y, deg = 1, kernel = gaussK)
      fit <- locpoly(t, y, bandwidth = h, drv = 0, range.x = range(t), degree = 1, gridsize = grid, truncate = FALSE)
      rgm <- cbind(rgm, fit$y)
      h_vec <- c(h_vec, h)
   }
   cbind(rgm, fit$x) %>% list(., h_vec) %>% return()
}



####derive gradient vector of LLE####
f_llrd <- function(d, h_vec, grid = 201L) {
   data <- d[, - ncol(d)]
   t_vec <- d[, ncol(d)]
   grid <- seq(min(t_vec), max(t_vec), length = grid)
   f_wt <- function(x_vec, x, h) {
     (
      1 / h * (2 * pi) ^ (- 1 / 2) * exp(
      - (x_vec - x) ^ 2 / 2 / h ^ 2
     )
     ) %>% diag() %>% return()
   }
   f_xt <- function(x_vec, x) {
      rep(1, length(x_vec)) %>% cbind(., x_vec - x) %>% return()
   }
   f_dxt <- function(x_vec) {
      n <- length(x_vec)
      cbind(numeric(n), rep(-1, n))
   }
   f_dwt <- function(x_vec, x, h) {
      (
         1 / h ^ 3 / sqrt(2 * pi) * (x_vec - x) * exp(
            - (x_vec - x) ^ 2 / 2 / h ^ 2
         )
      ) %>% diag() %>% return()
   }
   ####phi_i#####
   phi_mat <- NULL
   for(i in 1 : ncol(data)) {
      y <- data[, i]
      h <- h_vec[i]
      phi_d <- NULL
      ####t runs on I####
      for(t in grid) {
         xt <- f_xt(t_vec, t)
         wt <- f_wt(t_vec, t, h)
         dxt <- f_dxt(t_vec)
         dwt <- f_dwt(t_vec, t, h)
         f_hat <- solve(t(xt) %*% wt %*% xt) %*% (
         - (
            t(dxt) %*% wt %*% xt + t(xt) %*% dwt %*% xt + t(xt) %*% wt %*% dxt
         ) %*%  solve(t(xt) %*% wt %*% xt) %*% t(xt) %*% wt + t(dxt) %*% wt + t(xt) %*% dwt
         ) %*% y
         f_hat <- f_hat[1]
         phi_d <- c(phi_d, f_hat)
      }
      phi_mat <- cbind(phi_mat, phi_d)
   }
   return(cbind(phi_mat, grid))
}


###derive r_hat###
f_r <- function(d) {
   t <- d[, ncol(d)]
   data <- d[, - ncol(d)]
   f_hat <- NULL
   for(i in 1 : ncol(data)) {
      y <- data[, i]
      f_hat_i <- f_locpoly_vec(t, y, t)
      f_hat <- cbind(f_hat, f_hat_i)
   }
   tr_data <- var(data) %>% diag() %>% sum()
   tr_f_hat <- var(f_hat) %>% diag() %>% sum()
   trvar <- tr_data - tr_f_hat
   if (trvar > 0) {
      dim <- ncol(data)
   r <- (dim + 1) / (dim - 1) * (tr_data - tr_f_hat)
   return(sqrt(r))
   } else {
      (
         f_hat - data
      ) ^ 2 %>% rowSums() %>% sqrt() %>% max() %>% return()
   }
}


####calculate r_n with Algorithm 2 (2d and 3d)
f_rn_cent <- function(d, dim, h_vec, grid = 201L, alpha = 0.05, p = 2, flag = "MIN") {
   n <- nrow(d)
   t <- d[, dim + 1]
   h <- switch(flag, 
      "AM" = mean(h_vec) ^ p,
      "GM" = prod(h_vec) ^ (p / length(h_vec)),
      "MIN" = min(h_vec) ^ p,
      "MAX" = max(h_vec) ^ p
   )
   r_hat <- f_r(d)
   t_density <- density(t, kernel = "gaussian", n = grid, from = min(t), to = max(t))$y
   eta_hat <- t_density %>% min()
   c_hat <- (rk / (dim + 1) / eta_hat * qchisq(p = alpha, df = dim - 1, lower.tail = F)) %>% sqrt()
   return(c_hat * r_hat / sqrt(n * h))
}
