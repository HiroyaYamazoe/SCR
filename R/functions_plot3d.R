#' S.C.R for 3d Data
#' 
#' Construct S.C.R for 3d Data with Algorithms in the paper.
#' @param data Matrix or data.frame with three columns.
#' @param t Vector whose length is the same as the number of rows in the data. It represents the covariate. This is not always necessary.
#' @param div Numeric. The number of vertices in the cross-section of the region.
#' @param alpha Numeric. 1 - confidence coefficient.
#' @param h_hat "MAX", "AM", "GM" or "MIN". This parameter determines the choice of h in the Algorithms.
#' @param cent Boolean. TRUE, Algorithm2; FALSE, Algorithm1.
#' @param p Numeric. This value is used in Algorithm2. The recommend value is written in the paper.
#' @param grid_points Numeric. the number of evaluate points.
#' @param point_size Numeric. Dot size of the data points in 3d plot.
#' @param point_color String. Dot color of the data points in 3d plot.
#' @param curve_size Numeric. Width of the estimated curve in 3d plot.
#' @param curve_color String. Color of the estimated curve in 3d plot.
#' @param scr_alpha Numeric from 0 to 1. Transparency of the SCR in 3d plot.
#' @param scr_color String. The color of the SCR in 3d plot.
#' @param xlabel String. The label of x-axis.
#' @param ylabel String. The label of y-axis.
#' @param zlabel String. The label of z-axis.
#' @param aspect Boolean. If TRUE (recommended), the figure is drawn to match the actual aspect.
#' 
#' @return The figure of SCR and the radius of the cross-section of the region.
#' 
#' @import dplyr
#' @import KernSmooth
#' @import locpol
#' @import rgl
#' @export
#' 
scr3d <- function(data, t = NULL, div = 16, alpha = 0.05, h_hat = "AM", cent = TRUE, p = 2, grid_points = 201L, point_size = 1, point_color = "black", curve_size = 4, curve_color = "black", scr_alpha = 0.4, scr_color = "#afafaf", xlabel = "y1", ylabel = "y2", zlabel = "y3", aspect = TRUE) {
   if (is.null(t)) {
      t <- prcomp(data, center = TRUE)
      t <- (t$x[, 1]) %>% scale() %>% pnorm()
   }
   d <- cbind(data, t)
   llr_h <- f_llr(d)
   llr <- llr_h[[1]]
   h_vec <- llr_h[[2]]
   llr_d <- f_llrd(d, h_vec, grid = grid_points)
   n_mat <- f_n_3d(llr_d)
   if (cent == TRUE) {
      rn <- f_rn_cent(d, dim = 3, h_vec = h_vec, grid = grid_points, alpha = alpha, p = p, flag = h_hat)
   } else {
      rn <- f_rn_noncent_3d(d, llr_d, h_vec, flag = h_hat, grid = grid_points, alpha = alpha)
   }
   end_a <- f_end(llr_d, n_mat, div, llr, rn, 0)
   end_b <- f_end(llr_d, n_mat, div, llr, rn, 1)
   mes <- f_mesh(n_mat, div, llr, rn)
   mat <- f_connect(mes, end_a, end_b, div)
   ind <- f_ind(mat, div)
   mat <- mat %>% cbind(., rep(1, nrow(.))) %>% t()
   mat_alt <- cbind(c(0, 1, 0, 0), c(1, 0, 0, 0), c(0, 0, 0, 1), c(0, 0, 1, 0))
   ind[, ncol(ind) / div * (1 : div)] <- mat_alt %*% ind[, ncol(ind) / div * (1 : div)]
   qmes <- qmesh3d(vertices = mat, indices = ind)
   plot3d(data, xlab = xlabel, ylab = ylabel, zlab = zlabel, col = point_color, lwd = point_size)
   plot3d(llr, type = "l", col = curve_color, lwd = curve_size, add = TRUE)
   shade3d(qmes, alpha = scr_alpha, col = scr_color)
   if (aspect) {
      aspect3d("iso")
   }
   return(rn)
}

f_3d_all <- function(data, div = 8, alpha = 0.05, flag = "MIN", cent = TRUE, p = 2, grid = 201L, beta = 0) {
   ####dat:$y_1,y_2,y_3,t####
   x <- data
   if (beta == 1) {
      s <- prcomp(x, center = TRUE)
      x <- (s$x[, 1]) %>% scale() %>% pnorm() %>% cbind(x, .)
   }
   ####llr:y_1,y_2,y_3,t####
   llr_h <- f_llr(x)
   llr <- llr_h[[1]]
   h_vec <- llr_h[[2]]
   llrd <- f_llrd(x, h_vec, grid = grid)
   n_mat <- f_n_3d(llrd)
   if (cent == TRUE) {
      r_n <- f_rn_cent(x, dim = 3, h_vec = h_vec, grid = grid, alpha = alpha, p = p, flag = flag)
   } else {
      r_n <- f_rn_noncent_3d(x, llrd, h_vec, flag = flag, grid = grid, alpha)
   }
   end_a <- f_end(llrd, n_mat, div, llr, r_n, 0)
   end_b <- f_end(llrd, n_mat, div, llr, r_n, 1)
   mes <- f_mesh(n_mat, div, llr, r_n)
   mat <- f_connect(mes, end_a, end_b, div)
   ind <- f_ind(mat, div)
   mat <- mat %>% cbind(., rep(1, nrow(.))) %>% t()
   ####Handling of endpoints####
   mat_alt <- cbind(c(0, 1, 0, 0), c(1, 0, 0, 0), c(0, 0, 0, 1), c(0, 0, 1, 0))
   ind[, ncol(ind) / div * (1 : div)] <- mat_alt %*% ind[, ncol(ind) / div * (1 : div)]
   return(list(mat, ind, r_n))
}

f_rn_noncent_3d <- function(d, drv_mat, h_vec, flag, grid = 201L, alpha = 0.05) {
   dim <- 3
   n <- nrow(d)
   h <- switch(flag,
      "AM" = mean(h_vec),
      "GM" = prod(h_vec) ^ (1 / length(h_vec)),
      "MIN" = min(h_vec),
      "MAX" = max(h_vec)
   )
   r_hat <- f_r(d)
   t <- d[, dim + 1]
   u <- seq(min(t), max(t), length = 201)
   t_density <- density(t, kernel = "gaussian", n = 201L, from = min(t), to = max(t))$y
   eta_hat <- t_density %>% min()
   n_mat <- f_n_3d(drv_mat)
   n1 <- n_mat[, 1 : dim]
   n2 <- n_mat[, (dim + 1) : ncol(n_mat)]
   ####phi_2####
   f_hat <- NULL
   coef_vec <- c(2, 6, 12, 20)
   for(i in 1 : dim) {
      y <- d[, i]
      fit <- lm(y ~ poly(t, degree = 5, raw = TRUE))
      fit <- coefficients(fit) %>% as.vector() %>% .[3 : 6]
      f_hat <- cbind(f_hat, cbind(rep(1, length(u)), u, u ^2, u^3) %*% (fit * coef_vec))
   }
   n_max <- ((n1 * f_hat) ^ 2 + (n2 * f_hat) ^ 2) %>% rowSums()
   n_max <- (n_max * t_density ^ 2) %>% max()
   g_max <- (n * h ^ 5 * (dim + 1) ^ 2) / (4 * r_hat ^ 4 * rk ^ 2) * n_max
   chi <- qchisq(p = alpha, df = dim - 1, ncp = g_max, lower.tail = F)
   c_hat <- (rk / (dim + 1) / eta_hat * chi) %>% sqrt()
   r_n <- c_hat * r_hat / sqrt(n * h)
   return(r_n)
}


####pic up two orthogonal vectors
f_n_3d <- function(d_mat) {
   data <- d_mat[, -4]
   n1_mat <- cbind(data[, 2], - data[, 1], numeric(nrow(data)))
   n2_mat <- cbind(
      data[, 1] * data[, 2] * data[, 3] / (data[, 1] ^ 2 + data[, 2] ^ 2),
      data[, 3] * (
         1 - data[, 1] ^ 2 / (data[, 1] ^ 2 + data[, 2] ^ 2)
      ),
      - data[, 2]
   )
   n1_mat <- n1_mat / sqrt(apply(n1_mat ^ 2, 1, sum))
   n2_mat <- n2_mat / sqrt(apply(n2_mat ^ 2, 1, sum))
   cbind(n1_mat, n2_mat) %>% return()
}



####create mesh from two orthogonal vectors####
f_mesh <- function(n_mat, n = 8, llr, r) {
   theta <- 2 * pi * (1 : n) / n
   n1 <- n_mat[, 1 : 3]
   n2 <- n_mat[, 4 : 6]
   nn <- NULL
   for(i in 1 : length(theta)) {
      n <- n1 * cos(theta[i]) + n2 * sin(theta[i])
      n <- n * r + llr[, - ncol(llr)]
      nn <- rbind(nn, n)
   }
   return(nn)
}


###Handling of endpoints####
###ab:choice of endpoint###
f_end <- function(d_mat, n_mat, n = 8, llr, r, ab) {
   theta <- 2 * pi * (1 : n) / n
   m <- nrow(n_mat)
   n1 <- n_mat[c(1, m), 1 : 3]
   n2 <- n_mat[c(1, m), 4 : 6]
   if (ab == 0) {
      ###a###
      ####1.create all points on circle. 2.stretch it for the tip####
      ####create n phi_a####
      a <- llr[1, - ncol(llr)] %>% rep(., n) %>% matrix(ncol = n) %>% t()
      ####create n phi_d_a and standardize####
      a_d <- d_mat[1, -4] / sqrt(sum(d_mat[1, -4] ^ 2))
      a_d <- a_d %>% rep(., n) %>% matrix(ncol = n) %>% t()
      #n: the number of points to the tip and create r * 1/n , ... r*n/n
      r_vec <- r * (n : 1) / n
      nna <- NULL
      na <- NULL
      ####derive circle points####
      for(arg in theta) {
         na <- rbind(na, n1[1, ] * cos(arg) + n2[1, ] * sin(arg))
      }
      ####stretch the circle points for the tip####
      for(i in r_vec) {
         naa <- sqrt(r ^ 2 - i ^ 2) * na + a - a_d * i
         nna <- rbind(nna, naa)
      }
      ####sort it####
      nna_alt <- NULL
      for(i in 1 : n) {
         n_vec <- n * (0 : (n - 1)) + i
         nna_alt <- rbind(nna_alt, nna[n_vec, ])
      }
      return(nna_alt)
   } else if (ab == 1) {
      ###b###
      b <- llr[m, - ncol(llr)] %>% rep(., n) %>% matrix(ncol = n) %>% t()
      b_d <- d_mat[m, -4] / sqrt(sum(d_mat[m, - 4] ^ 2))
      b_d <- b_d %>% rep(., n) %>% matrix(ncol = n) %>% t()
      r_vec <- r * (1 : n) / n
      nnb <- NULL
      nb <- NULL
      for(arg in theta) {
         nb <- rbind(nb, n1[2, ] * cos(arg) + n2[2, ] * sin(arg))
      }
      for(i in r_vec) {
         nbb <- sqrt(r ^ 2 - i ^ 2) * nb + b + b_d * i
         nnb <- rbind(nnb, nbb)
      }
      ####sor it####
      nnb_alt <- NULL
      for(i in 1 : n) {
         n_vec <- n * (0 : (n - 1)) + i
         nnb_alt <- rbind(nnb_alt, nnb[n_vec, ])
      }
      return(nnb_alt)
   }
}



####given div * n times 3 matrix, create indices####
####it depends on the order of the matrix
f_ind <- function(mat, div = 8) {
   n <- nrow(mat) / div
   a <- 1 : (n - 1)
   b <- rbind(a, a + 1, a + n + 1, a + n)
   c <- NULL
   for(i in 1 : (div - 1)) {
      c <- cbind(c, b + n * (i - 1))
   }
   d <- cbind(c, rbind(a, a + 1, a + (div - 1) * n + 1, a + (div - 1) * n))
   return(d)
}

f_connect <- function(mj, ma, mb, div = 8) {
   n <- nrow(mj) / div
   mat <- NULL
   for(i in 1 : div) {
      mat <- rbind(mat, ma[(i - 1) * div + (1 : div), ], mj[(i - 1) * n + (1 : n), ], mb[(i - 1) * div + (1 : div), ])
   }
   return(mat)
}
   