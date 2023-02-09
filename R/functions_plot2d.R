#' SCR for 2d Data
#'
#' Construct SCR for 2d Data with Algorithms in Yamazoe and Naito (2023).
#' @param data Matrix or data.frame with two columns.
#' @param t Vector whose length is the same as the number of rows in the data. It represents the covariate. This is not always necessary.
#' @param alpha Numeric. 1 - confidence coefficient.
#' @param h_coef Numeric. The bandwidth used in local linear regression for phi_i is h_coef * h_ROT.
#' @param h_hat "MAX", "AM", "GM" or "MIN". This parameter determines the choice of h in the Algorithms.
#' @param cent Boolean. TRUE, Algorithm2; FALSE, Algorithm1.
#' @param p Numeric. This value is used in Algorithm2. The recommend value is written in the paper.
#' @param grid_points Numeric. the number of evaluate points.
#' @param point_size Numeric. Dot size of the data points in 2d plot.
#' @param point_color String. Dot color of the data points in 2d plot.
#' @param point_shape Numeric. Shape of the data points in 2d plot.
#' @param curve_size Numeric. Width of the estimated curve in 2d plot.
#' @param curve_color String. Color of the estimated curve in 2d plot.
#' @param curve_linetype String. Line type of the estimated curve in 2d plot.
#' @param scr_alpha Numeric from 0 to 1. Transparency of the SCR in 2d plot.
#' @param scr_color String. The inner color of the SCR in 2d plot.
#' @param scr_fill String. The boundary color of the SCR. It should be the same as the inner color.
#' @param xlabel String. The label of x-axis.
#' @param ylabel String. The label of y-axis.
#' @param xlimit Vector from the minimum value to the maximum value of x-axis to illustrate.
#' @param ylimit Vector from the minimum value to the maximum value of y-axis to illustrate.
#' @param plot_title String. The title of the figure.
#' @param aspect Boolean. If TRUE (recommended), the figure is drawn to match the actual aspect.
#' 
#' @return The figure of SCR, and list of the radius of the cross-section of the region and the bandwidth vector in local linear regression for phi_1 to phi_d.
#' 
#' @import ggplot2
#' @import KernSmooth
#' @import locpol
#' @import dplyr
#' @export 
#' 
scr2d <- function(data, t = NULL, alpha = 0.05, h_coef = 1, h_hat = "AM", cent = TRUE, p = 4 / 3, grid_points = 201L, point_size = 3, point_color = "black", point_shape = 1, curve_size = 1.5, curve_color = "black", curve_linetype = "solid", scr_alpha = 1, scr_color = "#dfdfdf", scr_fill = "#dfdfdf", xlabel = "y1", ylabel = "y2", xlimit = NULL, ylimit = NULL, plot_title = "", aspect = TRUE) {
   if (ncol(data) != 2) {
      stop("the data must be two dimensional.")
   }
   if (is.null(t)) {
      t <- prcomp(data, scale = TRUE)
      t <- (t$x[, 1]) %>% scale() %>% pnorm()
   }
   d <- cbind(data, t)
   llr_h <- f_llr(d, grid = grid_points, h_coef = h_coef)
   llr <- llr_h[[1]]
   h_vec <- llr_h[[2]]
   llr_d <- f_llrd(d, h_vec, grid = grid_points)
   n_mat <- f_nlu(llr_d, llr)
   if (cent == TRUE) {
      rn <- f_rn_cent(d, 2, h_vec, grid = grid_points, alpha = alpha, p = p, flag = h_hat)
   } else {
      rn <- f_rn_noncent_2d(d, llr_d, h_vec, flag = h_hat, grid = grid_points, alpha = alpha)
   }
   poi <- f_end_2d(llr_d, n_mat, llr, 16, rn)
   df <- data.frame(data)
   colnames(df) <- c("Y1", "Y2")
   llr <- data.frame(llr)
   colnames(llr) <- c("Y1", "Y2", "t")
   poi <- data.frame(poi)
   colnames(poi) <- c("Y1", "Y2")
   p <- ggplot(NULL) + geom_polygon(data = poi, aes(x = Y1, y = Y2), alpha = scr_alpha, col = scr_color, fill = scr_fill) + 
   geom_point(data = df, aes(x = Y1, y = Y2),shape = point_shape, size = point_size, color = point_color) + 
   geom_path(data = llr, aes(x = Y1, y = Y2), linewidth = curve_size, color = curve_color, linetype = curve_linetype) +
   mytheme + xlab(xlabel) + ylab(ylabel) + ggtitle(plot_title)
   if (!is.null(xlimit)) {
      p <- p + xlim(xlimit)
   }
   if (!is.null(ylimit)) {
      p <- p + ylim(ylimit)
   }
   if (aspect) {
      p <- p + coord_equal()
   }
   plot(p)
   return(list(r_n = rn, h_vector = h_vec))
}


f_end_2d <- function(d_mat, n_mat, llr, n = 16, r) {
   theta <-  r * (1 : n) / n
   a_d <- - d_mat[1, - ncol(d_mat)] / sqrt(sum(d_mat[1, -ncol(d_mat)] ^ 2))
   n_a_u <- n_mat[1, 1 : 2]
   n_a_l <- n_mat[1, 3 : 4]
   nau_vec <- NULL
   nal_vec <- NULL
   for (arg in theta) {
      nau <-  arg * a_d + sqrt(r ^ 2 - arg ^ 2) * n_a_u
      nal <- arg * a_d + sqrt(r ^ 2 - arg ^ 2) * n_a_l
      nau_vec <- rbind(nau_vec, nau)
      nal_vec <- rbind(nal_vec, nal)
   }
   llr_1 <- llr[1, -ncol(llr)] %>% rep(., n) %>% matrix(nrow = 2) %>% t()
   nau_vec <- nau_vec + llr_1
   nal_vec <- nal_vec + llr_1
   b_d <- d_mat[nrow(d_mat), - ncol(d_mat)] / sqrt(sum(d_mat[nrow(d_mat), - ncol(d_mat)] ^ 2))
   n_b_u <- n_mat[nrow(n_mat), 1 : 2]
   n_b_l <- n_mat[nrow(n_mat), 3 : 4]
   nbu_vec <- NULL
   nbl_vec <- NULL
   for (arg in theta) {
      nbu <- arg * b_d + sqrt(r ^ 2 - arg ^ 2) * n_b_u
      nbl <- arg * b_d + sqrt(r ^ 2 - arg ^ 2) * n_b_l
      nbu_vec <- rbind(nbu_vec, nbu)
      nbl_vec <- rbind(nbl_vec, nbl)
   }
   llr_n <- llr[nrow(llr), - ncol(llr)] %>% rep(., n) %>% matrix(nrow = 2) %>% t()
   nbu_vec <- nbu_vec + llr_n
   nbl_vec <- nbl_vec + llr_n
   rbind(nal_vec,
      apply(nau_vec, 2, rev),
      n_mat[1 : nrow(n_mat), 1 : 2] * r + llr[, - ncol(llr)],
      nbu_vec,
      apply(nbl_vec, 2, rev),
      apply(n_mat[1 : nrow(n_mat), 3 : 4] * r + llr[, - ncol(llr)], 2, rev)
   ) %>% return()
}

f_nlu <- function(d_mat, llr) {
   data <- d_mat[, -3]
   nl <- cbind(data[, 2], - data[, 1])
   nl <- nl / sqrt(apply(nl ^ 2, 1, sum))
   nu <- - nl
   return(cbind(nu, nl))
}

#####pic up an orthogonal vector (2d)#####
f_n_2d <- function(d_mat) {
   data <- d_mat[, - ncol(d_mat)]
   n1_mat <- cbind(data[, 2], - data[, 1])
   n1_mat <- n1_mat / sqrt(apply(n1_mat ^ 2, 1, sum))
   return(n1_mat)
}

####calculate r_n with Algorithm 1 (2d only)
f_rn_noncent_2d <- function(d, drv_mat, h_vec, flag, grid = 201L, alpha = 0.05) {
   dim <- 2
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
   n1 <- f_n_2d(drv_mat)
   f_hat <- NULL
   coef_vec <- c(2, 6, 12, 20)
   for(i in 1 : dim) {
      y <- d[, i]
      fit <- lm(y ~ poly(t, degree = 5, raw = TRUE))
      fit <- coefficients(fit) %>% as.vector(.) %>% .[3 : 6]
      fit_i <- cbind(rep(1, length(u)), u, u ^2, u^3) %*% (fit * coef_vec)
      f_hat <- cbind(f_hat, fit_i)
   }
   n_max <- (n1 * f_hat) ^ 2 %>% rowSums()
   n_max <- (t_density * n_max) %>% max()
   g_max <- (1 / r_hat ^ 4) * ((dim + 1) ^ 2 * n * h ^ 5) / (4 * rk ^ 2) * n_max
   chi <- qchisq(p = alpha, df = dim - 1, ncp = g_max, lower.tail = F)
   c_hat <- (rk / (dim + 1) / eta_hat * chi) %>% sqrt()
   r_n <- c_hat * r_hat / sqrt(n * h)
   return(r_n)
}

mytheme <- theme_bw() + theme(axis.title = element_text(size = rel(2.0)), axis.text = element_text(size = rel(2.0)))