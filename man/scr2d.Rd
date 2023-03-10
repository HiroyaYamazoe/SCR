% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_plot2d.R
\name{scr2d}
\alias{scr2d}
\title{SCR for 2d Data}
\usage{
scr2d(
  data,
  t = NULL,
  alpha = 0.05,
  h_coef = 1,
  h_hat = "AM",
  cent = TRUE,
  p = 4/3,
  grid_points = 201L,
  point_size = 3,
  point_color = "black",
  point_shape = 1,
  curve_size = 1.5,
  curve_color = "black",
  curve_linetype = "solid",
  scr_alpha = 1,
  scr_color = "#dfdfdf",
  scr_fill = "#dfdfdf",
  xlabel = "y1",
  ylabel = "y2",
  xlimit = NULL,
  ylimit = NULL,
  plot_title = "",
  aspect = TRUE
)
}
\arguments{
\item{data}{Matrix or data.frame with two columns.}

\item{t}{Vector whose length is the same as the number of rows in the data. It represents the covariate. This is not always necessary.}

\item{alpha}{Numeric. 1 - confidence coefficient.}

\item{h_coef}{Numeric. The bandwidth used in local linear regression for phi_i is h_coef * h_ROT.}

\item{h_hat}{"MAX", "AM", "GM" or "MIN". This parameter determines the choice of h in the Algorithms.}

\item{cent}{Boolean. TRUE, Algorithm2; FALSE, Algorithm1.}

\item{p}{Numeric. This value is used in Algorithm2. The recommend value is written in the paper.}

\item{grid_points}{Numeric. the number of evaluate points.}

\item{point_size}{Numeric. Dot size of the data points in 2d plot.}

\item{point_color}{String. Dot color of the data points in 2d plot.}

\item{point_shape}{Numeric. Shape of the data points in 2d plot.}

\item{curve_size}{Numeric. Width of the estimated curve in 2d plot.}

\item{curve_color}{String. Color of the estimated curve in 2d plot.}

\item{curve_linetype}{String. Line type of the estimated curve in 2d plot.}

\item{scr_alpha}{Numeric from 0 to 1. Transparency of the SCR in 2d plot.}

\item{scr_color}{String. The inner color of the SCR in 2d plot.}

\item{scr_fill}{String. The boundary color of the SCR. It should be the same as the inner color.}

\item{xlabel}{String. The label of x-axis.}

\item{ylabel}{String. The label of y-axis.}

\item{xlimit}{Vector from the minimum value to the maximum value of x-axis to illustrate.}

\item{ylimit}{Vector from the minimum value to the maximum value of y-axis to illustrate.}

\item{plot_title}{String. The title of the figure.}

\item{aspect}{Boolean. If TRUE (recommended), the figure is drawn to match the actual aspect.}
}
\value{
The figure of SCR, and list of the radius of the cross-section of the region and the bandwidth vector in local linear regression for phi_1 to phi_d.
}
\description{
Construct SCR for 2d Data with Algorithms in Yamazoe and Naito (2023).
}
