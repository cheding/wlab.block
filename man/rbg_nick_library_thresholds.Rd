% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rbg_nick_library_thresholds.R
\name{rbg_nick_library_thresholds}
\alias{rbg_nick_library_thresholds}
\title{A function to determine mutation background in nick library based on frequency}
\usage{
rbg_nick_library_thresholds(
  input_files_thre,
  input_dir,
  output_dir = NULL,
  block = c("1", "2", "3"),
  min_obs = 200,
  min_phenotypes = NULL
)
}
\arguments{
\item{input_files_thre}{libraries used for filtering}

\item{input_dir}{rbg_nick_library output directory}

\item{output_dir}{output path (default=input_dir)}

\item{block}{block information of input_files_thre (default=c('1','2','3'))}

\item{min_obs}{minimum permitted mutation reads}

\item{min_phenotypes}{minimum permitted frequency in different library}
}
\value{
Nothing
}
\description{
A function to determine mutation background in nick library based on frequency
}
