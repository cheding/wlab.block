% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rbg_synthetic_library.R
\name{rbg_synthetic_library}
\alias{rbg_synthetic_library}
\title{A Function to remove false mutations based on mutation backgrounds in synthetic library}
\usage{
rbg_synthetic_library(
  input_dir_list,
  output_dir = "./",
  position,
  codon = "NNK"
)
}
\arguments{
\item{input_dir_list}{list contain path to "variant_data_merge.tsv"}

\item{output_dir}{output path (default="./")}

\item{position}{list contain mutation background position}

\item{codon}{code type (default="NNK")}
}
\value{
Nothing
}
\description{
A Function to remove false mutations based on mutation backgrounds in synthetic library
}
