% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_coenoclines.r
\name{plot_coenoclines}
\alias{plot_coenoclines}
\title{plot standardized dispersal kernel}
\usage{
plot_coenoclines(sim_result = NULL, landscape = NULL, Ef = NULL,
  trait_Ef = NULL, trait_Ef_sd = NULL, ...)
}
\arguments{
\item{sim_result}{output from a simulation, default is NULL}

\item{landscape}{a landscape object, default is NULL}

\item{Ef}{numeric vector of environmental filter values for sites in a simulation result, default is NULL}

\item{trait_Ef}{numeric vector of niche positions of species in metacommunity, default is NULL}

\item{trait_Ef_sd}{numeric vector of niche breadths, default is NULL}
}
\description{
plot standardized dispersal kernel
}