# Jake Yeung
# Date of Creation: 2019-09-14
# File: ~/projects/scchic/scripts/scripts_analysis/revisions/probabilistic_celltype_assignment.R
# Probabilistic celltype assignment notebook

#' ---
#' title: Automatically assigning celltypes to scChIC-seq data from bulk ChIP-seq
#' author: Jake Yeung
#' path: "~/projects/scchic/scripts/scripts_analysis/revisions/probabilistic_celltype_assignment.R"
#' date: 2019-09-14
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---


#' ## Introduction
#' 
#' I mentioned before that I think the scChIC-seq data (and probably scRNA-seq) of a single cell contains a lot more information than people think, 
#' because we’re able to do things like unmix doubly stained cells on a cell-by-cell basis, if we define the right likelihood model.
#' 
#' The nice thing about the bone marrow is that there is publicly available chip-seq data from which to compare (namingly, Lara-Astiaso, Science 2014 from Amit lab). 
#' 
#' A naive way of assigning celltypes is to make a 2D summary of your data, create clusters, generate pseudobulk from clusters, and find which pseudobulk corresponds to which celltype, 
#' given a publicly available dataset (we did this by just calculating correlations, between pseudobulk and public data).
#' 
#' The output of this was less than satisfactory, mostly because it's difficult to compare across many correlations 
#' (interpretation of correlation nonlinear and complicated: correlation of 0.7 versus 0.6 is not the same as correlation of 0.9 versus 0.8. Difference between 0.9 and 0.8 is much bigger than 0.7 versus 0.6.).
#' 
#' It would be nice to systematically calculate the likelihood that the single-cell data was generated from a particular ChIP-seq bulk sample, then compare across many bulk samples to 
#' calculate the probability that a cell came from a particular sample. 
#' 
#' The main take home message is that there is a surprisingly large amount of information coming from scChIC-seq data of an individual cell such that you can do celltype calling without creating
#' pseudobulk (as long as you have good collection of annotated bulk data). I am pretty sure this also applies to scRNA-seq. 
#' 


#' ## The generative model of scChIC-seq data
#' 
#' We assume the single scChIC-seq data is being generated from a multinomial process (e.g., drawing balls in an urn). 
#' In the balls in and urn analogy, the multinomial process is parametrized by a vector $\vec{p}$ of length $|p| = K$ representing the probability of drawing a ball from each of the K colors.
#' In the scChIC-seq data, drawing a ball is sequencing a UMI, and the K colors represent the K different genomic regions from which to draw a UMI. 
#' 
#' The multinomial likelihood $L$ of seeing a vector $\vec{y}$ of UMI counts across the K genomic bins in a cell (total UMI counts in the cell, $N = \sum_i{y}), is:
#' $$L = P(\vec{y} | \vec{p}, N) \propto \prod_{k=1}^K{\left( p_k \right) ^ {y_k}}$$
#' 
#' If you have M annotated celltypes from which to select, then the game is to parametrize M probability vectors $\vec{p_1}, \vec{p_2}, ..., \vec{p_M}$ such that you can take your single-cell observation $\vec{y}$ and
#' calculate the likelihood for each of the probability vectors, i.e., calculate $L_1, L_2, ..., L_M$ and then select the model with the highest likelihood. We can therefore use $\vec{L}$ to calculate the probability
#' that a cell came a celltype.
#' 


#' ## Analysis and results
#' 
#' Let's load our packages

library(JFuncs)
library(dplyr)
library(ggplot2)
library(data.table)
library(topicmodels)
library(tidytext)
library(umap)
library(ggrepel)

library(tidyr)

library(hash)
library(igraph)


#' ## Conclusion and discussion


