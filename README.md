# EPI-clone
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/veltenlab/EPI-clone/blob/main/LICENSE)

## Introduction

This repository comprises code to reproduce the analysis of the [EPI-clone paper](https://doi.org/10.1038/s41586-025-09041-8). It uses the Seurat object from [figshare](https://doi.org/10.6084/m9.figshare.24204750.v1) to generate the plots from the paper. We refer to the [figures](figures) folder for an explanation on the datasets and the corresponding code. Additionally, utility functions are available in the [scripts](scripts) folder. The EPI-clone function within the [scripts](scripts) folder comprises a functionality to identify clones from scDNAm data in the form of a Seurat object.

## Installation

EPI-clone is a stand-alone R script that required the following dependencies:

* ggplot2 > v3.4.1
* Seurat > v4.3.0
* ROCR > v1.0
* fossil > v0.4.0
* reshape2 > v1.4.4.
* ComplexHeatmap > v2.14.0
* pheatmap > v1.0.12
* GenomicRanges > v1.50.2
* viridis > v0.6.2
* RnBeads > v2.20.0

It has been tested on R v4.2.2 on Ubuntu 22.04.2.

## Expected run time

For the example data set (13,000 cells), EPI-clone requires less than 5 minutes to finish.

## Citation

Scherer, M., Singh, I., Braun, M.M. et al. Clonal tracing with somatic epimutations reveals dynamics of blood ageing. Nature (2025). https://doi.org/10.1038/s41586-025-09041-8

## Contact

For questions and requests, please contact [Michael Scherer](mailto:michael.scherer@dkfz.de).
