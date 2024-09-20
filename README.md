# Diagnostic Evidence GAuge of Single cells (DEGAS): A flexible deep transfer learning framework for prioritizing cells in relation to disease
### https://doi.org/10.1186/s13073-022-01012-2
### Package development by: <br>Zhi Huang (https://github.com/huangzhii)<br>Travis S. Johnson (https://github.com/tsteelejohnson91)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

![DEGAS](figures/DEGAS.png "DEGAS")

## Installation
* Step 1 Install python3 and pip3 if not previously installed 
  * https://wiki.python.org/moin/BeginnersGuide/Download
  * https://www.python.org/downloads/
* Step 2 Install python packages from terminal
```bash
pip3 install tensorflow
pip3 install functools
pip3 install numpy
pip3 install math
```
* Step 3 Install devtools in R
```R
install.packages("devtools")
```
* Step 4 Install DEGAS in R
```R
library(devtools)
install_github("tsteelejohnson91/DEGAS")
```
* Step 5 Install packages useful for downstream analysis in R
```R
install.packages("Rtsne")
install.packages("ggplot2")
```
## Prerequisites

### OS
* OSX
* Linux

### Python packages
* tensorflow
* functools
* numpy
* math

### R
* Rtsne
* ggplot2

## Configurations tested

### Mac CPU
* R (4.0.1), Python (3.8.2), TensorFlow (2.3.1), NumPy (1.18.5), functools (3.8.2), math (3.8.2)
* R (4.1.0), Python (3.8.3), TensorFlow (2.5.0), NumPy (1.19.5), functools (3.8.3), math (3.8.3)
* R (3.5.1), Python (3.6.0), TensorFlow (2.3.1), NumPy (1.17.4), functools (3.6.0), math (3.6.0)

### Linux GPU
* R (3.4.4), Python (anaconda 3.6.5), TensorFlow GPU-enabled (1.9.0), NumPy (1.14.3), functools (anaconda 3.6.5), math (anaconda 3.6.5)

### R Package Versions
* Rtsne: 0.15, 
* ggplot2: 3.3.5, 3.3.0, 3.2.1
* See Session Info in the Simulation, GBM, AD, and MM examples for more version details
  * e.g. See bottom of https://github.com/tsteelejohnson91/DEGAS/blob/master/MM_example/MM_example.md

### Python package versions
* TensorFlow: 2.3.1, 2.5.0, 2.3.1, 1.9.0 (GPU)
* Numpy: 1.18.5, 1.19.5, 1.17.4, 1.14.3
* functools: 3.8.2, 3.8.3, 3.6.0, anaconda 3.6.5
* math: 3.8.2, 3.8.3, 3.6.0, anaconda 3.6.5

## DEGAS documentation and examples
* Documentation (https://github.com/tsteelejohnson91/DEGAS/blob/master/DEGAS_documentation.md)
* Simulation example (https://github.com/tsteelejohnson91/DEGAS/blob/master/Sim_example/Sim_example.md)
* GBM example (https://github.com/tsteelejohnson91/DEGAS/blob/master/GBM_example/GBM_example.md)
* AD example (https://github.com/tsteelejohnson91/DEGAS/blob/master/AD_example/AD_example.md)
* MM example (https://github.com/tsteelejohnson91/DEGAS/blob/master/MM_example/MM_example.md)
* BRCA spatial transcriptomics example (https://github.com/tsteelejohnson91/DEGAS/blob/master/BRCA_ST_example/BRCA_ST_example.md)

## Coming Soon
Diagnostic Evidence GAuge of Spatial transcriptomics (DEGAS)

### Usage tracking
Please click this link to help us track our usage.
https://tsteelejohnson91.github.io/DEGAS/




