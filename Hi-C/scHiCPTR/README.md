# scHiCPTR

## Introduction
scHiCPTR is an unsupervised pseudotime inference pipeline through dual graph refinement for single cell Hi-C data. It reconciles pseudotime inference in the case of circular and bifurcating topology.

## Installation

**OS**
- ubuntu 18.04

**Required Python Packages**

Make sure all the packages listed in the *requirements.txt* are installed.

- Python>=3.6
- numpy>=1.18.0
- pandas
- igraph>=0.9.8
- matplotlib>=3.3.2



**Install from Github**

First, download the folder *scHiCPTR*.

Second, install the package with following command:

```
$ cd scHiCPTR
$ pip install .
```


## Usage
**Input data format**

- Sparse matrices: For each cell, each intra-chromosomal contact matrice is stored in a .txt file. The file contain three columns separated by tab, i.e. sparse format of the contact matrix . The name of the file need to be in the format of 'cell_{id}_chr{nm}.txt'. {id} is the number of cells (like 1, 2, 3,...) and {nm} is the number of chromosome (like 1, 2, ..., X).

**Pseudotime inference**
- See examples/cellcycle.py

**Visualization**
- See examples/cellcycle.py
