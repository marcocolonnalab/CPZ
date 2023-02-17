# Transcriptomic atlas and interaction networks of brain cells in mouse CNS demyelination and remyelination


This repository contains the code used for snRNA-seq analysis by Hou et al. (2023) of the mouse cuprizone model.
[![DOI](https://zenodo.org/badge/602279096.svg)](https://zenodo.org/badge/latestdoi/602279096)

## Code
Included are the codes necessary to replicate the analyses.

 - `cell type clustering`: clustering analyses of all cells and subclustering of each cell type, using `Seurat`.
 - `nichenet`: ligand-receptor interaction analyses between other brain cells and microglia, using `nichenet`.
 
## Data
Raw and processed data are available at the Gene Expression Omnibus (GEO) database under accession number GSE204770. 



## Acknowledgements
This work was supoorted by the NIH (RF1 AG051485, R21 AG059176, and RF1 AG059082) (M. Colonna), Cure Alzheimer's Fund (M. Colonna).

We would like to thank Vincent Peng for advice on bioinformatic analyses.
