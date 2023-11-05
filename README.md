# An isoform-resolution transcriptomic atlas of colorectal cancer from long-read single-cell sequencing
This repository contains data analysis and demo code for the manuscript 

Li,Zhongxiao, et al. "An isoform-resolution transcriptomic atlas of colorectal cancer from long-read single-cell sequencing" \[[bioRxiv](https://www.biorxiv.org/content/10.1101/2023.04.21.536771v3)\]

## Prerequisites
The data analysis code in this repo was developed with the following Python dependencies:
- Python 3.8.13
- Numpy 1.22.4
- Scipy 1.7.0
- Pandas 1.3.2
- h5py 3.7.0
- matplotlib 3.5.2
- seaborn 0.11.2
- anndata 0.8.0
- scanpy 1.9.1

The code was developed with the following R dependencies:
- R 4.0.5
- anndata 0.7.5.6
- Seurat 4.0.1
- monocle3 1.0.1

If you are using a package of a different version, you need to modify the code if there is an API change.
## Part 1: 10x Illumina data analysis
<div align="center">
  <img src="images/10x_Illumina.png" width="350" height="225">
</div>

This part demonstrates how to use the processed 10x Illumina part of the atlas. Relevant analyses are in `10x-Illumina.ipynb`.

## Part 2: 10x PacBio data analysis
<div align="center">
  <img src="images/10x_PacBio.png" width="400" height="175">
</div>

This part demonstrates how to use the processed 10x PacBio part of the atlas. Relevant analyses are in `10x-PacBio.ipynb`.

## Part 3: Epithelial lineage analysis
<div align="center">
  <img src="images/lineage_analysis.png" width="375" height="250">
</div>

This part demonstrates the identified epithelial cell differentiation lineages. Relevant analyses are in `epithelial_lineage_analysis.ipynb`.

## Part 4: Neoantigen selection
<div align="center">
  <img src="images/neoantigen.png" width="275" height="300">
</div>

This part demonstrates how to obtain tumor-recurrent neoepitopes from a novel transcript and how to optimize a panel of neoepitopes with high binding affinities to patients' alleles. Relevant analyses are in `neoantigen_candidates.ipynb`.