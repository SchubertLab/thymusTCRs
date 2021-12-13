# Re-Investigation of Thymic development with advanced T cell Analysis tools
This Repo contains the student project for the lecture "Computational methods for single-cell biology" (Winter Semester 21/22)

## Introduction
Single-cell sequencing of T cell repertoires allows us to study their development in the thymus and further provides us insights into the response of the adaptive immune system towards diseases. While  transcriptomic data reveals the cell state, the T cell receptor (TCR) sequence indicates the cell’s functionality. Even though a joint examination of both modalities should therefore provide a more holistic view, studies on repertoires containing both data sources are often analysed for each modality, individually. Here, we will use advanced analysis tools, including our recently developed deep learning model for multiview embedding of T cells mvTCR, to re-investigate a dataset from human cells during thymic development.


Goals:    Analysis of scRNA and TCR data on a multi-modal dataset
          Showcase application of mvTCR (Deep Learning model for multiview embedding of T cells)


Data set:    ~10,000 scRNA and TCR sequenced T cells during human thymus development with donor, age, and cell-type annotation


## Supervisor:    
- Felix Drost (felix.drost@helmholtz-muenchen.de)
- Benjamin Schubert (benjamin.schubert@helmholtz-muenchen.de)

## Students:
- Coffey, Aubrey (aubray.coffey@tum.de
- Luo, Chengshi (chengshi.luo@tum.de)
- Pascu, Petru (petru.pascu@tum.de

## Methods:   
- Deep autoencoding networks
- Single-cell analysis (Scanpy)
- TCR analysis (Scirpy, IEDB tools)

## Literature:  
- Park, J.E., R.A. Botting, C. Dom ́ınguez Conde, D.M. Popescu, M. Lavaert, D.J. Kunz, I. Goh, E. Stephenson, R. Ragazzini, E. Tuck, et al. 2020. A cell atlas of human thymic development defines T cell repertoire formation. Science. 367.
- An Y, Drost F, Theis F, Schubert B, Lotfollahi M. Jointly Learning T-cell Receptor and Transcriptomic Information to Decipher the Immune Response. bioRxiv (2021).
- Sturm, G. Tamas, GS, ..., Finotello, F. Scirpy: A Scanpy extension for analyzing single-cell T-cell receptor sequencing data. Bioinformatics (2020).
- Wolf, F., Angerer, P. & Theis, F. SCANPY: large-scale single-cell gene expression data analysis. Genome Biol 19, 15 (2018).
- Luecken, M., Theis, F. Current best practices in single-cell RNA-seq analysis: a tutorial. Mol Syst Biol 15 (2019).

## Setup:
```
conda create --name thymusTCRs python=3.8
conda activate thymusTCRs
conda install nb_conda_kernels
pip install -r requirements.txt
```
