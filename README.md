# CCIP

## Overview

CCIP is a  machine learning method for predicting CTCF-mediated chromatin loops with transitivity.

CTCF-mediated chromatin loops underlie the formation of topological associating domains (TADs) and serve as the structural basis for transcriptional regulation. However, the formation mechanism of these loops remains unclear,and the genome-wide mapping of these loops is costly and difficult. 

Motivated by recent process on the formation mechanism of CTCF-mediated loops, we studied the possibility of making use of transitivity-related information of interacting CTCF anchors to predict CTCF loops computationally. In this context,transitivity arises when two CTCF anchors interact with a same third anchor by the loop extrusion mechanism and bring themselves close to each other spatially to form an indirect loop. 

We proposed an accurate and efficient two-stage random-forest-based machine learning method, CCIP (**C**TCF-mediated **C**hromatin **I**nteraction **P**rediction), to predict CTCF-mediated chromatin loops. Our two-stage learning approach makes it possible for us to train a prediction model by taking advantage of transitivity-related information as well as functional genome data and genomic data. 

## Required Package

CCIP could be installed in a linux-like system. The CCIP requires the following dependencies. We recommend to use [Anaconda python distribution](https://www.anaconda.com/what-is-anaconda/) for installation of the below packages.

1. Python (tested 3.6.10)
2. numpy (tested 1.18.1)
3. pandas (tested 0.24.2)
4. matplotlib (tested 3.1.1)
5. networkx (tested 2.4)
6. scikit-learn (tested 0.22.1)
7. joblib (tested 0.14.1)
8. bedtools (tested 2.29.2) 

## Installation

Download CCIP by

```shell
git clone https://github.com/gaolabXDU/CCIP
```

## Usage

Required data for GM12878, HeLa-S3, K562 and MCF-7 cell line are available in the CCIP/data directory. CTCF motif and CTCF age data are common for these cell lines while CTCF ChIP-seq and RAD21 ChIP-seq data are specific for these cell lines. you can use the shell script CCIP/code/test.sh to test the software. There are three main scripts for generating samples (generate_pairs.py), extracting features (generate_features.py) and training the predicting models (rf_graph.py).

The script generate_pairs.py is used to generate positive and negative samples for CCIP.

Output: pair_all_balance.csv

```
Usage: generate_pairs.py [options]

Options:
  -h|--help:            show this help message and exit
  -o|--output_path:     Path for output
  -c|--ctcf_file:       CTCF ChIP-seq data
  -m|--ctcf_motif_file: CTCF motif occurence data
  -p|--chia_pet_file:   CTCF ChIA-PET datahao
```

The script generate_features.py is used to extract features for samples generated from last step.

Output: samples.csv

```
Usage: generate_features.py [options]

Options:
  -h|--help:            show this help message and exit
  -o|--output_path:     Path for output
  -c|--ctcf_file:       CTCF ChIP-seq data
  -m|--ctcf_motif_file: CTCF motif occurence data
  -p|--chia_pet_file:   CTCF ChIA-PET data
  -r|--rad21_file:      RAD21 ChIP-seq data
  -a|--age_file:        CTCF age data
```

The script rf_graph.py is used to train the model and do ten fold cross validation.

Output: rf_base.model, rf_graph.model, cross_val_predict_prob.npy

```
Usage: rf_graph.py [options]

Options:
  -h|--help:            show this help message and exit
  -i|--input_file:      Samples for training the model
  -o|--output_path:     Output path for store the training results
```

The script rf_graph_chrom_cv.py is the cross chromosome validation version of rf_graph.py.

Output: rf_base.model, rf_graph.model, cross_val_predict_prob.npy

```
Usage: rf_graph_chrom_cv.py [options]

Options:
  -h|--help:            show this help message and exit
  -i|--input_file:      Samples for training the model
  -o|--output_path:     Output path for store the training results
```

The script rf_graph_test.py is used to predict the samples from one cell type using the trained model from another cell type.

Output: %s_%s_ccip_prob.npy (predicted probability of each sample)

```
Usage: rf_graph_test.py [options]

Options:
  -h|--help：            show this help message and exit
  -o|--output_path：     Output path for storing the testing results
  -m|--model_path：      Model file for predicting
  -s|--sample_file：     Sample file for predicting
  -M|--model_cell：      Model cell for predicting
  -S|--sample_cell：     Sample cell for predicting
```

For example, in CCIP/code/, we can run these scripts for GM12878 cell line:

```bash
python generate_pairs.py -c ../data/GM12878/CTCF_peak.bed\
                         -m ../data/GM12878/fimo.csv\
                         -p ../data/GM12878/gm12878_ctcf.interactions.intra.bedpe\
                         -o ../data/GM12878/output
python generate_features.py -c ../data/GM12878/CTCF_peak.bed\
                            -m ../data/GM12878/fimo.csv\
                            -p ../data/GM12878/gm12878_ctcf.interactions.intra.bedpe\
                            -o ../data/GM12878/output\
                            -r ../data/GM12878/rad21.narrowPeak\
                            -a ../data/GM12878/CTCF_age.bed
python rf_graph.py -i ../data/GM12878/output/sample.csv -o ../data/GM12878/output/model
```
