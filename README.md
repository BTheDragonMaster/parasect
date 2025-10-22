# PARASECT

Welcome to PARASECT: Predictive Algorithm for Resolving A-domain Specificity featurising Enzyme and Compound in Tandem. Detect NRPS AMP-binding domains from an amino acid sequence and predict their substrate specificity profile.

## Web application

You can find a live version of the web application [here](https://paras.bioinformatics.nl/).

## Database

Browse the data that PARAS and PARASECT were trained on [here](https://paras.bioinformatics.nl/query_database).

## Data submission

Do you have new datapoints that you think PARAS/PARASECT could benefit from in future versions? Submit your data [here](https://paras.bioinformatics.nl/data_annotation).

## Trained models

The trained models for PARAS and PARASECT can be found on Zenodo [here](https://zenodo.org/records/17224548).

## Command line installation

To install PARAS/PARASECT on the command line, run:

```angular2html
conda create -n paras python=3.9
conda activate paras

pip install paras
conda install -c bioconda hmmer
conda install -c bioconda hmmer2
conda install -c bioconda muscle==3.8.1551
```

For usage instructions, see our [wiki](https://github.com/BTheDragonMaster/parasect/wiki).
Note that the command line tool will download the models from zenodo upon the first run.