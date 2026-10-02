# PARASECT

Welcome to PARASECT: Predictive Algorithm for Resolving A-domain Specificity featurising Enzyme and Compound in Tandem. Detect NRPS AMP-binding domains from an amino acid sequence and predict their substrate specificity profile.

## Web application

You can find a live version of the web application [here](https://paras.bioinformatics.nl/).

## Database

Browse the data that PARAS and PARASECT were trained on [here](https://paras.bioinformatics.nl/query_database).

## License

The PARAS/PARASECT source code is licensed under the GNU Affero General Public License v3.0 (AGPL-3.0), see [LICENSE.txt](LICENSE.txt).

The PARAS/PARASECT database (the substrate-specificity training data, including `src/parasect/data/parasect.db`, `src/parasect/data/database_files/`, `src/parasect/data/model_metadata.txt`, and `condensed/data/dataset_paras.tsv`) is licensed separately under a Creative Commons Attribution 4.0 International (CC BY 4.0) license, see [LICENSE-DATABASE.txt](LICENSE-DATABASE.txt).

## Data submission

Do you have new datapoints that you think PARAS/PARASECT could benefit from in future versions? Submit your data [here](https://paras.bioinformatics.nl/data_annotation).

## Trained models

The trained models for PARAS and PARASECT can be found on Zenodo [here](https://zenodo.org/records/18682178).

## Command line installation

To install PARAS/PARASECT on the command line, run:

```bash
conda create -n paras python=3.9
conda activate paras

pip install paras
conda install -c bioconda hmmer
conda install -c bioconda hmmer2
conda install -c bioconda muscle==3.8.1551
```

For usage instructions, see our [wiki](https://github.com/BTheDragonMaster/parasect/wiki).
Note that the command line tool will download the models from zenodo upon the first run.