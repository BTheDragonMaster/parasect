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

The trained models for PARAS and PARASECT can be found on Zenodo [here](https://zenodo.org/records/23262276).

PARASECT models must be used with the fingerprint file from the same release:
`model.parasect.gz` with `fingerprints.txt`, and `bacterial_model.parasect.gz` with
`fingerprints_bacterial.txt`. The order of fingerprint columns is part of the model.
The CLI downloads and checks both files together in the model directory; legacy caches
without pairing metadata are refreshed automatically. If a different scikit-learn version
requires local retraining, the new fingerprints are saved alongside that trained model.

Python callers must pass the matching file explicitly to `run_parasect()`:

```python
from parasect.api import run_parasect
from parasect.core.helpers import prepare_model, load_parasect_bundle
from parasect.core.models import ModelType

model_path = prepare_model(ModelType.PARASECT, "models")
model, fingerprints_file = load_parasect_bundle(model_path)
# input_text contains FASTA data; temp_dir must be an existing directory.
results = run_parasect(input_text, "fasta", temp_dir, model,
                      fingerprints_file=fingerprints_file)
```

For the bacterial model, use `ModelType.PARASECT_BACTERIAL` and pass
`bacterial_only=True` to `run_parasect()`. Callers supplying their own classifier must
supply its matching fingerprint file; prediction no longer falls back to package data.

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
