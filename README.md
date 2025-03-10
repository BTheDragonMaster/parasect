# PARAS and PARASECT

Welcome to PARAS and PARASECT: Predictive Algorithm for Resolving A-domain Specificity (featurising Enzyme and Compound in Tandem). Detect NRPS AMP-binding domains from an amino acid sequence and predict their substrate specificity profile.

You can find our preprint on [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.01.08.631717v1).

## Information on branches

The current `master` branch contains all code that is needed to reproduce our findings in the manuscript.

The `webapp` branch contains minimal code from the `master` branch that is needed to run the web application.

You can find a live version of the web application [here](https://paras.bioinformatics.nl/).

## Trained models

The trained models for PARAS and PARASECT can be found on Zenodo [here](https://zenodo.org/records/13165500).

## Summary

PARAS and PARASECT are adenylation domain (A-domain) selectivity predictors that can be used to predict the peptide scaffold of NRPS systems. From a .gbk or .fasta file, they first detect all adenylation domains, extract their active site signatures, and use these signatures to predict which substrate is selected by each A-domain. Both PARAS and PARASECT are random forest models, implemented in scikit-learn.

### PARAS

PARAS makes predictions based on a 34 amino acid active site signature of an A-domain alone. From this signature, each of the 1000 decision trees in the random forest predicts a substrate. An overall prediction is reached through majority consensus, and the proportion of trees agreeing with that consensus is returned as a confidence score. E.g., if 560 trees predict that the substrate is alanine, and the remaining 440 trees predict glycine, PARAS will return alanine as the predicted substrate with a confidence score of 0.56. <br />

PARAS can be run in two modes: normal mode, which is restricted to predictions of 34 commonly occurring substrates for which enough training data was available for model validation; and all-substrate mode, which can predict any of the 252 substrates in our dataset. Model choice depends on the specific system that is being studied.

### PARASECT

PARASECT takes two inputs to make predictions: the 34 amino acid active site signature of an A-domain, and the SMILES string of a substrate. Based on these two inputs, it predicts the likelihood that the A-domain and the substrate interact. By default, PARASECT makes 34 such predictions, one for each of the commonly occurring substrates that the model was trained on, and returns the substrate with the greatest interaction probability as the predicted selectivity. Predictions for other/additional substrates can also be provided by uploading the substrate(s) in SMILES format.<br />

PARASECT confidence scores reflect how many of the 1000 trees in the forest predict interaction between the A-domain and the substrate. E.g., if 910 trees predict that an A-domain interactions with alanine, and 850 trees predict that an A-domain interactions with glycine, PARASECT will return alanine as the top predicted substrate with a confidence score of 0.91. Note that, while PARAS confidence scores always add up to 1 in total across all substrates, PARASECT confidence scores do not (in the above example, alanine would have a confidence score of 0.91, and glycine a confidence score of 0.85).<br />

If desired, both PARAS and PARASECT can be prompted to return more than 1 substrate prediction, sorted by their confidence score in descending order. This can give a broader overview of the types of substrates that are recognised, especially for A-domains for which no confident prediction can be given.<br />

To find more information on installing, running, and re-training PARAS and PARASECT, please navigate the pages listed below.

[Installation](https://github.com/BTheDragonMaster/parasect/wiki/Installation)<br />
[Running PARAS and PARASECT](https://github.com/BTheDragonMaster/parasect/wiki/Running-PARAS-and-PARASECT)<br />
[Retraining](https://github.com/BTheDragonMaster/parasect/wiki/Retraining)<br />
[Interpreting output](https://github.com/BTheDragonMaster/parasect/wiki/Interpreting-output)<br />
