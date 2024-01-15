from setuptools import setup, find_packages

VERSION = '0.0.1'
DESCRIPTION = 'PARASECT: Predictive Algorithm for Resolving A-domain Specificity featurising Enzyme and Compound in Tandem'
LONG_DESCRIPTION = 'Detect NRPS AMP-binding domains from an aa sequence and predict their substrate specificity profile'

setup(
    name="parasect",
    version=VERSION,
    author="Barbara Terlouw",
    author_email="barbara.terlouw@wur.nl",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    packages=find_packages(),
    package_data={"": ["*.hmm*",
                       'Apositions*',
                       "*.fasta",
                       "*.txt",
                       "*.tsv",
                       "*model.paras",
                       "*model_all_substrates.paras",
                       "*model.parasect",
                       "*model_onehot.paras",
                       "*model_onehot.parasect"
                       ]},
    scripts=['bin/parasect', 'bin/paras'])
