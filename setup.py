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
                       "*class_sequence_hmm.paras",
                       "*class_sequence_hmm_all.paras",
                       "*class_sequence_hmm.parasect",
                       "*class_sequence_hmm_onehot.paras",
                       "*class_sequence_hmm_onehot.parasect"
                       ]},
    scripts=['bin/parasect', 'bin/paras'])
