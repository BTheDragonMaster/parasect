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
    package_data={"": ["AMP-binding_full.hmm*",
                       'Apositions*',
                       "*.fasta",
                       "*.txt",
                       "*.tsv",
                       "*.classifier"]},
    scripts=['bin/parasect', 'bin/paras'])