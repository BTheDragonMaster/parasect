# Train PARAS for common substrates

This directory contains a highly condensed version of the code base for training the PARAS model for common substrates.

## Train the PARAS model for 34 common substrates, and perform a simple test:
```bash
python3 ./paras_condensed.py --data-dir <path-to-data-dir> --out-dir <path-to-out-dir>
```

## Submit a job to the PARAS web application:
```bash
python3 ./paras_outlink.py
```