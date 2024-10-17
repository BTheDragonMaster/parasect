# Training PARAS and PARASECT models used in the web application

This directory contains the training code for the PARAS and PARASECT models used in the web application.

## PARAS (common subsrrates)

```bash
python3 ./train_paras_common.py --data-dir <path-to-data-dir> --out-dir <path-to-out-dir>
```

## PARAS (all subsrrates)

```bash
python3 ./train_paras_all.py --data-dir <path-to-data-dir> --out-dir <path-to-out-dir>
```

## PARASECT

```bash
python3 ./train_parasect.py --data-dir <path-to-data-dir> --out-dir <path-to-out-dir>
```