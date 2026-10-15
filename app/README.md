# Webapp

## Models

The server needs the trained PARAS/PARASECT models in `app/models`. You don't place them
there by hand as `docker-compose` runs a one-shot `paras-models` init service that downloads
them from Zenodo, checks each md5, expands them and exits; `paras-server` won't start until it
has finished. Each PARASECT model is downloaded with its matching fingerprint file from the
same record. Complete, verified bundles already on disk are reused without network access.

For local (non-Docker) development, do the same thing by hand once:

```bash
python src/server/fetch_models.py
```

That's ~214 MB of download which expands to ~4 GB on disk, because the server memory-maps
the models (`mmap=True` in `routes/submit.py`) and joblib ignores `mmap_mode` on a compressed
file. `PARAS (all substrates)` alone accounts for ~3.5 GB of that.

Configure it with environment variables (set on `paras-models` in `docker-compose.yml`):

| Variable | Default | What it does |
|---|---|---|
| `PARAS_MODELS` | `all` | Which models to fetch, comma-separated, using the same keys the API takes: `parasAllSubstrates`, `parasCommonSubstrates`, `parasect`, `parasectBacterial`. Trim it to shrink the download. |
| `ZENODO_RECORD` | `23262276` | Which Zenodo release to pin to, including PARASECT fingerprint files. |
| `MODEL_DIR` | `app/models` | Where the files go. |
| `PARAS_MODELS_KEEP_EXISTING` | unset | Keep untracked PARAS model files. Untracked PARASECT files are rejected because their matching fingerprints cannot be verified. |

`PARAS_MODELS` can't drop `parasAllSubstrates` while `routes/data_annotation.py` and
`routes/submit_domain.py` request it unconditionally.

### Model version and scikit-learn go together

A random forest is a pickle, so the model file and the scikit-learn that reads it are one unit,
not two independent choices. scikit-learn added a field to the decision-tree node dtype in 1.3,
so a model pickled before that fails to load outright on anything newer. This is not a warning:

```
ValueError: node array from the pickle has an incompatible dtype
```

Three pins keep that pairing honest, and they move together:

| Where | Pin | Why |
|---|---|---|
| `docker-compose.yml` -> `ZENODO_RECORD` | `23262276` | Models and matching fingerprint files for scikit-learn 1.8.0. |
| `server-requirements.txt` | `scikit-learn==1.8.0` | The deployment runs exactly what the models were written with. |
| `server-environment.yml` | `python=3.11` | scikit-learn 1.8.0 requires Python >= 3.11. On 3.9 the ceiling below is unreachable and pip quietly settles for 1.6.1. |

`pyproject.toml` keeps a range (`scikit-learn>=1.6,<=1.8.0`) rather than an exact pin, because
it is the library rather than the deployment. The ceiling is a prediction-quality finding; the
floor exists so the range can never resolve to a pre-1.3 scikit-learn that cannot read the models
at all. Python 3.9 users get 1.6, which is verified against the current models.

To make a mismatch impossible to ship quietly, `fetch_models.py` loads every requested model
before exiting, and fails if one won't open. Since `MultiModelLoader` loads lazily, without that
check a broken pairing stays invisible until a user submits a prediction, long after the server
has come up and reported healthy. Set `PARAS_MODELS_VERIFY=0` to skip the load check;
the PARASECT model/fingerprint checksum checks still run.

### PARASECT fingerprints belong to the model

The init service installs `fingerprints.txt` beside `model.parasect`, and
`fingerprints_bacterial.txt` beside `bacterial_model.parasect`. Both files are fetched from
the selected Zenodo record, staged before installation, and recorded together in
`.fetch_manifest.json`. The server passes the matching path explicitly for both standard
and uploaded substrates. It never uses the package's bundled fingerprints for predictions.

The manifest includes SHA-256 checksums for the expanded model and fingerprint file. This
detects a changed column order even when the number of features and file size stay the same.
The lazy loader checks the pair before loading, and predictions check that the fingerprint
file still matches the cached model. Job provenance records the fingerprint filename and checksum.

On the first deployment of this change, existing PARASECT models with older manifests are
downloaded again with their sidecars, even if they already came from record `23262276`.
Records that lack fingerprint files cannot be used for PARASECT. Stop backend workers before
refreshing the files and restart them afterward; an already-loaded classifier is cached in memory.

To update the backend while retaining Redis results, run from `app/` after pulling the changes:

```bash
docker compose build paras-server
docker compose stop paras-server
docker compose run --rm --no-deps paras-models &&
docker compose up -d --no-deps --force-recreate paras-server
```

Confirm that the init logs show the intended record and successful model verification.
Existing stored predictions are unchanged and need a fresh submission to use the corrected pair.

The CLI pins the same record, in `ZENODO_RECORD` in `parasect/core/helpers.py`. The two are
separate constants because the webapp can be re-pinned per deployment without touching the core
package, but they should normally match, and both should move together with
`src/parasect/data/model_metadata.txt`.

## Run locally for development

### Server

Create a conda environment and install hmmer, hmmer2, and MUSCLE v3.8.1551 from bioconda (the
exact versions the domain-extraction pipeline needs) plus the Python dependencies:

```bash
conda create -n paras python=3.9
conda activate paras
conda install -c conda-forge -c bioconda hmmer hmmer2 "muscle=3.8.1551"
pip install -r server-requirements.txt
```

(`server-environment.yml` sets up the same thing, but only works inside `Dockerfile.server`, as
it points pip at an absolute `/app/server-requirements.txt`, which only exists inside that
container image. So don't `conda env create -f` it directly for local dev.)

`server-requirements.txt` only lists the webapp's own dependencies (Flask, gunicorn, etc.). It
does not install the core `paras`/`parasect` package itself. Install your local checkout in
editable mode so the server runs against whatever's actually on this branch:

```bash
pip install -e ..
```

You'll also need:
* Redis (e.g. `brew install redis` on macOS, or run it via Docker: `docker run --rm -p 6379:6379 redis:7-alpine`)

Prediction jobs are stored in Redis rather than in-process, so a worker process (or a whole extra
server replica) can pick up a job and any other worker can serve its results. The server needs
Redis running and reachable to work at all, even for local dev with a single process. By default
it looks for Redis at `redis://localhost:6379/0`; override with the `REDIS_URL` env var if yours
runs elsewhere.

Create an `.env` file in `app/` (`run_server.sh` loads it) with:

```
GITHUB_TOKEN=
TURNSTILE_SECRET_KEY=
REACT_APP_TURNSTILE_SITE_KEY=
```

These are only needed for the data-annotation submission flow (Cloudflare Turnstile verification
and opening the resulting GitHub issue). Everything else works fine with them left blank. The
last one is read by the client (see below), not the server.

### Client

First install NPM package manager and Node.js on your device.

Then install client side dependencies with NPM from `src/client/package.json`:

```bash
cd src/client
npm install
```

The network graph page (sigma.js/WebGL) needs a browser with WebGL support; everything else works
in any modern browser.

### Run 

Run the server in one terminal:

```bash
cd src/server
bash run_server.sh
```

This starts the Flask dev server on port 5001 (`FLASK_ENV=development`, auto-reload on).

Run the client in another terminal:

```bash 
cd src/client
npm start
```

This starts the React dev server on port 3000 and proxies `/api/*` requests to port 5001 (see the
`proxy` field in `src/client/package.json`).

Visit `https://localhost:3000/` in your browser to view the app.

## Run with Docker

Run the following script to build and run the app in a Docker container:

```bash
docker-compose up --build --force-recreate --remove-orphans -d
```

The app will be available at `https://localhost:4010/`. Run it from `app/` as shown. The compose
file's build context for `paras-server` is the repo root (`..`), not `app/`, so it can install
the `src/parasect/` checkout you actually have on disk (whatever branch/commit that is) into the
image, rather than pulling `parasect` from GitHub's default branch.

This also starts a `redis` container for prediction-job state (see above) and the `paras-server`
container waits for it to report healthy before starting. Gunicorn worker/thread counts are set
via the `GUNICORN_WORKERS`/`GUNICORN_THREADS` environment variables in `docker-compose.yml`
(default: 2 workers, 4 threads each).
