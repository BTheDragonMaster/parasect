# Webapp

First make sure all three PARAS and PARASECT models are present in the `app/models` folder.
You can download the models from [Zenodo](https://zenodo.org/records/13165500).

## Run locally for development

### Server 

Create a local environment with conda and install server side dependencies with pip from `src/server/requirements.txt`:

```bash
conda create -n paras python=3.9
conda activate paras
pip install -r src/server/requirements.txt
```

Install the following dependencies on your machine:
* hmmer2
* muscle (v.3.8.1551)
* Redis (e.g. `brew install redis` on macOS, or run it via Docker: `docker run --rm -p 6379:6379 redis:7-alpine`)

Prediction jobs are stored in Redis rather than in-process, so a worker process (or a whole extra
server replica) can pick up a job and any other worker can serve its results. The server needs
Redis running and reachable to work at all, even for local dev with a single process. By default
it looks for Redis at `redis://localhost:6379/0`; override with the `REDIS_URL` env var if yours
runs elsewhere.

### Client

First install NPM package manager and Node.js on your device.

Then install client side dependencies with NPM from `src/client/package.json`:

```bash
cd src/client
npm install
```

### Run 

Run the server in one terminal:

```bash
bash run_server.sh
```

Run the client in another terminal:

```bash 
cd src/client
npm start
```

Visit `https://localhost:3000/` in your browser to view the app.

## Run with Docker

Run the following script to build and runt he app in a Docker container:

```bash
docker-compose -p paras up --build --force-recreate --remove-orphans -d
```

The app will be available at `https://localhost:4010/`.

This also starts a `redis` container for prediction-job state (see above) and the `paras-server`
container waits for it to report healthy before starting. Gunicorn worker/thread counts are set
via the `GUNICORN_WORKERS`/`GUNICORN_THREADS` environment variables in `docker-compose.yml`
(default: 2 workers, 4 threads each); because job state now lives in Redis rather than in the Flask
process, it's safe to raise worker count above 1.