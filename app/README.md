# Webapp

## Run locally for development

### Server 

Create a local environment with conda and install server side dependencies with pip from `src/server/requirements.txt`:

```bash
conda create -n paras python=3.9
conda activate paras
pip install -r src/server/requirements.txt
```

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
python3 ./app/server/api.py
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
docker-compose up --build --force-recreate -d
```

The app will be available at `https://localhost:4010/`.