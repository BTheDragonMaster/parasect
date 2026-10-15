# Deploy to dev, then promote to production

This setup runs both sites on **one Docker server**, using local images. It needs
Git, Python 3.9+, and Docker Compose v2 with `up --wait` support. Run every command
below on that server from the repository root. GitHub pushes still only run tests;
deployment is an explicit command.

| Environment | Compose project | Default host port | Hostname |
|---|---|---|---|
| Dev | `paras-dev` | 4011 | `paras-dev.bioinformatics.nl` |
| Production | `paras` | 4010 | `paras.bioinformatics.nl` |

## One-time configuration

Commit the deployment changes and fetch/check out that commit on the server first.
The helper builds committed source with `git archive`; local edits, ignored files,
and untracked files are excluded. It does not switch branches or change your checkout.

```bash
cp app/deployment/dev.env.example app/.env.dev
cp app/deployment/prod.env.example app/.env.prod
chmod 600 app/.env.dev app/.env.prod
```

Edit both files. Copy the existing production secrets from `app/.env` into
`app/.env.prod`. Both files are gitignored. Keep dev's `GITHUB_TOKEN` empty unless
you intentionally want annotation tests to open real issues in the PARASECT repo.
Those submissions will fail without a token; predictions do not need these secrets.

Use the same `REACT_APP_TURNSTILE_SITE_KEY` in both files and authorize both hostnames
on that Cloudflare Turnstile widget. Use its corresponding `TURNSTILE_SECRET_KEY` in
both environments. The public site key is compiled into the frontend image; the
helper rejects a different key during promotion. Both keys can be empty if annotation
submission is not needed. Both public environments run Flask with debug disabled.

Have the server administrator verify DNS, TLS certificates, and the reverse proxy:

```text
https://paras.bioinformatics.nl      -> http://127.0.0.1:4010
https://paras-dev.bioinformatics.nl  -> http://127.0.0.1:4011
```

The ports are configurable in the environment files. `127.0.0.1` assumes the reverse
proxy runs directly on the Docker host. If the proxy is remote or containerized,
set `PARAS_BIND_ADDRESS` to an address it can reach and use that address in its upstream.
The domain alone does not select a Compose deployment: the proxy must route the two
hostnames to different ports. DNS/proxy/TLS configuration is outside this repo.

The helper keeps releases, deployment history, and downloaded models in `app/.deploy/`
(gitignored). Preserve this directory and the Docker images. It creates model paths
per environment and Zenodo record, so dev cannot replace production models and a
model-version rollback can reuse its previous files. Allow roughly 4 GB of model
storage per environment per model version. The model initializer uses UID 1000;
the deployment account should have that UID, or grant that UID write access to the
model directories. Docker Desktop handles host/container file ownership differently.

For storage outside the checkout, pass the **same** `--state-dir /srv/paras-deploy`
on every command. `--config-dir /srv/paras-config` can similarly select an external
directory containing `.env.dev` and `.env.prod`. Put these options before the subcommand.

## Everyday deployment

Fetch changes, then deploy a specific commit or branch to dev:

```bash
git fetch origin
python3 app/deploy.py dev origin/dev
```

You can use `HEAD`, a feature branch, a tag, or a commit ID instead of `origin/dev`.
This does not require adopting a new branch scheme. Only deploy commits that include
the environment-aware Compose configuration.

The helper builds `paras-server:<full-commit>` and `paras-client:<full-commit>` once,
saves the Compose file, model pins, public site key, and image IDs, and starts dev.
It waits for the model initialization, backend/Redis health, and frontend proxy health.
Existing releases reuse their saved images and model pins; use a new commit for a new
build, including changes to build-time configuration or model selection.

Open `https://paras-dev.bioinformatics.nl` and check a representative prediction,
the results page, database search, and whatever changed. Also check `/api/submit_quick`
if applicable: its redirect should stay on the dev hostname. `/health` checks the
frontend-to-backend connection and Redis; `/api/version` includes the deployed commit.
Health checks do not replace functional testing.

After testing, copy the full commit ID printed by the helper:

```bash
python3 app/deploy.py promote FULL_COMMIT_ID
```

Calling `promote` is the approval step. It checks that the release previously started
successfully on dev, verifies that the local image IDs still match, then deploys those
images with production's secrets and ports. It uses the saved Compose file and model
pins, with `--no-build --pull never`. The app images are never rebuilt or pulled during
promotion. Redis is pulled only if its image is missing on the server.

Omitting the commit ID promotes the most recent healthy dev deployment. Using an
explicit ID is useful when dev has moved on since your testing. Never overwrite the
release image tags or prune images that you may need for rollback.

## Status, failures, and rollback

```bash
python3 app/deploy.py status
python3 app/deploy.py rollback
```

Status shows recorded outcomes, including failed attempts; it is not a live health
check. Rollback restores the previous production release recorded by the helper,
including its Compose file, images, and model pin. To restore another older release,
use `promote FULL_COMMIT_ID`. A failed startup is not automatically rolled back;
inspect `docker logs paras-server` / `docker logs paras-client` (or the `paras-dev-`
equivalents), fix configuration or explicitly roll back.

Releases do not snapshot secrets or other environment-specific runtime settings.
If a failure came from changing these, restore the environment file too. If a build
failed before it saved a release and left an image tag behind, the next attempt
refuses to overwrite it. Remove only that unused partial-build image, or build a
new commit. The helper serializes deployments that use the same state directory.

## First production migration

The production project and container names remain `paras` / `paras-redis`,
`paras-server`, `paras-client`, and `paras-models-init`. This avoids a second stack
claiming port 4010 and avoids renaming/recreating Redis solely for this migration.
Verify the existing server actually uses these names before the first promotion.

Save the existing production commit/configuration and images before that first
promotion: the helper has no rollback history for deployments it did not create.
The new model directory is populated on first deployment. Production updates briefly
restart the app, and in-flight prediction threads can be interrupted. Redis holds
completed results only in memory: do not use `docker compose down`, `--force-recreate`
for the whole stack, or remove/recreate Redis to deploy an app update.

After migration, use this helper for deployments. Bare `docker compose up --build`
still supports the old local workflow, but uses `latest` and does not track releases.

If the two hostnames are on different servers, this local-image workflow needs an
image registry or explicit image and release-metadata transfer before promotion.
It does not deploy over SSH or change either public server from your laptop.
