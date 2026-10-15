#!/usr/bin/env python3
"""Deploy to dev, then promote the same local Docker images. Python stdlib only."""

import argparse
import fcntl
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile


APP = Path(__file__).resolve().parent
PROJECTS = {"dev": "paras-dev", "prod": "paras"}


class DeploymentError(Exception):
    pass


def run(*args, capture=False, **kwargs):
    result = subprocess.run(
        [str(arg) for arg in args], check=True, text=True,
        stdout=subprocess.PIPE if capture else None, **kwargs,
    )
    return result.stdout.strip() if capture else None


def read_json(path):
    return json.loads(path.read_text())


def write_json(path, value):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


class Deployment:
    def __init__(self, repo, config_dir, state_dir):
        self.repo = repo.resolve()
        self.config_dir = config_dir.resolve()
        self.state_dir = state_dir.resolve()

    def state(self, environment):
        path = self.state_dir / f"{environment}.json"
        return read_json(path) if path.exists() else {}

    def compose(self, file, environment, commit, models=None):
        config = self.config_dir / f".env.{environment}"
        if not config.is_file():
            raise DeploymentError(f"Missing {config}; copy the deployment/*.env.example first.")
        # Shell exports must not silently override a selected environment file.
        env = {
            key: value for key, value in os.environ.items()
            if not key.startswith(("COMPOSE_", "PARAS_", "GUNICORN_"))
            and key not in {"ZENODO_RECORD", "REACT_APP_TURNSTILE_SITE_KEY"}
        }
        env.update({
            "COMPOSE_PROJECT_NAME": PROJECTS[environment],
            "PARAS_IMAGE_TAG": commit,
            "PARAS_ENV_FILE": str(config),
            "PARAS_MODELS_DIR": str(
                self.state_dir / "models" / environment /
                (models["ZENODO_RECORD"] if models else "pending")
            ),
        })
        if models:
            env.update(models)
        command = [
            "docker", "compose", "--project-name", PROJECTS[environment],
            "--env-file", str(config), "--file", str(file),
        ]
        return command, env

    def compose_run(self, file, environment, commit, *args, models=None, capture=False):
        command, env = self.compose(file, environment, commit, models)
        return run(*command, *args, env=env, capture=capture)

    def config(self, file, environment, commit, models=None):
        # Do not persist or print the resolved config: it contains runtime secrets.
        config = json.loads(self.compose_run(
            file, environment, commit, "config", "--format", "json",
            models=models, capture=True,
        ))
        services = config["services"]
        expected = {
            "redis": "redis", "paras-models": "models-init",
            "paras-server": "server", "paras-client": "client",
        }
        for service, suffix in expected.items():
            if services[service].get("container_name") != f"{PROJECTS[environment]}-{suffix}":
                raise DeploymentError("This commit predates the dev/prod Compose configuration.")
        for service, image in (("paras-models", "paras-server"), ("paras-server", "paras-server"),
                               ("paras-client", "paras-client")):
            if services[service]["image"] != f"{image}:{commit}":
                raise DeploymentError("Compose must use the selected release's image tags.")
        return config

    def image_id(self, image):
        return run("docker", "image", "inspect", "--format", "{{.Id}}", image, capture=True)

    def verify_images(self, manifest):
        for image, expected in manifest["images"].items():
            if self.image_id(image) != expected:
                raise DeploymentError(f"Image {image} changed since it was built; refusing deployment.")

    def release(self, commit):
        if not re.fullmatch(r"[0-9a-f]{40}", commit):
            raise DeploymentError("Use the full commit ID printed by the dev deployment.")
        directory = self.state_dir / "releases" / commit
        if not (directory / "release.json").is_file():
            raise DeploymentError(f"Release {commit} has not been built on this server.")
        return directory, read_json(directory / "release.json")

    def deploy_dev(self, ref):
        commit = run("git", "-C", self.repo, "rev-parse", "--verify", "--end-of-options",
                     f"{ref}^{{commit}}", capture=True)
        directory = self.state_dir / "releases" / commit
        if not (directory / "release.json").exists():
            with tempfile.TemporaryDirectory(prefix="paras-build-") as temporary:
                source = Path(temporary) / "source"
                source.mkdir()
                archive = Path(temporary) / "source.tar"
                run("git", "-C", self.repo, "archive", "--format=tar", "-o", archive, commit)
                run("tar", "-xf", archive, "-C", source)
                file = source / "app" / "docker-compose.yml"
                config = self.config(file, "dev", commit)
                services = config["services"]
                model_env = services["paras-models"]["environment"]
                models = {key: model_env[key] for key in ("ZENODO_RECORD", "PARAS_MODELS")}
                if not re.fullmatch(r"[0-9]+", models["ZENODO_RECORD"]):
                    raise DeploymentError("ZENODO_RECORD must be a numeric Zenodo record ID.")
                images = [services[name]["image"] for name in ("paras-server", "paras-client")]
                for image in images:
                    existing = subprocess.run(
                        ["docker", "image", "inspect", image],
                        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                    )
                    if existing.returncode == 0:
                        raise DeploymentError(
                            f"Unrecorded image {image} already exists. Do not overwrite release tags; "
                            "use a new commit, or remove an unused image left by a failed build."
                        )
                print(f"Building committed source {commit} (working-tree edits are excluded).", flush=True)
                self.compose_run(file, "dev", commit, "build", "paras-server", "paras-client")
                manifest = {
                    "commit": commit, "models": models,
                    "site_key": services["paras-client"]["build"]["args"]["REACT_APP_TURNSTILE_SITE_KEY"] or "",
                    "images": {image: self.image_id(image) for image in images},
                    "dev_healthy": False,
                }
                directory.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(file, directory / "compose.yml")
                write_json(directory / "release.json", manifest)
        self.activate("dev", commit)
        directory, manifest = self.release(commit)
        manifest["dev_healthy"] = True
        write_json(directory / "release.json", manifest)
        command = ["python3", "app/deploy.py"]
        if self.config_dir != APP:
            command.extend(["--config-dir", str(self.config_dir)])
        if self.state_dir != APP / ".deploy":
            command.extend(["--state-dir", str(self.state_dir)])
        command.extend(["promote", commit])
        print(f"Dev is healthy at release {commit}. Test the website, then run: {shlex.join(command)}")

    def activate(self, environment, commit):
        directory, manifest = self.release(commit)
        file = directory / "compose.yml"
        self.verify_images(manifest)
        models = manifest["models"]
        config = self.config(file, environment, commit, models)
        services = config["services"]
        site_key = services["paras-client"]["build"]["args"]["REACT_APP_TURNSTILE_SITE_KEY"] or ""
        if site_key != manifest["site_key"]:
            raise DeploymentError("Turnstile site key differs from the tested image. Use the same widget in both environments.")
        other = "prod" if environment == "dev" else "dev"
        if (self.config_dir / f".env.{other}").exists():
            other_config = self.config(file, other, commit, models)
            port = services["paras-client"]["ports"][0]["published"]
            other_port = other_config["services"]["paras-client"]["ports"][0]["published"]
            if port == other_port:
                raise DeploymentError("Dev and production must use different PARAS_CLIENT_PORT values.")
        # Create bind mounts as the deployment user, not implicitly as root via
        # the Docker daemon. The model initializer runs as container UID 1000.
        (self.state_dir / "models" / environment / models["ZENODO_RECORD"]).mkdir(
            parents=True, exist_ok=True,
        )
        # Redis is infrastructure, not rebuilt with an application release. Pull it
        # only on first installation, preserving an existing Redis container/jobs.
        redis = services["redis"]["image"]
        inspected = subprocess.run(["docker", "image", "inspect", redis],
                                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if inspected.returncode:
            run("docker", "pull", redis)
        old = self.state(environment)
        previous = old.get("release") if old.get("status") == "healthy" else old.get("previous")
        if previous == commit:
            previous = old.get("previous")
        if previous == commit:
            previous = None
        state = {"release": commit, "previous": previous, "status": "deploying"}
        path = self.state_dir / f"{environment}.json"
        write_json(path, state)
        try:
            # Build paths in the saved Compose file are deliberately unused.
            # --force-recreate and down would also discard Redis job results.
            self.compose_run(file, environment, commit, "up", "--detach", "--no-build",
                             "--pull", "never", "--wait", "--wait-timeout", "900", models=models)
        except (subprocess.CalledProcessError, KeyboardInterrupt):
            state["status"] = "failed"
            write_json(path, state)
            raise
        state["status"] = "healthy"
        write_json(path, state)

    def promote(self, commit=None):
        if commit is None:
            dev = self.state("dev")
            if dev.get("status") != "healthy":
                raise DeploymentError("There is no healthy current dev release to promote.")
            commit = dev["release"]
        _, manifest = self.release(commit)
        if not manifest["dev_healthy"]:
            raise DeploymentError("This release has not completed a healthy dev deployment.")
        self.activate("prod", commit)
        print(f"Production is healthy at release {commit}.")

    def rollback(self):
        commit = self.state("prod").get("previous")
        if not commit:
            raise DeploymentError("No previous production release recorded by this helper.")
        self.promote(commit)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config-dir", type=Path, default=APP, help="Directory containing .env.dev and .env.prod")
    parser.add_argument("--state-dir", type=Path, default=APP / ".deploy", help="Persistent releases and model storage")
    commands = parser.add_subparsers(dest="action", required=True)
    dev = commands.add_parser("dev", help="Build a committed Git ref once and deploy it to dev")
    dev.add_argument("ref", nargs="?", default="HEAD")
    promote = commands.add_parser("promote", help="Promote a dev release without rebuilding; invocation is approval")
    promote.add_argument("release", nargs="?", help="Full commit ID; defaults to the last healthy dev release")
    commands.add_parser("rollback", help="Restore the previous recorded production release")
    commands.add_parser("status", help="Show recorded deployment status (not a live health check)")
    args = parser.parse_args()
    deploy = Deployment(APP.parent, args.config_dir, args.state_dir)
    try:
        deploy.state_dir.mkdir(parents=True, exist_ok=True)
        with (deploy.state_dir / "deploy.lock").open("w") as lock:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                raise DeploymentError("Another deployment is already running.")
            if args.action == "status":
                print(json.dumps({env: deploy.state(env) for env in PROJECTS}, indent=2))
            else:
                run("docker", "info", "--format", "{{.ServerVersion}}", capture=True)
                if args.action == "dev":
                    deploy.deploy_dev(args.ref)
                elif args.action == "promote":
                    deploy.promote(args.release)
                else:
                    deploy.rollback()
    except (DeploymentError, subprocess.CalledProcessError, OSError) as error:
        print(f"Deployment failed: {error}", file=sys.stderr)
        print("Check status and container logs. After a failed production update, use rollback if a previous release exists.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
