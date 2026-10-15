"""Exercise release orchestration with real Git/Compose and a simulated Docker engine.

No containers or models are needed. Run: python3 -m unittest discover -s tests/deployment -v
"""

import importlib.util
import contextlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("paras_deploy", ROOT / "app" / "deploy.py")
deploy = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(deploy)
REAL_RUN = subprocess.run


def has_compose():
    if not shutil.which("docker"):
        return False
    return REAL_RUN(["docker", "compose", "version"], stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL).returncode == 0


class DockerEngine:
    """Run Compose's real config parser; simulate only engine-dependent commands."""

    def __init__(self):
        self.images = {"redis:7-alpine": "sha256:redis"}
        self.commands = []
        self.fail_start = False
        self.build_sources = []

    def run(self, args, **kwargs):
        if args[0] != "docker":
            return REAL_RUN(args, **kwargs)
        self.commands.append(args)
        code, output = 0, ""
        if args[1:3] == ["image", "inspect"]:
            image = args[-1]
            code = 0 if image in self.images else 1
            output = self.images.get(image, "")
        elif args[1] == "compose":
            if "config" in args:
                return REAL_RUN(args, **kwargs)
            if "build" in args:
                commit = kwargs["env"]["PARAS_IMAGE_TAG"]
                file = Path(args[args.index("--file") + 1])
                self.build_sources.append((file.parent.parent / "payload.txt").read_text())
                for image in ("paras-server", "paras-client"):
                    self.images[f"{image}:{commit}"] = f"sha256:{image}-{commit}"
            elif "up" in args:
                code = 1 if self.fail_start else 0
            else:
                raise AssertionError(f"Unexpected Compose command: {args}")
        else:
            raise AssertionError(f"Unexpected Docker command: {args}")
        if code and kwargs.get("check"):
            raise subprocess.CalledProcessError(code, args)
        return subprocess.CompletedProcess(args, code, output, "")


@unittest.skipUnless(has_compose(), "Docker Compose v2 CLI is needed (no running engine required)")
class DeploymentTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.repo = self.root / "repo"
        (self.repo / "app").mkdir(parents=True)
        shutil.copyfile(ROOT / "app" / "docker-compose.yml", self.repo / "app" / "docker-compose.yml")
        (self.repo / "payload.txt").write_text("committed version")
        self.git("init", "-q")
        self.git("add", ".")
        self.git("-c", "user.name=Deployment Test", "-c", "user.email=test@example.invalid",
                 "-c", "commit.gpgsign=false", "commit", "-qm", "First release")
        self.commit = self.git("rev-parse", "HEAD")
        self.config_dir = self.root / "config"
        self.config_dir.mkdir()
        for env in ("dev", "prod"):
            shutil.copyfile(ROOT / "app" / "deployment" / f"{env}.env.example",
                            self.config_dir / f".env.{env}")
        self.state_dir = self.root / "state"
        self.state_dir.mkdir()
        self.deployment = deploy.Deployment(self.repo, self.config_dir, self.state_dir)
        self.engine = DockerEngine()
        self.patch = patch.object(deploy.subprocess, "run", self.engine.run)
        self.patch.start()
        self.addCleanup(self.patch.stop)
        output = contextlib.redirect_stdout(io.StringIO())
        output.__enter__()
        self.addCleanup(output.__exit__, None, None, None)

    def git(self, *args):
        return REAL_RUN(["git", "-C", str(self.repo), *args], check=True, text=True,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout.strip()

    def build(self):
        self.deployment.deploy_dev("HEAD")
        return self.deployment.release(self.git("rev-parse", "HEAD"))

    def test_promotion_reuses_committed_images_and_frozen_model_pin(self):
        (self.repo / "payload.txt").write_text("uncommitted edit must not ship")
        directory, manifest = self.build()
        self.assertEqual(self.engine.build_sources, ["committed version"])
        self.assertTrue(manifest["dev_healthy"])
        # Later checkout/config edits cannot change the tested release recipe.
        (self.repo / "app" / "docker-compose.yml").write_text("invalid: compose")
        with (self.config_dir / ".env.prod").open("a") as file:
            file.write("\nZENODO_RECORD=999\nPARAS_MODELS=parasect\nGITHUB_TOKEN=private-test-value\n")
        self.engine.commands.clear()
        self.deployment.promote()
        self.assertFalse(any("build" in command for command in self.engine.commands))
        command = next(command for command in self.engine.commands if "up" in command)
        self.assertIn("--no-build", command)
        self.assertEqual(command[command.index("--pull") + 1], "never")
        config = self.deployment.config(directory / "compose.yml", "prod", self.commit, manifest["models"])
        services = config["services"]
        self.assertEqual(services["paras-models"]["environment"]["ZENODO_RECORD"], "18682178")
        self.assertEqual(services["paras-server"]["environment"]["GITHUB_TOKEN"], "private-test-value")
        self.assertNotIn("private-test-value", (directory / "release.json").read_text())
        self.assertEqual(self.deployment.state("prod")["release"], self.commit)

    def test_isolation_and_shell_exports_cannot_override_selected_environment(self):
        file = self.repo / "app" / "docker-compose.yml"
        with patch.dict(os.environ, {"COMPOSE_PROJECT_NAME": "wrong", "PARAS_CLIENT_PORT": "9999",
                                     "PARAS_MODELS_DIR": "/shared", "PARAS_IMAGE_TAG": "latest"}):
            dev = self.deployment.config(file, "dev", self.commit)
            prod = self.deployment.config(file, "prod", self.commit)
        for environment, config, project, port in (("dev", dev, "paras-dev", "4011"),
                                                    ("prod", prod, "paras", "4010")):
            services = config["services"]
            self.assertEqual(services["redis"]["container_name"], f"{project}-redis")
            self.assertEqual(str(services["paras-client"]["ports"][0]["published"]), port)
            self.assertEqual(services["paras-server"]["environment"]["FLASK_ENV"], "production")
            source = services["paras-server"]["volumes"][0]["source"]
            self.assertIn(f"/models/{environment}/", source)
        self.assertNotEqual(dev["networks"], prod["networks"])

    def test_refuses_retagged_image_before_changing_production(self):
        self.build()
        self.engine.images[f"paras-server:{self.commit}"] = "sha256:unexpected"
        with self.assertRaisesRegex(deploy.DeploymentError, "changed since"):
            self.deployment.promote()
        self.assertEqual(self.deployment.state("prod"), {})

    def test_failed_dev_cannot_be_promoted_and_retry_does_not_rebuild(self):
        self.engine.fail_start = True
        with self.assertRaises(subprocess.CalledProcessError):
            self.build()
        self.assertEqual(self.deployment.state("dev")["status"], "failed")
        with self.assertRaisesRegex(deploy.DeploymentError, "healthy dev"):
            self.deployment.promote(self.commit)
        self.engine.fail_start = False
        self.build()
        self.assertEqual(len(self.engine.build_sources), 1)

    def test_site_key_mismatch_and_port_collision_are_rejected(self):
        self.build()
        with (self.config_dir / ".env.prod").open("a") as file:
            file.write("\nREACT_APP_TURNSTILE_SITE_KEY=different\n")
        with self.assertRaisesRegex(deploy.DeploymentError, "site key differs"):
            self.deployment.promote()
        with (self.config_dir / ".env.prod").open("a") as file:
            file.write("\nREACT_APP_TURNSTILE_SITE_KEY=\nPARAS_CLIENT_PORT=4011\n")
        with self.assertRaisesRegex(deploy.DeploymentError, "different PARAS_CLIENT_PORT"):
            self.deployment.promote()
        self.assertEqual(self.deployment.state("prod"), {})

    def test_rollback_after_failed_production_update_restores_previous_release(self):
        self.build()
        self.deployment.promote()
        first = self.commit
        (self.repo / "payload.txt").write_text("second release")
        self.git("add", ".")
        self.git("-c", "user.name=Deployment Test", "-c", "user.email=test@example.invalid",
                 "-c", "commit.gpgsign=false", "commit", "-qm", "Second release")
        self.build()
        self.engine.fail_start = True
        with self.assertRaises(subprocess.CalledProcessError):
            self.deployment.promote()
        self.assertEqual(self.deployment.state("prod")["previous"], first)
        self.engine.fail_start = False
        self.engine.commands.clear()
        self.deployment.rollback()
        self.assertEqual(self.deployment.state("prod")["release"], first)
        self.assertEqual(self.deployment.state("prod")["status"], "healthy")
        self.assertNotEqual(self.deployment.state("prod")["previous"], first)
        self.assertFalse(any("build" in command for command in self.engine.commands))

    def test_unrecorded_images_are_not_overwritten(self):
        self.engine.images[f"paras-server:{self.commit}"] = "sha256:other-build"
        with self.assertRaisesRegex(deploy.DeploymentError, "Unrecorded image"):
            self.build()
        self.assertEqual(self.engine.build_sources, [])

    def test_commit_before_environment_support_is_rejected(self):
        file = self.repo / "app" / "docker-compose.yml"
        file.write_text(file.read_text().replace("${COMPOSE_PROJECT_NAME:-paras}-", "paras-"))
        with self.assertRaisesRegex(deploy.DeploymentError, "predates"):
            self.deployment.config(file, "dev", self.commit)


if __name__ == "__main__":
    unittest.main()
