# -*- coding: utf-8 -*-

"""Constants used throughout the server package."""

import os
import shutil
import time

MODEL_DIR_LOCAL = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))),
    "models",
)
MODEL_DIR = os.getenv("MODEL_DIR", MODEL_DIR_LOCAL)
TEMP_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "temp")

DB_PATH_LOCAL = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "parasect.db",
)
DB_PATH = os.getenv("SQLITE_PATH", DB_PATH_LOCAL)


STALE_JOB_TEMP_DIR_SECONDS = 2 * 60 * 60  # 2 hours; normal jobs finish in seconds to minutes


def _sweep_stale_job_temp_dirs(max_age_seconds: int = STALE_JOB_TEMP_DIR_SECONDS) -> None:
    """Remove per-job temp subdirectories older than ``max_age_seconds``.

    Job temp dirs are always cleaned up in a ``finally`` block when a job
    finishes normally (success or a handled failure), but if a worker is
    killed abruptly mid-job (OOM, SIGKILL, a native crash in an hmmer/muscle
    subprocess) that ``finally`` never runs and the directory is orphaned.
    Normal jobs finish in seconds to minutes, so anything this old is
    orphaned, not just slow.

    Called from :func:`job_temp_dir`, i.e. opportunistically whenever a new
    job is about to start, rather than from a dedicated background thread:
    it's a cheap listdir + a few stat calls (there should rarely be more
    than a handful of entries), and piggybacking on work that's already
    about to happen means there's no separate scheduler/thread lifecycle to
    reason about.

    :param max_age_seconds: Age, in seconds, past which a job temp
        directory is considered orphaned.
    """
    if not os.path.isdir(TEMP_DIR):
        return
    now = time.time()
    for name in os.listdir(TEMP_DIR):
        path = os.path.join(TEMP_DIR, name)
        if os.path.isdir(path) and (now - os.path.getmtime(path)) > max_age_seconds:
            shutil.rmtree(path, ignore_errors=True)


def job_temp_dir(job_id: str) -> str:
    """Create (if needed) and return a temp directory scoped to a single job.

    The prediction pipeline (HMM/MUSCLE calls in ``parasect.core``) writes
    intermediate files under fixed names within whatever temp directory it is
    given. Handing every job its own subdirectory of ``TEMP_DIR`` keeps
    concurrent jobs (multiple gunicorn threads/workers processing requests
    at the same time) from overwriting each other's intermediate files.

    Also sweeps stale, orphaned job temp dirs left behind by past crashes
    (see :func:`_sweep_stale_job_temp_dirs`) before creating this one.

    :param job_id: Job ID; used as the subdirectory name.
    :type job_id: str
    :return: Path to the job-specific temp directory.
    :rtype: str
    """
    _sweep_stale_job_temp_dirs()

    path = os.path.join(TEMP_DIR, job_id)
    os.makedirs(path, exist_ok=True)
    return path


def cleanup_job_temp_dir(job_id: str) -> None:
    """Remove a job's temp directory and everything in it, if it exists.

    :param job_id: Job ID.
    :type job_id: str
    """
    shutil.rmtree(os.path.join(TEMP_DIR, job_id), ignore_errors=True)