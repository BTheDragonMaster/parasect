# -*- coding: utf-8 -*-

from __future__ import annotations

import json
import os
from typing import Any

import redis

JOB_TTL_SECONDS = 7 * 24 * 60 * 60  # 7 days, same retention as the old sweep

REDIS_URL = os.getenv("REDIS_URL", "redis://localhost:6379/0")

_client = redis.Redis.from_url(REDIS_URL, decode_responses=True)


def _key(job_id: str) -> str:
    return f"job:{job_id}"


def set_job(job_id: str, job: dict[str, Any]) -> None:
    """Create or fully overwrite a job's state.

    :param job_id: Job ID.
    :param job: Job state, e.g. status/message/results/timestamp.
    """
    _client.set(_key(job_id), json.dumps(job), ex=JOB_TTL_SECONDS)


def get_job(job_id: str) -> dict[str, Any] | None:
    """Retrieve a job's state.

    :param job_id: Job ID.
    :return: Job state, or None if the job does not exist (or has expired).
    """
    raw = _client.get(_key(job_id))
    if raw is None:
        return None
    return json.loads(raw)


def update_job(job_id: str, **fields: Any) -> None:
    """Update one or more fields of an existing job's state.

    Creates the job if it does not already exist, so this is also safe to
    use for the initial write.

    :param job_id: Job ID.
    :param fields: Fields to set/overwrite on the job.
    """
    job = get_job(job_id) or {}
    job.update(fields)
    set_job(job_id, job)


def ping() -> bool:
    """Check whether Redis is reachable.

    :return: True if Redis responded, False otherwise.
    """
    try:
        return bool(_client.ping())
    except redis.RedisError:
        return False
