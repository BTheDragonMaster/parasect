# -*- coding: utf-8 -*-

"""Common utilities for the server."""

from __future__ import annotations

from enum import Enum, auto
from typing import Any


class Status(Enum):
    """The status of a response.

    :cvar Success: The response was successful.
    :cvar Warning: The response was successful, but with warnings.
    :cvar Failure: The response was not successful.
    """

    Success = auto()
    Warning = auto()
    Failure = auto()
    Pending = auto()

    def __str__(self) -> str:
        """Return the string representation of the status.

        :return: The string representation of the status.
        :rtype: str
        """
        return self.name.lower()


class ResponseData:
    """A response object for the server."""

    def __init__(
        self,
        status: Status,
        payload: dict[str, Any] | None = None,
        message: str | None = None,
    ) -> None:
        """Create a new response object.

        :param status: The status of the response.
        :type status: Status
        :param payload: The payload of the response.
        :type payload: dict[str, Any] | None
        :param message: The message of the response.
        :type message: str | None
        """
        self.status = status
        self.payload = payload if payload is not None else dict()
        self.message = message if message is not None else "no message provided"

    def to_dict(self) -> dict[str, Any]:
        """Return the response as a dictionary.

        :return: The response as a dictionary.
        :rtype: dict[str, Any]
        """
        return dict(status=str(self.status), payload=self.payload, message=self.message)
