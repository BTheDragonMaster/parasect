# -*- coding: utf-8 -*-

"""Common utilities for the server."""

from enum import Enum, auto
from typing import Any, Dict, Optional


class Status(Enum):
    """The status of a response.

    :cvar Success: The response was successful.
    :cvar Warning: The response was successful, but with warnings.
    :cvar Failure: The response was not successful.
    """

    Success = auto()
    Warning = auto()
    Failure = auto()

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
        payload: Optional[Dict[str, Any]] = None,
        message: Optional[str] = None,
    ) -> None:
        """Create a new response object.

        :param status: The status of the response.
        :type status: Status
        :param payload: The payload of the response.
        :type payload: Optional[Dict[str, Any]]
        :param message: The message of the response.
        :type message: Optional[str]
        """
        self.status = status
        self.payload = payload if payload is not None else dict()
        self.message = message if message is not None else "no message provided"

    def to_dict(self) -> Dict[str, Any]:
        """Return the response as a dictionary.

        :return: The response as a dictionary.
        :rtype: Dict[str, Any]
        """
        return dict(status=str(self.status), payload=self.payload, message=self.message)
