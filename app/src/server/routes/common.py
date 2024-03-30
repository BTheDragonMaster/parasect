import typing as ty
from enum import Enum, auto

class Status(Enum):
    """
    The status of a response.

    :cvar Success: The response was successful.
    :cvar Warning: The response was successful, but with warnings.
    :cvar Failure: The response was not successful.
    """
    Success = auto()
    Warning = auto()
    Failure = auto()

    def __str__(self) -> str:
        return self.name.lower()
    
class ResponseData:
    def __init__(
        self, 
        status: Status, 
        payload: ty.Optional[dict] = None, 
        message: ty.Optional[str] = None 
    ) -> None:
        """
        Create a new response object.
        
        :param Status status: The status of the response.
        :param dict payload: The payload of the response.
        :param str message: The message of the response.
        """
        self.status = status
        self.payload = payload if payload is not None else dict()
        self.message = message if message is not None else "No message provided."

    def to_dict(self) -> ty.Dict[str, ty.Any]:
        """
        Return the response as a dictionary.
        
        :return: The response as a dictionary.
        :rtype: ty.Dict[str, ty.Any]
        """
        return dict(
            status=str(self.status),
            payload=self.payload,
            message=self.message
        )