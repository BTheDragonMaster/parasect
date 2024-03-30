from flask import Blueprint, Response, request 
import time

from .common import Status, ResponseData

blueprint_submit_parasect = Blueprint("submit_parasect", __name__)
@blueprint_submit_parasect.route("/api/submit_parasect", methods=["POST"])
def submit_parasect() -> Response:
    """
    Submit settings for prediction with Parasect model.
    
    :return: Response
    """
    data = request.get_json()

    print(data)

    # Wait to simulate processing.
    time.sleep(2)

    msg = "Not implemented."
    return ResponseData(Warning.Success, message=msg).to_dict()