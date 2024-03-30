from flask import Blueprint, Response, request 
import time

from .common import Status, ResponseData

blueprint_submit_paras = Blueprint("submit_paras", __name__)
@blueprint_submit_paras.route("/api/submit_paras", methods=["POST"])
def submit_paras() -> Response:
    """
    Submit settings for prediction with Paras model.
    
    :return: Response
    """
    data = request.get_json()

    print(data)

    # Wait to simulate processing.
    time.sleep(2)

    msg = "Settings submitted successfully."
    return ResponseData(Status.Success, message=msg).to_dict()