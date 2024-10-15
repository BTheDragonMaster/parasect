import json

import requests

url = "http://localhost:4000/api/submit_signature"

data = {
    "data": {
        "submissions": [
            {
                "protein_name": "dptA",
                "domain_start": 100,
                "domain_end": 200,
                "signature": "YOURMOM",
                "extended_signature": "LDQIFDVFVSEMSLIVGGEVNAYGPTETTVEATA",
            },
            {
                "protein_name": "dptA",
                "domain_start": 0,
                "domain_end": 100,
                "signature": "YOURMOM",
                "extended_signature": "LDQIFDVFVGEMSLIVGGEVNAYGPTETTVEATA",
            },
            {
                "protein_name": "dptB",
                "domain_start": 500,
                "domain_end": 600,
                "signature": "YOURMOM",
                "extended_signature": "LDQIFDVFVSEMSLIVGGEVNAYGPTETTVEATA",
            },
        ],
    }
}

headers = {"Content-Type": "application/json"}

response = requests.post(url, data=json.dumps(data), headers=headers).json()

job_id = response["payload"]["jobId"]

print(job_id)
