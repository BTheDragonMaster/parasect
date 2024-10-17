import json
import requests

url = "https://paras.bioinformatics.nl/api/submit_signature"

data = {
    "data": {
        "submissions": [
            {
                "protein_name": "SomeGeneName",
                "domain_start": 0,
                "domain_end": 1000,
                "extended_signature": "LAKAFDAFVAEGILISAGEVNAYGPTEVTVCATQ",
            },
            # add more submissions here for multiple domains ...
            # you can add multiple domains with the same protein_name
        ],
    }
}

headers = {"Content-Type": "application/json"}

# request will immediately return jobId
response = requests.post(url, data=json.dumps(data), headers=headers).json()

print(response) # response contains jobId

# in order to view the results, you can open the following link in the browser 
# with the jobId:
# https://paras.bioinformatics.nl/results/<jobId>
# replace <jobId> with the job ID

