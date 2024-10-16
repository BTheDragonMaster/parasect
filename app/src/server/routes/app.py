# -*- coding: utf-8 -*-

"""Module for defining the Flask app."""

import atexit
import time

from flask import Flask
from apscheduler.schedulers.background import BackgroundScheduler
from apscheduler.triggers.interval import IntervalTrigger


app = Flask(__name__)
app.config["JOB_RESULTS"] = dict()


def clear_job_results() -> None:
    """Clear the job results."""
    for job in list(app.config["JOB_RESULTS"]):

        # find out how long the job result has been stored
        timestamp = app.config["JOB_RESULTS"][job]["timestamp"]
        current_time = int(time.time())
        time_stored = current_time - timestamp

        # if stored longer than 7 days, delete the job result
        if time_stored >= 604800:  # 7 days in seconds
            del app.config["JOB_RESULTS"][job]


# set up a scheduler to run the clear_job_results function every week
scheduler = BackgroundScheduler()
scheduler.start()
scheduler.add_job(
    func=clear_job_results,
    trigger=IntervalTrigger(days=1),
    id="clear_job_results_job",
    name="clear job results every week",
    replace_existing=True,
)


# shut down the scheduler when exiting the app
atexit.register(lambda: scheduler.shutdown())


