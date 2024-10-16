# -*- coding: utf-8 -*-

"""Module for defining the Flask app."""

from flask import Flask
from apscheduler.schedulers.background import BackgroundScheduler
from apscheduler.triggers.interval import IntervalTrigger
import atexit


app = Flask(__name__)
app.config["JOB_RESULTS"] = dict()


def clear_job_results() -> None:
    """Clear the job results."""
    app.config["JOB_RESULTS"].clear()


# set up a scheduler to run the clear_job_results function every week
scheduler = BackgroundScheduler()
scheduler.start()
scheduler.add_job(
    func=clear_job_results,
    trigger=IntervalTrigger(weeks=1),
    id="clear_job_results_job",
    name="clear job results every week",
    replace_existing=True,
)


# shut down the scheduler when exiting the app
atexit.register(lambda: scheduler.shutdown())


