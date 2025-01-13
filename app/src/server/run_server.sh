#!/bin/bash

# Set environment variables
export FLASK_APP=api.py
export FLASK_ENV=development
export FLASK_DEBUG=1

# Run the Flask application
flask run