#!/bin/bash

# Load environment variables from .env in project root
set -a
source "$(dirname "$0")/../../.env"
set +a

# Set environment variables
export FLASK_APP=api.py
export FLASK_ENV=development
export FLASK_DEBUG=1

# Run the Flask application
flask run --port=5001