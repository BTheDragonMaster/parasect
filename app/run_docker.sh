#!/bin/bash

# Check if at least two arguments are provided.
if [ $# -lt 2 ]; then
  echo "Usage: $0 <SERVER_PORT> <CLIENT_PORT> [ -d ]"
  exit 1
fi

# Set the environment variables.
export SERVER_PORT=$1
export CLIENT_PORT=$2

# Check if the optional '-d' flag is provided as the third argument.
if [ $# -eq 3 ] && [ "$3" == "-d" ]; then
  # Run Docker Compose with the updated environment variables in detached mode.
  SERVER_PORT=$1 CLIENT_PORT=$2 docker-compose up --build --force-recreate -d
else
  # Run Docker Compose with the updated environment variables in the foreground.
  SERVER_PORT=$1 CLIENT_PORT=$2 docker-compose up --build --force-recreate
fi