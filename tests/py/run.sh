#!/bin/bash

PYTHONPATH=../..:../../build:$PYTHONPATH python3 -m unittest "$@"
