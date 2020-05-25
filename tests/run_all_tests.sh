#!/bin/bash

set -e

echo -e "\033[33mRunning Python unit tests...\033[39m"
(cd py && ./run.sh)

echo -e "\033[33mRunning C++ unit tests...\033[39m"
(cd ../build && ./libepidemics_unittests)
