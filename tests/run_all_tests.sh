#!/bin/bash

set -e

echo -e "\e[33mRunning Python unit tests...\e[39m"
(cd py && ./run.sh)

echo -e "\e[33mRunning C++ unit tests...\e[39m"
(cd ../build && ./libepidemics_unittests)
