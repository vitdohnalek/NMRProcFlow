#!/bin/bash
# =============================================================================
# Run unit tests inside the NMRProcFlow Docker container.
# Uses base-R only — no testthat dependency needed.
#
# Usage:
#   bash tests/docker_run_tests.sh             # run all tests
#   bash tests/docker_run_tests.sh test_utils  # run a specific test file
#
# Prerequisites:
#   - Docker image nmrprocflow/nmrprocflow:latest must exist
#   - Run from the project root directory
# =============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
IMAGE="nmrprocflow/nmrprocflow:latest"

if [ -n "$1" ]; then
    echo "Running tests/$1 in Docker..."
    sudo docker run --rm \
        -v "${PROJECT_DIR}:/project" \
        -w /project \
        "$IMAGE" \
        Rscript tests/run_tests.R "tests/${1}"
else
    echo "Running all tests in Docker..."
    sudo docker run --rm \
        -v "${PROJECT_DIR}:/project" \
        -w /project \
        "$IMAGE" \
        Rscript tests/run_tests.R
fi
