#!/bin/bash

set -o nounset
set -o errexit
set -o xtrace

MAIN_DIR={{cookiecutter.projects_base}}/{{cookiecutter.project_name}}
DATA_DIR=${MAIN_DIR}/data
