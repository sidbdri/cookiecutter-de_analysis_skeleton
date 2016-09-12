#!/bin/bash

export WORKON_HOME={{cookiecutter.virtualenv_home}}
export PROJECT_HOME={{cookiecutter.projects_base}}
source /usr/local/bin/virtualenvwrapper.sh

mkproject -f {{cookiecutter.project_name}}

pip install multiqc
