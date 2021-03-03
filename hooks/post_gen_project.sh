#!/usr/bin/env bash

cat <<EOT > rerun_cookiecutter.py
from cookiecutter.main import cookiecutter
import os

home = os.path.expanduser("~")
parameters={{ cookiecutter | jsonify }}

cookiecutter('git@github.com:sidbdri/cookiecutter-de_analysis_skeleton.git',
             no_input=True,
             checkout='master',
             overwrite_if_exists=False,
             output_dir=os.path.join(home, "projects"),
             extra_context=parameters
             )
EOT
