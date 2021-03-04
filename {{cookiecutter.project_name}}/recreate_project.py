#!/usr/bin/env python

from cookiecutter.main import cookiecutter
import os

HOME = os.path.expanduser("~")
BRANCH='master'
OVERWRITE_IF_EXISTS=True
OUT=os.path.join(HOME, "projects")

PARAMETERS={{ cookiecutter | jsonify }}

cookiecutter('git@github.com:sidbdri/cookiecutter-de_analysis_skeleton.git',
             no_input=True,
             checkout=BRANCH,
             overwrite_if_exists=OVERWRITE_IF_EXISTS,
             output_dir=OUT,
             extra_context=PARAMETERS
             )