
# one might have used/installed R packages in the project that are not in the docker image,
# for example, one installed a package called 'packageA' when developing the project,
# this package will be install under the project directory ./renv/library/,
# but will not automatically be recorded in the renv.lock file.
#
# Thus before commit project to github for a docker run, one needs to make sure that
# any R package that is used in the project is recorded in the renv.lock file.
# This can be done as below:


# check if any package is installed in the project library but not in the renv.lock file
R -e 'renv::status()'

# sync the lock file if there is any package that is not in the lock file
R -e 'renv::snapshot()'

# comfirm that the lock file is updated, below commander show
# "No issues found -- the project is in a consistent state."
R -e 'renv::status()'

# make sure the renv related files are included in the commit, these include below 4 files:
# renv.lock, renv/activate.R, renv/settings.json, .Rprofile