# renv

[renv](https://rstudio.github.io/renv/) is a R runtime package version management tool. It is used to manage R package
dependencies. It is similar to python's virtualenv. It is used to create a snapshot of the R package dependencies and
store it in a lock file. This lock file can be used to recreate the same environment in another machine. It is useful
for reproducibility of the R package environment.

We have setup renv in docker image 1.7.0, thus we now also set it up in cookiecutter, so all new project will use renv
to manage R package dependencies. A global renv cache as well as a copy of renv settings is setup on sidb at
/srv/data/renv, which store a snapshot of the R
packages we used on sidb.

Two cache folder are used in this setup, a global cache and a user cache. For more details RE multiple cache, please
see [here](https://rstudio.github.io/renv/articles/package-install.html):

- global cache
  The global cache was build by Xin from one of the
  docker [image](https://github.com/sidbdri/Docker/blob/7a16fcb20c47b6ce40585858d1005fe11a9c46c6/sidb_R/build_R_env_cache.sh)
  R env and is at /srv/data/renv/cache, which user will
  have read access but cannot edit.
- user cache
  The user cache is at ~/.cache/R/renv/cache, which is used to store the R packages you used in your R env, which is not
  available in the global cache.
-

In order to setup renv in your account, you will need to do the below steps **ONLY ONCE**:

```
# login to sidb and in a bash terminal do below:
export RENV_PATHS_CACHE="~/.cache/R/renv/cache:/srv/data/renv/cache"
mkdir /tmp/renvsetup_${USER} && cd /tmp/renvsetup_${USER}
R -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'));"
R -e "renv::settings\$snapshot.type('all');renv::init(bioconductor = TRUE, restart = TRUE);"
rm -rf /tmp/renvsetup_${USER}
# you should see some packages listed from below command
ls  ~/.cache/R/renv/cache/*/*/*/
```

The above steps will snapshot all the R packages you have installed in your sidb R env and store it in your local cache,
if the package is not available in the global cache.

In order to have renv manage your R package dependencies in your project. you will need to do the below steps **ONLY
ONCE**:

```
echo export RENV_PATHS_CACHE="~/.cache/R/renv/cache:/srv/data/renv/cache" >> ~/.bashrc
## if you use zsh, also uncomment and run below
# echo export RENV_PATHS_CACHE="~/.cache/R/renv/cache:/srv/data/renv/cache" >> ~/.zshrc
```

This will create the multiple cache setup for renv so that when you run R, it knows where to look for the R packages.

