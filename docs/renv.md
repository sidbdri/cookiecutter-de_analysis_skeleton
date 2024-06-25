# renv

[renv](https://rstudio.github.io/renv/) is a R runtime package version management tool. It is used to manage R package
dependencies. It is similar to python's virtualenv. It is used to create a snapshot of the R package dependencies and
store it in a lock file. This lock file can be used to recreate the same environment in another machine. It is useful
for reproducibility of the R package environment.

We have setup renv in docker image 1.7.0, thus we now also set it up in cookiecutter, so all new project will use renv
to manage R package dependencies. A global renv cache as well as a copy of renv settings is setup on sidb at
/srv/data/renv, which store a snapshot of the R packages we used on sidb.

Two cache folder are used in this setup, a global cache and a user cache. For more details RE multiple cache, please
see [here](https://rstudio.github.io/renv/articles/package-install.html):

- global cache
  The global cache was build by Xin at /srv/data/renv/cache_sidb, which user will have read access but cannot edit.
- user cache
  The user cache is at ~/.cache/R/renv/cache, which is used to store the R packages you used in your R env, which is not
  available in the global cache.

In order to setup renv in your account, you will need to do the below steps **ONLY ONCE**:

```
# login to sidb and in a bash terminal do below:
export RENV_PATHS_CACHE="~/.cache/R/renv/cache:/srv/data/renv/cache_sidb"
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
echo export RENV_PATHS_CACHE="~/.cache/R/renv/cache:/srv/data/renv/cache_sidb" >> ~/.bashrc
## if you use zsh, also uncomment and run below
# echo export RENV_PATHS_CACHE="~/.cache/R/renv/cache:/srv/data/renv/cache_sidb" >> ~/.zshrc
```

This will create the multiple cache setup for renv so that when you run R so it knows where to look for the R packages.

# renv for Rstudio

**this is for Owen**
Rstudio somehow did not pick up the above bash environment variable, so it will be explicitly set in the Rprofile.site

```bash
cat /opt/R-4.2.0/etc/Rprofile.site
# if setting not there, add below line
echo "Sys.setenv(RENV_PATHS_CACHE='~/.cache/R/renv/cache:/srv/data/renv/cache_sidb')" | sudo tee -a /opt/R-4.2.0/etc/Rprofile.site
# might also want to make sure proxy is set there is needed
echo "Sys.setenv(http_proxy='http://wwwcache.ed.ac.uk:3128/')" | sudo tee -a /opt/R-4.2.0/etc/Rprofile.site
echo "Sys.setenv(https_proxy='http://wwwcache.ed.ac.uk:3128/')" | sudo tee -a /opt/R-4.2.0/etc/Rprofile.site
```

# renv usage

For every new project created,

- a .Rprofile file will be created in the project root folder, which will tell R to use renv for version control.
- a renv.lock file will be created in the project root folder, which will store the snapshot of the R packages used in
  the project.

The first time you run R in the project folder, you will need to run:
```R
renv::restore(prompt = FALSE)
```
This will:
- Detect the renv settings files in the root folder.
- Identify any libraries listed in the .lock file that are missing from your project library folder.

For the missing packages, renv will:

- Link the packages found in the local/global cache to the project folder.
- Install the packages not found in any cache into the local cache folder, and then link them to the project folder.

Code has been added to the [project_setup.sh] to do this, so if you run the project_setup.sh, it will do the above steps
automatically.

At the end of the project development, you will want to update the .lock file if any dependency changes are made to the project.
You can do this by using:
```R
renv::snapshot()
```
The .lock file should be committed to the project repository so that the same environment can be recreated on another machine.


The below command give you some example of manually interact with renv.

- tell renv to link/install all the packages listed in the lock file into the project library folder:
```R
renv::restore(prompt = FALSE)
```

- When you installed new packages or change any package version in the project, before you commit the code to Github, you
will need to update the lock file by running:

```R
renv::snapshot()
```

- You can manually inspect any change of your project library compared to your .lock file by running:

```R
renv::status()
```

- If you want to exclude some package from getting into the lock file:

```R
renv::restore(exclude = c("package1", "package2"))
```