#!/usr/bin/env bash

# If this is the first time the script has been run, We find out the current
# master hash of transcript-utils and sidbdri-utils and put them in the
# config.sh, if they are not set. We also set up the directory as a git
# repository.
echo "##########################"
echo "# Running project setup. This script runs only once after cookiecutter init the project and will be self deleted afterwards."


function findHashFromBranchName {
    org=$1
    repo_name=$2
    commit_hash=$3
    echo $(curl --silent -H "Accept: application/vnd.github.VERSION.sha" \
    https://api.github.com/repos/${org}/${repo_name}/commits/${commit_hash})
}

if grep -Fq "unknown_hash" config.sh; then
    echo "Replace the unknown hash with the master branch hash"
    sed -i "s/unknown_hash_transcript-utils/$(findHashFromBranchName "sidbdri" "transcript-utils" "master")/" config.sh
    sed -i "s/unknown_hash_sidbdri-utils/$(findHashFromBranchName "sidbdri" "sidbdri-utils" "master")/" config.sh
    {% if cookiecutter.sargasso == "yes" %}
    sed -i "s/unknown_hash_sargasso/$(findHashFromBranchName "biomedicalinformaticsgroup" "Sargasso" "master")/" config.sh
    {% endif %}

fi

# Generate species specific de script
echo "Create species-specific R script"
{% for s in cookiecutter.species.split(' ') %}
cp diff_expr.R  diff_expr_{{ s }}.R
sed 's/unknown_species/{{ s }}/' ./meta_data.R > meta_data_{{ s }}.R
sed 's/unknown_species/{{ s }}/' ./diff_expr.R > diff_expr_{{ s }}.R
sed 's/unknown_species/{{ s }}/' ./diff_expr_tx.R > diff_expr_tx_{{ s }}.R
sed 's/unknown_species/{{ s }}/' ./rMATS.R > rMATS_{{ s }}.R
{% endfor %}

# remove the template
rm ./diff_expr.R ./meta_data.R ./rMATS.R  ./diff_expr_tx.R

echo "Init git"
git init
mv gitignore .gitignore

## delete myself so it only run again
echo "Delete setup script"
rm -- "$0"
#echo "Rename self to project_setup.sh.bk"
#mv "$0" "$0.bk"

echo "# project setup finish."
echo "##########################"
