#!/bin/bash

if [ $# -eq 0 ]
  then
    echo "No arguments supplied, exiting"
	exit
fi

git commit -a -m "$1"
git push
repos=($(awk -F= '{print $1}' rep-folders.txt))
echo ${repos[*]}
cd ..
for repo in ${repos[@]}; do
	cd ${repo}
	pwd
	git add -A
	git commit -a -m "$1"
	git push
	cd ..
done
