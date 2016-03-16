#!/bin/bash

source /root/.bashrc

if [ $# -eq 0 ]
  then
    echo "No arguments supplied, exiting"
	exit
fi

git pull
git commit -a -m "$1"
git push
repos=($(awk -F= '{print $1}' rep-folders.txt))
echo ${repos[*]}
cd ..
for repo in ${repos[@]}; do
	cd ${repo}
	pwd
	git pull
	git add -A
	git commit -a -m "$1"
	git lfs push origin master
	git push
	cd ..
done
