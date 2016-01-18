#!/bin/bash

readarray repos < rep-folders.txt
echo ${repos[*]}
cd ..
for repo in ${repos[@]}; do
	cd ${repo}
	pwd
	git pull
	cd ..
done
