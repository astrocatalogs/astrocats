#! /bin/bash

PATH=/opt/local/bin:/usr/local/bin:$PATH ; export PATH
LD_LIBRARY_PATH=/usr/local/lib:/opt/local/lib ; export LD_LIBRARY_PATH

cd /var/www/html/sne/sne/scripts
./import.py -u
./make-catalog.py
stamp=$(date +"%Y-%m-%d %k:%M")
./commit-and-push-repos.sh "Auto update: $stamp"
