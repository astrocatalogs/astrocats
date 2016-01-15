#!/bin/bash
repos=('sne-pre-1990' 'sne-1990-1999' 'sne-2000-2004' 'sne-2005-2009' 'sne-2010-2014' 'sne-2015-2019')
cd ..
for repo in ${repos[@]}; do
	cd ${repo}
	pwd
	git add -u
	git commit -a -m '$1' && git push
	cd ..
done
