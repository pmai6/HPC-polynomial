#! /bin/bash
git config --global credential.helper cache
git config --global credential.helper 'cache --timeout=3600'
git add .

echo 'Enter the commit message:'
read commitMessage

git commit -m "$commitMessage"

#echo 'Enter the name of the branch:'
#read branch

git push -u origin master

#read
