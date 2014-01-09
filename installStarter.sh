#!/bin/bash
# install starterator as python package
# 	using setuptools
# add desktop file to right place

# ignore above
# make .starterator if not exists
# NO get tarball of code/all files NO
# branch bzr repository
# run starterator.sh?
sudo apt-get update
sudo apt-get install ncbi-blast+ python-pip git
sudo pip install PyPDF2

starterator_dir = ${HOME}/Applications/Starterator
#make starterator directory
mkdir $starterator_dir
cd $starterator_dir
git clone https://
# dropbox link

cd $starterator_dir
./starterator