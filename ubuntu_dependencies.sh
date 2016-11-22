#!/bin/sh
echo  "Installing system dependencies via apt-get.."
apt-get -y install python-pip
apt-get -y install build-essential
apt-get -y install libxml2-dev
apt-get -y install libxslt-dev
apt-get -y install python-dev
apt-get -y install python-gtk2
apt-get -y install python-scipy
apt-get -y install python-matplotlib
apt-get -y install graphviz
apt-get -y install libhdf5-dev
apt-get -y install python-h5py
echo "Installing python dependencies via pip"
pip install -r requirements.txt
