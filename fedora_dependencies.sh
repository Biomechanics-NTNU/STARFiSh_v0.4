#!/bin/sh
echo "Installation of system dependencies via dnf.."
dnf -y install gcc
dnf -y install gcc-gfortran
dnf -y install gcc-c++
dnf -y install python-devel
dnf -y install pygtk2-devel
dnf -y install python2-matplotlib
dnf -y install graphviz
dnf -y install python-scipy
dnf -y install python-h5py
dnf -y install python-lxml
dnf -y install python-psutil
echo "Installing python dependencies via pip"
pip install -r requirements.txt
