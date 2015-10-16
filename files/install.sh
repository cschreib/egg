#!/bin/bash

# This script will download and install Ifni and its direct dependencies
# without requiring any interaction from you. Before running this script,
# you still need to install the indirect dependencies, i.e. cfitsio and
# libwcs, that are required to use the phy++ library.
#
# The script sequence is:
#
# 1) Download the phy++ library and install it on your system
# 2) Download Ifni and build it
# 3) Download the filter response curve library
#
# To use this script, simply copy the script to a directory where you keep
# your programs, make it executable, and run it. Example:
#
# cd ~/programs/
# cp ~/Downloads/install.sh .
# chmod +x install.sh
# ./install.sh
# rm install.sh
#
# You may have to edit the script if you do not use the system-default
# locations to install your libraries and codes. See below for more
# information.
#

# -----------------------------------------
# Configurable options
# -----------------------------------------
#
# Modify if needed.

# PHYPP_ROOT_DIR: the location where the phy++ library will be installed.
# Leave empty to install the library in the system default folders,
# or manually specify a directory if you do not want it to be installed
# system-wise or if you do not have root access on your computer.
# This has to be an absolute path.
#
# If you are installing in a non-system directory, set PHYPP_SYSTEM_DIR=0
# so that the script doesn't ask you for a root password.
#
# See the INSTALL file in the phypp-master directory or the doc.pdf file
# in the ifni-master directory if you get into trouble.

# PHYPP_ROOT_DIR="/opt/local"
# PHYPP_SYSTEM_DIR=0
PHYPP_ROOT_DIR=""
PHYPP_SYSTEM_DIR=1

# CFITSIO_ROOT_DIR: the location where the cfitsio library is installed
# on your computer, including headers and the compiled library.

# CFITSIO_ROOT_DIR="/usr/local/share/cfitsio"
CFITSIO_ROOT_DIR=""

# WCSLIB_ROOT_DIR: the location where the cfitsio library is installed
# on your computer, including headers and the compiled library.

# WCSLIB_ROOT_DIR="/usr/local/share/wcslib"
WCSLIB_ROOT_DIR=""

# -----------------------------------------
# The phy++ library
# -----------------------------------------

# Download and extract phy++
wget https://github.com/cschreib/phypp/archive/master.tar.gz
tar -xvzf master.tar.gz && rm master.tar.gz

# Configure it
mkdir -p phypp-master/build && cd phypp-master/build
cmake ../ -DCFITSIO_ROOT_DIR=$CFITSIO_ROOT_DIR \
          -DWCSLIB_ROOT_DIR=$WCSLIB_ROOT_DIR \
          -DCMAKE_INSTALL_PREFIX=$PHYPP_ROOT_DIR

# Install it
if [ $PHYPP_SYSTEM_DIR -eq 1 ]; then
    sudo make install
else
    make install
fi
source ~/.phypprc

# -----------------------------------------
# The filter library
# -----------------------------------------

# Download the filters
cd ../../
wget https://github.com/cschreib/filter-db/archive/master.tar.gz
tar -xvzf master.tar.gz && rm master.tar.gz
IFNI_FILTERS_PATH=`pwd`/filter-db-master

# -----------------------------------------
# Ifni
# -----------------------------------------

# Download Ifni
wget https://github.com/cschreib/ifni/archive/master.tar.gz
tar -xvzf master.tar.gz && rm master.tar.gz

# Configure it
mkdir -p ifni-master/build && cd ifni-master/build
cmake ../ -DPHYPP_ROOT_DIR=$PHYPP_ROOT_DIR

# Build it
make install
cd ../../

IFNI_PATH=`pwd`/ifni-master/bin

if [ $? -ne 0 ]; then
    echo ""
    echo ""
    echo "Oops, there was an error in the installation process."
    echo "Make sure that all the dependencies are properly installed"
    echo "and that your compiler is supported by phy++."
    echo ""
else
    echo ""
    echo ""
    echo "   ------------------------------------"
    echo "   Ifni has been successfuly installed!"
    echo "   ------------------------------------"
    echo ""
    echo "To finish the installation, just add the following two"
    echo "lines to your startup script (~/.bashrc and equivalent):"
    echo ""
    echo "export IFNI_PATH=\"$IFNI_PATH\""
    echo "export IFNI_FILTERS_PATH=\"$IFNI_FILTERS_PATH\""
    echo ""
fi
