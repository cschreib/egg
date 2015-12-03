#!/bin/bash

# This script will download and install EGG and its direct dependencies
# without requiring any interaction from you. Before running this script,
# you still need to install the indirect dependencies, i.e. cfitsio and
# libwcs, that are required to use the phy++ library.
#
# The script sequence is:
#
# 1) Download the phy++ library and install it on your system
# 2) Download EGG and build it
#
# To use this script, just make it executable, and run it. Your current
# directory does not matter and will not be modified. Example:
#
# chmod +x install.sh
# ./install.sh
#
# You may have to edit the script if you don't want to (or cannot) use
# the system-default locations to install your libraries and codes. See
# below for more information.


# -----------------------------------------
# Configurable options
# -----------------------------------------
#
# Modify if needed.

# INSTALL_ROOT_DIR: the location where the phy++ library and EGG will be
# installed. Leave empty to install them in the system default folders.
# You should manually specify this directory only if you do not want it
# to be installed system-wise or if you do not have root access on your
# computer. In any case this has to be an absolute path.
#
# See the INSTALL file in the phypp-master directory or the doc.pdf file
# in the egg-master directory if you get into trouble.

# Example:
# INSTALL_ROOT_DIR="/opt/local"
# Default: (system default)
INSTALL_ROOT_DIR=""
PHYPP_VERSION="master"
EGG_VERSION="1.0rc1"

# CFITSIO_ROOT_DIR: the location where the cfitsio library is installed
# on your computer, including headers (in the "include" subdirectory)
# and the compiled library (in the "lib" subdirectory).

# Example:
# CFITSIO_ROOT_DIR="/opt/local"
# Default: (system default)
CFITSIO_ROOT_DIR=""

# WCSLIB_ROOT_DIR: the location where the cfitsio library is installed
# on your computer, including headers (in the "include" subdirectory)
# and the compiled library (in the "lib" subdirectory).

# Example:
# WCSLIB_ROOT_DIR="/opt/local"
# Default: (system default)
WCSLIB_ROOT_DIR=""


# -----------------------------------------
# Prepare installation
# -----------------------------------------

function abort {
    echo ""
    echo ""
    echo "Oops, there was an error in the installation process."
    echo "Make sure that all the dependencies are properly installed"
    echo "and that your compiler is supported by phy++."
    echo ""
    exit 1
}

trap 'abort' 0
set -e

cd `mktemp -d 2>/dev/null || mktemp -d -t 'egg-tmp-dir'`

# -----------------------------------------
# The phy++ library
# -----------------------------------------

# Download and extract phy++
wget https://github.com/cschreib/phypp/archive/$PHYPP_VERSION.tar.gz \
    --no-check-certificate -O $PHYPP_VERSION.tar.gz
tar -xvzf $PHYPP_VERSION.tar.gz && rm $PHYPP_VERSION.tar.gz

if [ -n "$CFITSIO_ROOT_DIR" ]; then
    DCFITSIO_ROOT_DIR="-DCFITSIO_ROOT_DIR=$CFITSIO_ROOT_DIR"
fi
if [ -n "$WCSLIB_ROOT_DIR" ]; then
    DWCSLIB_ROOT_DIR="-DWCSLIB_ROOT_DIR=$WCSLIB_ROOT_DIR"
fi
if [ -n "$INSTALL_ROOT_DIR" ]; then
    DINSTALL_ROOT_DIR="-DCMAKE_INSTALL_PREFIX=$INSTALL_ROOT_DIR"
    DPHYPP_ROOT_DIR="-DPHYPP_ROOT_DIR=$INSTALL_ROOT_DIR"
fi

# Configure it
mkdir -p phypp-$PHYPP_VERSION/build && cd phypp-$PHYPP_VERSION/build
cmake ../ $DCFITSIO_ROOT_DIR $DWCSLIB_ROOT_DIR $DINSTALL_ROOT_DIR

# Extract install dir from CMake to check if we need sudo
if [ -z "$INSTALL_ROOT_DIR" ]; then
    INSTALL_ROOT_DIR=`cat CMakeCache.txt | grep -Po "(?<=CMAKE_INSTALL_PREFIX:PATH=).*$"`
fi

mkdir -p $INSTALL_ROOT_DIR

# Build and install phy++
make

if [ -w "$INSTALL_ROOT_DIR" ]; then
    make install
else
    sudo make install
fi
cd ../../


# -----------------------------------------
# EGG
# -----------------------------------------

# Download EGG
wget https://github.com/cschreib/egg/archive/$EGG_VERSION.tar.gz \
    --no-check-certificate -O $EGG_VERSION.tar.gz
tar -xvzf $EGG_VERSION.tar.gz && rm $EGG_VERSION.tar.gz

# Configure it
mkdir -p egg-$EGG_VERSION/build && cd egg-$EGG_VERSION/build
cmake ../ $DCFITSIO_ROOT_DIR $DWCSLIB_ROOT_DIR $DPHYPP_ROOT_DIR $DINSTALL_ROOT_DIR

# Build and install EGG
make
if [ -w "$INSTALL_ROOT_DIR" ]; then
    make install
else
    sudo make install
fi
cd ../../


# -----------------------------------------
# End of install, you made it!
# -----------------------------------------

trap : 0

echo ""
echo ""
echo "   -----------------------------------"
echo "   EGG has been successfuly installed!"
echo "   -----------------------------------"
echo ""
