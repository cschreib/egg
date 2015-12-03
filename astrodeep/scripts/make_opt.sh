# Description
# -----------
#
# This script will create the SkyMaker catalogs and configuration files then launch
# SkyMaker to build the UV-to-NIR images.
#

VERSION=1.0
CATBASE=egg-ad-0p25deg

CATDIR=../catalogs/v$VERSION
MAPDIR=../maps/v$VERSION
CATALOG=$CATDIR/$CATBASE.fits

SKYDIR=$CATDIR/skymaker

I2SKY_OPTIONS="verbose size_cap=0.3"

mkdir -p $MAPDIR

# Abort immediately on error
set -e

echo "-----------------------------"
echo "      Hubble WFC3 F160W      "
echo "-----------------------------"

egg-2skymaker cat=$CATALOG $I2SKY_OPTIONS band=hst-f160w \
    out=$SKYDIR/$CATBASE-f160w.cat \
    img_dir=$MAPDIR \
    template=goodss-hst-f160w.conf

for SKYCAT in $SKYDIR/$CATBASE-f160w*.cat; do
    sky $SKYCAT -c $(dirname $SKYCAT)/$(basename $SKYCAT .cat)-sky.conf
    rm $MAPDIR/$(basename $SKYCAT .cat)-sci.list
done
