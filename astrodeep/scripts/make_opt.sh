# Description
# -----------
#
# This script will create the SkyMaker catalogs and configuration files then launch
# SkyMaker to build the UV-to-NIR images.
#

VERSION=1.0
CATBASE=ifni-ad-0p25deg

PSFDIR=../psfs
CATDIR=../catalogs/v$VERSION
MAPDIR=../maps/v$VERSION
CATALOG=$CATDIR/$CATBASE.fits

SKYDIR=$CATDIR/skymaker

I2SKY_OPTIONS="verbose"

mkdir -p $MAPDIR

echo "-----------------------------"
echo "      Hubble WFC3 F160W      "
echo "-----------------------------"

../../bin/ifni-2skymaker $CATALOG $I2SKY_OPTIONS band=f160w \
    out=$SKYDIR/$CATBASE-f160w \
    img_dir=$MAPDIR \
    template=templates/hubble-f160w.conf

for SKYCAT in $SKYDIR/$CATBASE-f160w*.cat; do
    sky $SKYCAT -c $(dirname $SKYCAT)/$(basename $SKYCAT .cat)-sky.conf
    rm $MAPDIR/f160w/$(basename $SKYCAT .cat)-sci.list
done
