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
    SKYCONF=$(dirname $SKYCAT)/$(basename $SKYCAT .cat)-sky.conf
    sky $SKYCAT -c $SKYCONF
    egg-postskymaker conf=$SKYCONF
    rm $MAPDIR/$(basename $SKYCAT .cat)-sci.list
done


echo "-----------------------------"
echo "      Spitzer IRAC ch2       "
echo "-----------------------------"

egg-2skymaker cat=$CATALOG $I2SKY_OPTIONS band=spitzer-irac2 \
    out=$SKYDIR/$CATBASE-irac2.cat \
    img_dir=$MAPDIR \
    template=goodss-spitzer-irac2.conf

for SKYCAT in $SKYDIR/$CATBASE-irac2*.cat; do
    SKYCONF=$(dirname $SKYCAT)/$(basename $SKYCAT .cat)-sky.conf
    sky $SKYCAT -c $SKYCONF
    egg-postskymaker conf=$SKYCONF background=-0.0029
    rm $MAPDIR/$(basename $SKYCAT .cat)-sci.list
done
