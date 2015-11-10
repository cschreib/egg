# Description
# -----------
#
# This script will create the Spitzer MIPS and Herschel images
# corresponding to a given mock catalog, using the same observing
# conditions as in the GOODS-South field (same image depths).
#


VERSION=1.0
CATBASE=ifni-ad-0p25deg

PSFDIR=../psfs
MAPDIR=../maps/v$VERSION
CATDIR=../catalogs/v$VERSION
CATALOG=$CATDIR/$CATBASE.fits

NOISE_OPTIONS="verbose seed=42"
MAP_OPTIONS="verbose"


echo "-----------------------------"
echo "      Spitzer IRS 16um       "
echo "-----------------------------"
../../bin/ifni-gennoise $CATALOG $NOISE_OPTIONS \
    out=$MAPDIR/$CATBASE-spitzer-irs1 psf=$PSFDIR/spitzer-irs1.fits \
    aspix=0.9 rms=1.6e-2

../../bin/ifni-genmap   $CATALOG $MAP_OPTIONS \
    out=$MAPDIR/$CATBASE-spitzer-irs1-sci.fits psf=$PSFDIR/spitzer-irs1.fits \
    noise_map=$MAPDIR/$CATBASE-spitzer-irs1-noise.fits band=irs1 \
    flux_factor=21.8835


echo "-----------------------------"
echo "      Spitzer MIPS 24um      "
echo "-----------------------------"
../../bin/ifni-gennoise $CATALOG $NOISE_OPTIONS \
    out=$MAPDIR/$CATBASE-spitzer-m1 psf=$PSFDIR/spitzer-m1.fits \
    aspix=1.2 rms=3.4e-2

../../bin/ifni-genmap   $CATALOG $MAP_OPTIONS \
    out=$MAPDIR/$CATBASE-spitzer-m1-sci.fits psf=$PSFDIR/spitzer-m1.fits \
    noise_map=$MAPDIR/$CATBASE-spitzer-m1-noise.fits band=m1 \
    flux_factor=7.78838


echo "-----------------------------"
echo "     Herschel PACS 70um      "
echo "-----------------------------"
../../bin/ifni-gennoise $CATALOG $NOISE_OPTIONS \
    out=$MAPDIR/$CATBASE-herschel-p1 psf=$PSFDIR/herschel-p1.fits \
    aspix=1.2 rms=1.95e-5

../../bin/ifni-genmap   $CATALOG $MAP_OPTIONS \
    out=$MAPDIR/$CATBASE-herschel-p1-sci.fits psf=$PSFDIR/herschel-p1.fits \
    noise_map=$MAPDIR/$CATBASE-herschel-p1-noise.fits band=p1 \
    flux_factor=1.452e6


echo "-----------------------------"
echo "     Herschel PACS 100um     "
echo "-----------------------------"
../../bin/ifni-gennoise $CATALOG $NOISE_OPTIONS \
    out=$MAPDIR/$CATBASE-herschel-p2 psf=$PSFDIR/herschel-p2.fits \
    aspix=1.2 rms=8.23e-6

../../bin/ifni-genmap   $CATALOG $MAP_OPTIONS \
    out=$MAPDIR/$CATBASE-herschel-p2-sci.fits psf=$PSFDIR/herschel-p2.fits \
    noise_map=$MAPDIR/$CATBASE-herschel-p2-noise.fits band=p2 \
    flux_factor=1.672e6


echo "-----------------------------"
echo "     Herschel PACS 160um     "
echo "-----------------------------"
../../bin/ifni-gennoise $CATALOG $NOISE_OPTIONS \
    out=$MAPDIR/$CATBASE-herschel-p3 psf=$PSFDIR/herschel-p3.fits \
    aspix=2.4 rms=1.68e-5

../../bin/ifni-genmap   $CATALOG $MAP_OPTIONS \
    out=$MAPDIR/$CATBASE-herschel-p3-sci.fits psf=$PSFDIR/herschel-p3.fits \
    noise_map=$MAPDIR/$CATBASE-herschel-p3-noise.fits band=p3 \
    flux_factor=1.613e6


echo "-----------------------------"
echo "    Herschel SPIRE 250um     "
echo "-----------------------------"
../../bin/ifni-gennoise $CATALOG $NOISE_OPTIONS \
    out=$MAPDIR/$CATBASE-herschel-s1 psf=$PSFDIR/herschel-s1.fits \
    aspix=3.6 rms=1.66e-3

../../bin/ifni-genmap   $CATALOG $MAP_OPTIONS \
    out=$MAPDIR/$CATBASE-herschel-s1-sci.fits psf=$PSFDIR/herschel-s1.fits \
    noise_map=$MAPDIR/$CATBASE-herschel-s1-noise.fits band=s1 \
    flux_factor=1e6


echo "-----------------------------"
echo "    Herschel SPIRE 350um     "
echo "-----------------------------"
../../bin/ifni-gennoise $CATALOG $NOISE_OPTIONS \
    out=$MAPDIR/$CATBASE-herschel-s2 psf=$PSFDIR/herschel-s2.fits \
    aspix=4.8 rms=1.66e-3

../../bin/ifni-genmap   $CATALOG $MAP_OPTIONS \
    out=$MAPDIR/$CATBASE-herschel-s2-sci.fits psf=$PSFDIR/herschel-s2.fits \
    noise_map=$MAPDIR/$CATBASE-herschel-s2-noise.fits band=s2 \
    flux_factor=1e6


echo "-----------------------------"
echo "    Herschel SPIRE 500um     "
echo "-----------------------------"
../../bin/ifni-gennoise $CATALOG $NOISE_OPTIONS \
    out=$MAPDIR/$CATBASE-herschel-s3 psf=$PSFDIR/herschel-s3.fits \
    aspix=7.2 rms=1.88e-3

../../bin/ifni-genmap   $CATALOG $MAP_OPTIONS \
    out=$MAPDIR/$CATBASE-herschel-s3-sci.fits psf=$PSFDIR/herschel-s3.fits \
    noise_map=$MAPDIR/$CATBASE-herschel-s3-noise.fits band=s3 \
    flux_factor=1e6
