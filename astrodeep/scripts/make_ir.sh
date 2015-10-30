# Description
# -----------
#
# This script will create the Spitzer MIPS and Herschel images
# corresponding to a given mock catalog, using the same observing
# conditions as in the GOODS-South field (same image depths).
#


VERSION=1.0
CATBASE=ifni-ad-1deg

CATALOG=../catalogs/v$VERSION/$CATBASE.fits
MAPDIR=../maps/v$VERSION
PSFDIR=../psfs

NOISE_OPTIONS="verbose seed=42"
MAP_OPTIONS="verbose"


echo "-----------------------------"
echo "      Spitzer IRS 16um       "
echo "-----------------------------"
../../bin/ifni-gennoise $CATALOG $NOISE_OPTIONS \
    out=$MAPDIR/$CATBASE-spitzer-irs1 psf=$PSFDIR/spitzer-irs1.fits \
    aspix=0.9 rms=1.6e-2

../../bin/ifni-genmap   $CATALOG $MAP_OPTIONS \
    out=$MAPDIR/$CATBASE-spitzer-irs1_sci.fits psf=$PSFDIR/spitzer-irs1.fits \
    noise_map=$MAPDIR/spitzer-irs1_noise.fits band=irs1 \
    flux_factor=21.8835


echo "-----------------------------"
echo "      Spitzer MIPS 24um      "
echo "-----------------------------"
../../bin/ifni-gennoise $CATALOG $NOISE_OPTIONS \
    out=$MAPDIR/$CATBASE-spitzer-m1 psf=$PSFDIR/spitzer-m1.fits \
    aspix=1.2 rms=3.4e-2

../../bin/ifni-genmap   $CATALOG $MAP_OPTIONS \
    out=$MAPDIR/spitzer-m1_sci.fits psf=$PSFDIR/spitzer-m1.fits \
    noise_map=$MAPDIR/spitzer-m1_noise.fits band=m1 \
    flux_factor=7.78838


echo "-----------------------------"
echo "     Herschel PACS 70um      "
echo "-----------------------------"
../../bin/ifni-gennoise $CATALOG $NOISE_OPTIONS \
    out=$MAPDIR/$CATBASE-herschel-p1 psf=$PSFDIR/herschel-p1.fits \
    aspix=1.2 rms=1.95e-5

../../bin/ifni-genmap   $CATALOG $MAP_OPTIONS \
    out=$MAPDIR/herschel-p1_sci.fits psf=$PSFDIR/herschel-p1.fits \
    noise_map=$MAPDIR/herschel-p1_noise.fits band=p1 \
    flux_factor=1.452e6


echo "-----------------------------"
echo "     Herschel PACS 100um     "
echo "-----------------------------"
../../bin/ifni-gennoise $CATALOG $NOISE_OPTIONS \
    out=$MAPDIR/$CATBASE-herschel-p2 psf=$PSFDIR/herschel-p2.fits \
    aspix=1.2 rms=8.23e-6

../../bin/ifni-genmap   $CATALOG $MAP_OPTIONS \
    out=$MAPDIR/herschel-p2_sci.fits psf=$PSFDIR/herschel-p2.fits \
    noise_map=$MAPDIR/herschel-p2_noise.fits band=p2 \
    flux_factor=1.672e6


echo "-----------------------------"
echo "     Herschel PACS 160um     "
echo "-----------------------------"
../../bin/ifni-gennoise $CATALOG $NOISE_OPTIONS \
    out=$MAPDIR/$CATBASE-herschel-p3 psf=$PSFDIR/herschel-p3.fits \
    aspix=2.4 rms=1.68e-5

../../bin/ifni-genmap   $CATALOG $MAP_OPTIONS \
    out=$MAPDIR/herschel-p3_sci.fits psf=$PSFDIR/herschel-p3.fits \
    noise_map=$MAPDIR/herschel-p3_noise.fits band=p3 \
    flux_factor=1.613e6


echo "-----------------------------"
echo "    Herschel SPIRE 250um     "
echo "-----------------------------"
../../bin/ifni-gennoise $CATALOG $NOISE_OPTIONS \
    out=$MAPDIR/$CATBASE-herschel-s1 psf=$PSFDIR/herschel-s1.fits \
    aspix=3.6 rms=9.73e-4

../../bin/ifni-genmap   $CATALOG $MAP_OPTIONS \
    out=$MAPDIR/herschel-s1_sci.fits psf=$PSFDIR/herschel-s1.fits \
    noise_map=$MAPDIR/herschel-s1_noise.fits band=s1 \
    flux_factor=1e6


echo "-----------------------------"
echo "    Herschel SPIRE 350um     "
echo "-----------------------------"
../../bin/ifni-gennoise $CATALOG $NOISE_OPTIONS \
    out=$MAPDIR/$CATBASE-herschel-s2 psf=$PSFDIR/herschel-s2.fits \
    aspix=4.8 rms=9.39e-4

../../bin/ifni-genmap   $CATALOG $MAP_OPTIONS \
    out=$MAPDIR/herschel-s2_sci.fits psf=$PSFDIR/herschel-s2.fits \
    noise_map=$MAPDIR/herschel-s2_noise.fits band=s2 \
    flux_factor=1e6


echo "-----------------------------"
echo "    Herschel SPIRE 500um     "
echo "-----------------------------"
../../bin/ifni-gennoise $CATALOG $NOISE_OPTIONS \
    out=$MAPDIR/herschel-s3 psf=$PSFDIR/herschel-s3.fits \
    aspix=7.2 rms=1.14e-3

../../bin/ifni-genmap   $CATALOG $MAP_OPTIONS \
    out=$MAPDIR/herschel-s3_sci.fits psf=$PSFDIR/herschel-s3.fits \
    noise_map=$MAPDIR/herschel-s3_noise.fits band=s3 \
    flux_factor=1e6
