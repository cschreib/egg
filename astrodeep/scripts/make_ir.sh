# Description
# -----------
#
# This script will create the Spitzer MIPS and Herschel images
# corresponding to a given mock catalog, using the same observing
# conditions as in the GOODS-South field (same image depths).
#


VERSION=1.0
CATBASE=egg-ad-0p25deg

MAPDIR=../maps/v$VERSION
CATDIR=../catalogs/v$VERSION
CATALOG=$CATDIR/$CATBASE.fits

NOISE_OPTIONS="verbose seed=42"
MAP_OPTIONS="verbose"

# Abort immediately on error
set -e

echo "-----------------------------"
echo "      Spitzer IRS 16um       "
echo "-----------------------------"
egg-gennoise $CATALOG $NOISE_OPTIONS \
    out=$MAPDIR/$CATBASE-spitzer-irs16 psf=spitzer-irs16.fits \
    aspix=0.9 rms=1.6e-2

egg-genmap   $CATALOG $MAP_OPTIONS \
    out=$MAPDIR/$CATBASE-spitzer-irs16-sci.fits psf=spitzer-irs16.fits \
    noise_map=$MAPDIR/$CATBASE-spitzer-irs16-noise.fits band=spitzer-irs16 \
    flux_factor=21.8835


echo "-----------------------------"
echo "      Spitzer MIPS 24um      "
echo "-----------------------------"
egg-gennoise $CATALOG $NOISE_OPTIONS \
    out=$MAPDIR/$CATBASE-spitzer-mips24 psf=spitzer-mips24.fits \
    aspix=1.2 rms=3.4e-2

egg-genmap   $CATALOG $MAP_OPTIONS \
    out=$MAPDIR/$CATBASE-spitzer-mips24-sci.fits psf=spitzer-mips24.fits \
    noise_map=$MAPDIR/$CATBASE-spitzer-mips24-noise.fits band=spitzer-mips24 \
    flux_factor=7.78838


echo "-----------------------------"
echo "     Herschel PACS 70um      "
echo "-----------------------------"
egg-gennoise $CATALOG $NOISE_OPTIONS \
    out=$MAPDIR/$CATBASE-herschel-pacs70 psf=herschel-pacs70.fits \
    aspix=1.2 rms=1.95e-5

egg-genmap   $CATALOG $MAP_OPTIONS \
    out=$MAPDIR/$CATBASE-herschel-pacs70-sci.fits psf=herschel-pacs70.fits \
    noise_map=$MAPDIR/$CATBASE-herschel-pacs70-noise.fits band=herschel-pacs70 \
    flux_factor=1.452e6


echo "-----------------------------"
echo "     Herschel PACS 100um     "
echo "-----------------------------"
egg-gennoise $CATALOG $NOISE_OPTIONS \
    out=$MAPDIR/$CATBASE-herschel-pacs100 psf=herschel-pacs100.fits \
    aspix=1.2 rms=8.23e-6

egg-genmap   $CATALOG $MAP_OPTIONS \
    out=$MAPDIR/$CATBASE-herschel-pacs100-sci.fits psf=herschel-pacs100.fits \
    noise_map=$MAPDIR/$CATBASE-herschel-pacs100-noise.fits band=herschel-pacs100 \
    flux_factor=1.672e6


echo "-----------------------------"
echo "     Herschel PACS 160um     "
echo "-----------------------------"
egg-gennoise $CATALOG $NOISE_OPTIONS \
    out=$MAPDIR/$CATBASE-herschel-pacs160 psf=herschel-pacs160.fits \
    aspix=2.4 rms=1.68e-5

egg-genmap   $CATALOG $MAP_OPTIONS \
    out=$MAPDIR/$CATBASE-herschel-pacs160-sci.fits psf=herschel-pacs160.fits \
    noise_map=$MAPDIR/$CATBASE-herschel-pacs160-noise.fits band=herschel-pacs160 \
    flux_factor=1.613e6


echo "-----------------------------"
echo "    Herschel SPIRE 250um     "
echo "-----------------------------"
egg-gennoise $CATALOG $NOISE_OPTIONS \
    out=$MAPDIR/$CATBASE-herschel-spire250 psf=herschel-spire250.fits \
    aspix=3.6 rms=1.66e-3

egg-genmap   $CATALOG $MAP_OPTIONS \
    out=$MAPDIR/$CATBASE-herschel-spire250-sci.fits psf=herschel-spire250.fits \
    noise_map=$MAPDIR/$CATBASE-herschel-spire250-noise.fits band=herschel-spire250 \
    flux_factor=1e6


echo "-----------------------------"
echo "    Herschel SPIRE 350um     "
echo "-----------------------------"
egg-gennoise $CATALOG $NOISE_OPTIONS \
    out=$MAPDIR/$CATBASE-herschel-spire350 psf=herschel-spire350.fits \
    aspix=4.8 rms=1.66e-3

egg-genmap   $CATALOG $MAP_OPTIONS \
    out=$MAPDIR/$CATBASE-herschel-spire350-sci.fits psf=herschel-spire350.fits \
    noise_map=$MAPDIR/$CATBASE-herschel-spire350-noise.fits band=herschel-spire350 \
    flux_factor=1e6


echo "-----------------------------"
echo "    Herschel SPIRE 500um     "
echo "-----------------------------"
egg-gennoise $CATALOG $NOISE_OPTIONS \
    out=$MAPDIR/$CATBASE-herschel-spire500 psf=herschel-spire500.fits \
    aspix=7.2 rms=1.88e-3

egg-genmap   $CATALOG $MAP_OPTIONS \
    out=$MAPDIR/$CATBASE-herschel-spire500-sci.fits psf=herschel-spire500.fits \
    noise_map=$MAPDIR/$CATBASE-herschel-spire500-noise.fits band=herschel-spire500 \
    flux_factor=1e6
