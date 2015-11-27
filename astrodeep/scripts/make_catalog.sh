# Description
# -----------
#
# This script will create a mock catalog similar to the one released
# by the Astrodeep collaboration.
#

VERSION=1.0
CATBASE=egg-ad-0p25deg

CAT_OPTIONS="verbose seed=42 maglim=28 selection_band=hst-f160w area=0.25"

egg-gencat $CAT_OPTIONS out=../catalogs/v$VERSION/$CATBASE.fits

