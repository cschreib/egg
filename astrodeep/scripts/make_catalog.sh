# Description
# -----------
#
# This script will create a mock catalog similar to the one released
# by the Astrodeep collaboration.
#

VERSION=1.0
CATBASE=ifni-ad-1deg

CAT_OPTIONS="verbose seed=42"

../../bin/ifni-gencat $CAT_OPTIONS out=../catalogs/v$VERSION$/$CATBASE.fits \
     maglim=29 area=1
