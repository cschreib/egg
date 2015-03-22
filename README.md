# gencat
gencat is a tool to generate fake galaxy catalogs with realistic positions, morphologies and fluxes from the far-ultraviolet to the far-infrared. The generated catalogs are stored in binary FITS tables (column oriented). Another program, make_skymaker, is used to convert the generated catalog into ASCII tables suitable for ingestion by skymaker[1] to generate realistic images.

This set of tools can be used to test source extraction codes, or to evaluate the reliability of any map-based science (stacking, dropout identification, ...).

# Installation
You must have the phy++ library[2] installed on your machine to compile gencat.
CMake[3] is used as a build system to handle dependency checks and compilation in a cross-platform way.

# Usage
See the output of 'gencat help'.
For further and more detailed help, see the documentation in the 'doc/' folder.

[1]: http://www.astromatic.net/software/skymaker
[2]: http://cschreib.github.io/phypp/
[3]: http://www.cmake.org/
