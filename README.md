# Ifni
Ifni is a tool to generate fake galaxy catalogs with realistic positions, morphologies and fluxes from the far-ultraviolet to the far-infrared. The generated catalogs are stored in binary FITS tables (column oriented). Another program, ifni-2skymaker, is used to convert the generated catalog into ASCII tables suitable for ingestion by [SkyMaker] to generate realistic high resolution images.

This set of tools can be used to test source extraction codes, or to evaluate the reliability of any map-based science (stacking, dropout identification, ...).

# Installation
You must have the [phy++] library installed on your machine to compile Ifni.
[CMake] is used as a build system to handle dependency checks and compilation in a cross-platform way.

One you have installed all the dependencies, create yourself a directory called 'build' within this directory. Navigate to the 'build' directory with your terminal, and call:

    cmake ../
    make install

The binaries will be generated in the 'bin' folder. More detailed information can be found in the script 'install.sh' or in the documentation ('doc/' folder).

# Usage
See the output of 'ifni help' after installing the program.
For further and more detailed help, see the documentation in the 'doc/' folder.

[SkyMaker]: http://www.astromatic.net/software/skymaker
[phy++]: http://cschreib.github.io/phypp/
[CMake]: http://www.cmake.org/
