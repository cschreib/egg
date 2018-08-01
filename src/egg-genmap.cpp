#include <phypp.hpp>

const std::string egg_share_dir = file::directorize(EGG_SHARE_DIR);

void print_help();

int phypp_main(int argc, char* argv[]) {
    // Name of the FITS table containing the mock catalog
    std::string cat_file;
    // Name of the output FITS file in which to store the map
    std::string out;
    // Name of the FITS file of the PSF
    std::string psf_file;
    // Name of the pre-computed noise map
    std::string noise_map;
    // Name of the band from which the fluxes come from
    std::string band;

    // Conversion from uJy*PSF to map units
    double flux_factor = 1.0;
    // Filter the final map with beam smoothing
    bool beam_smoothed = false;
    // Beam smoothing uses the PSF by default, but you can choose an arbitrary Gaussian
    double smooth_fwhm = dnan;
    // Put the mean of the map to zero in post-processing
    bool zero_mean = false;
    // Do not put the sources on the map (only use the noise)
    bool no_source = false;
    // Do not apply sub-pixel corrections to the PSF of each source
    bool no_subpixel = false;
    // Extend the PSF to the full cutout
    bool extpsf = false;
    // Output the map in single precision instead of double precision
    bool cdouble = false;
    // Display the list of default available PSFs
    bool list_psfs = false;
    // Print some text on the console.
    bool verbose = false;

    if (argc < 2) {
        print_help();
        return 0;
    }

    read_args(argc, argv, arg_list(
        out, name(psf_file, "psf"), noise_map, band, beam_smoothed, smooth_fwhm,
        verbose, no_source, no_subpixel, flux_factor, zero_mean, list_psfs,
        name(extpsf, "extend_psf"), name(cdouble, "double"), name(cat_file, "cat")
    ));

    auto display_psf_list = [&]() {
        note("available default PSFs:");
        vec1s tmps = file::list_files(egg_share_dir+"psfs/*.fits");
        inplace_sort(tmps);
        for (auto t : tmps) {
            print(" - ", t);
        }
    };

    if (list_psfs) {
        display_psf_list();
        return 0;
    }

    if (cat_file.empty() || out.empty() || psf_file.empty() || noise_map.empty()) {
        print_help();
        return 0;
    }

    if (!file::exists(psf_file)) {
        psf_file = egg_share_dir+"psfs/"+psf_file;
    }
    if (!file::exists(psf_file)) {
        error("the PSF file '"+psf_file+"' could not be found");
        display_psf_list();
        return 1;
    }

    file::mkdir(file::get_directory(out));

    // Read the input catalog
    struct {
        vec1d ra, dec;

        vec2f flux;
        vec1s bands;
    } cat;

    uint_t idb;
    if (ends_with(cat_file, ".fits")) {
        fits::input_table table(cat_file);
        table.read_columns(ftable(cat.ra, cat.dec));
        table.read_column(fits::dim_promote, "flux", cat.flux);

        if (cat.flux.dims[1] > 1) {
            table.read_column("bands", cat.bands);
        } else {
            table.read_column(fits::missing, "bands", cat.bands);
        }

        table.close();

        // Find the band we are looking for
        if (cat.bands.empty() || (cat.bands.size() == 1 && band.empty())) {
            idb = 0;
        } else {
            vec1u ids = where(cat.bands == band);
            if (ids.empty()) {
                warning("no band named '", band, "' in this catalog");
                return 1;
            } else if (ids.size() > 1) {
                warning("multiple bands matching '", band, "'");
                return 1;
            }

            idb = ids[0];
        }
    } else {
        ascii::read_table(cat_file, cat.ra, cat.dec, ascii::columns(1,cat.flux));

        idb = 0;
    }

    // Read the PSF
    vec2d psf;
    fits::read(psf_file, psf);

    // Center the PSF
    int_t hsize = std::max(psf.dims[0], psf.dims[1])/2;
    if (hsize > 40) hsize = 40;
    vec1i mid = mult_ids(psf, max_id(psf));
    psf = subregion(psf, {mid[0]-hsize, mid[1]-hsize, mid[0]+hsize, mid[1]+hsize});

    // Compute FWHM
    double fwhm;
    vec2d r2 = generate_img(psf.dims, [&](double ix, double iy) {
        return sqr(ix - hsize) + sqr(iy - hsize);
    });

    {
        vec1u idok = where(psf > 0.2*max(psf) && r2 > 0);
        double res = mean(log(psf[idok]/max(psf))/r2[idok]);
        fwhm = 2*sqrt(-log(2)/res);

        if (verbose) note("FWHM=", fwhm, " pixels");
    }

    // If asked, fit the central part of the provided PSF with a 2D gaussian and use
    // it to fill in the pixels where the input PSF is zero.
    if (extpsf) {
        vec1u idz = where(psf == 0.0);
        psf[idz] = max(psf)*exp(-log(2)*r2[idz]/sqr(fwhm/2.0));
    }

    // Read noise map
    vec2d img;
    fits::header hdr;
    fits::read(noise_map, img, hdr);

    // Read astrometry
    astro::wcs w(hdr);
    vec1d tx, ty;
    astro::ad2xy(w, cat.ra, cat.dec, tx, ty);
    tx -= 1; ty -= 1;
    vec1i x = round(tx), y = round(ty);

    // Select sources to put on the map
    vec1b sel = x-hsize < int_t(img.dims[1]) && x+hsize >= 0 &&
                y-hsize < int_t(img.dims[0]) && y+hsize >= 0 &&
                is_finite(cat.flux(_,idb));

    vec1u gid = where(sel);

    if (verbose) note(gid.size(), " sources fall on this map");

    // Put the sources on the map
    vec2d tpsf = psf;
    if (!no_source) {
        if (verbose) note("placing sources...");

        auto pg = progress_start(gid.size());
        for (uint_t i : gid) {
            if (!no_subpixel) {
                tpsf = translate(psf, ty[i] - y[i], tx[i] - x[i]);
            }

            double flx = cat.flux(i,idb);

            vec1u idm, idp;
            subregion(img, {y[i]-hsize, x[i]-hsize, y[i]+hsize, x[i]+hsize}, idm, idp);
            img.safe[idm] += (flx/flux_factor)*tpsf.safe[idp];

            if (verbose) progress(pg, 127);
        }
    }

    if (beam_smoothed) {
        if (verbose) note("beam smoothing...");

        // Create smoothing kernel
        vec2d kernel;
        if (is_finite(smooth_fwhm)) {
            kernel = gaussian_profile(psf.dims, smooth_fwhm/2.355);
        } else {
            kernel = psf;
        }

        // Remove zero elements
        vec1u idz = where(kernel < 1e-5*max(abs(kernel)));
        kernel[idz] = 0.0;

        // Normalize to unity to preserve flux
        vec2d cpsf = convolve2d(psf, kernel);
        kernel /= cpsf(hsize,hsize);

        // Temporarilly remove NaN pixels from the map
        vec2b bad = !is_finite(img);
        img[where(bad)] = 0.0;

        // Convolve
        img = convolve2d(img, kernel);

        // Re-introduce NaN pixels, inflating them to reach one FWHM further
        // to avoid having pixels in the map which are were convolved with only
        // half their surroundings.
        bad = mask_inflate(bad, ceil(fwhm));
        img[where(bad)] = dnan;
    }

    if (zero_mean) {
        vec1u idg = where(is_finite(img));
        img -= mean(img[idg]);
    }

    // Save the simulated map
    if (verbose) note("saving map to disk...");

    if (cdouble) {
        fits::write(out, img, hdr);
    } else {
        vec2f fimg = img;
        fits::write(out, fimg, hdr);
    }

    if (verbose) note("done.");

    return 0;
}

void print_help() {
    using namespace format;

    auto argdoc = [](const std::string& name, const std::string& type,
        const std::string& desc) {

        std::string header = " - "+name+" "+type;
        print(header);

        std::string indent = "    ";
        vec1s w = wrap(indent+desc, 80, indent);

        for (auto& s : w) {
            print(s);
        }
    };

    print("egg-genmap v1.0rc1");
    print("usage: egg-genmap cat=... out=... psf=... noise_map=... band=... [options]\n");

    print("List of mandatory parameters (no default value):");
    argdoc("cat", "[string]", "name of the file containing the mock catalog");
    argdoc("out", "[string]", "name of the file into which the map will be saved");
    argdoc("psf", "[string]", "name of the file containing the PSF. Note: the performance "
        "of this program depends greatly on the dimensions of this PSF image. Avoid using "
        "large PSF images unless you have very bright sources in the catalog. Note that, "
        "if the provided file does not exists, the program will also search in the list "
        "of the default PSFs provided with EGG ("+egg_share_dir+
        "/psfs) for a file with the same name, and use it if it exists.");
    argdoc("noise_map", "[string]", "name of the file containing the noise map (e.g, "
        "created by egg-gennoise)");
    argdoc("band", "[string]", "name of the photometric band in the input catalog from "
        "which the fluxes should be read");
    print("");

    print("List of options:");
    argdoc("flux_factor", "[double]", "conversion factor from uJy to map units. Note that "
        "this factor should depend on the PSF. If the PSF is normalized to unit integral, "
        "then the conversion factor is from uJy/pixel to map units. If the PSF is "
        "normalized to unit peak flux, then the conversion factor is from uJy/beam to map "
        "units (default: none)");
    argdoc("beam_smoothed", "[flag]", "after the sources are placed on the map, convolve "
        "the map with the PSF (or an aribtrary Gaussian, see 'smooth_fwhm'). This is the "
        "standard procedure in the reduction of submillimeter images (default: false).");
    argdoc("smooth_fwhm", "[double]", "if 'beam_smoothed' is set, the default is to use "
        "the PSF itself to perform the convolution. This option allows choosing any other "
        "aribtrary Gaussian instead by specifying the FWHM manually (default: use PSF)");
    argdoc("extend_psf", "[flag]", "in case the PSF is truncated beyond a certain radius, "
        "enabling this option will extrapolate the available PSF data to fill the missing "
        "pixels (default: false). A Gaussian profile is fitted to the core, and used for "
        "the extrapolation. Note however that the PSF will never be larger than the "
        "dimensions of the provided PSF image.");
    argdoc("no_subpixel", "[flag]", "do not perform sub-pixel interpolations when placing "
        "galaxies on the map (default: false). With this option activated, the map "
        "making will be quite faster, at the expense of correctness.");
    argdoc("no_source", "[flag]", "do not put sources on the map, just keep the noise and "
        "apply the pixel conversion and beam smoothing (default: false)");
    argdoc("verbose", "[flag]", "print additional information in the standard output while "
        "the program is running (default: false)");
    argdoc("list_psfs", "[flag]", "display a list of all the default PSFs provided with EGG "
        "and exit");
    print("");
}
