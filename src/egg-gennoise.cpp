#include <phypp.hpp>

const std::string egg_share_dir = file::directorize(EGG_SHARE_DIR);

void print_help();

int phypp_main(int argc, char* argv[]) {
    // Random seed to use to generate the noise
    uint_t tseed = 42;
    // Name of the FITS table containing the mock catalog
    std::string cat_file;
    // Name of the output FITS file "basename" in which to store the map
    // If out="file", then the map will be saved in "file-noise.fits"
    std::string out;
    // Name of the FITS file of the PSF
    std::string psf_file;
    // Name of the file from which to copy the astrometry
    std::string astro;
    // If not the above, specify the pixel scale
    double aspix = dnan;
    // The amplitude of the noise
    double rms = 1.0;

    // Filter the final map with beam smoothing
    bool beam_smoothed = false;
    // Beam smoothing uses the PSF by default, but you can choose an arbitrary Gaussian
    double smooth_fwhm = dnan;
    // Flag all the border pixels that are not well covered by NaN
    bool clip_borders = false;
    // Also build a coverage map ("file-cov.fits")
    bool make_cov = false;
    // Also build an error map ("file-err.fits")
    bool make_err = false;
    // Save maps with double precision
    bool tdouble = false;
    // Display the list of default available PSFs
    bool list_psfs = false;
    // Print some text on the console.
    bool verbose = false;
    // Save the distance map (debug purposes, "file-dist.fits")
    bool save_dist = false;

    if (argc < 2) {
        print_help();
        return 0;
    }

    read_args(argc, argv, arg_list(
        name(tseed, "seed"), out, name(psf_file, "psf"), astro, aspix, rms, beam_smoothed,
        smooth_fwhm, list_psfs, verbose, make_err, make_cov, clip_borders, save_dist,
        name(cat_file, "cat"), name(tdouble, "double")
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

    if (cat_file.empty() || out.empty() || psf_file.empty() || (astro.empty() && !is_finite(aspix))) {
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

    std::string out_base = file::remove_extension(out);
    if (end_with(out_base, "-noise")) {
        out_base = erase_end(out_base, "-noise");
    }

    auto seed = make_seed(tseed);

    struct {
        vec1d ra, dec;
    } cat;

    if (end_with(cat_file, ".fits")) {
        fits::read_table(cat_file, ftable(cat.ra, cat.dec));
    } else {
        ascii::read_table(cat_file, ascii::find_skip(cat_file), cat.ra, cat.dec);
    }

    // Read the PSF
    if (verbose) note("reading and normalizing PSF...");

    vec2d psf;
    fits::read(psf_file, psf);

    // Center the PSF
    int_t hsize = std::max(psf.dims[0], psf.dims[1])/2;
    if (hsize > 40) hsize = 40;
    vec1i mid = mult_ids(psf, max_id(psf));
    psf = subregion(psf, {mid[0]-hsize, mid[1]-hsize, mid[0]+hsize, mid[1]+hsize});

    // Compute PSF FWHM
    double fwhm; {
        vec2d r2 = generate_img(psf.dims, [&](double ix, double iy) {
            return sqr(ix - hsize) + sqr(iy - hsize);
        });

        vec1u idok = where(psf > 0.2*max(psf) && r2 > 0);
        double res = mean(log(psf[idok]/max(psf))/r2[idok]);
        fwhm = 2*sqrt(-log(2)/res);

        if (verbose) note("FWHM=", fwhm, " pixels");
    }

    // Read astrometry
    if (verbose) note("configure astrometry...");

    vec2d noise;
    fits::header hdr; {
        if (astro.empty()) {
            // From parameters

            // Try a first naive projection where the first pixel is defined by the
            // most extreme galaxy in the corner of the field
            fits::make_wcs_header_params wcs_params;
            wcs_params.pixel_scale = aspix;
            wcs_params.pixel_ref_x = 1.0; wcs_params.pixel_ref_y = 1.0;
            wcs_params.sky_ref_ra = max(cat.ra); wcs_params.sky_ref_dec = min(cat.dec);

            if (!fits::make_wcs_header(wcs_params, hdr)) {
                error("could not make WCS header");
                return 1;
            }

            vec1d x, y;
            fits::ad2xy(fits::wcs(hdr), cat.ra, cat.dec, x, y);

            // Now, because of the spherical projection, some galaxies can have negative pixel
            // coordinates even though we picked the most extreme one in RA and Dec.
            // Figure out how negative this can be and adjust the reference pixel accordingly.
            // Also add an inset.
            double dx = 5*fwhm - (min(x) - 1);
            double dy = 5*fwhm - (min(y) - 1);
            wcs_params.pixel_ref_x += dx;
            wcs_params.pixel_ref_y += dy;

            if (!fits::make_wcs_header(wcs_params, hdr)) {
                error("could not make WCS header");
                return 1;
            }

            fits::ad2xy(fits::wcs(hdr), cat.ra, cat.dec, x, y);

            uint_t npx = ceil(max(x) + 5*fwhm);
            uint_t npy = ceil(max(y) + 5*fwhm);

            if (verbose) note("map will be ", npy, " x ", npx);

            noise.resize(npy, npx);
        } else {
            // From an existing FITS file
            fits::read(astro, noise, hdr);

            if (verbose) note("map will be ", noise.dims[0], " x ", noise.dims[1]);
        }
    }

    noise[_] = fnan;
    vec1u idin;

    if (clip_borders) {
        if (verbose) note("clipping field borders...");
        fits::wcs wcs(hdr);

        // Build convex hull around catalog sources
        vec1d sx, sy;
        fits::ad2xy(wcs, cat.ra, cat.dec, sx, sy);
        auto hull = build_convex_hull(sx, sy);

        // Compute distance of each pixel to the hull
        vec2d x = replicate(uindgen(noise.dims[1]), noise.dims[0]);
        vec2d y = transpose(replicate(uindgen(noise.dims[0]), noise.dims[1]));
        vec2d dist = convex_hull_distance(x, y, hull);

        if (save_dist) {
            fits::write(out_base+"-dist.fits", dist, hdr);
        }

        // Only use pixels that are more than 3 FWHM away from the border of the catalog
        idin = where(dist > 1.5*fwhm);
    } else {
        idin = uindgen(noise.size());
    }

    if (beam_smoothed) {
        if (verbose) note("estimate unfiltered noise level...");

        // Adjust the provided RMS to apply to the beam-smoothed map
        // Note: the map is not beam-smoothed yet, this will be done in 'egg-genmap'

        // Create the smoothing kernel
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

        // Create a random patch of noise and convolve it to see how it is affected
        vec2d rnd = randomn(seed, psf.dims[0]*4, psf.dims[1]*4);
        double orig_rms = stddev(rnd); // should be very close to 1
        rnd = convolve2d(rnd, kernel);
        double new_rms = stddev(rnd);

        // Correct the requested RMS so that it will apply to the final, beam-smoothed
        // map, but not the raw map with independent pixels
        rms *= orig_rms/new_rms;
    }

    if (verbose) note("generate noise...");
    noise[idin] = rms*randomn(seed, idin.size());

    // Save noise map
    if (verbose) note("write map to disk...");
    if (tdouble) {
        fits::write(out, noise, hdr);
    } else {
        fits::write(out, vec2f{noise}, hdr);
    }


    if (make_err) {
        // Build and save error map
        noise[idin] = rms;
        if (tdouble) {
            fits::write(out_base+"-err.fits", noise, hdr);
        } else {
            fits::write(out_base+"-err.fits", vec2f{noise}, hdr);
        }
    }

    if (make_cov) {
        // Build and save coverage map
        noise[idin] = 1.0;
        if (tdouble) {
            fits::write(out_base+"-cov.fits", noise, hdr);
        } else {
            fits::write(out_base+"-cov.fits", vec2f{noise}, hdr);
        }
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

    print("egg-gennoise v1.0rc1");
    print("usage: egg-gennoise cat=... out=... psf=... astro=... aspix=... rms=... [options]\n");

    print("List of mandatory parameters (no default value):");
    argdoc("cat", "[string]", "name of the file containing the mock catalog");
    argdoc("out", "[string]", "name of the file into which the noise map will be saved.");
    argdoc("psf", "[string]", "name of the file containing the PSF. Note that, "
        "if the provided file does not exists, the program will also search in the list "
        "of the default PSFs provided with EGG ("+egg_share_dir+
        "/psfs) for a file with the same name, and use it if it exists.");
    argdoc("astro", "[string]", "name of an existing FITS file from which the astrometry "
        "should be copied. The alternative is to specify the pixel scale (see below).");
    argdoc("aspix", "[double, arcsec/pixel]", "if 'astro' is not provided, one can also "
        "define the astrometry solution simply by giving the pixel scale. In this case "
        "the WCS data will be set so that the resulting image encompasses the whole "
        "input catalog.");
    argdoc("rms", "[double, map unit]", "specify the target pixel RMS of the image "
        "(pure Gaussian noise), as measured by the standard deviation of pixel values "
        "inside an aperture. This value will apply to the *final* map, i.e., after all "
        "post-processing (including in particular beam smoothing).");
    print("");

    print("List of options:");
    argdoc("beam_smoothed", "[flag]", "after the sources are placed on the map, convolve "
        "the map with the PSF (or an aribtrary Gaussian, see 'smooth_fwhm'). This is the "
        "standard procedure in the reduction of submillimeter images (default: false).");
    argdoc("smooth_fwhm", "[double]", "if 'beam_smoothed' is set, the default is to use "
        "the PSF itself to perform the convolution. This option allows choosing any other "
        "aribtrary Gaussian instead by specifying the FWHM manually (default: use PSF)");
    argdoc("clip_borders", "[flag]", "if set, pixels that are too close to the catalog "
        "border are flagged with NaN, since they lack contribution from nearby sources. "
        "This option ensures that the map statistics is uniform, and that no edge effects "
        "are present, but some galaxies will be clipped out.");
    argdoc("make_cov", "[flag]", "if set, save a coverage map (\"[out]-cov.fits\").");
    argdoc("make_err", "[flag]", "if set, save an RMS map (\"[out]-err.fits\").");
    argdoc("double", "[flag]", "write maps on the disk with double precision (default: "
        "use single precision)");
    argdoc("verbose", "[flag]", "print additional information in the standard output while "
        "the program is running (default: false)");
    argdoc("list_psfs", "[flag]", "display a list of all the default PSFs provided with EGG "
        "and exit");
    print("");
}
