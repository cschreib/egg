#include <phypp.hpp>

void print_help();

// SkyMaker pixel type (see PIXTYPE in skymaker/src/fits/fitscat.h)
using SkyPixType = float;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_help();
        return 1;
    }

    // Read command line arguments
    std::string band;
    std::string out_file;
    double maglim = dnan;
    double aspix = 0.06;
    double inset = 5.0;
    bool npix = false;
    bool verbose = false;
    float size_cap = fnan;
    std::string save_pixpos;

    read_args(argc-1, argv+1, arg_list(
        band, name(out_file, "out"), maglim, aspix, npix, size_cap, save_pixpos, verbose
    ));

    // Read input catalog
    if (verbose) {
        note("reading catalog...");
    }

    struct {
        vec1u id;
        vec1d ra, dec;

        vec1f disk_angle,  disk_radius,  disk_ratio;
        vec1f bulge_angle, bulge_radius, bulge_ratio;

        vec2f flux_disk, flux_bulge;
        vec1s bands;
    } cat;

    fits::read_table(argv[1], ftable(
        cat.id, cat.ra, cat.dec,
        cat.disk_angle, cat.disk_radius, cat.disk_ratio,
        cat.bulge_angle, cat.bulge_radius, cat.bulge_ratio,
        cat.flux_disk, cat.flux_bulge, cat.bands
    ));

    vec1u idb = where(cat.bands == band);
    if (idb.empty()) {
        error("no band named '", band, "' in this catalog");
        return 1;
    } else if (idb.size() > 1) {
        error("multiple bands named '", band, "' in this catalog");
        return 1;
    }

    if (out_file.empty()) {
        out_file = file::remove_extension(argv[1])+"-"+band+".cat";
    }

    file::mkdir(file::get_directory(out_file));

    uint_t b = idb[0];

    // Build the quantities expected by SkyMaker
    if (verbose) {
        note("convert parameters for ingestion by SkyMaker...");
    }

    vec1f mag = uJy2mag(cat.flux_disk(_,b)+cat.flux_bulge(_,b));
    vec1f bt = cat.flux_bulge(_,b)/(cat.flux_disk(_,b)+cat.flux_bulge(_,b));

    // Convert effective radius to scale length, which is what SkyMaker expects
    // (NB: only for the disk, the bulge size has to be the effective radius)
    cat.disk_radius /= 1.678;

    // Try first to make a single, big image
    if (verbose) {
        note("define image size and pixel coordinates...");
    }

    // Roughly estimate galaxy size to check for border effects [pixels]
    vec1d gsize = 10*max(cat.disk_radius*1.678, cat.bulge_radius)/aspix;

    // Try a first naive projection where the first pixel is defined by the
    // most extreme galaxy in the corner of the field
    fits::make_wcs_header_params wcs_params;
    wcs_params.pixel_scale = aspix;
    wcs_params.pixel_ref_x = 1.0; wcs_params.pixel_ref_y = 1.0;
    wcs_params.sky_ref_ra = max(cat.ra); wcs_params.sky_ref_dec = min(cat.dec);

    fits::header hdr;
    if (!fits::make_wcs_header(wcs_params, hdr)) {
        error("could not make WCS header");
        return 1;
    }

    vec1d x, y;
    fits::wcs wcs_big(hdr);
    fits::ad2xy(wcs_big, cat.ra, cat.dec, x, y);

    // Now, because of the spherical projection, some galaxies can have negative pixel
    // coordinates even though we picked the most extreme one in RA and Dec.
    // Figure out how negative this can be and adjust the reference pixel accordingly.
    // Also add the inset that was requested by the user, and take into account the
    // extent of each galaxy
    double dx = inset/aspix - (min(x - gsize) - 1);
    double dy = inset/aspix - (min(y - gsize) - 1);
    wcs_params.pixel_ref_x += dx;
    wcs_params.pixel_ref_y += dy;

    if (!fits::make_wcs_header(wcs_params, hdr)) {
        error("could not make WCS header");
        return 1;
    }

    wcs_big = fits::wcs(hdr);
    fits::ad2xy(wcs_big, cat.ra, cat.dec, x, y);

    // Compute expected full image size in GB
    uint_t nx = ceil(max(x + gsize) + inset/aspix);
    uint_t ny = ceil(max(y + gsize) + inset/aspix);

    // Select only the galaxies with valid of bright enough fluxes
    vec1u ids;
    if (is_finite(maglim)) {
        ids = where(mag < maglim);
        if (ids.empty()) {
            error("no source brighter than ", maglim);
            note("minimum magnitude is ", min(mag));
            return 1;
        }
    } else {
        ids = uindgen(x.size());
        vec1u idbad = where(!is_finite(mag));
        mag[idbad] = 99;
    }

    double img_size = double(nx)*ny*sizeof(SkyPixType)/pow(1024.0, 3);

    // Inline function to write a catalog to disk
    auto write_catalog = [&](std::string ofile, const vec1u& oids){
        // Write catalog
        file::write_table_hdr(ofile, 16,
            {"type", "x", "y", "mag", "bt",
            "bulge_radius", "bulge_ratio", "bulge_angle",
            "disk_radius", "disk_ratio", "disk_angle"},
            replicate(200u, oids.size()), x[oids], y[oids],
            mag[oids], bt[oids],
            cat.bulge_radius[oids], cat.bulge_ratio[oids], cat.bulge_angle[oids],
            cat.disk_radius[oids], cat.disk_ratio[oids], cat.disk_angle[oids]
        );

        // Write FITS header
        {
            std::ofstream ohdr(file::remove_extension(out_file)+"-hdr.txt");
            vec1s shdr = cut(hdr, 80);
            for (auto& s : shdr) {
                ohdr << s << "\n";
            }
        }
    };

    vec1s tile(x.size());

    if (is_finite(size_cap) && img_size > size_cap) {
        // The whole image would be too big, we need to split it into smaller sections
        // of equal size.
        if (verbose) {
            note("full image would be ", nx, " x ", ny, " (", img_size, " GB)");
        }

        // Estimate the number of such sections
        const double fudge = 0.95; // make sure we are below the maximum size
        double nsec = 1.0/(fudge*size_cap/img_size);

        // Estimate the number of X and Y sections
        double ar = nx/double(ny);
        double tnsx = sqrt(ar*nsec);
        double tnsy = nsec/tnsx;

        uint_t nsx, nsy;
        if (fabs(round(tnsx) - tnsx) < fabs(round(tnsy) - tnsy)) {
            nsx = round(tnsx);
            nsy = ceil(nsec/nsx);
        } else {
            nsy = round(tnsy);
            nsx = ceil(nsec/nsy);
        }

        nx = ceil(nx/double(nsx));
        ny = ceil(ny/double(nsy));

        if (npix || verbose) {
            double sec_size = double(nx)*ny*sizeof(SkyPixType)/pow(1024.0, 3);
            note("will require ", nsx, "x", nsy, " sections (", sec_size, " GB each)");
            note("section dimensions: ", nx, ",", ny);
        }

        if (verbose) {
            note("write catalogs...");
        }

        std::string out_base = file::remove_extension(out_file);

        tile.resize(x.size());

        // Keep "old" coordinates from the whole image
        vec1d ox = x;
        vec1d oy = y;
        // And these will be the new coordinates of each galaxy
        vec1d tx = x;
        vec1d ty = y;

        auto pg = progress_start(nsx*nsy);
        for (uint_t iy : range(nsy))
        for (uint_t ix : range(nsx)) {
            // Pick galaxies that fall in this section
            // NB: pixel coordinates in FITS standard are in [1,N] and '1' is the center
            // of the first pixel (which therefore covers 0.5 to 1.5)
            vec1u idi = ids[where(ox[ids]+gsize[ids] >= nx*ix+0.5 && ox[ids]-gsize[ids] <= nx*(ix+1)+0.5 &&
                                  oy[ids]+gsize[ids] >= ny*iy+0.5 && oy[ids]-gsize[ids] <= ny*(iy+1)+0.5)];

            // This is the "strict" list, i.e., only galaxies whose center lies in this tile
            vec1u idis = where(ox[idi] >= nx*ix+0.5 && ox[idi] < nx*(ix+1)+0.5 &&
                               oy[idi] >= ny*iy+0.5 && oy[idi] < ny*(iy+1)+0.5);

            // Convert pixel coordinates to local frame
            wcs_params.pixel_ref_x = 1-nx*double(ix);
            wcs_params.pixel_ref_y = 1-ny*double(iy);

            if (!fits::make_wcs_header(wcs_params, hdr)) {
                error("could not make WCS header");
                return 1;
            }

            vec1d ttx, tty;
            fits::ad2xy(fits::wcs(hdr), cat.ra[idi], cat.dec[idi], ttx, tty);
            x[idi]        = ttx;       y[idi]        = tty;
            tx[idi[idis]] = ttx[idis]; ty[idi[idis]] = tty[idis];

            // Write the catalog for this section
            std::string six = align_right(strn(ix+1), strn(nsx).size(), '0');
            std::string siy = align_right(strn(iy+1), strn(nsy).size(), '0');
            out_file = out_base+"-"+six+"-"+siy+".cat";
            write_catalog(out_file, idi);

            tile[idi[idis]] = file::remove_extension(file::get_basename(out_file));

            if (verbose) progress(pg);
        }

        std::swap(tx, x);
        std::swap(ty, y);
    } else {
        // We will make a single big image
        if (npix || verbose) {
            note("image dimensions: ", nx, ",", ny, " (", img_size, " GB)");
        }

        if (verbose) {
            note("write catalog...");
        }

        write_catalog(out_file, ids);

        tile = replicate(file::remove_extension(file::get_basename(out_file)), x.size());
    }

    if (!save_pixpos.empty()) {
        note("write pixel positions...");
        file::mkdir(file::get_directory(save_pixpos));
        fits::write_table(save_pixpos, ftable(cat.id, x, y, tile));
    }

    if (verbose) {
        note("done.");
    }

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

    print("ifni-2skymaker v1.0");
    print("usage: ifni-2skymaker cat.fits band=... [options]\n");

    print("List of options:");
    argdoc("band", "[string]", "name of the photometric band for which a SkyMaker "
        "input catalg will be created (required)");
    argdoc("out", "[string]", "name of the output catalog (default: [catalog]-[band].cat)");
    argdoc("maglim", "[float]", "galaxies fainter than this magnitude will not be "
        "present in the output catalog (default: none)");
    argdoc("aspix", "[double, arcsec/pixel]", "pixel size of the SkyMaker image "
        "(default: 0.06\")");
    argdoc("inset", "[double, arcsec]", "add an empty border around the input catalog "
        "to prevent galaxies from being truncated by the edge of the image "
        "(default: 5\")");
    argdoc("npix", "[flag]", "ask the program to estimate the optimal image dimensions "
        "that should be put into the SkyMaker configuration file (in pixels)");
    argdoc("size_cap", "[float, GB]", "split the survey in multiple tiles of equal area "
        "that will each be rendered into a separate image, so that the size of each "
        "such image is less than 'size_cap' (expressed in GigaBytes, default: none). NB: "
        "SkyMaker requires the whole image to fit inside your RAM memory, so this value "
        "should be no larger than the total amount of RAM memory available on your "
        "computer.");
    argdoc("save_pixpos", "[string]", "specify the name of the file into which the pixel "
        "coordinates and the name of the image of each galaxy will be written (default: "
        "not saved). This is most useful when 'size_cap' is used, as this catalog "
        "provides and easy way to figure out in which of the tile each galaxy is located.");
    argdoc("verbose", "[flag]", "display information about the progress of the program "
        "in the terminal");
    print("");
}
