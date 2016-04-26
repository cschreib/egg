#include <phypp.hpp>

const std::string egg_share_dir = file::directorize(EGG_SHARE_DIR);

void print_help();

// SkyMaker pixel type (see PIXTYPE in skymaker/src/fits/fitscat.h)
using SkyPixType = float;

int phypp_main(int argc, char* argv[]) {
    if (argc < 2) {
        print_help();
        return 1;
    }

    // Read command line arguments
    std::string band;
    std::string cat_file;
    std::string out_file;
    std::string img_dir;
    std::string template_file;
    double maglim = dnan;
    double inset = 0.0;
    bool strict_clip = false;
    bool verbose = false;
    bool list_templates = false;
    float size_cap = fnan;
    std::string save_pixpos;

    read_args(argc, argv, arg_list(
        band, name(cat_file, "cat"), name(out_file, "out"), img_dir, maglim, size_cap,
        save_pixpos, verbose, strict_clip, name(template_file, "template"), list_templates
    ));

    if (!img_dir.empty()) {
        img_dir = file::directorize(img_dir);
    }

    auto print_template_list = [&]() {
        note("available default templates:");
        vec1s tmps = file::list_files(egg_share_dir+"skymaker-templates/*.conf");
        inplace_sort(tmps);
        for (auto t : tmps) {
            print(" - ", t);
        }
    };

    if (list_templates) {
        print_template_list();
        return 0;
    }

    // Check for missing mandatory parameters
    bool bad = false;
    if (cat_file.empty()) {
        error("missing input catalog file (cat=...)");
        bad = true;
    }
    if (band.empty()) {
        error("missing band ID (band=...)");
        bad = true;
    }
    if (template_file.empty()) {
        error("missing SkyMaker template configuration file (template=...)");
        bad = true;
    } else {
        if (!file::exists(template_file)) {
            template_file = egg_share_dir+"skymaker-templates/"+template_file;
        }
        if (!file::exists(template_file)) {
            error("the template file '"+template_file+"' could not be found");
            print_template_list();
            bad = true;
        }
    }

    if (bad) return 1;

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
        vec1f lambda;
    } cat;

    fits::read_table(cat_file, ftable(
        cat.id, cat.ra, cat.dec,
        cat.disk_angle, cat.disk_radius, cat.disk_ratio,
        cat.bulge_angle, cat.bulge_radius, cat.bulge_ratio,
        cat.flux_disk, cat.flux_bulge, cat.bands, cat.lambda
    ));

    if (verbose) {
        note("found ", cat.id.size(), " galaxies");
    }

    vec1u idb = where(cat.bands == band);
    if (idb.empty()) {
        error("no band named '", band, "' in this catalog");
        return 1;
    } else if (idb.size() > 1) {
        error("multiple bands named '", band, "' in this catalog");
        return 1;
    }

    uint_t b = idb[0];

    // Read template SkyMaker configuration file
    double aspix = fnan;
    vec1s sky_param, sky_value; {
        std::ifstream file(template_file);

        uint_t l = 0;
        std::string line;
        while (std::getline(file, line)) {
            ++l;

            line = trim(line);
            if (line.empty() || line[0] == '#') continue;

            auto pos = line.find_first_of(" \t");
            if (pos == line.npos) {
                note("reading ", template_file, ":", l);
                error("ill formed line, expected \"PARAMETER   VALUE\"");
                return 1;
            }

            std::string param = line.substr(0, pos);
            std::string value = trim(line.substr(pos));

            if (param == "LISTCOORD_TYPE" || param == "IMAGE_NAME" || param == "IMAGE_SIZE" ||
                param == "IMAGE_HEADER") {
                // We will set these ourselves
                continue;
            }

            sky_param.push_back(param);
            sky_value.push_back(value);

            if (param == "PIXEL_SIZE" && !from_string(value, aspix)) {
                note("reading ", template_file, ":", l);
                error("could not read pixel size value '", value, "'");
                return 1;
            }
        }
    }

    if (!is_finite(aspix)) {
        error("missing 'PIXEL_SIZE' parameter from SkyMaker configuration file");
        return 1;
    }

    // Make sure that the PSF file is given with an absolute path
    if (count(sky_param == "PSF_NAME") != 0) {
        uint_t i = where_first(sky_param == "PSF_NAME");
        if (!file::is_absolute_path(sky_value[i])) {
            sky_value[i] = file::get_directory(template_file)+sky_value[i];
        }
    }

    // If missing, add some parameters we can provide
    if (count(sky_param == "LISTCOORD_TYPE") == 0) {
        sky_param.push_back("LISTCOORD_TYPE");
        sky_value.push_back("PIXEL");
    }

    if (count(sky_param == "WAVELENGTH") == 0) {
        sky_param.push_back("WAVELENGTH");
        sky_value.push_back(strn(cat.lambda[b]));

    }

    if (count(sky_param == "IMAGE_TYPE") == 0) {
        sky_param.push_back("IMAGE_TYPE");
        sky_value.push_back("SKY");
    }

    if (out_file.empty()) {
        out_file = file::remove_extension(cat_file)+"-"+band+".cat";
    }

    // Now start the real job
    file::mkdir(file::get_directory(out_file));

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
    vec1d gsize;
    if (strict_clip) {
        gsize = replicate(0.0, cat.ra.size());
    } else {
        gsize = 10*max(cat.disk_radius*1.678, cat.bulge_radius)/aspix;
    }

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

        if (verbose) {
            note("selecting ", ids.size(), "/", cat.id.size(), " galaxies");
        }
    } else {
        ids = uindgen(x.size());
        vec1u idbad = where(!is_finite(mag));
        mag[idbad] = 99;
    }

    double img_size = double(nx)*ny*sizeof(SkyPixType)/pow(1024.0, 3);

    // Inline function to write a catalog to disk
    auto write_catalog = [&](std::string ofile, const vec1u& oids, uint_t sx, uint_t sy){
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

        std::string obase = file::remove_extension(ofile);

        // Write FITS header
        std::string hdr_file = obase+"-hdr.txt";
        std::ofstream ohdr(hdr_file);
        vec1s shdr = cut(hdr, 80);
        for (auto& s : shdr) {
            ohdr << s << "\n";
        }

        // Write SkyMaker configuration file
        std::ofstream oconf(obase+"-sky.conf");
        std::string img_name = file::get_basename(obase);

        // Copy the parameters from the template
        vec1s tpar = sky_param;
        vec1s tval = sky_value;

        // Add our own
        tpar.push_back("IMAGE_NAME");
        tval.push_back(img_dir+img_name+"-sci.fits");
        tpar.push_back("IMAGE_SIZE");
        tval.push_back(strn(sx)+","+strn(sy));
        tpar.push_back("IMAGE_HEADER");
        tval.push_back(hdr_file);

        // Write to file
        tpar = align_left(tpar, max(length(tpar)));
        for (uint_t i : range(tpar)) {
            oconf << tpar[i] << " " << tval[i] << "\n";
        }
    };

    vec1u tilex(x.size());
    vec1u tiley(x.size());

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

        if (verbose) {
            double sec_size = double(nx)*ny*sizeof(SkyPixType)/pow(1024.0, 3);
            note("will require ", nsx, "x", nsy, " sections (", sec_size, " GB each)");
            note("section dimensions: ", nx, ",", ny);
        }

        if (verbose) {
            note("write catalogs...");
        }

        tilex.resize(x.size());
        tiley.resize(x.size());

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
            vec1u idi = ids[where(ox[ids]+gsize[ids] >= nx*ix+0.5     &&
                                  ox[ids]-gsize[ids] <= nx*(ix+1)+0.5 &&
                                  oy[ids]+gsize[ids] >= ny*iy+0.5     &&
                                  oy[ids]-gsize[ids] <= ny*(iy+1)+0.5)];

            // This is the "strict" list, i.e., only galaxies whose center lies in this tile
            vec1u idis = where(ox[idi] >= nx*ix+0.5 && ox[idi] < nx*(ix+1)+0.5 &&
                               oy[idi] >= ny*iy+0.5 && oy[idi] < ny*(iy+1)+0.5);

            // Convert pixel coordinates to local frame
            wcs_params.pixel_ref_x = dx+1-nx*double(ix);
            wcs_params.pixel_ref_y = dy+1-ny*double(iy);

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

            auto spl = file::split_extension(out_file);
            write_catalog(spl.first+"-"+six+"-"+siy+spl.second, idi, nx, ny);

            tilex[idi[idis]] = ix+1;
            tiley[idi[idis]] = iy+1;

            if (verbose) progress(pg);
        }

        std::swap(tx, x);
        std::swap(ty, y);
    } else {
        // We will make a single big image
        if (verbose) {
            note("image dimensions: ", nx, ",", ny, " (", img_size, " GB)");
        }

        if (verbose) {
            note("write catalog...");
        }

        write_catalog(out_file, ids, nx, ny);

        tilex = replicate(1u, x.size());
        tiley = replicate(1u, x.size());
    }

    if (!save_pixpos.empty()) {
        note("write pixel positions...");
        file::mkdir(file::get_directory(save_pixpos));
        fits::write_table(save_pixpos, ftable(cat.id, x, y, tilex, tiley));
    }

    if (verbose) {
        note("done.");
    }

    // TODO: fix SkyMaker and border objects

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

    print("egg-2skymaker v1.0rc1");
    print("usage: egg-2skymaker cat=... band=... template=... [options]\n");

    print("List of mandatory parameters (no default):");
    argdoc("cat", "[string]", "path to the mock catalog (FITS table)");
    argdoc("band", "[string]", "name of the photometric band for which a SkyMaker "
        "input catalg will be created");
    argdoc("template", "[string]", "path to a template SkyMaker configuration file. "
        "This template file must contain the 'PIXEL_SIZE' parameter. If they are missing, "
        "this program will add 'WAVELENGTH', 'IMAGE_TYPE' and 'LISTCOORD_TYPE' (which has "
        "to be 'PIXELS'). It will also take care of supplying 'IMAGE_NAME', 'IMAGE_SIZE' "
        "and 'IMAGE_HEADER'. Note that, if the provided file does not exists, the "
        "program will also search in the list of the pre-built templates provided with "
        "EGG ("+egg_share_dir+"/skymaker-templates) for a file with the same name, and "
        "use it if it exists.");

    print("");
    print("List of options:");
    argdoc("out", "[string]", "name of the output catalog (default: [catalog]-[band].cat)");
    argdoc("maglim", "[float]", "galaxies fainter than this magnitude will not be "
        "present in the output catalog (default: none)");
    argdoc("inset", "[double, arcsec]", "add an empty border around the input catalog "
        "to prevent galaxies from being truncated by the edge of the image "
        "(default: 0\"). Note that even if this option is set to zero, an automatic inset "
        "is still added by the program, taking into account the extents of each galaxy "
        "that is close to the border. See 'strict_clip'.");
    argdoc("strict_clip", "[flag]", "if this flag is set, the boundaries of the image will "
        "be computed simply from the coordinates of the galaxies in the input catalog. By "
        "default, the program also takes into account their sizes, to make sure that no "
        "galaxy is severely truncated close to the edge of the image.");
    argdoc("size_cap", "[float, GB]", "split the survey in multiple tiles of equal area "
        "that will each be rendered into a separate image, so that the size of each "
        "such image is less than 'size_cap' (expressed in GigaBytes, default: none). NB: "
        "SkyMaker requires the whole image to fit inside your RAM memory, so this value "
        "should be no larger than the total amount of RAM memory available on your "
        "computer.");
    argdoc("save_pixpos", "[string]", "save the pixel coordinates and the name of the "
        "tile of each galaxy in this file (default: not saved). This is most useful when "
        "'size_cap' is used, as this catalog provides and easy way to figure out in which "
        "of the tile each galaxy is located.");
    argdoc("verbose", "[flag]", "display information about the progress of the program "
        "in the terminal");
    argdoc("list_templates", "[flag]", "print the list of available pre-built templates "
        "and exit");

    print("");
}
