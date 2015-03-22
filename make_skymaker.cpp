#include <phypp.hpp>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print("usage: make_skymaker gencat.fits band=... maglim=... out=...");
        return 1;
    }

    std::string band;
    std::string out_file;
    double maglim = dnan;
    double aspix = 0.06;

    read_args(argc-1, argv+1, arg_list(band, name(out_file, "out"), maglim, aspix));

    struct {
        vec1d ra, dec;

        vec1f disk_angle,  disk_radius,  disk_ratio;
        vec1f bulge_angle, bulge_radius, bulge_ratio;

        vec2f flux_disk, flux_bulge;
        vec1s bands;
    } cat;

    fits::read_table(argv[1], ftable(
        cat.ra, cat.dec,
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

    uint_t b = idb[0];

    vec1f mag = uJy2mag(cat.flux_disk(_,b)+cat.flux_bulge(_,b));
    vec1f bt = cat.flux_bulge(_,b)/(cat.flux_disk(_,b)+cat.flux_bulge(_,b));

    fits::make_wcs_header_params wcs_params;
    wcs_params.pixel_scale = aspix;
    wcs_params.pixel_ref_x = 0.0; wcs_params.pixel_ref_y = 0.0;
    wcs_params.sky_ref_ra = max(cat.ra); wcs_params.sky_ref_dec = min(cat.dec);
    move_ra_dec(wcs_params.sky_ref_ra, wcs_params.sky_ref_dec, +5.0, -5.0);

    fits::header hdr;
    if (!fits::make_wcs_header(wcs_params, hdr)) {
        error("could not make WCS header");
        return 1;
    }

    vec1d x, y;
    fits::ad2xy(fits::wcs(hdr), cat.ra, cat.dec, x, y);

    vec1u ids;
    if (finite(maglim)) {
        ids = where(mag < maglim);
        if (ids.empty()) {
            error("no source brighter than ", maglim);
            note("maximum magnitude is ", max(mag));
            return 1;
        }
    } else {
        ids = uindgen(x.size());
        vec1u idbad = where(!finite(mag));
        mag[idbad] = 99;
    }

    if (out_file.empty()) {
        out_file = erase_end(argv[1], ".fits")+"-"+band+".cat";
    }

    file::mkdir(file::get_directory(out_file));

    file::write_table_hdr(out_file, 16,
        {"type", "x", "y", "mag", "bt",
        "bulge_radius", "bulge_ratio", "bulge_angle",
        "disk_radius", "disk_ratio", "disk_angle"},
        replicate(200u, ids.size()), x[ids], y[ids],
        mag[ids], bt[ids],
        cat.bulge_radius[ids], cat.bulge_ratio[ids], cat.bulge_angle[ids],
        cat.disk_radius[ids], cat.disk_ratio[ids], cat.disk_angle[ids]
    );

    // Write FITS header to a file
    {
        std::ofstream ohdr(erase_end(out_file, ".cat")+"-hdr.txt");
        vec1s shdr = cut(hdr, 80);
        for (auto& s : shdr) {
            ohdr << s << "\n";
        }
    }

    return 0;
}
