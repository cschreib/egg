#include <vif.hpp>
#include "egg-utils.hpp"

using namespace vif;
using namespace vif::astro;

const std::string egg_share_dir = file::directorize(EGG_SHARE_DIR);
const std::string filters_dir_env = file::directorize(system_var("EGG_FILTERS_PATH", ""));
const std::string filters_dir = filters_dir_env.empty() ?
    file::directorize(FILTER_DB_DIR) : filters_dir_env;

void print_help();

int vif_main(int argc, char* argv[]) {
    // Initialization
    // --------------

    // Define the survey position and boundaries
    // Warning: positions too close to the angular poles should be avoided.
    double ra0 = 53.558750;   // in degrees
    double dec0 = -27.176001; // in degrees
    double area = dnan;       // in degrees^2

    // Define the stellar mass range and binning
    double mmin = dnan;
    double mmax = 12.0;
    double bin_dm = 0.05;

    // Define instead a fixed magnitude cut
    double maglim = dnan;
    std::string selection_band;

    // Define the redshift range
    double zmin = 0.05;
    double zmax = 10.5;
    // Define the redshift binning
    double bin_dz = 0.05;
    double max_dz = 0.5;
    double min_dz = 0.001;

    // Clustering parameters
    double clust_r0 = 0.05;          // clustering outer truncation radius in degree
    double clust_r1 = 0.0006;        // clustering inner truncation radius in degree
    uint_t clust_eta = 5;            // number of sub haloes per halo
    double clust_lambda = 6.0;       // radius shrinking factor of subhaloes
    double clust_fclust_him = 0.4;   // fraction of clustered vs uniform positions for
                                     // high-mass galaxies
    double clust_fclust_lom = 0.2;   // fraction of clustered vs uniform positions for
                                     // low-mass galaxies
    double clust_fclust_mlim = 10.5; // threshold between high and low mass for the above
    double clust_urnd_mlim = 8.0;    // threshold in mass below which there is no clustering

    // Dispersion of the Main Sequence in SFR, at fixed Mstar
    // Changing this value will alter the quality of the simulation.
    // Only do this if you know what you are doing.
    double ms_disp = 0.3;

    // Disable some components of the simulation
    bool no_pos = false;         // do not generate positions
    bool no_clust = false;       // do not generate clustering
    bool no_flux = false;        // do not generate fluxes
    bool no_stellar = false;     // do not generate fluxes from stellar origin
    bool no_nebular = false;     // do not generate fluxes from nebular origin
    bool no_dust = false;        // do not generate fluxes from dust origin
    bool no_random = false;      // disable all randomization of parameters
                                 // (just passive, M*, z, and position + pos angle)

    // Absorption form intergalactic medium (IGM)
    // Valid values: none, constant, madau95, inoue14
    std::string igm = "inoue14";
    std::string naive_igm; // [deprecated]


    // Debug - testing
    bool no_passive_lir = false; // do not generate dust for passive galaxies
    bool magdis_tdust = false;   // use Magdis+12 Tdust evolution (lower)


    // Save the full spectrum of each galaxy to a file
    // Warning: the current implementation will consume a lot of RAM memory
    bool save_sed = false;
    std::string seds_file = "";

    // Mass function file
    std::string mass_func_file = egg_share_dir+"mass_func_candels.fits";

    // SED libraries
    std::string ir_lib_file = egg_share_dir+"ir_lib_cs17.fits";
    std::string opt_lib_file;
    std::string opt_lib_imf = "salpeter";
    double opt_lib_min_step = 0.001; // um

    // Filter library
    std::string filter_db_file = filters_dir+"db.dat";
    bool filter_flambda = false; // set to true to get equivalent of FAST FILTER_FORMAT=1
    bool filter_photons = false; // set to true to get equivalent of FAST FILTER_FORMAT=1
    bool filter_normalize = true;
    bool trim_filters = false;

    // Output file
    std::string out_file = "";

    // Existing catalog to start from (empty: none)
    std::string input_cat_file = "";

    // Seed for random number generation
    uint_t tseed = 42;

    // Cosmological parameters
    // Warning: even if these are wrong or inaccurate, they should not be changed
    // since the whole simulation relies on these parameters. The output observables,
    // i.e., fluxes and morphologies, do not assume any cosmology, and can be used
    // regardless of your preference. This is not the case, however, for the other
    // intermediate quantities like Mstar or SFR.
    std::string tcosmo = "std";

    // Display information about the simulation process in the terminal
    bool verbose = false;

    // Display some help
    bool help = false;
    std::string list_bands;

    // Photometric bands for which to generate observed fluxes
    vec1s bands = {"vimos-u", "hst-f435w", "hst-f606w", "hst-f775w",
        "hst-f814w", "hst-f850lp", "hst-f105w", "hst-f125w", "hst-f140w",
        "hst-f160w", "hawki-Ks", "spitzer-irac1", "spitzer-irac2", "spitzer-irac3",
        "spitzer-irac4", "spitzer-irs16", "spitzer-mips24", "herschel-pacs70",
        "herschel-pacs100", "herschel-pacs160", "herschel-spire250", "herschel-spire350",
        "herschel-spire500"};

    // Photometric bands for which to generate absolute magnitudes
    vec1s rfbands;

    // Read command line arguments
    // ---------------------------
    read_args(argc, argv, arg_list(
        ra0, dec0, area, mmin, maglim, zmin, zmax, name(bin_dz, "dz"), min_dz, max_dz,
        name(bin_dm, "dm"), ms_disp,
        no_pos, no_clust, no_flux, no_stellar, no_dust, no_passive_lir, no_random, no_nebular,
        save_sed, name(mass_func_file, "mass_func"),
        name(ir_lib_file, "ir_lib"), name(opt_lib_file, "opt_lib"), opt_lib_imf,
        name(out_file, "out"), name(filter_db_file, "filter_db"),
        verbose, name(tseed, "seed"), name(tcosmo, "cosmo"),
        name(input_cat_file, "input_cat"), selection_band, bands, rfbands, help, list_bands,
        clust_r0, clust_r1, clust_lambda, clust_eta, clust_fclust_mlim, clust_fclust_lom,
        clust_fclust_him, clust_urnd_mlim, magdis_tdust, igm, naive_igm, filter_photons,
        filter_flambda, trim_filters, filter_normalize, opt_lib_min_step
    ));

    if (help) {
        print_help();
        return 0;
    }

    if (!list_bands.empty()) {
        if (list_bands == "1") {
            print("List of available bands:");
        } else {
            print("List of available bands (filter: '", list_bands, "'):");
        }

        if (!file::exists(filter_db_file)) {
            error("could not find filter database file '", filter_db_file, "'");
            return 1;
        } else {
            vec1s fils; vec1f rlam, width;
            auto filter_db = read_filter_db(filter_db_file);
            for (auto fil : filter_db) {
                // Apply the regex provided by the user to only disply some of the filters
                if (list_bands != "1" && !regex_match(fil.first, list_bands)) {
                    continue;
                }

                // Read filter
                vec1d lam, res;
                fits::read_table(fil.second, ftable(lam, res));

                if (lam.size() != res.size()) {
                    warning("mismatch in wavelength and response for filter ", fil.first);
                    continue;
                }

                if (max(res) <= 0) {
                    warning("zero transmission for filter ", fil.first);
                    continue;
                }

                if (res.empty()) {
                    warning("filter ", fil.first, " is empty");
                    continue;
                }

                // Compute reference wavelength
                double trlam = integrate(lam, lam*res);

                // Compute bandwidth
                vec1u idg = where(res > 0.5*max(res));
                double twidth = lam[idg.back()] - lam[idg.front()];

                fils.push_back(fil.first);
                rlam.push_back(trlam);
                width.push_back(twidth);
            }

            // Sort the filter list by instrument & reference wavelength
            vec1u ids = indgen(fils.size());
            inplace_sort(ids, [&](uint_t i1, uint_t i2) {
                std::string inst1, inst2;
                inst1 = split(fils[i1], "-")[0];
                inst2 = split(fils[i2], "-")[0];
                if (inst1 != inst2) return inst1 < inst2;
                else                return rlam[i1] < rlam[i2];
            });

            // Print
            vec1s srlam = to_string_vector(rlam)+" um, ";
            vec1s swidth = to_string_vector(width);
            vec1s item = " - "+align_left(fils, max(length(fils)))+
                "  ref-lam = "+align_left(srlam, max(length(srlam)))+
                "FWHM = "+swidth+" um";
            for (uint_t i : ids) {
                print(item[i]);
            }
            return 0;
        }
    }

    if (!is_finite(area) && input_cat_file.empty()) {
        error("please specify either the area of the mock catalog (area=...), or provide");
        error("an input catalog with galaxy redshifts, masses, UVJ flag and positions "
            "(input_cat=...)");
        note("type 'egg-gencat help' for more information");
        return 1;
    }

    if (!(is_finite(maglim) && !selection_band.empty()) && !is_finite(mmin) && input_cat_file.empty()) {
        error("please specify either the minimum stellar mass (mmin=...), or a magnitude "
            "limit (maglim=..., selection_band=...), or provide an input catalog with "
            "galaxy redshifts, masses, UVJ flag and positions (input_cat=...)");
        note("type 'egg-gencat help' for more information");
        return 1;
    }

    if (zmin <= 0.0) {
        error("minimum redshift must be > 0 (zmin=...)");
        return 1;
    }

    if (min_dz > max_dz) {
        error("min_dz must be larger than max_dz (got ", min_dz, " and ", max_dz, ")");
        return 1;
    }

    if (no_dust && no_stellar) {
        no_flux = true;
    }

    if (out_file.empty()) {
        out_file = "egg-"+today()+".fits";
    }

    if (verbose) {
        note("will save catalog in '", out_file, "'");
    }

    if (save_sed && seds_file.empty()) {
        seds_file = file::remove_extension(out_file)+"-seds.dat";

        if (verbose) {
            note("will save spectra in '", seds_file, "'");
        }
    }

    if (!naive_igm.empty()) {
        warning("the 'naive_igm' option is now deprecated");
        warning("please use the 'igm=...' option instead");

        if (naive_igm == "1") {
            igm = "constant";
        } else {
            igm = "madau95";
        }
    }

    if (igm != "none" && igm != "constant" && igm != "madau95" && igm != "inoue14") {
        error("unknown IGM prescription '", igm, "'");
        error("possible values are: 'none', 'constant', 'madau95', or 'inoue14'");
        return 1;
    } else {
        if (verbose) {
            note("using IGM prescription '", igm, "'");
        }
    }

    if (opt_lib_file.empty()) {
        if (igm == "constant") {
            opt_lib_file = egg_share_dir+"opt_lib_fast.fits";
        } else {
            opt_lib_file = egg_share_dir+"opt_lib_fast_noigm.fits";
        }
    }

    // Initialize random seed and cosmology
    // ------------------------------------
    auto seed = make_seed(tseed);
    auto cosmo = get_cosmo(tcosmo);

    // Initizalize filters
    // -------------------

    if (verbose) {
        note("initializing filters...");
    }

    // Get filters
    auto filter_db = read_filter_db(filter_db_file);

    vec<1,filter_t> filters;
    if (!get_filters(filter_db, bands, filters)) {
        return 1;
    }

    vec<1,filter_t> rffilters;
    if (!get_filters(filter_db, rfbands, rffilters)) {
        return 1;
    }

    // Perform adjustments on the filter response curves, if asked
    auto adjust_filter = vectorize_lambda([&](filter_t& fil) {
        // Truncate
        if (trim_filters) {
            vec1u idg = where(fil.res/max(fil.res) > 1e-3);
            vif_check(!idg.empty(), "filter has no usable data");
            fil.lam = fil.lam[idg[0]-_-idg[-1]];
            fil.res = fil.res[idg[0]-_-idg[-1]];
        }

        // Apply filter definition
        if (filter_flambda) {
            // Filter is defined such that it must integrate f_lambda and not f_nu
            // f_lambda*r(lambda) ~ f_nu*[r(lambda)/lambda^2]
            fil.res /= sqr(fil.lam);
        }
        if (filter_photons) {
            // Filter is defined such that it integrates photons and not energy
            // n(lambda)*r(lambda) ~ f(lambda)*[r(lambda)*lambda]
            fil.res *= fil.lam;
        }

        // Re-normalize filter
        double norm = integrate(fil.lam, fil.res);
        fil.rlam = integrate(fil.lam, fil.lam*fil.res)/norm;
        if (filter_normalize) {
            fil.res /= norm;
        }

        return fil;
    });

    if (trim_filters || filter_flambda || filter_photons) {
        adjust_filter(filters);
        adjust_filter(rffilters);
    }

    // Sort the bands by increasing wavelength
    auto filter_lam = vectorize_lambda([](const filter_t& f) { return f.rlam; });
    vec1f lambda = filter_lam(filters);
    vec1u ids = sort(lambda);
    bands = bands[ids];
    lambda = lambda[ids];
    filters = filters[ids];

    vec1f rflambda = filter_lam(rffilters);
    ids = sort(rflambda);
    rfbands = rfbands[ids];
    rflambda = rflambda[ids];
    rffilters = rffilters[ids];


    // Initialize SED libraries
    // ------------------------

    if (verbose) {
        note("initializing SED libraries...");
    }

    // The IR library. When using the generic code for IR libraries, each template must
    // be normalized to unit LIR. If a more subtle calibration needs to be used, then
    // a custom code has to be written.
    struct {
        vec2d lam, sed;
    } ir_lib;

    struct {
        vec2d lam, dust, pah;
        vec1d tdust, lir_dust, lir_pah, l8_dust, l8_pah;
    } ir_lib_cs17;

    uint_t nirsed;
    if (file::get_basename(ir_lib_file) == "ir_lib_cs17.fits") {
        // My library, calibrated in unit Mdust
        fits::read_table(ir_lib_file, ftable(
            ir_lib_cs17.lam, ir_lib_cs17.dust, ir_lib_cs17.pah,
            ir_lib_cs17.tdust, ir_lib_cs17.lir_dust, ir_lib_cs17.lir_pah,
            ir_lib_cs17.l8_dust, ir_lib_cs17.l8_pah
        ));

        nirsed = ir_lib_cs17.lam.dims[0];
    } else {
        // Generic library
        auto cols = fits::read_table_columns(ir_lib_file);
        for (auto& c : cols) {
            if (c.name == "SED") {
                if (c.dims.size() == 1) {
                    struct {
                        vec1d lam, sed;
                    } tmp;

                    fits::read_table(ir_lib_file, ftable(tmp.lam, tmp.sed));

                    ir_lib.lam.resize(1, tmp.lam.size());
                    ir_lib.sed.resize(1, tmp.lam.size());

                    ir_lib.lam(0,_) = tmp.lam;
                    ir_lib.sed(0,_) = tmp.sed;
                } else if (c.dims.size() == 2) {
                    fits::read_table(ir_lib_file, ftable(ir_lib.lam, ir_lib.sed));
                } else {
                    error("IR library must contain either 1D or 2D columns LAM and SED");
                    return 1;
                }

                break;
            }
        }

        if (ir_lib.sed.empty() || ir_lib.lam.empty()) {
            error("missing LAM and/or SED columns in the IR library file");
            return 1;
        }

        nirsed = ir_lib.sed.dims[0];
    }

    // The optical library. Each template must be normalized to solar luminosity
    // per unit stellar mass (in solar mass).
    // The library is binned in redshift, U-V and V-J colors.
    struct {
        vec3f lam, sed;
        vec2b use;
        vec2f av;
        vec2f bvj, buv;
    } opt_lib;

    fits::read_table(opt_lib_file, ftable(
        opt_lib.lam, opt_lib.sed, opt_lib.use, opt_lib.av, opt_lib.bvj, opt_lib.buv
    ));

    // Adjust IMF if necessary
    if (opt_lib_imf == "chabrier" || opt_lib_imf == "kroupa") {
        opt_lib.sed *= e10(-0.2);
    } else if (opt_lib_imf != "salpeter") {
        error("unknown IMF '", opt_lib_imf, "'");
        return 1;
    }

    // Adjust resolution if necessary
    if (opt_lib_min_step > 0) {
        if (verbose) {
            note("rebinning optical library...");
        }

        // Make a copy of library
        vec1u idf = mult_ids(opt_lib.use, where_first(opt_lib.use));
        vec1d olam = opt_lib.lam(idf[0],idf[1],_);
        vec3d osed = opt_lib.sed;

        // Compute resolution of current grid
        uint_t onpt = olam.size();
        vec1d dlam(onpt);
        dlam[0] = olam[1] - olam[0];
        dlam[1-_-(onpt-2)] = 0.5*(olam[2-_] - olam[0-_-(onpt-3)]);
        dlam[onpt-1] = olam[onpt-1] - olam[onpt-2];

        // Find range where resolution is too high
        uint_t id1 = where_first(dlam < opt_lib_min_step);
        uint_t id2 = where_last(dlam < opt_lib_min_step);

        if (id1 != npos) {
            // Define new grid
            vec1d nlam;

            // Copy first elements that are already at low resolution
            if (id1 != 0) {
                append(nlam, olam[_-(id1-1)]);
            }

            // Create low resolution grid
            uint_t nnew = floor((olam[id2] - olam[id1])/opt_lib_min_step) - 1;
            vec1d newlam = olam[id1] + 0.5*opt_lib_min_step + opt_lib_min_step*indgen<double>(nnew);
            append(nlam, newlam);
            uint_t nid2 = nlam.size()-1;

            // Copy last elements that are already at low resolution
            if (id2 != onpt-1) {
                append(nlam, olam[(id2+1)-_]);
            }

            if (verbose) {
                note("reduced SED samples from ", olam.size(), " to ", nlam.size());
            }

            // Regrid SEDs
            opt_lib.lam = replicate(0.0f, opt_lib.use.dims,nlam.size());
            opt_lib.sed = replicate(0.0f, opt_lib.use.dims,nlam.size());
            for (uint_t iuv : range(opt_lib.use.dims[0]))
            for (uint_t ivj : range(opt_lib.use.dims[1])) {
                if (!opt_lib.use(iuv,ivj)) continue;

                opt_lib.lam(iuv,ivj,_) = nlam;
                if (id1 != 0) {
                    opt_lib.sed(iuv,ivj,_-(id1-1)) = osed(iuv,ivj,_-(id1-1));
                }
                for (uint_t l : range(newlam)) {
                    opt_lib.sed(iuv,ivj,id1+l) = integrate(olam, osed(iuv,ivj,_),
                        newlam[l]-0.5*opt_lib_min_step, newlam[l]+0.5*opt_lib_min_step)/opt_lib_min_step;
                }
                if (id2 != onpt-1) {
                    opt_lib.sed(iuv,ivj,(nid2+1)-_) = osed(iuv,ivj,(id2+1)-_);
                }
            }
        }
    }

    // Make sure that it contains at least one valid SED
    if (count(opt_lib.use) == 0) {
        error("the provided optical SED library does not contain any valid SED");
        note("'use' was false for all the UVJ positions");
        return 1;
    }

    // Manual tuning of the M/L ratio as a function of z
    auto get_m2l_cor = vectorize_lambda([](double z) {
       return interpolate({0.15,0.15,0.0,0.0,-0.6}, {0.0,0.45,1.3,6.0,8.0}, z);
    });

    // Output catalog format
    // ---------------------

    struct {
        vec1u id;

        // Sky properties
        vec1d ra, dec;
        vec1b clustered;

        // Physical properties
        vec1f z, d, m, sfr, rsb;

        vec1f disk_angle,  disk_radius,  disk_ratio;
        vec1f bulge_angle, bulge_radius, bulge_ratio;
        vec1f bt, m_disk, m_bulge;
        vec1f av_disk, av_bulge;

        vec1f tdust, ir8, fpah, mdust;

        vec1b passive;

        // Flux properties
        vec1f rfuv_bulge, rfuv_disk, rfvj_bulge, rfvj_disk;
        vec1f sfrir, sfruv, irx;
        vec1f m2lcor, lir;
        vec1u ir_sed;
        vec1u opt_sed_bulge;
        vec1u opt_sed_disk;
        vec2f igm_abs;

        vec2f flux;
        vec2f flux_disk;
        vec2f flux_bulge;
        vec1s bands;
        vec1f lambda;

        vec2f rfmag;
        vec2f rfmag_disk;
        vec2f rfmag_bulge;
        vec1s rfbands;
        vec1f rflambda;

        vec1f vdisp;
        vec1f avlines_disk;
        vec1f avlines_bulge;
        vec2f line_lum;
        vec2f line_lum_disk;
        vec2f line_lum_bulge;
        vec1s lines;
        vec1f line_lambda;

        vec2f zb;

        std::string cmd;
    } out;

    out.lambda = lambda;
    out.bands = bands;
    out.rflambda = rflambda;
    out.rfbands = rfbands;

    // Generate redshift, masses, and positions
    // ----------------------------------------

    // Either generate everything from scratch, or provide
    // all these data in input, i.e., from an existing catalog.

    // First define the functions to generate UVJ colors to have a prior on colors
    // This is used for computing the stellar mass limit for a given magnitude cut.
    // It is also used later to actually generate the colors of the simulated galaxies.
    auto gen_blue_uvj = [&seed](const vec1f& m, const vec1f& z, vec1f& uv, vec1f& vj, bool norand) {
        // Calibrate "UVJ vector" from mass and redshift
        vec1f a0 = 0.58*erf(m-10) + 1.39;
        vec1f as = -0.34 + 0.3*max(m-10.35, 0.0);
        vec1f a = min(a0 + as*min(z, 3.3), 2.0);

        if (!norand) {
            vec1d rnd_amp = 0.2 + (0.25 + 0.12*clamp((z-0.5)/2.0, 0, 1))*
                max(1.0 - 2.0*abs(m - (10.3 + 0.4*erf(z - 1.5))), 0.0);
            a += rnd_amp*randomn(seed, m.size());
        }

        // Move in the UVJ diagram according to the UVJ vector
        double slope = 0.65;
        double theta = atan(slope);
        vj = 0.0  + a*cos(theta);
        uv = 0.45 + a*sin(theta);

        if (!norand) {
            vj += 0.15*randomn(seed, m.size());
            uv += 0.15*randomn(seed, m.size());
        }

        // Add an additional global color offset depending on redshift
        uv += 0.4*max((0.5 - z)/0.5, 0.0);
        vj += 0.2*max((0.5 - z)/0.5, 0.0);
    };

    auto gen_red_uvj = [&seed](const vec1f& m, const vec1f& z, vec1f& uv, vec1f& vj, bool norand) {
        double pvj = 1.25, puv = 1.85;
        vec1f mspread = 0.1*(m - 11.0);

        if (!norand) {
            mspread += 0.1*randomn(seed, m.size());
        }

        mspread = clamp(mspread, -0.1, 0.2);

        vj = pvj + mspread;
        uv = puv + mspread*0.88;

        if (!norand) {
            vj += 0.1*randomn(seed, m.size());
            uv += 0.1*randomn(seed, m.size());
        }

        // Add an additional global color offset depending on redshift
        uv += 0.4*max((0.5 - z)/0.5, 0.0);
        vj += 0.2*max((0.5 - z)/0.5, 0.0);
    };

    auto get_sed_uvj = [&opt_lib](const vec1f& uv, const vec1f& vj, vec1u& sed, vec1f& av, bool doprint) {
        uint_t snb = opt_lib.buv.dims[1];

        sed.resize(uv.size());
        av.resize(uv.size());
        sed[_] = npos;
        av[_] = fnan;

        auto pg = progress_start(snb*snb);
        for (uint_t u  : range(snb)) {
            vec1b ibu = in_bin_open(uv, opt_lib.buv, u);
            for (uint_t v  : range(snb)) {
                // Find galaxies which are in this bin
                vec1u idl = where(ibu && in_bin_open(vj, opt_lib.bvj, v));

                if (!idl.empty()) {
                    // Find the closest usable SED in the library
                    uint_t lu = u, lv = v;
                    if (!astar_find(opt_lib.use, lu, lv)) {
                        error("could not find good SEDs in the optical library!");
                        note("this should not have happened...");
                        note(u, ", ", v);
                        return false;
                    }

                    // Assign the SED to these galaxies
                    uint_t ised = flat_id(opt_lib.use, lu, lv);
                    sed[idl] = ised;
                    av[idl] = opt_lib.av[ised];
                }

                if (doprint) progress(pg);
            }
        }

        return true;
    };

    auto get_opt_sed = [&opt_lib](double m, uint_t ised, vec1f& orlam, vec1f& orsed) {
        vec1u did = mult_ids(opt_lib.use, ised);
        orlam = opt_lib.lam(did[0],did[1],_);
        orsed = e10(m)*opt_lib.sed(did[0],did[1],_);
    };

    if (!input_cat_file.empty()) {
        if (verbose) {
            note("reading input catalog '", input_cat_file, "'...");
        }

        // Load data from an existing catalog
        if (ends_with(input_cat_file, ".fits")) {
            fits::input_table tbl(input_cat_file);

            // Read main parameters
            tbl.read_columns(fits::narrow, ftable(
                out.ra, out.dec, out.z, out.m
            ));

            // See if there is an 'ID' column, else create our own
            if (!tbl.read_column("id", out.id)) {
                out.id = indgen(out.ra.size());
            }

            // See if there is either 'passive' or 'passive_prob'
            vec1f passive_prob;
            tbl.read_columns(fits::missing | fits::narrow, ftable(
                out.lir, out.tdust, out.ir8, out.passive, passive_prob
            ));

            if (out.passive.empty()) {
                if (passive_prob.empty()) {
                    error("you must provide the 'passive' flag in the input catalog, "
                        "or prodive the 'passive_prob' value");
                    return 1;
                }

                // Generate passive flag according to the given probability
                out.passive = random_coin(seed, passive_prob);
            }

            // See if there is data on the IR spectrum
            tbl.read_columns(fits::missing | fits::narrow, ftable(
                out.lir, out.tdust, out.ir8
            ));

            if (verbose) {
                vec1s found;
                if (!out.lir.empty())   found.push_back("LIR");
                if (!out.tdust.empty()) found.push_back("Tdust");
                if (!out.ir8.empty())   found.push_back("IR8");

                if (!found.empty()) {
                    note("found user provided values for "+collapse(found, ", ")+" in the input catalog");
                }
            }
        } else {
            ascii::read_table(input_cat_file, out.id, out.ra, out.dec, out.z, out.m, out.passive);
        }

        // Verify input
        uint_t nbad = count(out.z <= 0 || out.m <= 4 || out.m > 13);
        if (nbad != 0) {
            error("input catalog contains ", nbad, " unphysical galax", (nbad == 1 ? "y" : "ies"));
            note("must have z > 0 and 4 < m < 13");
            return 1;
        }

        if (verbose) {
            note("found ", out.z.size(), " galaxies");
        }

        out.d = lumdist(out.z, cosmo);
    } else {
        // Generate all from scratch

        // Initialize redshift bins
        // ------------------------
        if (verbose) {
            note("initializing redshift bins...");
        }

        if (bin_dz < 0.0) {
            error("dz cannot be negative");
            return 1;
        } else if (bin_dz >= 2.0) {
            error("dz cannot be >= 2.0");
            return 1;
        }

        // Build the redshift bins with logarithmic bins
        // To prevent too narrow bins with very few galaxies, we can impose a minimum
        // redshift width (this feature is only used when making the clustering slices)
        // or, conversely, a maximum redshift width to avoid generating too large bins.
        auto make_zbins = [](double zstart, double zend, double dz, double midz, double madz) {
            vec1d tzb;

            vif_check(!is_nan(midz) || !is_nan(midz) || midz < madz,
                "min_dz must be smaller than max_dz (got ", midz, " and ", madz, ")");

            // Generate redshifts incrementally using z(i+1) = z(i)*(1+dz) until we reach zend
            double z1 = zstart;
            while (z1 < zend) {
                tzb.push_back(z1);
                z1 *= 1.0 + dz;
                if (is_finite(midz) && z1 - tzb.back() < midz) {
                    z1 = tzb.back() + midz;
                } else if (is_finite(madz) && z1 - tzb.back() > madz) {
                    z1 = tzb.back() + madz;
                }
            }

            // Now determine wether we need an extra bin for zend, or if we can just enlarge
            // the last bin to include it.
            double min_end_dz = 0.5*zend*dz;
            if (is_finite(midz)) {
                min_end_dz = min(min_end_dz, midz);
            }

            if (tzb.size() > 1 && zend - tzb.back() < min_end_dz) {
                // Make the last bin a bit bigger to include zend
                tzb.back() = zend;
            } else {
                // Add zend to make the last bin
                tzb.push_back(zend);
            }

            return make_bins_from_edges(tzb);
        };

        vec2d zb = make_zbins(zmin, zmax, bin_dz, min_dz, max_dz);
        uint_t nz = zb.dims[1];

        if (verbose) {
            note(nz, " redshift slices");
            vec1f tdz = bin_width(zb);
            note("min dz: ", min(tdz), ", max dz: ", max(tdz));
        }

        // Compute the corresponding volume of the Universe
        vec2d vb = vuniverse(zb, cosmo);

        // Initialize mass functions
        // -------------------------
        if (verbose) {
            note("reading mass functions...");
        }

        struct {
            vec2f zb, mb;
            vec1f zx, mx;
            vec2d mz_active, mz_passive;
            vec1d z_active, z_passive;
            std::string imf;
        } mf;

        fits::read_table(mass_func_file, "zb", mf.zb, "mb", mf.mb, "active", mf.mz_active,
            "passive", mf.mz_passive);
        fits::read_table_loose(mass_func_file, "imf", mf.imf);

        // Correct for the IMF, if not Salpeter
        if (mf.imf.empty()) {
            warning("no IMF specification given in the mass function file, assuming Salpeter");
            mf.imf = "salpeter";
        } else {
            if (mf.imf == "salpeter") {
                // Nothing to do
            } else if (mf.imf == "chabrier") {
                mf.mb += 0.24;
            }
        }

        if (verbose) {
            note("found ", mf.zb.dims[1], " redshift bins and ", mf.mb.dims[1], " mass bins");
        }

        double mmin_strict = min(mf.mb);

        // Check that the provided mass functions cover all the parameter space we need
        if (max(mf.zb) < max(zb)) {
            error("mass functions do not cover redshifts > ", max(mf.zb));
            note("requested: ", max(zb));
            return 1;
        }
        if (min(mf.zb) > min(zb)) {
            error("mass functions do not cover redshifts < ", min(mf.zb));
            note("requested: ", min(zb));
            return 1;
        }

        // Compute evolving mass limit
        // ---------------------------
        vec1f z_mmin(nz);

        if (is_finite(maglim)) {
            if (verbose) {
                note("estimating redshift-dependent mass limit...");
            }

            filter_t filsel;
            if (!get_filter(filter_db, selection_band, filsel)) {
                return 1;
            }

            // A magnitude limit was requested.
            // Use the relations between mass and UVJ colors to have an estimate of the typical
            // M/L ratio (hence flux) of a galaxy of a given mass at any redshift, and use that
            // to estimate the typical minimum mass that is observed for a given flux limit.
            // We add a safety margin to this estimate to take into account the scatter in the
            // M/L ratios. Note that here we want to be complete at the requested magnitude limit,
            // which means we should go down in mass such that the brightest galaxies at that mass
            // are not seen.

            double flim = mag2uJy(maglim);

            auto pg = progress_start(nz);
            for (uint_t z : range(nz)) {
                double tz = bin_center(zb(_,z));
                double td = lumdist(tz, cosmo);

                // Generate M* - flux relation
                // First generate M* and z
                vec1f tmx = rgen(mmin_strict, mmax, 1000);
                vec1f tzx = replicate(tz, tmx.size());
                // Then UVJ colors
                vec1f tuvx, tvjx;
                gen_blue_uvj(tmx, tzx, tuvx, tvjx, true);
                // Then find the corresponding SED in the library
                vec1u tsx;
                vec1f tav;
                get_sed_uvj(tuvx, tvjx, tsx, tav, false);
                // And estimate the flux in the selection band
                double mlcor = get_m2l_cor(tz);
                vec1d tflx(tmx.dims);
                for (uint_t it : range(tmx)) {
                    vec1f rlam, rsed;
                    get_opt_sed(tmx[it]-mlcor, tsx[it], rlam, rsed);

                    vec1f lam = rlam*(1.0 + tz);
                    vec1f sed = lsun2uJy(tz, td, rlam, rsed);

                    tflx[it] = sed2flux(filsel, lam, sed);
                }

                // Minimum mass is obtained from the M* - flux relation
                // And we add a safety margin
                uint_t idf = where_first(tflx - flim > 0);
                if (idf != npos) {
                    z_mmin[z] = max(mmin_strict, tmx[idf] - 0.2);
                } else {
                    z_mmin[z] = mmax;
                }

                if (verbose) progress(pg);
            }

            mmin = min(z_mmin);

            if (verbose) {
                note("will generate masses from as low as ", mmin, ", up to ", mmax);
            }
        } else {
            // If no magnitude limit is requested, just use a fixed mass limit.
            z_mmin[_] = mmin;
        }

        {
            // Rebin the mass functions to our chosen redshift and mass bins
            vec1d old_zx = bin_center(mf.zb);
            vec1d old_mx = bin_center(mf.mb);

            mf.zb = zb;
            mf.zx = bin_center(mf.zb);

            mf.mb = make_bins(mmin, mmax, ceil((mmax - mmin)/bin_dm));
            mf.mx = bin_center(mf.mb);

            mf.mz_active  = e10(rebin(log10(mf.mz_active),  old_zx, old_mx, mf.zx, mf.mx));
            mf.mz_passive = e10(rebin(log10(mf.mz_passive), old_zx, old_mx, mf.zx, mf.mx));

            // Remove extrapolations at the borders of the bins
            vec1u idzb = where(mf.zx < old_zx.front());
            for (uint_t iz : idzb) mf.mz_active(iz,_)  = mf.mz_active(idzb.back()+1,_);
            for (uint_t iz : idzb) mf.mz_passive(iz,_) = mf.mz_passive(idzb.back()+1,_);

            idzb = where(mf.zx > old_zx.back());
            for (uint_t iz : idzb) mf.mz_active(iz,_)  = mf.mz_active(idzb.front()-1,_);
            for (uint_t iz : idzb) mf.mz_passive(iz,_) = mf.mz_passive(idzb.front()-1,_);

            // Build redshift functions by integrating over stellar mass
            vec2d tvb = vuniverse(mf.zb, cosmo);
            mf.z_active.resize(mf.zb.dims[1]);
            mf.z_passive.resize(mf.zb.dims[1]);

            double conv = area*sqr(dpi/180.0)/(4.0*dpi);

            for (uint_t z : range(mf.zb.dims[1])) {
                vec1u idm = where(mf.mx > z_mmin[z]);
                double dv_dz = bin_width(tvb(_,z))/bin_width(mf.zb(_,z));
                mf.z_active[z] = conv*integrate(mf.mb(_,idm), mf.mz_active(z,idm))*dv_dz;
                mf.z_passive[z] = conv*integrate(mf.mb(_,idm), mf.mz_passive(z,idm))*dv_dz;
            }
        }

        out.zb = zb;

        // Initialize sky positions
        // ------------------------

        // The shape of the generated survey can be changed here. The survey boundaries
        // are defined as a *convex* polygon whose vertices are given in RA and Dec
        // coordinates below.
        // The default is to generate a square of requested area around the reference position
        // given by ra0 and dec0.
        convex_hull<double> hull; {
            double dd = sqrt(area)/2.0;
            double dr = sqrt(area)/2.0/cos(dec0*dpi/180.0);

            vec1d hdec = {dec0-dd, dec0-dd, dec0+dd, dec0+dd};
            vec1d hra = {ra0-dr, ra0+dr, ra0+dr, ra0-dr};

            hull = build_convex_hull(hra, hdec);
        }

        // Generate redshifts
        // ------------------
        if (verbose) {
            note("generating redshifts...");
        }

        // Compute how many active and passive galaxies we will actually generate
        vec1u z_nactive(nz);
        vec1u z_npassive(nz);
        for (uint_t iz : range(nz)) {
            double dz = bin_width(zb(_,iz));

            double nact = mf.z_active[iz]*dz;
            z_nactive[iz] = floor(nact);
            nact -= z_nactive[iz];
            z_nactive[iz] += random_coin(seed, nact);

            double npas = mf.z_passive[iz]*dz;
            z_npassive[iz] = floor(npas);
            npas -= z_npassive[iz];
            z_npassive[iz] += random_coin(seed, npas);
        }

        uint_t nactive = total(z_nactive);
        uint_t npassive = total(z_npassive);
        uint_t ngal = nactive + npassive;

        if (verbose) {
            note("generated ", ngal, " galaxies");
        }

        auto random_pdf_bin = [&seed](const vec2d& b, const vec1d& x, const vec1d& p, uint_t n) {
            vec1d tx = {b(0,0)};
            append(tx, x);
            tx.push_back(b(1,x.size()-1));
            vec1d tp = {p[0]};
            append(tp, p);
            tp.push_back(p.back());
            return random_pdf(seed, tx, tp, n);
        };

        vec1f tz = random_pdf_bin(mf.zb, mf.zx, mf.z_active, nactive);
        append(out.z, tz);
        append(out.passive, replicate(false, nactive));

        tz = random_pdf_bin(mf.zb, mf.zx, mf.z_passive, npassive);
        append(out.z, tz);
        append(out.passive, replicate(true, npassive));

        out.id = indgen(out.z.size());
        out.d = lumdist(out.z, cosmo);

        // Generate masses
        // ---------------

        if (verbose) {
            note("generating masses...");
        }

        out.m.resize(ngal);

        for (uint_t iz : range(nz)) {
            vec1u idm = where(mf.mx > z_mmin[iz]);
            if (idm.empty()) continue;

            vec1b inz = in_bin_open(out.z, zb, iz);

            vec1u idl = where(inz && !out.passive);
            vec1f m = random_pdf_bin(mf.mb(_,idm), mf.mx[idm], mf.mz_active(iz,idm), idl.size());
            out.m[idl] = m;

            idl = where(inz && out.passive);
            m = random_pdf_bin(mf.mb(_,idm), mf.mx[idm], mf.mz_passive(iz,idm), idl.size());
            out.m[idl] = m;

            // Cleanup
            idl = where(inz and out.m < z_mmin[iz]);
            out.m[idl] = z_mmin[iz];
        }

        // Generate random positions
        // -------------------------

        if (!no_pos) {
            if (verbose) {
                note("generating sky positions...");
            }

            // Compute the hull area
            double fa = field_area_hull(hull);
            // Compute the maximum distance and field center
            double rmax = 0;
            double cra = 0, cdec = 0;
            for (uint_t i : range(hull))
            for (uint_t j : range(i+1, hull.size())) {
                double d = angdist(hull.x[i], hull.y[i], hull.x[j], hull.y[j])/3600.0;
                if (d > rmax) {
                    cra = 0.5*(hull.x[i] + hull.x[j]);
                    cdec = 0.5*(hull.y[i] + hull.y[j]);
                    rmax = d;
                }
            }

            rmax /= 2.0; // we want a radius, not a diameter

            if (fa < 0.2*dpi*sqr(clust_r0)) {
                note("requested area is too small, clustering is disabled");
                no_clust = true;
            } if (fa > 100.0) {
                error("the area of the field (", fa, " deg^2) is too large, this is not "
                    "yet supported");
                note("please choose an area that is lower than 100 deg^2");
                // reason: see below
                return 1;
            }

            // Note: To properly handle generating positions on the sphere, we will assume
            // that the field is relatively small (spanning less than 10 degree in any
            // direction) and generate the whole field as a planar surface around RA = 90
            // and Dec = 0. We rotate the generated coordinates to the requested
            // field center at the end.

            double ang = vernal_angle(90.0, cdec);

            // Rotate the convex hull to center at (90,0)
            hull.x += 90.0 - cra;
            angrot_vernal(hull.x, hull.y, -ang, hull.x, hull.y);

            // Rebuild the hull
            hull = build_convex_hull(hull.x, hull.y);

            // Compute bounding box of the rotated hull
            vec1d rra = {min(hull.x), max(hull.x)};
            vec1d rdec = {min(hull.y), max(hull.y)};

            if (no_clust) {
                // Uniformly random positions

                auto status = randpos_uniform_box(seed, ngal,
                    rra, rdec, out.ra, out.dec,
                    [&](double xra, double xdec) {
                        return in_convex_hull(xra, xdec, hull);
                    }
                );

                if (!status.success) {
                    error("generating positions");
                    error("failed: ", status.failure);
                    return 1;
                }

                out.clustered = replicate(false, ngal);

            } else {
                // Clustered positions

                // The initial radius
                double rstart = clust_r0;
                double rend = clust_r1;

                out.ra.resize(ngal);
                out.dec.resize(ngal);
                out.clustered.resize(ngal);

                // Clustering redshift slices (fixed to preserve clustering amplitude)
                vec2f czb = make_zbins(0.0, std::max(10.0, zmax), 0.1, 0.1, 0.5);
                uint_t ncz = czb.dims[1];

                auto pg = progress_start(ncz);
                for (uint_t iz : range(ncz)) {
                    // Generate clustered positions within this redshift window.
                    vec1u idz = where(in_bin_open(out.z, czb, iz));
                    uint_t z_ngal = idz.size();

                    if (z_ngal != 0) {
                        vec1d ra, dec;
                        ra.reserve(z_ngal);
                        dec.reserve(z_ngal);

                        // Select the number of simulated positions
                        // We impose a minimum number of 1000 to prevent issues with low
                        // densities and also generate a fraction more objects than we
                        // need.
                        uint_t nsim_min = 1000;
                        double fudge = 1.5;
                        uint_t nsim = fudge*(z_ngal > nsim_min ? z_ngal : nsim_min);

                        // Compute the expected source density
                        double rho = nsim/fa;

                        // Initialize the algorithm parameters
                        randpos_power_options opt;
                        opt.eta = clust_eta;
                        opt.lambda = clust_lambda;

                        // Compute how many deeper levels we need to go with respect to
                        // rstart to reach the local density
                        opt.levels = ceil(max(0, log(rho*dpi*sqr(rstart)))/log(opt.eta));

                        // Compute how many homogeneous starting positions we need
                        // above rstart
                        uint_t nstart = ceil(nsim/pow(opt.eta, opt.levels));

                        // and how many homogenous positions we need at rend
                        if (rend != 0) {
                            uint_t maxl = round(log(rstart/rend)/log(opt.lambda));
                            if (opt.levels > maxl) {
                                opt.eta_end = pow(opt.eta, opt.levels - maxl + 1);
                                opt.levels = maxl;
                            }
                        }

                        // Generate starting positions
                        vec1d sra, sdec;
                        if (nstart <= 1) {
                            // Just use the center of the field
                            sra = {90.0}; sdec = {0.0};
                        } else {
                            // Generate honogenous positions
                            auto status = randpos_uniform_box(seed, nstart,
                                rra, rdec, sra, sdec,
                                [&](double tra, double tdec) {
                                    return in_convex_hull(tra, tdec, hull);
                                }
                            );

                            if (!status.success) {
                                if (verbose) print("");
                                error("generating starting positions for redshift bin ",
                                    iz);
                                error("failed: ", status.failure);
                                return 1;
                            }
                        }

                        // Now generate the positions using the original Soneira & Pebbles
                        // algorithm
                        for (uint_t i : range(sra)) {
                            vec1d tra, tdec;
                            auto status = randpos_power_circle(seed, sra[i], sdec[i],
                                rstart, tra, tdec,
                                [&](double ttra, double ttdec, double threshold) {
                                    double dist = convex_hull_distance(
                                        ttra, ttdec, hull
                                    );
                                    return dist > -threshold;
                                },
                                opt
                            );

                            if (!status.success) {
                                if (verbose) print("");
                                error("generating clustered positions for redshift bin ",
                                    iz);
                                error("failed: ", status.failure);
                                return 1;
                            }

                            append(ra, tra);
                            append(dec, tdec);
                        }

                        // Now we have too many positions, just pick the amount we need
                        vec1u idg = shuffle(seed, indgen(ra.size()))[_-(z_ngal-1)];
                        ra = ra[idg]; dec = dec[idg];

                        vec1b cls = replicate(true, z_ngal);

                        // Each galaxy has a probability of being clustered
                        vec1d prnd = replicate(clust_fclust_lom, z_ngal);
                        // We treat massive and low mass galaxies differently,
                        // assuming more clustering for the most massive objects.
                        prnd[where(out.m[idz] > clust_fclust_mlim)] = clust_fclust_him;
                        // and no clustering at all below a certain mass
                        prnd[where(out.m[idz] < clust_urnd_mlim)] = 0;
                        // then reduce clustering at low redshifts (Béthermin+15)
                        {
                            vec1u idlz = where(out.z[idz] < 1.0);
                            prnd[idlz] *= max(0.0, out.z[idz[idlz]] - 0.5)/0.5;
                        }

                        prnd = 1-prnd;

                        vec1u idr = where(random_coin(seed, prnd));
                        uint_t nrnd = idr.size();

                        vec1d tra, tdec;
                        auto status = randpos_uniform_box(seed, nrnd, rra, rdec,
                            tra, tdec,
                            [&](double xra, double xdec) {
                                return in_convex_hull(xra, xdec, hull);
                            }
                        );

                        if (!status.success) {
                            if (verbose) print("");
                            error("generating non-clustered positions for redshift bin ",
                                iz);
                            error("failed: ", status.failure);
                            return 1;
                        }

                        ra[idr] = tra;
                        dec[idr] = tdec;
                        cls[idr] = false;

                        // Store the new positions
                        out.ra[idz] = ra;
                        out.dec[idz] = dec;
                        out.clustered[idz] = cls;
                    }

                    if (verbose) progress(pg);
                }
            }

            // Rotate back the generated coordinates
            angrot_vernal(out.ra, out.dec, ang, out.ra, out.dec);
            out.ra += cra - 90.0;
            normalize_coordinates(out.ra, out.dec);
        } else {
            out.ra = replicate(dnan, ngal);
            out.dec = replicate(dnan, ngal);
            out.clustered = replicate(false, ngal);
        }
    }

    // Get counts and IDs of the active and passive populations
    const vec1u ida = where(!out.passive);
    const vec1u idp = where(out.passive);
    const uint_t nactive = ida.size();
    const uint_t npassive = idp.size();
    const uint_t ngal = out.z.size();

    // Bulge to total mass ratio
    // Calibrated from Lang et al. (2014), assuming no redshift evolution
    out.bt.resize(ngal);

    out.bt[ida] = e10(0.27*(out.m[ida] - 10.0) - 0.7);
    out.bt[idp] = e10(0.1*(out.m[idp] - 10.0) - 0.3);

    // Add random scatter
    if (!no_random) {
        out.bt *= e10(0.2*randomn(seed, ngal));
    }

    out.bt = clamp(out.bt, 0.0, 1.0);

    out.m_bulge = out.m + log10(out.bt);
    vec1u idnf = where(!is_finite(out.m_bulge));
    out.m_bulge[idnf] = 0;

    out.m_disk = out.m + log10(1.0-out.bt);
    idnf = where(!is_finite(out.m_disk));
    out.m_disk[idnf] = 0;

    // Generate morphology
    // -------------------
    // Note: Most of this is calibrated using the single Sersic fits of
    // Arjen van der Wel in the CANDELS fields.

    if (verbose) {
        note("generating morphology...");
    }

    // Declare functions to assign sizes
    auto size_disk = vectorize_lambda([](double z, double m) {
        if (z < 1.7) {
            return e10(0.41 - 0.22*log10(1.0+z) + 0.2*(min(m, 10.6) - 9.35));
        } else {
            return e10(0.62 - 0.7*log10(1.0+z) + 0.2*(min(m, 10.6) - 9.35));
        }
    });

    auto size_bulge = vectorize_lambda([&size_disk](double z, double m, double bt) {
        if (z < 0.5) {
            return e10(0.78 - 0.6*log10(1.0+z) + 0.56*(m - 11.25));
        } else {
            return e10(0.90 - 1.3*log10(1.0+z) + 0.56*(m - 11.25));
        }
    });

    // Bulge and disk have the same position angle, which is completely random
    out.disk_angle = (randomu(seed, ngal) - 0.5)*180.0;
    out.bulge_angle = out.disk_angle;

    // Bulges: calibration from n>2.5 galaxies and M* > 10.5
    vec1f bulge_ratio_x =
        {0.0, 0.1, 0.15,  0.25,  0.35,  0.45,  0.55,   0.65,   0.75,   0.85,   0.95,  1.0};
    vec1f bulge_ratio_p =
        {0.0, 0.0,  165.0, 428.0, 773.0, 914.0, 1069.0, 1191.0, 1154.0, 1067.0, 639.0, 450.0};
    out.bulge_ratio = random_pdf(seed, bulge_ratio_x, bulge_ratio_p, ngal);

    out.bulge_radius = size_bulge(out.z, out.m, out.bt);

    if (!no_random) {
        out.bulge_radius *= e10(0.2*randomn(seed, ngal));
    }

    // Disks: calibration from n<1.5 galaxies and M* > 9.0
    vec1f disk_ratio_x =
        {0.0, 0.1, 0.15,  0.25,  0.35,   0.45,   0.55,   0.65,  0.75,  0.85,  0.95,  1.0};
    vec1f disk_ratio_p =
        {0.0, 0.0,  313.0, 900.0, 1143.0, 1127.0, 1059.0, 898.0, 775.0, 548.0, 318.0, 200.0};
    out.disk_ratio = random_pdf(seed, disk_ratio_x, disk_ratio_p, ngal);

    out.disk_radius = size_disk(out.z, out.m);

    if (!no_random) {
        out.disk_radius *= e10(0.17*randomn(seed, ngal));
    }

    // Convert proper sizes to angular sizes
    vec1d psize = propsize(out.z, cosmo);
    out.disk_radius /= psize;
    out.bulge_radius /= psize;

    // Generate SFR
    // ------------
    if (verbose) {
        note("generating SFR...");
    }

    out.sfr.resize(ngal);
    out.rsb.resize(ngal);

    vec1f sfrms;

    {
        // Active galaxies
        vec1f az = log10(1.0 + out.z);
        sfrms = out.m - 9.67 + 1.82*az - 0.38*sqr(max(0.0, out.m - 9.59 - 2.22*az));

        out.sfr[ida] = sfrms[ida];

        // Passive galaxies
        out.sfr[idp] = min(sfrms[idp], 0.5*(out.m[idp] - 11) + az[idp] - 0.6);

        if (!no_random) {
            // Add dispersion
            out.sfr[ida] += ms_disp*randomn(seed, nactive);
            out.sfr[idp] += 0.4*randomn(seed, npassive);

            // Add starbursts
            vec1u idsb = where(random_coin(seed, 0.033, nactive));
            out.sfr[ida[idsb]] += 0.72;
        }

        // Compute RSB and get final SFR in linear units
        out.rsb = out.sfr - sfrms;
        out.sfr = e10(out.sfr);
    }

    // Generate SFR_IR and SFR_UV
    // --------------------------
    {
        // IRX = SFR_IR/SFR_UV
        // Calibrated from Herschel stacks
        double am = 10.5;
        double airx = 1.2;
        double ss = 0.45, so = 0.35;
        vec1f slope = ss*min(out.z, 3.0) + so;
        out.irx = e10(slope*(out.m - am) + airx);

        if (!no_random) {
            out.irx *= e10(0.4*randomn(seed, ngal));
        }

        if (no_passive_lir) {
            out.irx[idp] = 0.0;
        }

        out.sfrir = out.sfr/(1.0 + 1.0/out.irx);
        out.sfruv = out.sfr/(1.0 + out.irx);
    }

    // Generate UVJ colors
    // -------------------
    out.rfvj_disk.resize(ngal);
    out.rfvj_bulge.resize(ngal);
    out.rfuv_disk.resize(ngal);
    out.rfuv_bulge.resize(ngal);

    // Disks are always blue for active galaxies and red for passive galaxies
    {
        vec1f tuv, tvj;
        gen_blue_uvj(out.m[ida], out.z[ida], tuv, tvj, no_random);
        out.rfuv_disk[ida] = tuv;
        out.rfvj_disk[ida] = tvj;

        gen_red_uvj(out.m[idp], out.z[idp], tuv, tvj, no_random);
        out.rfuv_disk[idp] = tuv;
        out.rfvj_disk[idp] = tvj;
    }

    // Bulges of bulge-dominated galaxies or passive galaxies are always red
    // Other bulges have a 50% chance of being red
    {
        vec1b red_bulge = out.bt > 0.6 || out.passive;
        if (!no_random) {
            red_bulge = red_bulge || random_coin(seed, 0.5, ngal);
        }

        vec1f tuv, tvj;
        vec1u idb = where(red_bulge);
        gen_red_uvj(out.m[idb], out.z[idb], tuv, tvj, no_random);
        out.rfuv_bulge[idb] = tuv;
        out.rfvj_bulge[idb] = tvj;

        idb = where(!red_bulge);
        gen_blue_uvj(out.m[idb], out.z[idb], tuv, tvj, no_random);
        out.rfuv_bulge[idb] = tuv;
        out.rfvj_bulge[idb] = tvj;
    }

    // Generate redshift-dependent mass-to-light ratio correction
    out.m2lcor = get_m2l_cor(out.z);

    // Assign optical SED
    // ------------------

if (!no_stellar) {
    if (verbose) {
        note("assigning optical SEDs...");
    }

    if (!get_sed_uvj(out.rfuv_bulge, out.rfvj_bulge, out.opt_sed_bulge, out.av_bulge, verbose)) {
        return 1;
    }

    if (!get_sed_uvj(out.rfuv_disk,  out.rfvj_disk,  out.opt_sed_disk, out.av_disk, verbose)) {
        return 1;
    }
} else {
    out.opt_sed_bulge = replicate(npos, ngal);
    out.opt_sed_disk = replicate(npos, ngal);
    out.av_bulge = replicate(fnan, ngal);
    out.av_disk = replicate(fnan, ngal);
}

    // Generate IR properties
    // ---------------------

    if (verbose) {
        note("generate IR properties...");
    }

    // Tdust and IR8 as observed in stacks and detections of Herschel galaxies
    vec1f olir = out.lir;
    vec1f oir8 = out.ir8;
    vec1f otdust = out.tdust;

    vec1f tdust_ms = 32.9 + 4.60*(out.z - 2.0) - 0.77;
    out.tdust = tdust_ms;

    if (magdis_tdust) {
        vec1u idz = where(out.z > 1);
        out.tdust[idz] = 2.0*(clamp(out.z[idz], 0.0, 2.0)-1.0) + 27.0;
    }

    // Starbursts are warmer
    out.tdust += 10.1*out.rsb;
    // Massive galaxies are colder (= downfall of SFE)
    out.tdust -= 1.5*max(0.0, 2.0 - out.z)*clamp(out.m - 10.7, 0.0, 1.0);

    out.ir8 = (4.08 + 3.29*clamp(out.z-1.0, 0.0, 1.0))*0.81;
    // Starburst have larger IR8
    out.ir8 *= e10(0.66*max(0.0, out.rsb));
    // Low-mass galaxies have larger IR8
    out.ir8 *= e10(-1.0*clamp(out.m - 10.0, -1, 0.0));

    if (!no_random) {
        // Add some random scatter
        out.tdust += 0.12*tdust_ms*randomn(seed, ngal);
        out.ir8 *= e10(0.18*randomn(seed, ngal));
    }

    out.lir = out.sfrir/1.72e-10;

    // Function to update a source's properties if its LIR is changed
    auto update_properties = [&out,&sfrms](uint_t id) {
        out.sfrir[id] = out.lir[id]*1.72e-10;
        out.sfr[id] = out.sfruv[id] + out.sfrir[id];
        out.irx[id] = out.sfrir[id]/out.sfruv[id];
        double orsb = out.rsb[id];
        out.rsb[id] = log10(out.sfr[id]) - sfrms[id];
        out.tdust[id] += 6.42*(out.rsb[id] - orsb);
        out.ir8[id] *= e10(0.40*(max(0.0, out.rsb[id]) - max(0.0, orsb)));
    };

    // Keep the values provided by the user in the input catalog (if any)
    if (!olir.empty()) {
        vec1u idg = where(is_finite(olir));

        // For each provided LIR, find the source at a similar redshift and M* with
        // the closest LIR, and give it the old LIR of the provided source.
        // This preserves as much as possible the overall LIR distribution
        for (uint_t i : idg) {
            vec1u idl = where(abs(out.m[i] - out.m) < 0.3 &&
                abs(out.z[i] - out.z)/(1.0 + out.z) < 0.1 &&
                !is_finite(olir) && out.passive[i] == out.passive);

            if (idl.empty()) {
                // No one? Try loosening the constraints a little
                idl = where(abs(out.m[i] - out.m) < 0.4 &&
                    abs(out.z[i] - out.z)/(1.0 + out.z) < 0.2 &&
                    !is_finite(olir) && out.passive[i] == out.passive);
            }

            if (idl.empty()) {
                // Found one, give it the old LIR and adjust all its other
                // properties to reflect this change
                uint_t mid = min_id(abs(log10(olir[i]/out.lir[idl])));
                out.lir[mid] = out.lir[i];
                update_properties(mid);
            }

            out.lir[i] = olir[i];
            update_properties(i);
        }
    }

    if (!oir8.empty()) {
        vec1u idg = where(is_finite(oir8));
        out.ir8[idg] = oir8[idg];
    }

    if (!otdust.empty()) {
        vec1u idg = where(is_finite(otdust));
        out.tdust[idg] = otdust[idg];
    }

    // Make sure the values are valid
    out.ir8 = clamp(out.ir8, 0.48, 27.5); // range allowed by IR library
    out.fpah = clamp(1.0/(1.0 - (331 - 691*out.ir8)/(193 - 6.98*out.ir8)), 0.0, 1.0);
    out.lir[where(!is_finite(out.lir))] = 0.0;


    // Assign IR SED
    // -------------

    if (verbose) {
        note("assigning IR SED...");
    }

    // Floating point SED index, we will arrange it into a integer later
    vec1d fir_sed;

    if (nirsed == 1) {
        fir_sed = replicate(0, nactive);
    } else if (file::get_basename(ir_lib_file) == "ir_lib_ce01.fits") {
        // The Chary & Elbaz 2001 library, redshift evolution calibrated from stacks
        vec1u tidsed = {26, 26, 40, 54, 53, 52, 52};
        vec1f tz = {0.57, 1.0, 1.5, 2.1, 2.9, 4.0, 6.0};
        fir_sed = interpolate(tidsed, tz, out.z);

        // Temperature offset as function of RSB (not calibrated, but see Magnelli+13)
        fir_sed += 15*clamp(out.rsb/ms_disp, -2.0, 2.0);
    } else if (file::get_basename(ir_lib_file) == "ir_lib_m12.fits") {
        // The Magdis et al. 2012 library, using their reported redshift evolution
        vec1f tz = {0.0125, 0.1625, 0.4625, 0.8125, 1.15, 1.525, 2.0, 2.635};
        fir_sed = interpolate(indgen<float>(nirsed), tz, out.z);

        // Temperature offset as function of RSB (not calibrated, but see Magnelli+13)
        fir_sed += clamp(out.rsb/ms_disp, -2.0, 2.0);
    } else if (file::get_basename(ir_lib_file) == "ir_lib_cs17.fits") {
        // My own library, using calibration from stacks and detections
        fir_sed = interpolate(indgen<float>(nirsed), ir_lib_cs17.tdust, out.tdust);
    } else {
        error("no calibration code available for the IR library '", ir_lib_file, "'");
        return 1;
    }

    out.ir_sed = clamp(round(fir_sed), 0, nirsed-1);

    if (file::get_basename(ir_lib_file) == "ir_lib_cs17.fits") {
        // My library

        // Get exact fPAH
        out.fpah = clamp(1.0/(1.0 - (ir_lib_cs17.lir_pah[out.ir_sed] - ir_lib_cs17.l8_pah[out.ir_sed]*out.ir8)/
            (ir_lib_cs17.lir_dust[out.ir_sed] - ir_lib_cs17.l8_dust[out.ir_sed]*out.ir8)), 0.0, 1.0);

        // SED is in unit Mdust, normalize it to our dust mass
        out.mdust = out.lir/(ir_lib_cs17.lir_dust[out.ir_sed]*(1.0 - out.fpah)
            + ir_lib_cs17.lir_pah[out.ir_sed]*out.fpah);
    } else {
        // Placeholders
        out.fpah = replicate(0.04, ngal);
        out.mdust = out.lir/1e3;
    }

    // Pre-compute IGM absorption
    // ---------------------------

    if (verbose) {
        note("generating IGM absorption");
    }

    out.igm_abs.resize(ngal, 3);
    for (uint_t i : range(ngal)) {
        double z = out.z[i];

        if (igm == "madau95") {
            // Madau+95 IGM (implementation from FAST)
            // http://adsabs.harvard.edu/abs/1995ApJ...441...18M

            {
                double l0 = 1050.0*(1.0 + z);
                double l1 = 1170.0*(1.0 + z);
                uint_t nstep = 100;
                vec1d tl = rgen(l0, l1, nstep);
                vec1d ptau = exp(-3.6e-3*pow(tl/1216.0, 3.46));
                out.igm_abs(i,2) = mean(ptau);
            }

            {
                double l0 = 920.0*(1.0 + z);
                double l1 = 1015.0*(1.0 + z);
                uint_t nstep = 100;
                vec1d tl = rgen(l0, l1, nstep);
                vec1d ptau = exp(-1.7e-3*pow(tl/1026.0, 3.46) - 1.2e-3*pow(tl/972.5, 3.46) -
                    9.3e-4*pow(tl/950.0, 3.46));
                out.igm_abs(i,1) = mean(ptau);
            }

            // No transmission below 912A
            out.igm_abs(i,0) = 0.0;
        } else if (igm == "inoue14") {
            // Inoue+14 IGM prescriptioms.
            // Probability distribution functions calibrated on Monte Carlo realizations
            // kindly provided by Akio Inoue (see Inoue+08 for method)
            // http://adsabs.harvard.edu/abs/2014MNRAS.442.1805I
            // http://adsabs.harvard.edu/abs/2008MNRAS.387.1681I

            static const vec2f igm_z1 = {
                {0.00000,-0.119247,-0.836013,0.188466,5.58061,1.80948,2.23333,
                 2.58904,3.05070,3.44932,3.64749,4.33884,4.80831},
                {4.71243,4.51667,5.13483,5.51325,6.02321,6.11072,3.96055,
                 6.16846,6.20577,4.11253,6.25897,6.30988,6.30475},
                {4.29922,4.43374,5.81153,6.07670,6.19915,4.01990,6.22758,
                 6.23932,4.18000,6.29191,6.30732,6.34660,4.29822}
            };
            static const vec2f igm_dz1 = {
                {0.00000,0.103763,1.53267,0.831935,0.101257,1.30343,1.23153,
                 1.18386,1.18305,1.21584,1.22993,1.25208,1.15790},
                {1.65837,1.78514,1.81298,1.76652,1.33753,1.19458,2.24961,
                 1.11023,1.07492,2.16499,1.03259,0.998597,1.00750},
                {1.83693,1.92457,1.44038,1.07901,0.980853,2.13450,1.00173,
                 0.996720,2.25818,0.979706,0.979425,0.962433,2.04183}
            };
            static const vec2f igm_z2 = {
                {0.00000,0.00000,8.03488,4.25279,0.989755,6.97760,0.00000,
                 0.00000,7.46371,5.79003,5.31787,4.35838,4.12159},
                {3.13640,4.21302,4.01985,3.91207,3.85061,3.91916,6.14354,
                 3.99942,4.05655,6.24326,4.14120,4.22585,4.28847},
                {5.85101,5.60883,3.88240,3.87634,3.98672,-9.52333,4.10508,
                 4.13700,6.26376,4.22116,4.24243,4.30912,9.86437}
            };
            static const vec2f igm_dz2 = {
                {0.00000,0.00000,0.209441,0.0120493,1.27333,0.0810017,0.00000,
                 0.00000,5.04822,3.71614,3.30993,2.40054,2.10543},
                {5.31391,5.23251,3.21628,2.77969,2.37366,2.28242,1.14531,
                 2.22511,2.19319,1.04437,2.15259,2.11453,2.11378},
                {9.48358,6.52196,2.79285,2.47621,2.31278,-1.12091,2.28194,
                 2.27401,0.988508,2.24063,2.23221,2.20610,0.233647}
            };

            static const vec1f quantiles = {
                0.0, 0.001, 0.023, 0.05, 0.16, 0.35, 0.50, 0.65, 0.84, 0.95, 0.977, 0.999, 1.0
            };

            for (uint_t l : range(igm_z1.dims[0])) {
                // Evaluate CDF at z
                vec1d trans_quant(quantiles.size());
                for (uint_t q : range(quantiles)) {
                    if (igm_dz1.safe(l,q) == 0.0) continue;
                    trans_quant.safe[q] = 0.5*(erf((igm_z1.safe(l,q) - z)/igm_dz1.safe(l,q)) + 1.0);

                    if (igm_dz2.safe(l,q) == 0.0) continue;
                    trans_quant.safe[q] *= 0.5*(erf((igm_z2.safe(l,q) - z)/igm_dz2.safe(l,q)) + 1.0);
                }

                if (no_random) {
                    // Use median transmission
                    out.igm_abs(i,l) = trans_quant[quantiles.size()/2];
                } else{
                    // Randomize IGM
                    out.igm_abs(i,l) = interpolate(trans_quant, quantiles, randomu(seed));
                }
            }
        }
    }

    // Generate emission lines
    // -----------------------

    if (!no_nebular) {
        if (verbose) {
            note("generating emission lines");
        }

        out.lines       = {"c2_157", "n2_205", "c1_609", "co10", "co21", "co32", "co43", "co54", "co65", "co76"};
        out.line_lambda = {157.71,   205.18,   609.14,   2600.9, 1300.4, 867.0,  650.27, 520.24, 433.57, 371.66};

        append(out.lines,       vec1s{"halpha", "hbeta", "hgamma", "hdelta", "n2_6583", "n2_6548", "o3_5007", "o3_4959", "o2_3727", "lyalpha"});
        append(out.line_lambda, vec1f{0.65628,  0.48613, 0.43405,  0.41017,  0.65835,   0.65480,   0.50068,   0.49589,   0.37274,   0.12157});

        out.line_lum_bulge.resize(ngal, out.lines.size());
        out.line_lum_disk.resize(ngal, out.lines.size());

        // TODO: compare to Francesco's paper

        // Pre-compute some more properties
        // Metallicity
        vec1d metal = replicate(9.07, ngal);
        vec1d mu32 = (out.m-0.2) - 0.32*(log10(out.sfr) - 0.2);
        vec1u idl = where(mu32 < 10.36);
        metal[idl] = 8.9 + 0.47*(mu32[idl] - 10);
        metal -= 8.69; // in solar metalicity
        // Broken FMR (Bethermin+16), tweaked to reproduce [OIII]/Hbeta @ z=2 (Dickey+16)
        idl = where(out.z > 1);
        metal[idl] -= 0.2*clamp(out.z[idl] - 1, 0.0, 1.0);
        metal = clamp(metal, -0.6, 1.0);
        // Gas-to-dust ratio
        vec1d gdr = 2.23 - metal;
        if (!no_random) gdr += 0.04*randomn(seed, ngal);
        // Gas masses
        vec1d mgas = out.mdust*e10(gdr); // HI + H2
        vec1d mh2 = mgas*0.3;
        if (!no_random) mh2 *= e10(0.2*randomn(seed, ngal));
        // Attenuation of lines (Pannella+15 for redshift dependence)
        out.avlines_disk = 0.5*(log10(out.sfr/out.sfruv)*0.95 + out.av_disk);
        out.avlines_disk *= interpolate(vec1d{1.7, 1.3, 1.0, 1.0}, vec1d{0,1,2,100}, out.z);
        if (!no_random) out.avlines_disk += 0.1*randomn(seed, ngal);
        out.avlines_bulge = out.av_bulge;
        if (!no_random) out.avlines_bulge += 0.1*randomn(seed, ngal);
        out.avlines_disk  = clamp(out.avlines_disk,  0.0, 6.0);
        out.avlines_bulge = clamp(out.avlines_bulge, 0.0, 6.0);
        // Ly-alpha escape fraction (Hayes+10)
        vec1d fescape_disk = 0.445*e10(-0.4*(out.av_disk/4.05)*17.8);
        // correction to avoid over-producing Ly-alpha at z~1-2
        fescape_disk *= e10(interpolate(vec1d{-1.2, -1.2, -0.8, 0, 0}, vec1d{0.5, 1.2, 2.2, 3.0, 4.0}, out.z));
        if (!no_random) fescape_disk *= e10(0.4*randomn(seed, ngal));
        fescape_disk = clamp(fescape_disk, 0.0, 1.0);
        // Velocity dispersion (Stott+16)
        out.vdisp = e10(0.12*(out.m - 10) + 1.78);
        if (!no_random) out.vdisp *= e10(0.3*randomn(seed, ngal));

        uint_t l;

        // Far IR lines

        l = where_first(out.lines == "c2_157");
        out.line_lum_bulge(_,l) = 0.0;
        out.line_lum_disk(_,l) = e10(1.017*(log10(out.sfr)-0.2) + 7.13); // deLooze+11
        if (!no_random) out.line_lum_disk(_,l) *= e10(0.27*randomn(seed, ngal));

        l = where_first(out.lines == "n2_205");
        out.line_lum_bulge(_,l) = 0.0;
        out.line_lum_disk(_,l) = e10(1.053*(log10(out.sfr)-0.2) + 5.31); // Zhao+13
        if (!no_random) out.line_lum_disk(_,l) *= e10(0.26*randomn(seed, ngal));

        l = where_first(out.lines == "c1_609");
        out.line_lum_bulge(_,l) = 0.0;
        out.line_lum_disk(_,l) = 0.00246*mh2; // Bothwell+17

        l = where_first(out.lines == "co10");
        out.line_lum_bulge(_,l) = 0.0;
        out.line_lum_disk(_,l) = 4.93e-5*mh2; // alphaCO = 4.6

        // CO ladder
        vec1d sled = {2.3, 3.7, 4.5, 7.5, 8.5, 8.3}; // Bothwell+13 SMGs
        for (uint_t dl = 0; dl < 6; ++dl) {
            out.line_lum_bulge(_,l+1+dl) = 0.0;
            out.line_lum_disk(_,l+1+dl) = out.line_lum_disk(_,l)*sled[dl];
        }

        // UV-optical lines

        l = where_first(out.lines == "halpha");
        out.line_lum_bulge(_,l) = 0.0;
        out.line_lum_disk(_,l) = e10(log10(out.sfr) + 7.519); // Kennicutt+98
        if (!no_random) out.line_lum_disk(_,l) *= e10(0.15*randomn(seed, ngal));

        // Case B recombination
        vec1d recomb = {0.35, 0.163, 0.0862};
        for (uint_t dl = 0; dl < 3; ++dl) {
            out.line_lum_bulge(_,l+1+dl) = 0.0;
            out.line_lum_disk(_,l+1+dl) = out.line_lum_disk(_,l)*recomb[dl];
        }

        // [NII]/Halpha, [OIII]/Hbeta calibrated vs Halpha, Hbeta and FMR metallicity in SDSS
        // plus [NII]/Halpha offset from MOSDEF

        auto poly3 = [](const vec1d& x, const vec1d& coef) {
            return coef[0] + x*coef[1] + sqr(x)*coef[2] + pow(x,3)*coef[3];
        };

        vec1d lo3mn2 = poly3(metal, {0.59518604, -2.7283054, -0.23229248, 2.0151786});
        if (!no_random) {
            lo3mn2 += poly3(metal, {0.32970659, 0.043645415, -1.2231136, -2.4836792})*randomn(seed, ngal);
        }
        vec1d lo3pn2 = poly3(lo3mn2, {-0.69337981, 0.66460936, -0.37040100, 0.047935498});
        if (!no_random) {
            lo3pn2 += poly3(lo3mn2, {0.086819227, -0.078283561, 0.093276637, -0.025069864})*randomn(seed,ngal);
        }

        vec1d lha = out.line_lum_disk(_,where_first(out.lines == "halpha"));
        vec1d n2offset = interpolate(vec1d{0.0, 0.0, 0.3, 0.3}, vec1d{0.0, 1.0, 2.3, 2.4}, out.z)
            *interpolate(vec1d{1.0, 1.0, 0.0, 0.0}, vec1d{9, 10.0, 10.5, 11.0}, out.m); // Shapley+15
        vec1d ln2 = lha*e10(0.5*(lo3pn2 - lo3mn2) + n2offset);

        l = where_first(out.lines == "n2_6583");
        out.line_lum_bulge(_,l) = 0.0;
        out.line_lum_disk(_,l) = 0.75*ln2;
        l = where_first(out.lines == "n2_6548");
        out.line_lum_bulge(_,l) = 0.0;
        out.line_lum_disk(_,l) = 0.25*ln2;

        vec1d lhb = out.line_lum_disk(_,where_first(out.lines == "hbeta"));
        vec1d lo3 = lhb*e10(0.5*(lo3pn2 + lo3mn2));

        l = where_first(out.lines == "o3_5007");
        out.line_lum_bulge(_,l) = 0.0;
        out.line_lum_disk(_,l) = 0.7*lo3;
        l = where_first(out.lines == "o3_4959");
        out.line_lum_bulge(_,l) = 0.0;
        out.line_lum_disk(_,l) = 0.3*lo3;

        // [OII] calibrated on Hbeta and [OIII]/Hbeta in SDSS

        l = where_first(out.lines == "o2_3727");
        vec1d o3hb = log10(lo3/lhb);
        out.line_lum_bulge(_,l) = 0.0;
        out.line_lum_disk(_,l) = lhb*e10(poly3(o3hb, {0.51869283, 0.29843257, -0.49610294, -0.11728400}));
        if (!no_random) {
            vec1d scatter = poly3(o3hb, {0.061403650, -0.0081335335, 0.078420795, -0.032656826});
            out.line_lum_disk(_,l) *= e10(scatter*randomn(seed, ngal));
        }

        // Ly-alpha based on Sobral+18

        l = where_first(out.lines == "lyalpha");
        out.line_lum_bulge(_,l) = 0.0;
        out.line_lum_disk(_,l) = 8.7*lha*fescape_disk;

        // Apply reddening and IGM

        auto calzetti2000 = vectorize_lambda([&](double l) {
            // http://adsabs.harvard.edu/abs/2000ApJ...533..682C

            const double iRv = 1.0/3.33;

            if (l <= 0.63) {
                l = (2.659*iRv)*(-2.156 + 1.509/l - 0.198*pow(l, -2) + 0.011*pow(l, -3)) + 1.0;
            } else {
                l = (2.659*iRv)*(-1.857 + 1.040/l) + 1.0;
            }

            if (l < 0) l = 0;

            return l;
        });

        vec1d dust_law = calzetti2000(out.line_lambda);
        dust_law[where_first(out.lines == "lyalpha")] = 0; // extinction already folded in

        for (uint_t i : range(ngal)) {
            vec1d igm_lines(out.lines.size());
            for (uint_t l : range(out.lines)) {
                if (out.line_lambda[l] < 0.0600) {
                    igm_lines[l] = 0.0;
                } else if (out.line_lambda[l] < 0.0912) {
                    igm_lines[l] = out.igm_abs(i,0);
                } else if (out.line_lambda[l] < 0.1026) {
                    igm_lines[l] = out.igm_abs(i,1);
                } else if (out.line_lambda[l] < 0.1216) {
                    igm_lines[l] = out.igm_abs(i,2);
                } else {
                    igm_lines[l] = 1.0;
                }
            }

            out.line_lum_disk(i,_) *= e10(-0.4*out.avlines_disk[i]*dust_law)*igm_lines;
            out.line_lum_bulge(i,_) *= e10(-0.4*out.avlines_bulge[i]*dust_law)*igm_lines;
        }

        // Compute total luminosities
        out.line_lum = out.line_lum_bulge + out.line_lum_disk;
    }

if (!no_flux) {

    // Compute fluxes
    // --------------

    out.flux_disk.resize(ngal, bands.size());
    out.flux_bulge.resize(ngal, bands.size());

    out.rfmag_disk.resize(ngal, rfbands.size());
    out.rfmag_bulge.resize(ngal, rfbands.size());

    if (verbose) {
        note("computing fluxes", (rffilters.empty() ? "" : " and absolute magnitudes"), "...");
    }

    auto get_ir_sed = [&](uint_t i, vec1f& irlam, vec1f& irsed) {
        uint_t s = out.ir_sed[i];

        if (file::get_basename(ir_lib_file) == "ir_lib_cs17.fits") {
            // My library

            // Build combined SED
            irlam = ir_lib_cs17.lam(s,_);
            irsed = ir_lib_cs17.dust(s,_)*(1.0 - out.fpah[i])
                + ir_lib_cs17.pah(s,_)*out.fpah[i];

            // SED is in unit Mdust, normalize it to our dust mass
            irsed *= out.mdust[i];
        } else {
            // Generic library, SEDs in units of LIR
            irlam = ir_lib.lam(s,_);
            irsed = out.lir[i]*ir_lib.sed(s,_);
        }
    };

    thread::worker sed_saver;
    uint_t nsed = 0;
    std::atomic<uint_t> nwritten(0);
    std::ofstream seds_data;
    vec1u save_sed_bulge_start, save_sed_bulge_nbyte;
    vec1u save_sed_disk_start,  save_sed_disk_nbyte;

    if (save_sed) {
        seds_data.open(seds_file, std::ios::binary);
        save_sed_bulge_start.resize(ngal);
        save_sed_bulge_nbyte.resize(ngal);
        save_sed_disk_start.resize(ngal);
        save_sed_disk_nbyte.resize(ngal);
    }

    enum class component {
        bulge, disk
    };

    auto add_emission_line = [&](const vec1f& lam, vec1f& sed, double l0, double lum, double vdisp) {
        const uint_t n = lam.size();
        double lw = l0*(vdisp/2.998e5);
        auto b = bounds(lam, l0 - 10*lw, l0 + 10*lw);
        if (b[1] == 0 || b[0] == n-1) return;
        if (b[0] == npos) b[0] = 0;
        if (b[1] == npos) b[1] = n-1;

        vec1f llow(b[1]-b[0]+1);
        vec1f lup(b[1]-b[0]+1);
        for (uint_t i : range(llow)) {
            uint_t l = b[0] + i;
            llow[i] = (l != 0 ? 0.5*(lam[l-1]+lam[l]) : lam[l] - 0.5*(lam[l+1]-lam[l]));
            lup[i] = (l != n-1 ? 0.5*(lam[l]+lam[l+1]) : lam[l] + 0.5*(lam[l]-lam[l-1]));
        }

        sed[b[0]-_-b[1]] += integrate_gauss(llow, lup, l0, lw, lum*l0);
    };

    auto get_flux = [&](const vec1f& m, const vec1u& optsed, const vec2f& llum,
        vec2f& flux, vec2f& rfmag, bool no_ir, component cmp) {

        auto pg1 = progress_start(ngal);
        for (uint_t i : range(ngal)) {
            // Build the full rest-frame SED
            vec1f rlam, rsed;
            if (no_ir) {
                // Just the stellar component
                get_opt_sed(m[i]-out.m2lcor[i], optsed[i], rlam, rsed);
            } else if (no_stellar) {
                // Just the dust component
                get_ir_sed(i, rlam, rsed);
            } else {
                // Both stellar and dust components
                vec1f orlam, orsed;
                vec1f irlam, irsed;
                get_opt_sed(m[i]-out.m2lcor[i], optsed[i], orlam, orsed);
                get_ir_sed(i, irlam, irsed);

                // Combine IR SED with optical SED
                merge_add(orlam, irlam, orsed, irsed, rlam, rsed);
            }

            // Apply IGM absorption to continuum
            if (igm != "none" && igm != "constant") {
                uint_t lz = lower_bound(rlam, 0.0600);
                uint_t l0 = lower_bound(rlam, 0.0912);
                uint_t l1 = lower_bound(rlam, 0.1026);
                uint_t l2 = lower_bound(rlam, 0.1216);

                for (uint_t l : range(lz)) {
                    rsed.safe[l] = 0.0;
                }
                for (uint_t l : range(lz, l0)) {
                    rsed.safe[l] *= out.igm_abs(i,0);
                }
                for (uint_t l : range(l0, l1)) {
                    rsed.safe[l] *= out.igm_abs(i,1);
                }
                for (uint_t l : range(l1, l2)) {
                    rsed.safe[l] *= out.igm_abs(i,2);
                }
            }

            // Add emission lines
            if (!no_nebular) {
                for (uint_t l : range(out.lines)) {
                    add_emission_line(rlam, rsed, out.line_lambda[l], llum(i,l), out.vdisp[i]);
                }
            }

            // Obtain rest-frame magnitudes
            if (!rffilters.empty()) {
                // Redshift the SED @ 10pc
                vec1f sed = lsun2uJy(0.0, 1e-5, rlam, rsed);

                for (uint_t b : range(rffilters)) {
                    double mag = uJy2mag(sed2flux(rffilters[b], rlam, sed));
                    rfmag(i,b) = (is_finite(mag) ? mag : +99);
                }
            }

            // Redshift the SED
            vec1f lam = rlam*(1.0 + out.z[i]);
            vec1f sed = lsun2uJy(out.z[i], out.d[i], rlam, rsed);

            if (save_sed) {
                ++nsed;
                sed_saver.push(
                    [&nwritten,&seds_data,&save_sed_bulge_start,&save_sed_bulge_nbyte,
                    &save_sed_disk_start,&save_sed_disk_nbyte,cmp,i,lam,sed]() {

                    auto sp = seds_data.tellp();

                    seds_data.write(
                        reinterpret_cast<const char*>(lam.data.data()),
                        lam.size()*sizeof(decltype(lam[0]))
                    );
                    seds_data.write(
                        reinterpret_cast<const char*>(sed.data.data()),
                        sed.size()*sizeof(decltype(sed[0]))
                    );

                    uint_t nval = seds_data.tellp()-sp;

                    if (cmp == component::bulge) {
                        save_sed_bulge_start[i] = sp;
                        save_sed_bulge_nbyte[i] = nval;
                    } else {
                        save_sed_disk_start[i] = sp;
                        save_sed_disk_nbyte[i] = nval;
                    }

                    ++nwritten;
                });
            }

            // Integrate the SED to get broadband fluxes
            for (uint_t b : range(filters)) {
                double flx = sed2flux(filters[b], lam, sed);
                flux(i,b) = (is_finite(flx) ? flx : 0.0);
            }

            if (verbose) progress(pg1, 127);
        }
    };

    // Get flux for the bulge
    if (!no_stellar) {
        get_flux(out.m_bulge, out.opt_sed_bulge, out.line_lum_bulge,
            out.flux_bulge, out.rfmag_bulge, true, component::bulge);
    }

    get_flux(out.m_disk, out.opt_sed_disk, out.line_lum_disk,
        out.flux_disk, out.rfmag_disk, no_dust, component::disk);

    if (save_sed) {
        if (verbose) {
            note("writing ", nsed, " SEDs to disk...");

            auto pg = progress_start(nsed);
            do {
                print_progress(pg, nwritten);
                thread::sleep_for(0.2);
            } while (nwritten != nsed);

            print_progress(pg, nsed);
        }

        sed_saver.wait();
        seds_data.close();

        fits::write_table(file::remove_extension(seds_file)+"-lookup.fits",
            "id", out.id,
            "bulge_start", save_sed_bulge_start, "bulge_nbyte", save_sed_bulge_nbyte,
            "disk_start", save_sed_disk_start, "disk_nbyte", save_sed_disk_nbyte,
            "elem_size", sizeof(float)
        );
    }

    out.flux = out.flux_disk + out.flux_bulge;
    out.rfmag = uJy2mag(mag2uJy(out.rfmag_disk) + mag2uJy(out.rfmag_bulge));
}

    // Save the catalog to a file
    // --------------------------

    if (verbose) {
        note("saving catalog...");
    }

    out.cmd = make_cmd(argc, argv);

    file::mkdir(file::get_directory(out_file));

    fits::write_table(out_file, ftable(
        out.id, out.ra, out.dec, out.clustered, out.z, out.d, out.m,
        out.sfr, out.rsb, out.m2lcor,
        out.disk_angle, out.disk_radius, out.disk_ratio,
        out.bulge_angle, out.bulge_radius, out.bulge_ratio,
        out.bt, out.m_disk, out.m_bulge, out.av_disk, out.av_bulge, out.passive,
        out.rfuv_bulge, out.rfuv_disk, out.rfvj_bulge, out.rfvj_disk,
        out.sfrir, out.sfruv, out.irx, out.lir, out.ir_sed,
        out.opt_sed_bulge, out.opt_sed_disk,
        out.tdust, out.ir8, out.fpah, out.mdust, out.igm_abs,
        out.zb, out.cmd
    ));

    if (!bands.empty()) {
        fits::update_table(out_file, ftable(
            out.flux, out.flux_disk, out.flux_bulge,
            out.bands, out.lambda
        ));
    }

    if (!out.lines.empty()) {
        fits::update_table(out_file, ftable(
            out.vdisp, out.avlines_disk, out.avlines_bulge,
            out.line_lum, out.line_lum_disk, out.line_lum_bulge,
            out.lines, out.line_lambda
        ));
    }

    if (!rfbands.empty()) {
        fits::update_table(out_file, ftable(
            out.rfmag, out.rfmag_disk, out.rfmag_bulge,
            out.rfbands, out.rflambda
        ));
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

    print("egg-gencat v1.4.0");
    print("usage: egg-gencat [options]\n");

    print("List of generic options:");
    argdoc("verbose", "[flag]", "print additional information in the standard output while "
        "the program is running (default: false)");
    argdoc("help", "[flag]", "print this text and exit");
    argdoc("list_bands", "[flag/string]", "print the list of available photometric bands "
        "and exit. If a string is given in argument, it is understood as a POSIX regular "
        "expression to filter the output list.");
    // Warning: this should not be changed, do not make it part of the public interface
    // argdoc("cosmo", "[string] set of cosmological parameters (possible values: "+
    //     collapse(cosmo_list(), ", ")+", default: std)");
    argdoc("seed", "[uint]", "number used to initialize the random number generator "
        "(default: 42)");
    argdoc("out", "[string]", "name of the output catalog "
        "(default: egg-"+today()+".fits)");
    print("");

    print("List of library related options:");
    argdoc("mass_func", "[string]", "FITS file containing the mass functions (default: "
        "mass_func_candels.fits)");
    argdoc("ir_lib", "[string]", "FITS file containing the IR SED library (default: "
        "ir_lib_ce01.fits)");
    argdoc("opt_lib", "[string]", "FITS file containing the optical SED library (default: "
        "opt_lib_fast.fits)");
    argdoc("filter_db", "[string]", "location of the filter database file (default: "+
        file::directorize(filters_dir)+"db.dat)");
    print("");

    print("List of component related options:");
    argdoc("input_cat", "[string]", "name of the file containing the input catalog. Must be "
        "either (a) a FITS table containing: RA & DEC (degrees), Z, M (log Msun, Salpeter "
        "IMF), PASSIVE (0 or 1), and an optional ID column, or (b) an ASCII table "
        "containing the columns ID, RA, DEC, Z, M and PASSIVE, in that specific order.");
    argdoc("igm", "[string]", "recipe to use for intergalactic medium (IGM) absorption. Possible "
        "values are: 'none' (no IGM), 'constant' (z=1 mean IGM for all galaxies), 'madau95' (IGM "
        "from Madau+95, as implemented in FAST), and 'inoue14' (randomized IGM). Default is "
        "'inoue14'.");
    argdoc("no_pos", "[flag]", "do not generate galaxy positions on the sky");
    argdoc("no_clust", "[flag]", "do not generate clustering in galaxy positions");
    argdoc("no_flux", "[flag]", "do not generate fluxes");
    argdoc("no_nebular", "[flag]", "do not generate emission lines");
    argdoc("no_random", "[flag]", "disable most randomization in the simulation, and use "
        "fully deterministic recipes");
    print("");

    print("List of sky position related options:");
    argdoc("ra0", "[double, degrees]", "sky position of the center of the field "
        "(right ascension, default: 53.558750)");
    argdoc("dec0", "[double, degrees]", "sky position of the center of the field "
        "(declination, default: -27.176001)");
    argdoc("area", "[double, square degrees]", "sky area occupied by the generated field "
        "(default: none)");
    print("");

    print("List of galaxy properties related options:");
    argdoc("mmin", "[double, log10 msun]", "minimum stellar mass generated (default: none)");
    argdoc("mmax", "[double, log10 msun]", "maximum stellar mass generated (default: 12)");
    argdoc("maglim", "[double, AB mag]", "maximum magnitude that will be generated "
        "(default: none). Note that, when set, this parameter overrides 'mmin'.");
    argdoc("selection_band", "[string]", "if 'maglim' is set, name of band in which the "
        "magnitude cut is applied (default: none)");
    argdoc("zmin", "[double]", "minimum redshift generated (default: 0.05)");
    argdoc("zmax", "[double]", "maximum redshift generated (default: 10.5)");
    argdoc("min_dz", "[double]", "minimum size of a redshift bin (default: 0.001)");
    argdoc("max_dz", "[double]", "maximum size of a redshift bin (default: 0.5)");
    argdoc("dz", "[double]", "size of a redshift bin, as a fraction of (1+z) "
        "(default: 0.05)");
    argdoc("dm", "[double, dex]", "size of a mass bin (default: 0.05)");
    argdoc("ms_disp", "[double, dex]", "scatter of the main sequence (default: 0.3)");
    print("");

    print("List of flux related options:");
    argdoc("no_dust", "[flag]", "SEDs will not contain the dust continuum and PAH emission");
    argdoc("no_stellar", "[flag]", "SEDs will not contain the stellar emission");
    argdoc("bands", "[string array]", "optical and IR bands for which to generate fluxes");
    argdoc("rfbands", "[string array]", "optical and IR bands for which to generate "
        "absolute magnitudes");
    argdoc("save_sed", "[flag]", "save the full SEDs of each simulated galaxies in "
        "'seds_file' in binary format. The header file '[seds_file]-header.fits' can be used to "
        "locate the SED of each galaxy and each component in the file (starting position and length "
        "are given in bytes). 'egg-getsed' can take care of this dirty work for you. Be warned: "
        "these SEDs will consume a lot of disk space (about 40KB/galaxy).");
    argdoc("seds_file", "[string]", "binary file in which the SEDs of each galaxy should be "
        "saved if 'save_sed' is activated (default: '[out]-seds.dat'). In addition, a \"header\" "
        "file will also be created in '[seds_file]-header.fits' to provide the position of each SED "
        "in the file.");
    print("");
}
