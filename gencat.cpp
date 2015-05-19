#include <phypp.hpp>

void print_help(std::string filter_db_file);

int main(int argc, char* argv[]) {
    // Initialization
    // --------------
    double ra0 = 53.558750;
    double dec0 = -27.176001;
    double area = 0.08;
    double mmin = 7.0;
    double mmax = 12.0;
    double maglim = dnan;
    double zmin = 0.01;
    double zmax = 9.0;
    double min_dz = 0.05;
    double bin_dz = 0.1;
    double bin_dm = 0.05;

    double ms_disp = 0.3;

    bool no_pos = false;
    bool no_clust = false;
    bool no_flux = false;
    bool no_stellar = false;
    bool no_dust = false;

    bool save_sed = false;

    std::string mass_func_file = "mass_func_candels.fits";
    std::string ir_lib_file = "ir_lib_cs15.fits";
    std::string opt_lib_file = "opt_lib_fast.fits";
    std::string out_file = "";
    std::string filter_db_file = data_dir+"fits/filter-db/db.dat";

    uint_t tseed = 42;
    std::string tcosmo = "std";
    bool verbose = false;
    bool help = false;

    std::string selection_band = "f160w";
    vec1s bands = {"vimos_u", "f435w", "f606w", "f775w", "f814w", "f850lp",
        "f105w", "f125w", "f140w", "f160w", "hawki_Ks", "i1", "i2", "i3", "i4",
        "irs1", "m1", "p1", "p2", "p3", "s1", "s2", "s3"};

    // Read command line arguments
    // ---------------------------
    read_args(argc, argv, arg_list(
        ra0, dec0, area, mmin, maglim, zmin, zmax, name(bin_dz, "dz"), min_dz,
        name(bin_dm, "dm"), ms_disp,
        no_pos, no_clust, no_flux, no_stellar, no_dust, save_sed,
        name(mass_func_file, "mass_func"),
        name(ir_lib_file, "ir_lib"), name(opt_lib_file, "opt_lib"),
        name(out_file, "out"), name(filter_db_file, "filter_db"),
        verbose, name(tseed, "seed"), name(tcosmo, "cosmo"),
        selection_band, bands
    ));

    if (help) {
        print_help(filter_db_file);
        return 0;
    }

    if (no_dust && no_stellar) {
        no_flux = true;
    }

    if (out_file.empty()) {
        out_file = "gencat-"+today()+".fits";
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

    auto filter_db = read_filter_db(filter_db_file);

    // Sort the bands by increasing wavelength
    vec1f lambda;
    for (auto& band : bands) {
        auto fil = get_filter(filter_db, band);
        lambda.push_back(fil.rlam);
    }

    vec1u ids = sort(lambda);
    bands = bands[ids];
    lambda = lambda[ids];

    auto filters = get_filters(filter_db, bands);

    // Initialize SED libraries
    // ------------------------

    // The IR library. When using the generic code for IR libraries, each template must
    // be normalized to unit LIR. If a more subtle calibration needs to be used, then
    // a custom code has to be written.
    struct {
        vec2d lam, sed;
    } ir_lib;

    struct {
        vec2d lam, dust, pah;
        vec1d tdust, lir_dust, lir_pah;
    } ir_lib_cs15;

    uint_t nirsed;
    if (ir_lib_file == "ir_lib_cs15.fits") {
        // My library, calibrated in unit Mdust
        fits::read_table(ir_lib_file, ftable(
            ir_lib_cs15.lam, ir_lib_cs15.dust, ir_lib_cs15.pah,
            ir_lib_cs15.tdust, ir_lib_cs15.lir_dust, ir_lib_cs15.lir_pah
        ));

        nirsed = ir_lib_cs15.lam.dims[0];
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

    // The optical library. Each template must be normalized to unit stellar mass.
    // The library is binned in redshift, U-V and V-J colors.
    struct {
        vec3f lam, sed;
        vec2b use;
        vec2f bvj, buv;
    } opt_lib;

    fits::read_table(opt_lib_file, ftable(
        opt_lib.lam, opt_lib.sed, opt_lib.use, opt_lib.bvj, opt_lib.buv
    ));

    // Make sure that it contains at least one valid SED
    if (count(opt_lib.use) == 0) {
        error("the provided optical SED library does not contain any valid SED");
        note("'use' was false for all the UVJ positions");
        return 1;
    }

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
    // To prevent too narrow bins with very little galaxies, we impose a minimum
    // redshift width.
    vec2d zb = [&](){
        vec1d tzb;
        double z1 = zmin;
        while (z1 < zmax) {
            tzb.push_back(z1);
            z1 *= (1.0 + bin_dz/2.0)/(1.0 - bin_dz/2.0);
            if (z1 - tzb.back() < min_dz) z1 = tzb.back() + min_dz;
        }

        if ((zmax - tzb.back())/(1.0 + zmax) > 0.5*bin_dz && zmax - tzb.back() > min_dz) {
            tzb.push_back(zmax);
        } else {
            tzb.back() = zmax;
        }

        return make_bins(tzb);
    }();

    if (verbose) {
        vec1f tdz = bin_width(zb);
        note("min dz: ", min(tdz), ", max dz: ", max(tdz));
    }

    // Compute the corresponding volume of the Universe
    vec2d vb = vuniverse(zb, cosmo);

    uint_t nz = zb.dims[1];

    if (verbose) {
        note(nz, " redshift slices");
    }

    // Compute evolving mass limit
    // ---------------------------
    vec1f z_mmin(nz);

    if (is_finite(maglim)) {
        if (verbose) {
            note("estimating redshift-dependend mass limit...");
        }

        if (count(bands == selection_band) == 0) {
            error("selection_band (", selection_band, ") is not found in bands");
            return 1;
        }

        // A magnitude limite was requested.
        // Use the optical library to estimate a rough mass completeness, and compute
        // the minimum mass needed to be 90% complete at the requested magnitude.
        // To do so, we compute the flux in the selection band of all the templates
        // in the library (within the corresponding redshift interval). Since the library
        // is in flux per unit mass, we can just divide the requested flux limit by each
        // computed flux to estimate the minimum detectable mass for each template. By
        // picking the 10% lowest value, we should be close to 90% completeness.

        double comp = 0.9;
        uint_t bsel = where(bands == selection_band)[0];
        double flim = mag2uJy(maglim);
        uint_t nb = opt_lib.buv.dims[1];

        for (uint_t z : range(nz)) {
            double tz = bin_center(zb(_,z));
            double td = lumdist(tz, cosmo);

            vec1f tm;
            for (uint_t iu : range(nb))
            for (uint_t iv : range(nb)) {
                if (!opt_lib.use(iu,iv)) continue;

                vec1f rlam = opt_lib.lam(iu,iv,_);
                vec1f rsed = opt_lib.sed(iu,iv,_);

                vec1f lam = rlam*(1.0 + tz);
                vec1f sed = lsun2uJy(tz, td, rlam, rsed);

                double funitmass = sed2flux(filters[bsel], lam, sed);
                double mlim = log10(flim/funitmass);

                tm.push_back(mlim);
            }

            z_mmin[z] = percentile(tm, 1.0 - comp);
        }

        mmin = min(z_mmin);

        if (verbose) {
            note("will generate masses from as low as ", mmin, ", up to ", mmax);
        }
    } else {
        // If no magnitude limit is requested, just use a fixed mass limit.
        z_mmin[_] = mmin;
    }

    // Initialize sky positions
    // ------------------------

    // The shape of the generated survey can be changed here. The survey boundaries
    // are defined as a *convex* polygon whose vertices are given in RA and Dec
    // coordinates below.
    // The default is to generate a square of requested area around the reference position
    // given by ra0 and dec0.
    vec1d hra, hdec; {
        double dd = sqrt(area)/2.0;
        double dr = sqrt(area)/2.0/cos(dec0*dpi/180.0);

        hdec = {dec0-dd, dec0-dd, dec0+dd, dec0+dd, dec0-dd};
        hra = {ra0-dr, ra0+dr, ra0+dr, ra0-dr, ra0-dr};
    }

    vec1u hull = uindgen(hra.size());

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
    if (min(mf.mb) > mmin) {
        error("mass functions do not cover masses < ", min(mf.mb));
        note("requested: ", mmin);
        return 1;
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

    // Output catalog format
    // ---------------------
    struct {
        vec1u id;

        // Sky properties
        vec1d ra, dec;

        // Physical properties
        vec1f z, d, m, sfr, rsb;

        vec1f disk_angle,  disk_radius,  disk_ratio;
        vec1f bulge_angle, bulge_radius, bulge_ratio;
        vec1f bt, m_disk, m_bulge;

        vec1f tdust, fpah, mdust;

        vec1b passive;

        // Flux properties
        vec1f rfuv_bulge, rfuv_disk, rfvj_bulge, rfvj_disk;
        vec1f sfrir, sfruv, irx;
        vec1f lir;
        vec1u ir_sed;
        vec1u opt_sed_bulge;
        vec1u opt_sed_disk;

        vec2f flux;
        vec2f flux_disk;
        vec2f flux_bulge;
        vec1s bands;
        vec1f lambda;

        vec2f zb;

        std::string cmd;
    } out;

    out.zb = zb;

    out.lambda = lambda;
    out.bands = bands;

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

    vec1u ida = uindgen(nactive);
    vec1u idp = uindgen(npassive) + nactive;

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

    out.id = uindgen(out.z.size());
    out.d = lumdist(out.z, cosmo);

    // Generate masses
    // ---------------
    if (verbose) {
        note("generating masses...");
    }

    out.m.resize(ngal);

    for (uint_t iz : range(nz)) {
        vec1b inz = in_bin_open(out.z, zb, iz);
        vec1u idm = where(mf.mx > z_mmin[iz]);

        vec1u idl = where(inz && !out.passive);
        vec1f m = random_pdf_bin(mf.mb(_,idm), mf.mx[idm], mf.mz_active(iz,idm), idl.size());
        out.m[idl] = m;

        idl = where(inz && out.passive);
        m = random_pdf_bin(mf.mb(_,idm), mf.mx[idm], mf.mz_passive(iz,idm), idl.size());
        out.m[idl] = m;
    }

    // Bulge to total mass ratio
    // Calibrated from Lang et al. (2014), assuming no redshift evolution
    out.bt = clamp(e10(0.27*(out.m[ida] - 10.0) - 0.7 + 0.2*randomn(seed, nactive)), 0.0, 1.0);
    vec1f pbt = clamp(e10(0.1*(out.m[idp] - 10.0) - 0.3 + 0.2*randomn(seed, npassive)), 0.0, 1.0);
    append(out.bt, pbt);
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

    // Bulge and disk have the same position angle, which is completely random
    out.disk_angle = (randomu(seed, ngal) - 0.5)*90.0;
    out.bulge_angle = out.disk_angle;

    // Calibration from n>2.5 galaxies and M* > 10.5
    vec1f bulge_ratio_x =
        {0.0, 0.05, 0.15,  0.25,  0.35,  0.45,  0.55,   0.65,   0.75,   0.85,   0.95,  1.0};
    vec1f bulge_ratio_p =
        {0.0, 0.0,  165.0, 428.0, 773.0, 914.0, 1069.0, 1191.0, 1154.0, 1067.0, 639.0, 450.0};
    out.bulge_ratio = random_pdf(seed, bulge_ratio_x, bulge_ratio_p, ngal);

    out.bulge_radius = 8.0*e10(-1.3*log10(1.0+out.z) + 0.56*(out.m - 11.25)
        + 0.2*randomn(seed, ngal));

    // Calibration from n<1.5 galaxies and M* > 9.0
    vec1f disk_ratio_x =
        {0.0, 0.05, 0.15,  0.25,  0.35,   0.45,   0.55,   0.65,  0.75,  0.85,  0.95,  1.0};
    vec1f disk_ratio_p =
        {0.0, 0.0,  313.0, 900.0, 1143.0, 1127.0, 1059.0, 898.0, 775.0, 548.0, 318.0, 200.0};
    out.disk_ratio = random_pdf(seed, disk_ratio_x, disk_ratio_p, ngal);

    vec1f fz = -0.3*log10(1.0+out.z);
    double zmid = 1.7;
    vec1u idhz = where(out.z > zmid);
    fz[idhz] = -0.7*log10((1.0 + out.z[idhz])/(1.0 + zmid)) - 0.3*log10(1.0 + zmid);
    out.disk_radius = 2.8*e10(fz + 0.2*(out.m - 9.35)
        + 0.25*randomn(seed, ngal));

    // Use similar size for bulges of disk-dominated galaxies
    vec1u idd = where(out.bt < 0.5);
    out.bulge_radius[idd] = 2.8*e10(fz[idd] + 0.2*(out.m[idd] - 9.35)
        + 0.25*randomn(seed, idd.size()));

    vec1d psize = propsize(out.z, cosmo);
    out.disk_radius /= psize;
    out.bulge_radius /= psize;

    // Generate SFR
    // ------------
    if (verbose) {
        note("generating SFR...");
    }

    {
        // Active galaxies
        vec1f az = log10(1.0 + out.z);
        vec1f sfr = out.m[ida] - 9.67 + 1.82*az[ida]
            - 0.38*sqr(max(0.0, out.m[ida] - 9.59 - 2.22*az[ida]));

        // Add dispersion
        vec1f rsb = ms_disp*randomn(seed, nactive);

        // Add starbursts
        vec1u idsb = where(random_coin(seed, 0.033, nactive));
        rsb[idsb] += 0.72;

        append(out.rsb, rsb);

        // Make final SFR
        sfr += rsb;
        append(out.sfr, e10(sfr));

        // Passive galaxies
        sfr = 0.5*(out.m[idp] - 11)
            + az[idp] - 0.6
            + 0.45*randomn(seed, npassive);
        append(out.rsb, replicate(0.0, npassive));
        append(out.sfr, e10(sfr));
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
        out.irx = e10(slope*(out.m - am) + airx + 0.4*randomn(seed, ngal));

        out.sfrir = out.sfr/(1.0 + 1.0/out.irx);
        out.sfruv = out.sfr/(1.0 + out.irx);
    }

    // Generate UVJ colors
    // -------------------
    out.rfvj_disk.resize(ngal);
    out.rfvj_bulge.resize(ngal);
    out.rfuv_disk.resize(ngal);
    out.rfuv_bulge.resize(ngal);

    auto gen_blue = [&seed](const vec1f& m, const vec1f& z, vec1f& uv, vec1f& vj) {
        // Calibrate "attenuation vector" from mass and redshift
        // See 'uvj_track.pro'
        vec1f a0 = 0.58*erf(m-10) + 1.39;
        vec1f as = -0.34 + 0.3*max(m-10.35, 0.0);
        vec1f a = min(a0 + as*z, 2.0) + 0.1*randomn(seed, m.size());

        // Move in the UVJ diagram according to "attenuation"
        double slope = 0.65;
        double theta = atan(slope);
        vj = 0.0  + a*cos(theta) + 0.12*randomn(seed, m.size());
        uv = 0.45 + a*sin(theta) + 0.12*randomn(seed, m.size());
    };

    auto gen_red = [&seed](const vec1f& m, vec1f& uv, vec1f& vj) {
        double pvj = 1.25, puv = 1.85;
        vec1f mspread = clamp(0.1*(m - 11.0) + 0.1*randomn(seed, m.size()), -0.1, 0.2);

        vj = pvj + mspread      + 0.1*randomn(seed, m.size());
        uv = puv + mspread*0.88 + 0.1*randomn(seed, m.size());
    };

    // Disks are always blue
    gen_blue(out.m, out.z, out.rfuv_disk, out.rfvj_disk);

    // Bulges of bulge-dominated galaxies are always red
    // Other bulges have a 50% chance of being red
    {
        vec1b red_bulge = out.bt > 0.6 || random_coin(seed, 0.5, ngal);

        vec1f tuv, tvj;
        vec1u idb = where(red_bulge);
        gen_red(out.m[idb], tuv, tvj);
        out.rfuv_bulge[idb] = tuv;
        out.rfvj_bulge[idb] = tvj;

        idb = where(!red_bulge);
        gen_blue(out.m[idb], out.z[idb], tuv, tvj);
        out.rfuv_bulge[idb] = tuv;
        out.rfvj_bulge[idb] = tvj;
    }

    // Assign optical SED
    // ------------------

if (!no_flux) {
    if (verbose) {
        note("assigning optical SEDs...");
    }

    auto get_sed = [&](const vec1f& uv, const vec1f& vj, vec1u& sed) {
        uint_t snb = opt_lib.buv.dims[1];

        sed.resize(uv.size());
        sed[_] = npos;

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
                }

                if (verbose) progress(pg);
            }
        }

        return true;
    };

    if (!get_sed(out.rfuv_bulge, out.rfvj_bulge, out.opt_sed_bulge)) {
        return 1;
    }

    if (!get_sed(out.rfuv_disk,  out.rfvj_disk,  out.opt_sed_disk)) {
        return 1;
    }
}

    // Generate IR properties
    // ---------------------

    if (verbose) {
        note("generate IR properties...");
    }

    // Tdust and f_PAH as observed in stacks and detections of Herschel galaxies
    out.tdust.resize(ngal);
    out.fpah.resize(ngal);

    out.tdust[ida] = min(24.6*pow(1.0+out.z[ida], 0.37), 29.7*pow(1.0+out.z[ida], 0.20))
        // Starbursts are warmer
        + 6.53*out.rsb[ida]
        // Add some random scatter
        + (4.05 + 1.43*min(out.z[ida], 2.5))*randomn(seed, ida.size());

    out.fpah[ida] = (0.04 + 0.035*(1.0-0.85*clamp(out.z[ida], 1.0, 2.0)))
        // Starburst have weaker PAH
        *e10(-0.34*out.rsb[ida])
        // Add some random scatter
        *e10(0.20*randomn(seed, ida.size()));

    out.fpah = clamp(out.fpah, 0.0, 1.0);

    out.lir = out.sfrir/1.72e-10;
    vec1u idb = where(!is_finite(out.lir));
    out.lir[idb] = 0.0;

    // Assign IR SED
    // -------------

    if (verbose) {
        note("assigning IR SED...");
    }

    // Floating point SED index, we will arrange it into a integer later
    vec1d fir_sed;

    if (nirsed == 1) {
        fir_sed = replicate(0, nactive);
    } else if (ir_lib_file == "ir_lib_ce01.fits") {
        // The Chary & Elbaz 2001 library, redshift evolution calibrated from stacks
        fir_sed = interpolate({26, 26, 40, 54, 53, 52, 52}, {0.57, 1.0, 1.5, 2.1, 2.9, 4.0, 6.0}, out.z[ida])
            // Temperature offset as function of RSB (not calibrated, but see Magnelli+13)
            + 15*clamp(out.rsb[ida]/ms_disp, -2.0, 2.0);
    } else if (ir_lib_file == "ir_lib_m12.fits") {
        // The Magdis et al. 2012 library, using their reported redshift evolution
        fir_sed = interpolate(findgen(nirsed), {0.0125, 0.1625, 0.4625, 0.8125, 1.15, 1.525, 2.0, 2.635}, out.z[ida])
            // Temperature offset as function of RSB (not calibrated, but see Magnelli+13)
            + clamp(out.rsb[ida]/ms_disp, -2.0, 2.0);
    } else if (ir_lib_file == "ir_lib_cs15.fits") {
        // My own library, using calibration from stacks and detections
        fir_sed = interpolate(findgen(nirsed), ir_lib_cs15.tdust, out.tdust[ida]);
    } else {
        error("no calibration code available for the IR library '", ir_lib_file, "'");
        return 1;
    }

    append(out.ir_sed, clamp(round(fir_sed), 0, nirsed-1));
    append(out.ir_sed, replicate(0u, npassive));


if (!no_flux) {

    // Compute fluxes
    // --------------

    out.flux_disk.resize(ngal, bands.size());
    out.flux_bulge.resize(ngal, bands.size());

    if (verbose) {
        note("computing fluxes...");
    }

    if (ir_lib_file == "ir_lib_cs15.fits") {
        out.mdust.resize(ngal);
    }

    vec2f slam, ssed;
    if (save_sed) {
        double l0 = 0.5*min(lambda);
        double l1 = 2.0*max(lambda);
        vec1f x = rgen_log(l0, l1, 800*min(1.0, log10(l1/l0)/4.6));
        slam = replicate(x, ngal);
        ssed.resize(slam.dims);
    }

    auto get_opt_sed = [&](double m, uint_t ised, vec1f& orlam, vec1f& orsed) {
        vec1u did = mult_ids(opt_lib.use, ised);
        orlam = opt_lib.lam(did[0],did[1],_);
        orsed = e10(m)*opt_lib.sed(did[0],did[1],_);
    };

    auto get_ir_sed = [&](uint_t i, vec1f& irlam, vec1f& irsed) {
        if (ir_lib_file == "ir_lib_cs15.fits") {
            // My library
            irlam = ir_lib_cs15.lam(out.ir_sed[i],_);
            irsed = ir_lib_cs15.dust(out.ir_sed[i],_)*(1.0 - out.fpah[i])
                + ir_lib_cs15.pah(out.ir_sed[i],_)*out.fpah[i];

            // SED is in unit Mdust
            out.mdust[i] = out.lir[i]/(ir_lib_cs15.lir_dust[out.ir_sed[i]]*(1.0 - out.fpah[i])
                + ir_lib_cs15.lir_pah[out.ir_sed[i]]*out.fpah[i]);
            irsed *= out.mdust[i];
        } else {
            // Generic library, SEDs in units of LIR
            irlam = ir_lib.lam(out.ir_sed[i],_);
            irsed = out.lir[i]*ir_lib.sed(out.ir_sed[i],_);
        }
    };

    auto get_flux = [&](const vec1f& m, const vec1u& optsed, vec2f& flux, bool no_ir) {
        auto pg1 = progress_start(ngal);
        for (uint_t i : range(ngal)) {
            vec1f rlam, rsed;
            if (no_ir) {
                // Just the stellar component
                get_opt_sed(m[i], optsed[i], rlam, rsed);
            } else if (no_stellar) {
                // Just the dust component
                get_ir_sed(i, rlam, rsed);
            } else {
                // Both stellar and dust components
                vec1f orlam, orsed;
                vec1f irlam, irsed;
                get_opt_sed(m[i], optsed[i], orlam, orsed);
                get_ir_sed(i, irlam, irsed);

                // Combine IR SED with optical SED
                merge_add(orlam, irlam, orsed, irsed, rlam, rsed);
            }

            // Redshift
            vec1f lam = rlam*(1.0 + out.z[i]);
            vec1f sed = lsun2uJy(out.z[i], out.d[i], rlam, rsed);

            if (save_sed) {
                ssed(i,_) = interpolate(sed, lam, slam(i,_));
            }

            // Integrate SED to get broadband fluxes
            for (uint_t b : range(filters)) {
                double flx = sed2flux(filters[b], lam, sed);
                flux(i,b) = (is_finite(flx) ? flx : 0.0);
            }

            if (verbose) progress(pg1, 127);
        }
    };

    // Get flux for the bulge
    if (!no_stellar) {
        get_flux(out.m_bulge, out.opt_sed_bulge, out.flux_bulge, true);
    }

    get_flux(out.m_disk,  out.opt_sed_disk,  out.flux_disk, no_dust);

    out.flux = out.flux_disk + out.flux_bulge;

    if (save_sed) {
        if (verbose) {
            note("saving SEDs...");
        }

        fits::write_table(file::remove_extension(out_file)+"_seds.fits",
            "lam", slam, "sed", ssed
        );

        slam.clear();
        ssed.clear();
    }
}

    // Generate random positions
    // -------------------------

if (!no_pos) {
    if (verbose) {
        note("generating sky positions...");
    }

    out.ra.resize(ngal);
    out.dec.resize(ngal);

    auto pg = progress_start(nz);
    for (uint_t im : range(2))
    for (uint_t iz : range(nz)) {
        // Per redshift bin, we treat massive and low mass galaxies separately.
        // Assuming more clustering for the most massive objects (im == 0).

        vec1u idz = where(in_bin_open(out.z, zb, iz) && (im == 0 ? out.m < 10.5 : out.m >= 10.5));
        uint_t z_ngal = idz.size();

        // Use clustering for only 40% of the sample (low mass) or 70% (high mass)
        double rnd_frac = (no_clust ? 1.0 : im == 0 ? 0.6 : 0.3);
        uint_t nrnd = std::min(z_ngal, uint_t(ceil(rnd_frac*z_ngal)));
        uint_t nclust = z_ngal - nrnd;

        randpos_power_options opt;
        opt.power = 0.5;    // power law of index 0.5
        opt.levels = 4;     // 4 levels in the Soneira & Pebbles
        opt.overload = 2.0; // generate twice more than needed, then pick half
        opt.nsrc = nclust;

        vec1d ra, dec;
        auto status = randpos_power(seed, hull, hra, hdec, ra, dec, opt);
        if (!status.success) {
            error("failed: ", status.failure);
            return 1;
        }

        randpos_uniform_options opt2;
        opt2.nsrc = nrnd;

        vec1d tra, tdec;
        status = randpos_uniform(seed, {min(hra), max(hra)}, {min(hdec), max(hdec)},
            [&](const vec1d& xra, const vec1d& xdec) {
            return in_convex_hull(xra, xdec, hull, hra, hdec);
        }, tra, tdec, opt2);

        if (!status.success) {
            error("failed: ", status.failure);
            return 1;
        }

        append(ra, tra);
        append(dec, tdec);

        vec1u rid = shuffle(seed, uindgen(z_ngal));

        out.ra[idz] = ra[rid];
        out.dec[idz] = dec[rid];

        if (verbose) progress(pg);
    }
}

    // Save the catalog to a file
    // --------------------------

    if (verbose) {
        note("saving catalog...");
    }

    out.cmd = make_cmd(argc, argv);

    file::mkdir(file::get_directory(out_file));

    fits::write_table(out_file, ftable(
        out.id, out.ra, out.dec, out.z, out.d, out.m, out.sfr, out.rsb,
        out.disk_angle, out.disk_radius, out.disk_ratio,
        out.bulge_angle, out.bulge_radius, out.bulge_ratio,
        out.bt, out.m_disk, out.m_bulge, out.passive,
        out.rfuv_bulge, out.rfuv_disk, out.rfvj_bulge, out.rfvj_disk,
        out.sfrir, out.sfruv, out.irx, out.lir, out.ir_sed,
        out.opt_sed_bulge, out.opt_sed_disk,
        out.flux, out.flux_disk, out.flux_bulge,
        out.tdust, out.fpah, out.mdust,
        out.bands, out.lambda, out.zb, out.cmd
    ));

    return 0;
}

void print_help(std::string filter_db_file) {
    using namespace format;

    print("gencat v1.0");
    paragraph("usage: gencat [options]");

    header("List of generic options:");
    bullet("verbose", "[flag] print additional information in the standard output while "
        "the program is running");
    bullet("help", "[flag] print this text and exit");
    bullet("cosmo", "[string] set of cosmological parameters (possible values: "+
        collapse(cosmo_list(), ", ")+", default: std)");
    bullet("seed", "[uint] number used to initialize the random number generator "
        "(default: 42)");
    bullet("out", "[string] number used to initialize the random number generator "
        "(default: gencat-"+today()+".fits)");
    print("");

    header("List of library related options:");
    bullet("mass_func", "[string] FITS file containing the mass functions (default: "
        "mass_func_candels.fits)");
    bullet("ir_lib", "[string] FITS file containing the IR SED library (default: "
        "ir_lib_ce01.fits)");
    bullet("opt_lib", "[string] FITS file containing the optical SED library (default: "
        "opt_lib_fast.fits)");
    bullet("filter_db", "[string] location of the filter database file (default: "+
        data_dir+"fits/filter-db/db.dat)");
    print("");

    header("List of component related options:");
    bullet("no_pos", "[flag] do not generate galaxy positions on the sky");
    bullet("no_clust", "[flag] do not generate clustering in galaxy positions");
    bullet("no_flux", "[flag] do not generate fluxes");
    print("");

    header("List of sky position related options:");
    bullet("ra0", "[double, degrees] sky position of the center of the field "
        "(right ascension, default: 53.558750)");
    bullet("dec0", "[double, degrees] sky position of the center of the field "
        "(declination, default: -27.176001)");
    bullet("area", "[double, square degrees] sky area occupied by the generated field "
        "(default: 0.08)");
    print("");

    header("List of galaxy properties related options:");
    bullet("mmin", "[double, log10 msun] minimum stellar mass generated (default: 7)");
    bullet("mmax", "[double, log10 msun] maximum stellar mass generated (default: 12)");
    bullet("maglim", "[double, AB mag] maximum magnitude that will be generated "
        "(default: none). Note that, when set, this parameter overrides 'mmin'.");
    bullet("zmin", "[double] minimum redshift generated (default: 0.01)");
    bullet("zmax", "[double] maximum redshift generated (default: 10)");
    bullet("min_dz", "[double] minimum size of a redshift bin (default: 0.05)");
    bullet("dz", "[double] size of a redshift bin, as a fraction of (1+z) "
        "(default: 0.1)");
    bullet("dm", "[double, dex] size of a mass bin (default: 0.05)");
    bullet("ms_disp", "[double, dex] scatter of the main sequence (default: 0.3)");
    print("");

    header("List of flux related options:");
    bullet("selection_band", "[string] if 'maglim' is set, name of band in which the "
        "magnitude cut is applied (default: none)");
    bullet("bands", "[string array] optical and IR bands for which to generate fluxes");
    bullet("save_sed", "[flag] save the full SEDs of each simulated galaxies in "
        "[out]_seds.fits (WARNING: will use a lot of memory, be sure to only use this "
        "flag for small simulations)");
    print("");

    header("List of available bands:");
    if (!file::exists(filter_db_file)) {
        warning("could not find filter database file '", filter_db_file, "'");
        note("when running gencat, please set the options 'no_flux'");
        note("and do not use the 'maglim' feature");
    } else {
        vec1s fils;
        auto filter_db = read_filter_db(filter_db_file);
        for (auto fil : filter_db) {
            fils.push_back(fil.first);
        }

        paragraph(collapse(fils, ", "));
    }
}
