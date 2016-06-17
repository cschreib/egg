#include <phypp.hpp>

vec1d schechter2(vec1d m, double mstar1, double index1, double phistar1,
    double mstar2, double index2, double phistar2) {
    vec1d tm1 = e10(m-mstar1);
    vec1d tm2 = e10(m-mstar2);
    return log(10.0)*(exp(-tm1)*phistar1*pow(tm1, 1+index1) + exp(-tm2)*phistar2*pow(tm2, 1+index2));
}

int phypp_main(int argc, char* argv[]) {
    double mmin = 4.0, mmax = 13.0;
    double zmax = 10.5;
    double dm = 0.05;
    std::string mf;

    read_args(argc, argv, arg_list(mmin, mmax, dm, zmax, mf));

    vec1f zl, zu;
    vec1f a_mstar1, a_mstar2, a_phistar1, a_phistar2, a_index1, a_index2;
    vec1f p_mstar1, p_mstar2, p_phistar1, p_phistar2, p_index1, p_index2;
    std::string imf;

    if (mf == "tomczak") {
        // Tomczak et al. (2013)

        zl = {0.2, 0.5,  0.75, 1.0,  1.25, 1.5, 2.0, 2.5};
        zu = {0.5, 0.75, 1.0,  1.25, 1.5,  2.0, 2.5, 3.0};

        // Active
        a_mstar1 = {10.59, 10.65, 10.56, 10.44, 10.69, 10.59, 10.58, 10.61};
        a_mstar2 = a_mstar1;
        a_phistar1 = {-2.67, -2.97, -2.81, -2.98, -3.04, -3.37, -4.30, -4.95};
        a_index1 = {-1.08, -0.97, -0.46, 0.53, -0.55, 0.75, 2.06, 2.36};
        a_phistar2 = {-4.46, -3.34, -3.36, -3.11, -3.59, -3.28, -3.28, -3.71};
        a_index2 = {-2.00, -1.58, -1.61, -1.44, -1.62, -1.47, -1.38, -1.67};

        // Passive
        p_mstar1 = {10.75, 10.68, 10.63, 10.63, 10.49, 10.77, 10.69, 9.95};
        p_mstar2 = p_mstar1;
        p_phistar1 = {-2.76, -2.67, -2.81, -3.03, -3.36, -3.41, -3.59, -4.22};
        p_index1 = {-0.47, -0.10, 0.04, 0.11, 0.85, -0.19, -0.37, -0.62};
        p_phistar2 = {-5.21, -4.29, -4.40, -4.80, -3.72, -3.91, -6.95, -4.51};
        p_index2 = {-1.97, -1.69, -1.51, -1.57, -0.54, -0.18, -3.07, 2.51};

        a_phistar1 = e10(a_phistar1);
        a_phistar2 = e10(a_phistar2);
        p_phistar1 = e10(p_phistar1);
        p_phistar2 = e10(p_phistar2);

        imf = "chabrier";
    } else {
        // My stellar mass functions in GS

        zl = {0.3, 0.7, 1.2, 1.9, 2.5, 3.5};
        zu = {0.7, 1.2, 1.9, 2.5, 3.5, 4.5};

        // Active
        a_mstar1 = {11.0000, 11.0000, 11.0000, 11.0000, 11.0000, 11.0000};
        a_mstar2 = {10.6418, 10.7292, 10.6717, 10.8404, 10.9443, 11.0000};
        a_phistar1 = {0.000898887, 0.000718160, 0.000465684, 0.000213874, 0.000212404, 4.44724e-05};
        a_index1 = {-1.40000, -1.40000, -1.50000, -1.57000, -1.60000, -2.00000};
        a_phistar2 = {8.30778e-05, 0.000404045, 0.000417749, 0.000406023, 9.06860e-05, 0.0};
        a_index2 = {0.500000, 0.500000, 0.500000, 0.00000, 0.500000, 0.500000};

        // Passive
        p_mstar1 = {11.0000, 11.0000, 11.0000, 11.0000, 11.0000, 11.0000};
        p_mstar2 = {11.0426, 10.8601, 10.8342, 11.0471, 10.9439, 11.0000};
        p_phistar1 = {7.77453e-05, 3.54586e-05, 2.29979e-05, 1.00000e-05, 0.00000, 0.00000};
        p_index1 = {-1.65000, -1.60000, -1.25000, -1.00000, -1.00000, -1.35000};
        p_phistar2 = {0.00154472, 0.00104263, 0.000624682, 0.000173119, 0.000122278, 3.00000e-05};
        p_index2 = {-0.481039, 0.0594024, 0.296244, -0.166611, -0.263124, -0.300000};

        imf = "salpeter";
    }

    // Note: z = 0 is obtained from Baldry et al. (2012)
    // This mass function is in Chabrier IMF, so we may need to convert that
    // to match the mass functions at other redshifts.
    // NB: the IMF can be anything, but it has to be consistent at all z...
    float b12_factor = 1.0;
    if (imf == "salpeter") {
        b12_factor = 1.0/1.8;
    }

    prepend(zu,         vec1f{min(zl)});
    prepend(zl,         vec1f{0.0});

    prepend(a_mstar1,   vec1f{10.72});
    prepend(a_mstar2,   vec1f{10.72});
    prepend(a_phistar1, vec1f{0.71e-3});
    prepend(a_phistar2, vec1f{0.0});
    prepend(a_index1,   vec1f{-1.45});
    prepend(a_index2,   vec1f{0.3});

    prepend(p_mstar1,   vec1f{10.72});
    prepend(p_mstar2,   vec1f{10.72});
    prepend(p_phistar1, vec1f{3.25e-3});
    prepend(p_phistar2, vec1f{0.08e-3});
    prepend(p_index1,   vec1f{-0.45});
    prepend(p_index2,   vec1f{-1.45});

    // Note: z > 4.0 is obtained by keeping the shape of the last redshift bin and
    //       decreasing phistar, following the total stellar mass density of
    //       Grazian et al. (2015)
    uint_t ilast = a_mstar1.size()-1;
    uint_t nhzb = ceil((zmax - max(zu))/0.5);
    double dhz = (zmax - max(zu))/nhzb;

    append(zl, findgen(nhzb)*dhz + max(zu));
    append(zu, findgen(nhzb)*dhz + max(zu) + dhz);

    auto g15_rhostar = vectorize_lambda([](double z){
        return e10(-0.43*z);
    });

    vec1d decrease = g15_rhostar(0.5*(zl+zu)[(ilast+1)-_])/g15_rhostar(0.5*(zl[ilast] + zu[ilast]));

    for (uint_t i : range(decrease)) {
        a_mstar1.push_back(a_mstar1[ilast]);
        a_mstar2.push_back(a_mstar2[ilast]);
        a_index1.push_back(a_index1[ilast]);
        a_index2.push_back(a_index2[ilast]);
        a_phistar1.push_back(decrease[i]*a_phistar1[ilast]);
        a_phistar2.push_back(decrease[i]*a_phistar2[ilast]);

        p_mstar1.push_back(p_mstar1[ilast]);
        p_mstar2.push_back(p_mstar2[ilast]);
        p_index1.push_back(p_index1[ilast]);
        p_index2.push_back(p_index2[ilast]);
        p_phistar1.push_back(decrease[i]*p_phistar1[ilast]);
        p_phistar2.push_back(decrease[i]*p_phistar2[ilast]);
    }

    struct {
        vec2f zb, mb;
        vec2d active, passive;
        std::string imf;
    } out;

    out.imf = imf;

    out.zb.resize(2, zl.size());
    out.zb(0,_) = zl;
    out.zb(1,_) = zu;
    out.mb = make_bins(rgen(mmin - dm, mmax + dm, ceil((mmax - mmin)/dm) + 2));

    uint_t nm = out.mb.dims[1];
    uint_t nz = out.zb.dims[1];

    out.active.resize(nz, nm);
    out.passive.resize(nz, nm);

    vec1d mx = 0.5*(out.mb(0,_) + out.mb(1,_));

    for (uint_t z : range(nz)) {
        double factor = 0.0;
        if (z == 0) {
            factor = log10(b12_factor);
        }

        out.active(z,_) = schechter2(factor + mx,
            a_mstar1[z], a_index1[z], a_phistar1[z], a_mstar2[z], a_index2[z], a_phistar2[z]);
        out.passive(z,_) = schechter2(factor + mx,
            p_mstar1[z], p_index1[z], p_phistar1[z], p_mstar2[z], p_index2[z], p_phistar2[z]);
    }

    fits::write_table("mass_func_candels.fits", ftable(
        out.zb, out.mb, out.active, out.passive, out.imf
    ));

    return 0;
}
