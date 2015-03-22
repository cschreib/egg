#include <phypp.hpp>

double schechter2(double m, double mstar, double index1, double phistar1,
    double index2, double phistar2) {

    double tm = e10(m-mstar);
    return log(10.0)*exp(-tm)*(phistar1*pow(tm, 1+index1) + phistar2*pow(tm, 1+index2));
}

vec1d schechter2(const vec1d& m, double mstar, double index1, double phistar1,
    double index2, double phistar2) {

    vec1d res(m.dims);
    for (uint_t i : range(m)) {
        res[i] = schechter2(m[i], mstar, index1, phistar1, index2, phistar2);
    }

    return res;
}

double schechter2m(double m, double mstar1, double index1, double phistar1,
    double mstar2, double index2, double phistar2) {

    double tm1 = e10(m-mstar1);
    double tm2 = e10(m-mstar2);
    return log(10.0)*(exp(-tm1)*phistar1*pow(tm1, 1+index1) + exp(-tm2)*phistar2*pow(tm2, 1+index2));
}

vec1d schechter2m(const vec1d& m, double mstar1, double index1, double phistar1,
    double mstar2, double index2, double phistar2) {

    vec1d res(m.dims);
    for (uint_t i : range(m)) {
        res[i] = schechter2m(m[i], mstar1, index1, phistar1, mstar2, index2, phistar2);
    }

    return res;
}

int main(int argc, char* argv[]) {
    double mmin = 4.0, mmax = 13.0;
    double dm = 0.05;

    read_args(argc, argv, arg_list(mmin, mmax, dm));

    vec1f zl = {0.0, 0.3, 0.7, 1.2, 1.9, 2.5, 3.5, 4.5, 5.5, 6.5};
    vec1f zu = {0.3, 0.7, 1.2, 1.9, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5};

    // Active
    // Note: first redshift bin is from Baldry et al. (2012) (converted to Salpeter)
    // Note: z=[0.3,4.5] is my MF in GS
    vec1f a_mstar    = {11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0};
    vec1f a_phistar1 = {0.00062, 0.000950642, 0.000713135, 0.000621699, 0.000349725, 0.000254334, 3.12181e-05};
    vec1f a_index1   = {-1.37000, -1.37000, -1.37000, -1.37000, -1.37000, -1.44088, -1.83184};
    vec1f a_phistar2 = {0.0, 7.18916e-11, 0.000157160, 7.48937e-05, 9.97429e-05, 4.53267e-05, 9.24631e-06};
    vec1f a_index2   = {0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000};

    // Note: z > 4.0 is obtained by keeping the shape of the z=[3.0,4.0] bin and
    //       decreasing phistar, more or less following Grazian et al. (2014)
    uint_t ilast = a_mstar.size()-1;
    vec1f decrease = {0.63, 0.32, 0.13};
    for (uint_t i : range(3)) {
        a_mstar.push_back(a_mstar[ilast]);
        a_index1.push_back(a_index1[ilast]);
        a_index2.push_back(a_index2[ilast]);
        a_phistar1.push_back(decrease[i]*a_phistar1[ilast]);
        a_phistar2.push_back(decrease[i]*a_phistar2[ilast]);
    }

    // Passive
    // Note: first redshift bin is from Baldry et al. (2012) (converted to Salpeter)
    vec1f p_mstar1 = {11.0000, 11.0000, 11.0000, 11.0000, 11.0000, 11.0000, 11.0000};
    vec1f p_mstar2 = {11.0000, 10.9133, 10.8555, 10.8487, 11.0127, 10.9438, 11.0000};
    vec1f p_phistar1 = {0.000168882, 0.000168882, 3.57714e-05, 1.70957e-05, 6.70509e-13, 2.19864e-13, 0.000};
    vec1f p_index1 = {-1.450000, -1.50000, -1.60000, -1.35000, -1.35000, -1.35000, -1.35000};
    vec1f p_phistar2 = {0.00325, 0.00167701, 0.00105043, 0.000632619, 0.000122571, 0.000122904, 0.0001};
    vec1f p_index2 = {0.08, -0.136292, 0.0702109, 0.261546, -0.334458, -0.262291, -0.30000};

    // Note: z > 4.0 is obtained by keeping the shape of the z=[3.0,4.0] bin and
    //       decreasing phistar, more or less following the trend of Grazian et al. (2014)
    ilast = p_mstar1.size()-1;
    for (uint_t i : range(3)) {
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
        std::string imf = "salpeter";
    } out;

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
        out.active(z,_)  = schechter2(mx,
            a_mstar[z], a_index1[z], a_phistar1[z], a_index2[z], a_phistar2[z]);
        out.passive(z,_) = schechter2m(mx,
            p_mstar1[z], p_index1[z], p_phistar1[z], p_mstar1[z], p_index2[z], p_phistar2[z]);
    }

    fits::write_table("mass_func_candels.fits", ftable(
        out.zb, out.mb, out.active, out.passive, out.imf
    ));

    return 0;
}
