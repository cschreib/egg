#ifndef EGG_UTILS_HPP
#define EGG_UTILS_HPP

// Function to combine two SED templates (covering different wavelength regions) into one
template<typename TX1, typename TX2, typename TY1, typename TY2, typename TX, typename TY>
void merge_add(const vec<1,TX1>& x1, const vec<1,TX2>& x2,
    const vec<1,TY1>& y1, const vec<1,TY2>& y2,
    vec<1,TX>& x, vec<1,TY>& y) {

    phypp_check(x1.dims == y1.dims, "incompatible dimensions between X1 and Y1 (",
        x1.dims, " vs. ", y1.dims, ")");
    phypp_check(x2.dims == y2.dims, "incompatible dimensions between X2 and Y2 (",
        x2.dims, " vs. ", y2.dims, ")");

    uint_t n1 = x1.size(), n2 = x2.size();
    x.clear(); x.reserve(n1+n2);
    y.clear(); y.reserve(n1+n2);

    uint_t i1 = 0, i2 = 0;
    while (i1 < n1 || i2 < n2) {
        if (i1 == n1) {
            x.push_back(x2.safe[i2]);
            y.push_back(y2.safe[i2]);
            ++i2;
        } else if (i2 == n2) {
            x.push_back(x1.safe[i1]);
            y.push_back(y1.safe[i1]);
            ++i1;
        } else {
            if (x1.safe[i1] < x2.safe[i2]) {
                x.push_back(x1.safe[i1]);

                if (i2 == 0) {
                    y.push_back(y1.safe[i1]);
                } else {
                    y.push_back(y1.safe[i1] + interpolate(
                        y2.safe[i2-1], y2.safe[i2], x2.safe[i2-1], x2.safe[i2], x1.safe[i1]
                    ));
                }

                ++i1;
            } else if (x1.safe[i1] > x2.safe[i2]) {
                x.push_back(x2.safe[i2]);

                if (i1 == 0) {
                    y.push_back(y2.safe[i2]);
                } else {
                    y.push_back(y2.safe[i2] + interpolate(
                        y1.safe[i1-1], y1.safe[i1], x1.safe[i1-1], x1.safe[i1], x2.safe[i2]
                    ));
                }

                ++i2;
            } else {
                x.push_back(x1.safe[i1]);
                y.push_back(y1.safe[i1] + y2.safe[i2]);

                ++i1;
                ++i2;
            }
        }
    }
}

#endif
