#include <phypp.hpp>

void print_help();

struct sky_conf {
    std::string filename;
    vec1s param, value;
};

bool read_sky_conf(const std::string& filename, sky_conf& conf) {
    conf.filename = filename;

    std::ifstream file(filename);

    uint_t l = 0;
    std::string line;
    while (std::getline(file, line)) {
        ++l;

        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        auto pos = line.find_first_of(" \t");
        if (pos == line.npos) {
            note("reading ", filename, ":", l);
            error("ill formed line, expected \"PARAMETER   VALUE\"");
            return false;
        }

        std::string param = line.substr(0, pos);
        std::string value = trim(line.substr(pos));

        conf.param.push_back(param);
        conf.value.push_back(value);
    }

    return true;
}

template<typename T>
bool sky_get_param(const sky_conf& conf, const std::string& param, T& value) {
    uint_t pid = where_first(conf.param == param);
    if (pid == npos) {
        note("reading ", conf.filename);
        error("missing ", param, " parameter");
        return false;
    } else {
        if (!from_string(conf.value[pid], value)) {
            note("reading ", conf.filename);
            error("could not read value of ", param, " into a ", pretty_type(value));
            return false;
        }

        return true;
    }
}

int phypp_main(int argc, char* argv[]) {
    if (argc < 2) {
        print_help();
        return 1;
    }

    // Read command line arguments
    std::string config_file;
    double background = 0.0;

    read_args(argc, argv, arg_list(
        name(config_file, "conf"), background
    ));

    // Read template SkyMaker configuration file
    sky_conf conf;
    if (!read_sky_conf(config_file, conf)) return 1;

    // Extract the data we need
    std::string img_file;
    double exposure = dnan;
    double magzp = dnan;
    double bgmag = dnan;
    double gain = dnan;
    double aspix = dnan;

    bool bad = false;
    if (!sky_get_param(conf, "IMAGE_NAME", img_file)) bad = true;
    if (!sky_get_param(conf, "EXPOSURE_TIME", exposure)) bad = true;
    if (!sky_get_param(conf, "MAG_ZEROPOINT", magzp)) bad = true;
    if (!sky_get_param(conf, "BACK_MAG", bgmag)) bad = true;
    if (!sky_get_param(conf, "GAIN", gain)) bad = true;
    if (!sky_get_param(conf, "PIXEL_SIZE", aspix)) bad = true;
    if (bad) return 1;

    // Now start the real job
    fits::image img_fits(img_file);

    vec2d img;
    img_fits.read(img);

    // ADU -> ADU/sec (/exposure)
    // NB: SkyMaker already takes into account the gain
    img /= exposure;

    // Background subtraction
    img += background - e10(0.4*(magzp - bgmag))*sqr(aspix);

    // Write back the post-processed image
    img_fits.update(img);

    // Set some common FITS keywords
    img_fits.write_keyword("BUNIT", "DN/sec", "Units of image data");
    img_fits.write_keyword("EXPTIME", 1.0, "[sec] Effective integration time per pixel");
    img_fits.write_keyword("GAIN", gain, "[e/DN] Instrumental gain conversion");
    img_fits.write_keyword("MAGZERO", magzp, "AB magnitude for 1 DN/sec");
    img_fits.write_keyword("FLUXCONV", e10(0.4*(23.9-magzp)),
        "[microJy per DN/sec] Flux conversion factor");

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

    print("egg-postskymaker v1.0rc1");
    print("usage: egg-postskymaker conf=...\n");

    print("List of mandatory parameters (no default):");
    argdoc("conf", "[string]", "path to the SkyMaker configuration file");

    print("");
}
