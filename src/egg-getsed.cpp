#include <phypp.hpp>

void print_help();

int main(int argc, char* argv[]) {
    uint_t id = npos;
    std::string component, seds, out;
    read_args(argc, argv, arg_list(seds, out, id, component));

    if (seds.empty() || id == npos || component.empty()) {
        print_help();
        return 0;
    }

    if (out.empty()) {
        out = file::remove_extension(seds)+"-"+component+"-"+strn(id)+".fits";
    } else {
        file::mkdir(file::get_directory(out));
    }

    vec1u tstart, tnbyte;
    fits::read_table(file::remove_extension(seds)+"-header.fits",
        component+"_start", tstart, component+"_nbyte", tnbyte
    );

    uint_t start = tstart[id], nbyte = tnbyte[id], npt = nbyte/sizeof(float);

    vec1f lambda(npt), flux(npt);

    std::ifstream file(seds);
    file.seekg(start);
    file.read(reinterpret_cast<char*>(lambda.data.data()), nbyte/2);
    file.read(reinterpret_cast<char*>(flux.data.data()), nbyte/2);

    fits::write_table(out, ftable(lambda, flux));

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

    print("egg-getsed v1.0rc1");
    print("usage: egg-getsed [options]\n");

    print("List of options:");
    argdoc("seds", "[string]", "file containing the SEDs (mandatory)");
    argdoc("id", "[uint]", "ID of the galaxy to extract (zero-based and mandatory) ");
    argdoc("component", "[string]", "name of the galaxy component to extract (disk or bulge, "
        "mandatory)");
    argdoc("out", "[string]", "FITS file in which the SED will be extracted (default: "
        "[out]-[component]-[id].fits)");
    print("");
}
