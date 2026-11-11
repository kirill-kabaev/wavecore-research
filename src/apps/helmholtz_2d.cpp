#include "wavecore/helmholtz/version.hpp"

#include <iostream>
#include <string_view>

namespace {

void print_help(std::ostream& output) {
    output << "Usage: helmholtz_2d [OPTION]\n"
              "\n"
              "WaveCore Research 2D Helmholtz solver scaffold.\n"
              "Numerical solution is not implemented yet.\n"
              "\n"
              "Options:\n"
              "  --help     Show this help message\n"
              "  --version  Show the program version\n";
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc == 1) {
        print_help(std::cout);
        return 0;
    }

    if (argc != 2) {
        std::cerr << "helmholtz_2d: expected zero or one option\n";
        return 2;
    }

    const std::string_view argument{argv[1]};

    if (argument == "--help") {
        print_help(std::cout);
        return 0;
    }

    if (argument == "--version") {
        std::cout << "helmholtz_2d " << wavecore::helmholtz::version() << '\n';
        return 0;
    }

    std::cerr << "helmholtz_2d: unknown option: " << argument << '\n';
    std::cerr << "Try 'helmholtz_2d --help' for usage information.\n";
    return 2;
}
