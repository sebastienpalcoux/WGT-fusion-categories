// Enumerate one parent from the very same frozen GMP recursion used to search it.
#define main frozen_search_entry_point
#include "exact.cpp"
#undef main

int main(int argc, char** argv) {
    try {
        if (argc != 6 || std::string(argv[1]) != "--manifest")
            throw std::invalid_argument("use --manifest R DEPTH Q fixed,prefix");
        int rank = std::stoi(argv[2]), depth = std::stoi(argv[3]), q = std::stoi(argv[4]);
        if (q < 1 || q >= rank) throw std::invalid_argument("multiplier outside rank");
        std::stringstream input(argv[5]);
        std::vector<Z> prefix;
        std::string token;
        while (std::getline(input,token,',')) {
            Z b(token,10);
            if (b < 1 || (!prefix.empty() && b < prefix.back()))
                throw std::invalid_argument("invalid denominator prefix");
            prefix.push_back(b);
        }
        if (prefix.empty() || prefix.size() >= static_cast<unsigned>(depth))
            throw std::invalid_argument("split depth must exceed fixed prefix length");
        Search(rank, depth, q, std::move(prefix)).run();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "INCOMPLETE " << e.what() << "\n";
        return 2;
    }
}
