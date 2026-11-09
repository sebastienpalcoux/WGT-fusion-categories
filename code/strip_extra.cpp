// Checks finite word-support exclusions, not dimensions. GCC or Clang:
// g++ -std=c++20 -O3 -march=native strip_extra.cpp -o strip_extra
// ./strip_extra              or: ./strip_extra PRIME MINIMUM
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <set>
#include <vector>

using Mask = unsigned __int128; // one bit per opposite pair, H < 128
constexpr int positions = 7;

int weight(Mask mask) {
    return std::popcount(std::uint64_t(mask)) + std::popcount(std::uint64_t(mask >> 64));
}
int first_pair(Mask mask) {
    if (auto low = std::uint64_t(mask)) return std::countr_zero(low);
    return 64 + std::countr_zero(std::uint64_t(mask >> 64));
}

struct WordCheck {
    int prime, pairs, minimum, overlap = 1;
    Mask full;
    std::vector<Mask> signed_bit, twice;
    std::vector<std::vector<Mask>> pair_sums;
    std::vector<std::vector<int>> least_third;
    std::array<int, positions> support{0, 1};
    std::uint64_t seeds = 0, self_passes = 0, survivors = 0;

    WordCheck(int s, int h) : prime(s), pairs((s - 1) / 2), minimum(h),
        full((Mask(1) << pairs) - 1), signed_bit(s), twice(pairs),
        pair_sums(pairs, std::vector<Mask>(pairs)), least_third(s, std::vector<int>(s)) {
        while (overlap * (overlap + 1) < minimum) ++overlap;
        std::vector<int> representative(pairs);
        for (int base = 2; base < prime; ++base) {
            std::fill(signed_bit.begin(), signed_bit.end(), Mask(0));
            int x = 1, j = 0;
            while (j < pairs && !signed_bit[x]) {
                representative[j] = x;
                signed_bit[x] = signed_bit[prime - x] = Mask(1) << j;
                x = x * base % prime;
                ++j;
            }
            if (j == pairs) break;
        }
        for (int i = 0; i < pairs; ++i) {
            int x = representative[i];
            twice[i] = signed_bit[2 * x % prime];
            for (int j = 0; j < pairs; ++j) {
                int y = representative[j];
                pair_sums[i][j] = signed_bit[(x + y) % prime]
                                | signed_bit[(x - y + prime) % prime];
            }
        }

        std::vector<int> inverse(prime), least(prime);
        for (int a = 1; a < prime; ++a)
            for (int b = 1; b < prime; ++b)
                if (a * b % prime == 1) inverse[a] = b;
        for (int x = 2; x < prime; ++x) {
            std::array<int, 3> triple{0, 1, x};
            least[x] = prime;
            do {
                auto [a, b, c] = triple;
                int image = (c - a + prime) * inverse[(b - a + prime) % prime] % prime;
                least[x] = std::min(least[x], image);
            } while (std::next_permutation(triple.begin(), triple.end()));
        }
        for (int y = 1; y < prime; ++y)
            for (int z = y + 1; z < prime; ++z)
                least_third[y][z] = least[z * inverse[y] % prime];
    }

    Mask dilate(Mask m, int shift) const {
        return ((m << shift) & full) | (m >> (pairs - shift));
    }
    Mask canonical(Mask m) const {
        Mask best = m;
        for (int u = 1; u < pairs; ++u) best = std::min(best, dilate(m, u));
        return best;
    }
    // Input describes {0,+/-x_i}. Add the pairs 2x_i and x_i+/-x_j.
    Mask differences(Mask input) const {
        Mask result = input;
        for (Mask rest = input; rest; rest &= rest - 1) {
            int i = first_pair(rest);
            result |= twice[i];
            for (Mask others = rest & (rest - 1); others; others &= others - 1)
                result |= pair_sums[i][first_pair(others)];
        }
        return result;
    }
    bool self_test(Mask m) const {
        for (int u = 1; u <= pairs / 2; ++u) {
            Mask common = m & dilate(m, u);
            if (weight(common) < overlap || weight(differences(common)) < minimum)
                return false;
        }
        return true;
    }
    bool potentially_feasible(Mask seed) const {
        std::vector<Mask> forced{seed};
        std::set<Mask> known{seed};
        for (std::size_t i = 0; i < forced.size(); ++i) {
            if (i > 40) return true; // Deriving no contradiction is inconclusive.
            for (std::size_t j = 0; j <= i; ++j) {
                for (int u = 0; u < pairs; ++u) {
                    Mask common = forced[i] & dilate(forced[j], u);
                    int size = weight(common);
                    if (size < overlap) return false;
                    if (size > 6) continue; // Omit these consequences.
                    Mask next = differences(common);
                    int count = weight(next);
                    if (count < minimum) return false;
                    if (2 * (1 + 2 * count) >= prime) continue; // Test at both parities.
                    next = canonical(next);
                    if (known.insert(next).second) {
                        if (!self_test(next)) return false;
                        forced.push_back(next);
                    }
                }
            }
        }
        return true;
    }

    void extend(int used, int next_point, Mask partial) {
        if (used == positions) {
            ++seeds;
            if (weight(partial) < minimum || !self_test(partial)) return;
            ++self_passes;
            if (potentially_feasible(canonical(partial))) ++survivors;
            return;
        }
        for (int x = next_point; x <= prime - (positions - used); ++x) {
            int third = used == 2 ? x : support[2];
            bool minimal_third = true;
            Mask updated = partial;
            for (int i = 0; i < used && minimal_third; ++i) {
                updated |= signed_bit[x - support[i]];
                for (int j = i + 1; j < used; ++j)
                    if (least_third[support[j] - support[i]][x - support[i]] < third) {
                        minimal_third = false;
                        break;
                    }
            }
            if (!minimal_third) continue;
            int remaining_edges = 21 - (used + 1) * used / 2;
            if (weight(updated) + remaining_edges < minimum) continue;
            support[used] = x;
            extend(used + 1, x + 1, updated);
        }
    }
    bool check_four_pair_intersections() const {
        std::uint64_t tested = 0, self_pass = 0, undecided = 0;
        for (int i = 1; i < pairs; ++i)
            for (int j = i + 1; j < pairs; ++j)
                for (int k = j + 1; k < pairs; ++k) {
                    ++tested;
                    Mask common = 1 | (Mask(1) << i) | (Mask(1) << j) | (Mask(1) << k);
                    Mask seed = differences(common);
                    if (weight(seed) < minimum || !self_test(seed)) continue;
                    ++self_pass;
                    if (potentially_feasible(canonical(seed))) ++undecided;
                }
        std::cout << "four-pair " << prime << " " << minimum << " " << tested
                  << " " << self_pass << " " << undecided << std::endl;
        return undecided == 0;
    }

    bool run() {
        extend(2, 2, signed_bit[1]);
        std::cout << prime << " " << minimum << " " << seeds << " "
                  << self_passes << " " << survivors << std::endl;
        return survivors == 0;
    }
};

int main(int argc, char** argv) {
    if (argc == 3) {
        int prime = std::atoi(argv[1]), minimum = std::atoi(argv[2]);
        if (prime < 5 || prime > 251 || prime % 2 == 0 || minimum < 1 || minimum > 21) return 2;
        for (int d = 3; d * d <= prime; d += 2) if (prime % d == 0) return 2;
        return WordCheck(prime, minimum).run() ? 0 : 1;
    }
    if (argc != 1) return 2;
    if (!WordCheck(127,17).check_four_pair_intersections()) return 1;
    std::array<std::array<int,2>, 9> checks{{
        {127,19}, {131,17}, {137,17}, {139,17}, {149,18}, {151,18}, {157,19}, {163,19}, {167,19},
    }};
    for (auto [s,h] : checks) if (!WordCheck(s,h).run()) return 1;
}
