// Enumerate the formal codegrees of an arbitrary integral braided category.
// Compile with g++ -O2 -std=c++17 formal_sieve.cpp -o formal_sieve.
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using Integer = std::uint64_t;
using Wide = unsigned __int128;
using Factors = std::vector<std::pair<Integer, unsigned>>;

std::vector<Integer> small_primes;
std::vector<Integer> prefix;
const Integer even_threshold = 29160000;
const Integer odd_threshold = 6087156075;

Integer checked(Wide value) {
    if (value > std::numeric_limits<Integer>::max())
        throw std::overflow_error("Increase the working precision");
    return static_cast<Integer>(value);
}

Factors factor(Integer value) {
    Factors result;
    for (Integer prime : small_primes) {
        if (prime > value / prime) break;
        unsigned exponent = 0;
        while (value % prime == 0) {
            value /= prime;
            ++exponent;
        }
        if (exponent) result.emplace_back(prime, exponent);
    }
    if (value > 1) {
        Integer last = small_primes.back();
        if (last <= value / last)
            throw std::overflow_error("Increase the prime sieve");
        result.emplace_back(value, 1);
    }
    return result;
}

bool prime_family(const Factors& left, const Factors& right) {
    std::map<Integer, unsigned> powers;
    for (const auto& [p, a] : left) powers[p] += a;
    for (const auto& [p, a] : right) powers[p] += a;
    std::vector<unsigned> sorted;
    for (const auto& [p, a] : powers) if (a >= 2) sorted.push_back(a);
    if (sorted.size() <= 2) return true;
    if (sorted.size() == 3 && (sorted[0] <= 5 || sorted[1] <= 4 || sorted[2] <= 3))
        return true;
    std::sort(sorted.begin(), sorted.end());
    if (sorted.size() == 3)
        return sorted[0] <= 4 && sorted[1] <= 5 && sorted[2] <= 7;
    if (sorted.size() == 4) {
        auto d = sorted[0], c = sorted[1], b = sorted[2], a = sorted[3];
        return b <= 3 || (d == 2 && c <= 3 &&
             (a / 2 + b / 2 <= 5 || (powers[2] >= 4 && a <= 9 && b <= 5)));
    }
    if (sorted.size() == 5)
        return sorted.back() <= 3 ||
            (sorted.back() == 4 && sorted[3] <= 3 && sorted[0] == 2);
    return false;
}

std::vector<Integer> divisors(const Factors& factors) {
    std::vector<Integer> result{1};
    for (const auto& [prime, exponent] : factors) {
        const auto count = result.size();
        Integer power = 1;
        for (unsigned i = 0; i < exponent; ++i) {
            power *= prime;
            for (std::size_t j = 0; j < count; ++j)
                result.push_back(result[j] * power);
        }
    }
    return result;
}

// If p/q=1/x+1/D, set k=D/x. Then k divides q and D=q(k+1)/p.
void last_two(Integer p, Integer q, Integer common_multiple, Integer cutoff) {
    if (Wide(q) * (Wide(q) + 1) < Wide(cutoff) * p) return;
    const bool starts_two = prefix.front() == 2;
    const bool has_three = std::find(prefix.begin(), prefix.end(), 3) != prefix.end();
    auto factors_q = factor(q);
    for (Integer k : divisors(factors_q)) {
        if ((k + 1) % p) continue;
        Integer multiplier = (k + 1) / p;
        Wide dimension = Wide(q) * multiplier;
        Wide penultimate = dimension / k;
        if (dimension < cutoff || dimension % common_multiple || penultimate < prefix.back())
            continue;
        if (starts_two && (dimension % 4 != 2 || penultimate % 2 == 0 ||
                          (has_three && penultimate == 7))) continue;
        if (prime_family(factors_q, factor(multiplier))) continue;
        std::cout << prefix.size() + 2;
        for (Integer x : prefix) std::cout << ' ' << x;
        std::cout << ' ' << checked(penultimate) << ' ' << checked(dimension) << '\n';
    }
}

void search(unsigned count, Integer p, Integer q, Integer common_multiple) {
    bool starts_two = prefix.front() == 2;
    Integer cutoff = starts_two ? 2 * odd_threshold : even_threshold;
    if (count == 2) {
        last_two(p, q, common_multiple, cutoff);
        return;
    }
    Integer first = std::max(prefix.back(), q / p + 1);
    Integer last = checked(Wide(count) * q / p);
    if (first > last) return;
    if (count == 3) {
        // The lax two-term upper bound is convex in x after
        // multiplication by the positive linear denominator p*x-q.
        auto small = [&](Integer x) {
            Integer b = checked(Wide(q) * x);
            return Wide(b) * (Wide(b) + 1) < Wide(cutoff) * (Wide(p) * x - q);
        };
        if (small(first) && small(last)) return;
    }
    bool has_three = std::find(prefix.begin(), prefix.end(), 3) != prefix.end();
    for (Integer x = first; x <= last; ++x) {
        if (starts_two && (x % 2 == 0 || (has_three && x == 7))) continue;
        Integer num = checked(Wide(p) * x - q);
        Integer den = checked(Wide(q) * x);
        Integer g = std::gcd(num, den);
        Integer next_lcm = checked(Wide(common_multiple / std::gcd(common_multiple, x)) * x);
        prefix.push_back(x);
        search(count - 1, num / g, den / g, next_lcm);
        prefix.pop_back();
    }
}

// Rank thirteen, odd dimension, first paired codegree three. The paired
// prefix 3,7,43 is handled by the separate square-denominator computation.
Integer odd_large = 0, odd_open = 0;
unsigned odd_rank = 13;
void odd_last(Integer p, Integer q, Integer common_multiple) {
    if (Wide(q) * (2 * Wide(q) + 1) < Wide(odd_threshold) * p) return;
    auto factors_q = factor(q);
    for (Integer k : divisors(factors_q)) {
        if ((2 * k + 1) % p) continue;
        Integer h = (2 * k + 1) / p;
        Wide D = Wide(q) * h;
        if (D < odd_threshold || D % common_multiple || D / k < prefix.back()) continue;
        Integer prime = prefix.front();
        if (D % (prime * prime) != prime || (D / k) % prime == 0) continue;
        ++odd_large;
        if (!prime_family(factors_q, factor(h))) {
            ++odd_open;
            std::cout << "uncovered dimension " << checked(D) << '\n';
        }
    }
}

void odd_search(unsigned pairs, Integer p, Integer q, Integer common_multiple) {
    if (pairs == 1) {
        odd_last(p, q, common_multiple);
        return;
    }
    Integer first = std::max(prefix.back(), 2 * q / p + 1) | 1;
    Integer last = checked(Wide(2 * pairs + 1) * q / p);
    if (first > last) return;
    if (pairs == 2) {
        auto small = [&](Integer x) {
            Integer b = checked(Wide(q) * x);
            return Wide(b) * (2 * Wide(b) + 1) <
                   Wide(odd_threshold) * (Wide(p) * x - 2 * q);
        };
        if (small(first) && small(last)) return;
    }
    for (Integer x = first; x <= last; x += 2) {
        Integer prime = prefix.front();
        if (x % prime == 0) continue;
        if (odd_rank == 13 && prefix.size() == 2 && prefix[1] == 7 && x == 43)
            continue;
        if (odd_rank == 15 && prime == 3 && prefix.size() == 1 && x < 17)
            continue;
        if (odd_rank == 15 && prime == 5 && prefix.size() == 2 && x < 13)
            continue;
        Integer num = checked(Wide(p) * x - 2 * q);
        Integer den = checked(Wide(q) * x);
        Integer g = std::gcd(num, den);
        Integer next_lcm = checked(Wide(common_multiple / std::gcd(common_multiple, x)) * x);
        prefix.push_back(x);
        odd_search(pairs - 1, num / g, den / g, next_lcm);
        prefix.pop_back();
    }
}

int main(int argc, char** argv) {
    std::vector<bool> composite(1000005);
    for (Integer i = 2; i < composite.size(); ++i) {
        if (composite[i]) continue;
        small_primes.push_back(i);
        for (Integer j = i * i; j < composite.size(); j += i) composite[j] = true;
    }
    if (argc == 2 && std::string(argv[1]) == "--odd13") {
        prefix = {3};
        odd_search(5, 1, 3, 3);
        std::cout << "Rank 13, paired minimum 3, outside 3,7,43: "
                  << odd_large << " above cutoff, " << odd_open << " uncovered.\n";
        return odd_open != 0;
    }
    if (argc == 2 && std::string(argv[1]) == "--odd15") {
        odd_rank = 15;
        prefix = {3};
        odd_search(6, 1, 3, 3);
        std::cout << "Rank 15, minimum 3, next at least 17: "
                  << odd_large << " above cutoff, " << odd_open << " uncovered.\n";
        if (odd_open) return 1;
        odd_large = 0;
        prefix = {5, 5};
        odd_search(5, 1, 5, 5);
        std::cout << "Rank 15, minimum 5, next after 5,5 at least 13: "
                  << odd_large << " above cutoff, " << odd_open << " uncovered.\n";
        return odd_open != 0;
    }
    if (argc != 1) throw std::invalid_argument("Use no option, --odd13, or --odd15");
    // Ranks 1 and 2 have numerical dimensions 1 and 2.
    for (unsigned rank = 3; rank <= 8; ++rank) {
        for (Integer x = 2; x <= rank; ++x) {
            prefix = {x};
            search(rank - 1, x - 1, x, x);
        }
    }
}
