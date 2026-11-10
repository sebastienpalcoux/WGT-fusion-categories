#pragma once
#include <algorithm>
#include <cstdint>
#include <vector>

// Minimum label count over each choice of phases outside the odd-valuation primes.
inline unsigned complementary_floor(std::uint64_t square_class,
                                    const std::vector<unsigned>& primes) {
    unsigned phases = 1;
    for (auto mask = square_class & ~std::uint64_t(1); mask; mask &= mask - 1)
        phases *= (primes[__builtin_ctzll(mask)] + 1) / 2;
    unsigned labels = (square_class & 1) ? 2 * phases : phases;
    for (auto mask = square_class & ~std::uint64_t(1); mask; mask &= mask - 1) {
        const unsigned choices = (primes[__builtin_ctzll(mask)] + 1) / 2;
        labels = std::max(labels, phases / choices * (choices + 1));
    }
    return labels;
}

template<class PacketState>
inline unsigned missing_joint_labels(const PacketState& state, std::uint64_t eligible,
        std::uint64_t square_class, std::uint64_t support, unsigned fiber_size,
        const std::vector<unsigned>& primes) {
    unsigned debt = 0;
    for (auto mask = support & ~square_class & ~std::uint64_t(1); mask; mask &= mask - 1) {
        const unsigned i = __builtin_ctzll(mask);
        const unsigned step = (primes[i] - 1) / 2;
        // Closed runs supply whole packets. The open run may continue in the suffix.
        unsigned supplied = state.packets[i] * step;
        if (eligible & (std::uint64_t(1) << i)) supplied += state.run % step;
        const unsigned required = 2 * step * fiber_size;
        if (required > supplied) debt += required - supplied;
    }
    return debt;
}
