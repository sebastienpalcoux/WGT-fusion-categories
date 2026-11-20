// A collision in this cache only sends an extra branch to the exact tail solver.
namespace tail_cache {
struct Entry {
    unsigned requests = 0;
    unsigned shift = 0;
    std::vector<U> bits;
    std::uint64_t hash(U energy) const {
        return (energy * UINT64_C(11400714819323198485)) >> shift;
    }
    void store(U energy) {
        auto b = hash(energy);
        bits[b >> 6] |= U(1) << (b & 63);
    }
    bool contains(U energy) const {
        auto b = hash(energy);
        return bits[b >> 6] & (U(1) << (b & 63));
    }
};

static bool possible(U L, Mask support, Wide E) {
    // With L<2^31, 1+d^2+e^2 is representable in uint64_t at all pairs.
    if (L >= (U(1) << 31)) return true;
    if (E > UINT64_MAX) return false;
    static std::unordered_map<U,Entry> cache;
    static std::size_t bytes = 0;
    if (cache.size() > 16000 || bytes > (16U << 20)) {
        cache.clear();
        bytes = 0;
    }
    auto &entry = cache[L];
    if (!entry.bits.empty()) return entry.contains(static_cast<U>(E));
    if (++entry.requests < 64) return true;
    const auto &ds = nonunit_divisors(L, support);
    if (ds.size()>384 || entry.requests < ds.size()*2) return true;
    const std::size_t pairs = ds.size()*(ds.size()+1)/2;
    unsigned log = 6;
    while ((U(1)<<log) < 16*pairs) ++log;
    entry.shift=64-log;
    entry.bits.assign(U(1)<<(log-6), 0);
    bytes += entry.bits.capacity()*sizeof(U);

    std::vector<U> squares;
    for (Mask a=(support|Mask(2)) & ~Mask(5); a; a &= a-1) {
        const U p=primes[__builtin_ctzll(a)];
        squares.push_back(p*p);
    }
    for (std::size_t j=0; j<ds.size(); ++j) {
        const U offset=1+ds[j]*ds[j];
        for (std::size_t i=0; i<=j; ++i) {
            U target=offset+ds[i]*ds[i];
            // Mark all N/k² with k supported on the admissible odd primes.
            // Omitting gcd(k,d e)=1 here is only a relaxation.
            auto mark=[&](auto &&self, U value, std::size_t first)->void {
                entry.store(value);
                for (std::size_t t=first; t<squares.size(); ++t)
                    if (value%squares[t]==0) self(self, value/squares[t], t);
            };
            mark(mark,target,0);
        }
    }
    return entry.contains(static_cast<U>(E));
}
}
