// Exact enumeration of ordinary simple perfect modular dimension obstructions.
// Compile with c++ -O3 -std=c++17 search.cpp -lgmpxx -lgmp -o search.
#include <gmpxx.h>
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using Z = mpz_class;
using Mask = std::uint64_t;
using Array = std::array<unsigned, 32>;

static Z sqrt_floor(const Z& value) {
    Z root;
    mpz_sqrt(root.get_mpz_t(), value.get_mpz_t());
    return root;
}
static Z gcd(const Z& a, const Z& b) {
    Z out;
    mpz_gcd(out.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
    return out;
}
static bool divides(unsigned p, const Z& x) {
    return mpz_divisible_ui_p(x.get_mpz_t(), p) != 0;
}
static bool divides(const Z& p, const Z& x) {
    return mpz_divisible_p(x.get_mpz_t(), p.get_mpz_t()) != 0;
}
static bool square(const Z& x) {
    return mpz_perfect_square_p(x.get_mpz_t()) != 0;
}

using Ternary = std::map<Z,unsigned>;

// Inductively take the least possible outgoing singleton count modulo 3.
// An outgoing singleton at z consumes two labels of denominator 2z.
static bool close_ternary(Ternary& pending, const Z& z, unsigned count,
                          const Z& next_b, unsigned slots) {
    if (!divides(3,z)) {
        auto it=pending.find(z);
        unsigned inherited=it==pending.end()?0:it->second;
        if (it!=pending.end()) pending.erase(it);
        if (count<inherited) return false;
        unsigned singletons=(count-inherited)%3;
        if (singletons) pending.emplace(2*z,2*singletons);
    }
    if (!pending.empty() && pending.begin()->first<next_b) return false;
    unsigned demand=0;
    for(const auto& item:pending) demand+=item.second;
    return demand<=slots;
}

struct Packets {
    // For each eligible prime, count labels, largest run, and whole
    // packets of (p-1)/2 equal labels, across all runs processed so far.
    Array count{}, maximum{}, packets{};
    unsigned run = 0;
};

class Search {
    int rank, q = 0, split_depth, fixed_q;
    std::vector<unsigned> primes;
    std::vector<Z> prefix, fixed;
    Mask qmask = 0;
    std::map<Z, std::vector<Z>> divisor_cache;
    Z nodes = 0, types = 0, shards = 0;

    static Mask bit(unsigned i) { return Mask(1) << i; }
    static int size(Mask m) { return __builtin_popcountll(m); }
    Mask mask_of(const Z& n) const {
        Mask m = 0;
        for (unsigned i = 0; i < primes.size(); ++i)
            if (divides(primes[i], n)) m |= bit(i);
        return m;
    }
    bool smooth_on(const Z& n, Mask allowed, Mask& mask) const {
        mask = 0;
        Z tail = n;
        while (allowed && tail > 1) {
            unsigned i = __builtin_ctzll(allowed);
            allowed &= allowed - 1;
            unsigned p = primes[i];
            if (!divides(p, tail)) continue;
            mask |= bit(i);
            do { tail /= p; } while (divides(p, tail));
        }
        return tail == 1;
    }
    template<class F>
    void smooth_range(const Z& lo, const Z& hi, Z x, Mask choices,
                      Mask used, F& accept) const {
        if (!choices) {
            if (x >= lo) accept(x, used);
            return;
        }
        unsigned i = 63 - __builtin_clzll(choices);
        choices ^= bit(i);
        unsigned p = primes[i];
        if (!choices && p == 2) {
            if (x < lo) {
                Z quotient = (lo-1)/x;
                mp_bitcnt_t shift = mpz_sizeinbase(quotient.get_mpz_t(), 2);
                x <<= shift;
                used |= bit(i);
            }
            while (x <= hi) {
                accept(x, used);
                x *= 2;
                used |= bit(i);
            }
            return;
        }
        smooth_range(lo, hi, x, choices, used, accept);
        x *= p;
        while (x <= hi) {
            smooth_range(lo, hi, x, choices, used | bit(i), accept);
            x *= p;
        }
    }
    Packets append(const Packets& old, bool equal_previous, Mask bmask) const {
        Packets st = old;
        st.run = equal_previous ? old.run+1 : 1;
        for (unsigned i = 1; i < primes.size(); ++i) {
            if (!((qmask | bmask) & bit(i))) continue;
            unsigned m = (primes[i]-1)/2;
            ++st.count[i];
            st.maximum[i] = std::max(st.maximum[i], st.run);
            st.packets[i] += st.run % m == 0;
        }
        return st;
    }
    bool can_complete(const Packets& st, Mask bmask, Mask support, unsigned h) const {
        for (unsigned i = 1; i < primes.size(); ++i) {
            if (!((qmask | support) & bit(i))) continue;
            unsigned m = (primes[i]-1)/2;
            unsigned current = ((qmask | bmask) & bit(i)) ? st.run : 0;
            unsigned max_run = std::max(st.maximum[i], h+current);
            if (max_run < m) return false;
            if (qmask & bit(i)) continue;
            if (st.count[i] + h < 2*m) return false;
            unsigned packet_cap = st.packets[i] + (h + current % m)/m;
            if (max_run < 2*m && packet_cap < 3) return false;
        }
        return true;
    }
    const std::vector<Z>& mixed_divisors(const Z& L, Mask support) {
        auto old = divisor_cache.find(L);
        if (old != divisor_cache.end()) return old->second;
        std::vector<Z> divisors{1};
        Z n = L;
        for (unsigned i = 0; i < primes.size(); ++i) {
            if (!(support & bit(i))) continue;
            unsigned p = primes[i];
            size_t len = divisors.size();
            Z power = 1;
            while (divides(p, n)) {
                n /= p;
                power *= p;
                for (size_t j = 0; j < len; ++j)
                    divisors.push_back(divisors[j]*power);
            }
        }
        divisors.erase(std::remove_if(divisors.begin(), divisors.end(),
            [&](const Z& d){return size(mask_of(d)) < 2;}), divisors.end());
        std::sort(divisors.begin(), divisors.end());
        if (divisor_cache.size() >= 8192) divisor_cache.clear();
        return divisor_cache.emplace(L, std::move(divisors)).first->second;
    }
    bool tail_covered(const Z& L, Mask possible) const {
        // Upper box at two slots. A new scale has no factors 2 or 5.
        if (size(possible) <= 2) return true;
        if (size(possible) != 3) return false;
        int position = 0;
        for (unsigned i = 0; i < primes.size(); ++i) {
            if (!(possible & bit(i))) continue;
            unsigned p = primes[i];
            if (p == 2 || p == 5) {
                int exponent = (q % p == 0);
                for (Z n=L; divides(p,n); n/=p) exponent += 2;
                if (exponent <= 5-position) return true;
            }
            ++position;
        }
        return false;
    }
    void finish_pair(const Z& L, const Z& k, const Z& d, const Z& e,
                     Mask support, const Packets& old) {
        if (gcd(k,d) != 1 || gcd(k,e) != 1) return;
        Mask km;
        if (!smooth_on(k, (qmask|support|mask_of(3)) & ~mask_of(10), km)) return;
        Z s = k*L, b = s/e, c = s/d;
        if (b < prefix.back()) return;
        Mask bm = mask_of(b), cm = mask_of(c), sm = support | km;
        Packets st = append(old, b == prefix.back(), bm);
        st = append(st, c == b, cm);
        if (!can_complete(st, cm, sm, 0)) return;
        std::vector<Z> degrees{1,d,e};
        for (auto i=prefix.rbegin(); i!=prefix.rend(); ++i) {
            Z x = s/(*i);
            if (size(mask_of(x)) < 2) return;
            degrees.push_back(x);
        }
        for (unsigned i=0; i<primes.size(); ++i) {
            if (!(sm & bit(i))) continue;
            int free = 0;
            for (size_t j=1; j<degrees.size(); ++j)
                if (!divides(primes[i],degrees[j])) ++free;
            if (free < 2) return;
        }
        ++types;
        std::cout << "[";
        for (size_t i=0; i<degrees.size(); ++i) std::cout << (i?",":"") << degrees[i];
        std::cout << "]\n";
    }
    void pair_tail(const Z& L, const Z& E, Mask support, const Packets& st) {
        Mask possible = qmask | support | mask_of(3);
        if (tail_covered(L, possible)) return;
        const auto& ds = mixed_divisors(L, support);
        if (ds.empty()) return;
        Z ll=L*L, lo = std::max(prefix.back(), Z(sqrt_floor(ll/E)+1));
        Z hi = sqrt_floor(73*ll/(36*E));
        if (lo > hi) return;
        if (hi-lo <= std::max<size_t>(32,2*ds.size())) {
            for (Z b=lo; b<=hi; ++b) {
                Mask bm;
                if (!smooth_on(b, possible, bm)) continue;
                Z g=gcd(L,b), k=b/g, e=L/g;
                if (size(mask_of(e)) < 2 || divides(2,k) || divides(5,k)) continue;
                Z d2=E*k*k-e*e-1;
                if (d2<36 || d2>e*e || !square(d2)) continue;
                Z d=sqrt_floor(d2);
                if (divides(d,L) && size(mask_of(d))>=2) finish_pair(L,k,d,e,support,st);
            }
            return;
        }
        // Square residues turn the Cartesian product of divisors into a join.
        std::vector<std::pair<Z,Z>> by_residue;
        for (const Z& e:ds) by_residue.emplace_back(e*e%E,e);
        std::sort(by_residue.begin(),by_residue.end());
        for (const Z& d:ds) {
            Z a=d*d+1, residue=(E-a%E)%E;
            auto i=std::lower_bound(by_residue.begin(),by_residue.end(),std::make_pair(residue,d));
            for (; i!=by_residue.end() && i->first==residue; ++i) {
                const Z& e=i->second;
                Z k2=(a+e*e)/E;
                if (k2>0 && square(k2)) finish_pair(L,sqrt_floor(k2),d,e,support,st);
            }
        }
    }
    void visit(const Z& L, const Z& E, Mask support, Mask bmask, const Packets& st, const Ternary& ternary) {
        ++nodes;
        // The sharper local block bounds cap the number of largest
        // labels at two for codegree eight. The rank-dependent codegree nine
        // ceilings preserve earlier enumeration certificates. At codegree nine 3 cannot
        // divide the degree M/3, and the prefix lcm divides M.
        if (rank >= 23) {
            const Z f = q * prefix.front() * prefix.front();
            bool all_maximum = prefix.front() == prefix.back();
            unsigned nine_cap = rank >= 27 ? 2 : (rank==26 ? 3 : 4);
            if (f == 9 && (divides(9,L) || (all_maximum && prefix.size() > nine_cap))) return;
            if (rank >= 26 && all_maximum) {
                unsigned f0 = f.get_ui();
                for (unsigned prime: primes)
                    if (prime > 2 && f0 % prime == 0 && f0 % (prime*prime) != 0 &&
                        prefix.size() > f0 - 1 - (prime-1)/2) return;
                if ((f0 == 16 || f0 == 25 || f0 == 27) &&
                    f0*(unsigned(rank)-f0) >= 25 && prefix.size() >= f0-1) return;
            }
            if (rank >= 24 && (f == 8 || f == 9)) {
                // Every possible new denominator prime outside qL is at most h+1.
                unsigned h = rank - 1 - static_cast<unsigned>(prefix.size());
                Mask possible = qmask | support;
                for (unsigned i = 0; i < primes.size() && primes[i] <= h+1; ++i)
                    possible |= bit(i);
                if (size(possible) < 4) return;
            }
            unsigned eight_cap = rank == 23 ? 5 : 2;
            if (f == 8 && all_maximum && prefix.size() > eight_cap) return;

        }
        // For ranks >=20, the Gauss-sum saturation lemma rules out f-1
        // equal maximal labels at f=8 or f=9.
        if (rank >= 20 && prefix.front() == prefix.back()) {
            Z f = q * prefix.front() * prefix.front();
            if ((f == 8 || f == 9) && prefix.size()+1 >= f) return;
        }
        unsigned h=rank-1-static_cast<unsigned>(prefix.size());
        if (!can_complete(st,bmask,support,h)) return;
        // At three remaining slots a new factor two in M/L has valuation
        // one and requires E odd (reduce E*k*k = 1+sum(d_i^2) mod 8).
        if (rank>=23 && (h==3 || h==4)) {
            Mask future=support|qmask|(h==3?Mask(3):Mask(7));
            if (size(future)<=2) return;
            if (h==3 && size(future)==3 && (future&1)) {
                unsigned e=(q%2==0)+2*mpz_scan1(L.get_mpz_t(),0)+(divides(2,E)?0:2);
                if(e<=5) return;
            }
        }
        if (split_depth && prefix.size()==static_cast<size_t>(split_depth)) {
            std::cout << q << "\t";
            for (size_t i=0;i<prefix.size();++i) std::cout << (i?",":"") << prefix[i];
            std::cout << "\n";
            ++shards;
            return;
        }
        if (h==2) { pair_tail(L,E,support,st); return; }
        Z ll=L*L, lo=std::max(prefix.back(),Z(sqrt_floor(ll/E)+1));
        Z hi=sqrt_floor((36*h+1)*ll/(36*E));
        if (lo>hi) return;
        Mask allowed=qmask|support;
        for (unsigned i=0;i<primes.size() && primes[i]<=h+1;++i) allowed|=bit(i);
        auto take=[&](const Z& b, Mask bm) {
            if (rank>=23 && q==1 && prefix.front()==3 && divides(9,b)) return;
            Ternary child;
            if (rank>=27 && q%3==0) {
                child=ternary;
                if(b!=prefix.back() && !close_ternary(child,prefix.back(),st.run,b,h)) return;
            }
            Z g=gcd(L,b), k=b/g, d=L/g;
            Packets next=append(st,b==prefix.back(),bm);
            prefix.push_back(b);
            visit(L*k,E*k*k-d*d,support|bm,bm,next,child);
            prefix.pop_back();
        };
        if (prefix.size()<fixed.size()) {
            const Z& b=fixed[prefix.size()];
            Mask bm;
            if (lo<=b && b<=hi && smooth_on(b,allowed,bm)) take(b,bm);
        } else if (hi-lo<=32) {
            for (Z b=lo;b<=hi;++b) {
                Mask bm;
                if (smooth_on(b,allowed,bm)) take(b,bm);
            }
        } else smooth_range(lo,hi,1,allowed,0,take);
    }

public:
    Search(int r, int depth, int q0, std::vector<Z> fixed_prefix)
        : rank(r), split_depth(depth), fixed_q(q0), fixed(std::move(fixed_prefix)) {
        // Packet counters are <= rank. 2*rank+1 <= 97 has 25 primes,
        // fitting our 32 counter slots and 64-bit masks.
        if (rank<6 || rank>48) throw std::invalid_argument("rank must lie in 6..48");
        if (depth && (depth<1 || depth>rank-3)) throw std::invalid_argument("invalid split depth");
        if (fixed.size()>static_cast<size_t>(rank-3)) throw std::invalid_argument("prefix too long");
        for (unsigned p=2;p<=unsigned(2*rank+1);++p) {
            bool prime=true;
            for (unsigned j=2;j*j<=p;++j) if (p%j==0) prime=false;
            if (prime) primes.push_back(p);
        }
    }
    void run() {
        for (q=1;q<rank;++q) {
            if (fixed_q && q!=fixed_q) continue;
            bool sqfree=true;
            for (unsigned p:primes) if (q%(p*p)==0) sqfree=false;
            if (!sqfree) continue;
            qmask=mask_of(q);
            for (unsigned b=1;q*b*b<unsigned(rank);++b) {
                if ((!fixed.empty() && fixed.front()!=b) || q*b*b<8) continue;
                if (b==1 && size(qmask)==1) continue; // prime formal codegree
                Mask bm=mask_of(b);
                prefix={Z(b)};
                visit(b, q*b*b-1, bm,bm,append(Packets{},false,bm),Ternary{});
            }
        }
        std::cerr << "COMPLETE " << (split_depth?"manifest":"search")
                  << " rank=" << rank << " nodes=" << nodes << " types=" << types
                  << " shards=" << shards << "\n";
    }
};

int main(int argc,char** argv) {
    try {
        if (argc!=3 && argc!=4 && argc!=5)
            throw std::invalid_argument("use --all R, --manifest R DEPTH, or --shard R Q b1,b2,...");
        std::string mode=argv[1];
        int rank=std::stoi(argv[2]);
        if (mode=="--all" && argc==3) Search(rank,0,0,{}).run();
        else if (mode=="--manifest" && argc==4) Search(rank,std::stoi(argv[3]),0,{}).run();
        else if (mode=="--shard" && argc==5) {
            int q=std::stoi(argv[3]);
            if (q<1 || q>=rank) throw std::invalid_argument("invalid square-free multiplier");
            std::vector<Z> fixed;
            std::stringstream stream(argv[4]);
            std::string token;
            while (std::getline(stream,token,',')) {
                Z b(token,10);
                if (b<1 || (!fixed.empty() && b<fixed.back())) throw std::invalid_argument("invalid prefix");
                fixed.push_back(b);
            }
            if (fixed.empty()) throw std::invalid_argument("empty prefix");
            Search(rank,0,q,fixed).run();
        } else throw std::invalid_argument("bad arguments");
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "INCOMPLETE " << e.what() << "\n";
        return 2;
    }
}
