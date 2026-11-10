// Squared Egyptian fractions for ordinary simple perfect integral modular types.
// Compile: c++ -O3 -std=c++17 ordinary_native.cpp -o search
// Run: ./search R 8 [Q [fixed,initial,denominators]] > candidates.jsonl
// Candidates are sorted internally by degree. The output order across types is
// traversal order. All range errors stop the run without a completion record.
// Compile with -DUNPRUNED_AUDIT to disable the Gauss and early family tests;
// the remaining filters then match the rank-20 GMP search in this archive.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <string>
#include <numeric>
#include <unordered_map>
#include <stdexcept>
#include <vector>
using U = std::uint64_t;
using Wide = __uint128_t;
using Mask = std::uint64_t;
using Counts = std::array<unsigned char,32>;
static std::vector<unsigned> primes;
static std::vector<U> prefix;
static U nodes, types;
static std::vector<U> fixed;
#include "joint_capacity.hpp"
static unsigned fiber_size = 1;
struct PacketState {
    Counts maximum{}, count{}, packets{}, doubles{};
    int run=0;
};
static PacketState append(const PacketState& old,bool equal,Mask eligible) {
    PacketState st=old;
    st.run=equal?old.run+1:1;
    for(Mask a=eligible & ~Mask(1);a;a&=a-1) {
        unsigned i=__builtin_ctzll(a),m=(primes[i]-1)/2;
        ++st.count[i]; st.maximum[i]=std::max<int>(st.maximum[i],st.run);
        st.packets[i]+=(st.run%m==0);
        st.doubles[i]+=(st.run%(2*m)==0);
    }
    return st;
}
static bool possible(const PacketState& st,Mask eligible,Mask qmask,Mask support,unsigned h) {
    if (missing_joint_labels(st,eligible,qmask,support,fiber_size,primes) > h) return false;
    for(Mask a=(support|qmask) & ~Mask(1);a;a&=a-1) {
        unsigned i=__builtin_ctzll(a), m=(primes[i]-1)/2;
        unsigned current=(eligible&(Mask(1)<<i))?st.run:0;
        unsigned mx=std::max<unsigned>(st.maximum[i],h+current);
        if(mx<m) return false;
        if((qmask&1) && !(qmask&(Mask(1)<<i))) {
            unsigned cap=st.packets[i]+(h+current%m)/m;
            unsigned dcap=st.doubles[i]+(h+current%(2*m))/(2*m);
            if(cap+dcap<6) return false;
        }
        // A nontrivial odd-primary phase uses at least three labels with
        // odd ternary valuation, or two with odd dyadic valuation.
        if (!(qmask&(Mask(1)<<i)) && (qmask&3)) {
            unsigned required=(qmask&2)?6:4;
            if(st.count[i]+h < required*m ||
               st.packets[i]+(h+current%m)/m < required) return false;
        }
        if(!(qmask&(Mask(1)<<i)) && (st.count[i]+h<2*m ||
           (mx<2*m && st.packets[i]+(h+current%m)/m<3))) return false;
    }
    return true;
}

// For ternary defect-one blocks, each run at b consists of triples,
 // inherited pairs from b/2, and singleton labels demanding two at 2b.
struct Ternary {
    std::array<U,40> b{};
    std::array<unsigned char,40> quantity{};
    unsigned first=0, last=0;
};
static bool end_ternary_run(Ternary& t, U z, unsigned count, U next_b, int h) {
    if (z%3) {
        unsigned inherited=0;
        if (t.first<t.last && t.b[t.first]==z) inherited=t.quantity[t.first++];
        if (count<inherited) return false;
        unsigned singles=(count-inherited)%3;
        if (singles) {
            if (z>UINT64_MAX/2) throw std::overflow_error("double denominator overflow");
            t.b[t.last]=2*z;
            t.quantity[t.last++]=2*singles;
        }
    }
    if (t.first<t.last && t.b[t.first]<next_b) return false;
    unsigned needed=0;
    for(unsigned i=t.first;i<t.last;++i) needed+=t.quantity[i];
    return needed<=unsigned(h);
}
// Odd-multiplicity defect-one runs need distinct ramified partner
// orbits, unless the final incomplete run acquires another label.
static bool dyadic_pairs(int h) {
    unsigned odd_runs=0, partners=0, potential=0;
    for(unsigned i=0;i<prefix.size();) {
        unsigned j=i+1;
        while(j<prefix.size() && prefix[i]==prefix[j]) ++j;
        unsigned n=j-i;
        U b=prefix[i]; unsigned height=__builtin_ctzll(b);
        if(height==0) {
            if(j<prefix.size() || !h) odd_runs+=n%2;
        } else if(height<=2) {
            unsigned size=1u<<(2*height);
            partners+=n/size;
            if(j==prefix.size() && n%size && n%size+unsigned(h)>=size) {
                unsigned spent=size-n%size;
                potential=1+(h-spent)/4;
            }
        }
        i=j;
    }
    return odd_runs<=partners+std::max<unsigned>(potential,h/4);
}
static int rank_limit, minimum_codegree = 8;
static U square_root(Wide x) {
    U y = static_cast<U>(std::sqrt(static_cast<long double>(x)));
    while (Wide(y)*y > x) --y;
    while (y < UINT64_MAX && Wide(y+1)*(y+1) <= x) ++y;
    return y;
}
static Mask prime_mask(Wide x) {
    Mask mask = 0;
    for (unsigned i=0; i<primes.size(); ++i) if (x%primes[i]==0)
        mask |= Mask(1)<<i;
    return mask;
}
static bool sieve(U x, Mask allowed, Mask &mask) {
    mask=0;
    while (allowed && x>1) {
        int i=__builtin_ctzll(allowed);
        U p=primes[i];
        allowed &= allowed-1;
        if(x%p) continue;
        mask |= Mask(1)<<i;
        do { x/=p; } while(x%p==0);
    }
    return x==1;
}
template<class Call>
static void smooth(U low,U high,U product,Mask active,Mask used,Call &call) {
    if(!active) {
        if(product>=low) call(product,used);
        return;
    }
    int i=63-__builtin_clzll(active);
    U p=primes[i];
    active ^= Mask(1)<<i;
    smooth(low, high, product,active,used,call);
    while(product<=high/p) {
        product*=p;
        smooth(low, high,product,active,used|(Mask(1)<<i),call);
    }
}
// Smooth-number tables are shared by all nodes with the same prime mask.
struct SmoothTable {
    U bound=0;
    std::vector<std::pair<U,Mask>> items;
};
static std::unordered_map<Mask,SmoothTable> smooth_cache;
template<class F>
static void enumerate_smooth(U lo,U hi,Mask allowed,F& callback) {
    if(hi-lo<=20) {
        for(U b=lo;b<=hi;++b) {Mask mask; if(sieve(b,allowed,mask)) callback(b,mask);}
        return;
    }
    // Restrict allocation on larger supports or very large ranges.
    if(__builtin_popcountll(allowed)>6 || hi>10000000000ULL) {
        smooth(lo,hi,1,allowed,0,callback);
        return;
    }
    auto& tab=smooth_cache[allowed];
    if(hi>tab.bound) {
        U next=std::max<U>(16384,tab.bound);
        while(next<hi) next*=4;
        std::vector<std::pair<U,Mask>> entries;
        auto collect=[&](U v, Mask mask){ entries.emplace_back(v,mask); };
        smooth(1,next,1,allowed,0,collect);
        std::sort(entries.begin(),entries.end());
        tab.items.swap(entries); tab.bound=next;
    }
    auto it=std::lower_bound(tab.items.begin(),tab.items.end(),std::make_pair(lo,Mask(0)));
    size_t pos=static_cast<size_t>(it-tab.items.begin());
    while(pos<tab.items.size() && tab.items[pos].first<=hi) {
        auto [value,mask]=tab.items[pos++];
        callback(value,mask);
    }
}
static void pair_tail(int q, Mask qmask,U L,Mask support,Wide E,
                      const PacketState &st) {
    Mask future=support|qmask|prime_mask(3);
    if(__builtin_popcountll(future)<3) return;
    // At the two-tail only the 2 and 5 heights are frozen.
    if(__builtin_popcountll(future)==3) {
        int pos=0;
        for(unsigned i=0;i<primes.size();++i) if(future&(Mask(1)<<i)) {
            U p=primes[i];
            if(p==2 || p==5) {
                int e=(q%p==0);
                U n=L;
                while(n%p==0){n/=p;e+=2;}
                if(e<=5-pos) return;
            }
            ++pos;
        }
    }
    // The terminal scale k is odd. Thus E is congruent modulo eight
    // to 1 plus two squares. Values 0, 4 and 7 are impossible.
    unsigned residue8=static_cast<unsigned>(E & 7);
    if(residue8==0 || residue8==4 || residue8==7) return;
    auto emit=[&](U k,U d,U e) {
            Wide root=Wide(k)*L;
            if(root/e < prefix.back()) return;
            if(std::gcd(k,d)!=1 || std::gcd(k,e)!=1) return;
            Wide b=root/e,c=root/d;
            Mask km; if(!sieve(k,future,km) || !(k%2) || !(k%5))return;
            Mask bm=prime_mask(b),cm=prime_mask(c),real_support=support|km;
            auto p1=append(st,b==prefix.back(),bm|qmask);
            auto p2=append(p1,c==b,cm|qmask);
            if(!possible(p2,cm|qmask,qmask,real_support,0)) return;
            {
                std::vector<Wide> degrees={1,d,e};
                for(auto jt=prefix.rbegin();jt!=prefix.rend();++jt) {
                    Wide x=root/(*jt);
                    if(__builtin_popcountll(prime_mask(x))<2) return;
                    degrees.push_back(x);
                }
                for(unsigned i=0;i<primes.size();++i) if(real_support&(Mask(1)<<i)){
                    int free=0;
                    int positive=0;
                    for(unsigned j=1;j<degrees.size();++j) {
                        if(degrees[j]%primes[i]) ++free;
                        if((root/degrees[j])%primes[i]==0) ++positive;
                    }
                    if(i>0 && !(qmask&(Mask(1)<<i)) && positive<int(primes[i])-1) return;
                    if(free<2) return;
                }
                if (++types == 0) throw std::overflow_error("type count exceeds 64 bits");
                std::cout<<"[";
                for(unsigned i=0;i<degrees.size();++i) {
                    Wide x=degrees[i];std::string digits;
                    do {digits.push_back(char(48+x%10)); x/=10;} while(x);
                    std::reverse(digits.begin(),digits.end()); std::cout<<(i?",":"")<<digits;
                }
                std::cout<<"]\n";
            }
    };
    if(L>(U(1)<<60)) throw std::overflow_error("tail interval exceeds bit bound");
    Wide ll=Wide(L)*L;
    U lo=std::max(prefix.back(),square_root(ll/E)+1);
    U hi=square_root((73*ll)/(36*E));
    if(lo>hi) return;
    // Enumerate the admissible values of b=M/e directly. The interval
    // is narrow multiplicatively and b is smooth on the known support.
    // This avoids building and sorting a divisor-residue table per node.
    const Mask allowed=support|qmask|Mask(2);
    auto test=[&](U b, Mask){
        U g=std::gcd(L,b), k=b/g, e=L/g;
        if((k&1)==0 || k%5==0 || e<6) return;
        Wide ee=Wide(e)*e;
        Wide mass=E*k*k; // <=(73/36)e^2 by the lower bound on b
        if(mass<ee+37 || mass>2*ee+1) return;
        Wide d2=mass-ee-1;
        U residue=d2&63;
        static constexpr U squares=0x202021202030213ULL;
        if(!(squares & (U(1)<<residue))) return;
        U d=square_root(d2);
        if(Wide(d)*d==d2 && L%d==0 &&
           __builtin_popcountll(prime_mask(d))>=2 && __builtin_popcountll(prime_mask(e))>=2)
            emit(k,d,e);
    };
    enumerate_smooth(lo,hi,allowed,test);

}
static void dfs(int q,Mask qmask,U L,Mask support,Wide E,int h,
                const PacketState &st,Mask bmask, const Ternary& tern) {
    if (++nodes == 0) throw std::overflow_error("node count exceeds 64 bits");
    if(!possible(st,bmask|qmask,qmask,support,h)) return;
    if(q%2==0 && !dyadic_pairs(h)) return;
#ifndef UNPRUNED_AUDIT
    // At f=8 or 9, f-1 copies of the maximal degree are impossible
    // at rank at least 20 by the Gauss-sum multiplicity lemma.
    U f=q*prefix.front()*prefix.front();
    // Keep the historical rank-21/22 certificates reproducible. At rank 23
    // the certified run permits five maxima of codegree eight, then checks
    // parity. Subsequent ranks apply the final two-label cap and four-prime condition.
    if (rank_limit >= 23) {
        bool all_maximum = prefix.front() == prefix.back();
        if (rank_limit >= 26 && all_maximum) {
            for (unsigned p:primes)
                if (p>2 && f%p==0 && f%(p*p)!=0 &&
                    prefix.size() > f-1-(p-1)/2) return;
            if ((f==16 || f==25 || f==27) && f*(rank_limit-f)>=25 &&
                prefix.size() >= f-1) return;
        }
        unsigned nine_cap = rank_limit >= 27 ? 2 : (rank_limit==26 ? 3 : 4);
        if (f == 9 && (L % 9 == 0 || (all_maximum && prefix.size() > nine_cap))) return;
        if (rank_limit >= 24 && (f == 8 || f == 9)) {
            Mask future = qmask | support;
            for (unsigned i=0; i<primes.size() && primes[i] <= unsigned(h+1); ++i)
                future |= Mask(1) << i;
            if (__builtin_popcountll(future) < 4) return;
        }
        unsigned eight_cap = rank_limit == 23 ? 5 : 2;
        if (f == 8 && all_maximum && prefix.size() > eight_cap) return;
    }
    if(rank_limit>=20 && (f==8 || f==9) && prefix.size()==f-1 && prefix.front()==prefix.back()) return;
    // At three slots, 2|k for k=M/L forces three odd tail degrees,
    // so E k^2 = 1+sum d^2 is 4 mod 8. Consequently v_2(k)=1 and E
    // is odd. Bound the dyadic exponent before using the BDP family.
    // The low bits in these masks correspond to the primes 2, 3 and 5.
    if(h==3 || h==4) {
        Mask future=support|qmask|(h==3?Mask(3):Mask(7));
        int num=__builtin_popcountll(future);
        if(num<=2) return;
        if(h==3 && num==3 && (future&1)) {
            int e=(q%2==0)+2*__builtin_ctzll(L)+((E&1)?2:0);
            if(e<=5) return;
        }
    }
#endif
    // With L <= 2^58 all interval numerators fit in 128 bits.
    // The residual feasibility check also bounds the denominator 36E.
    if (L>(U(1)<<58)) throw std::overflow_error("native bound, split or use GMP");
    Wide ll=Wide(L)*L;
    if(E>ll*(36*h+1)/36) return;
    if(h==2) {pair_tail(q,qmask,L,support,E,st); return;}
    U low=std::max(prefix.back(),square_root(ll/E)+1);
    U high=square_root((36*h+1)*ll/(Wide(36)*E));
    if(low>high) return;
    Mask allowed=qmask|support|Mask(1);
    const unsigned debt=missing_joint_labels(st,bmask|qmask,qmask,support,fiber_size,primes);
    const unsigned last_prime=1+(h-debt)/fiber_size;
    for(unsigned i=0;i<primes.size() && primes[i]<=last_prime;++i)
        allowed|=Mask(1)<<i;
    auto next=[&](U b,Mask bmask) {
        Ternary child=tern;
        if(rank_limit >= 27 && q%3==0 && b!=prefix.back() &&
           !end_ternary_run(child,prefix.back(),st.run,b,h)) return;
        if(rank_limit>=23 && q==1 && prefix.front()==3 && b%9==0) return;
        U g=std::gcd(b,L), k=b/g, d=L/g;
        if(b==0 || L>(UINT64_MAX/k)) throw std::overflow_error("prefix LCM too large");
        if(k>1 && E > (~Wide(0))/k/k) throw std::overflow_error("residual square too large");
        Wide Enew=Wide(E)*k*k-Wide(d)*d;

        auto next_st=append(st,b==prefix.back(),qmask|bmask);
        prefix.push_back(b);
        dfs(q,qmask,L*k,support|bmask,Enew,h-1,next_st,bmask,child);
        prefix.pop_back();
    };
    if(prefix.size()<fixed.size()) {
        U b=fixed[prefix.size()];Mask mask;
        if(low<=b && b<=high && sieve(b,allowed,mask)) next(b,mask);
    } else enumerate_smooth(low,high,allowed,next);
}
int main(int argc,char** argv) {
    if(argc<2 || argc>5) return 2;
    rank_limit=std::stoi(argv[1]);
    if(argc>=3) minimum_codegree=std::stoi(argv[2]);
    if(rank_limit<6 || rank_limit>40 || minimum_codegree<5 || minimum_codegree>8) return 2;
    for(int p=2;p<=2*rank_limit+1;++p) {
        bool prime=true; for(int i=2;i*i<=p;++i) if(p%i==0) prime=false;
        if(prime) primes.push_back(p);
    }
    if(argc>=4 && (std::stoi(argv[3])<1 || std::stoi(argv[3])>=rank_limit)) return 2;
    if(argc==5) {
        std::stringstream input(argv[4]);std::string token;
        while(std::getline(input,token,',')) fixed.push_back(std::stoull(token));
        if(fixed.empty() || fixed.front()==0 || !std::is_sorted(fixed.begin(),fixed.end()) || fixed.size()>unsigned(rank_limit-3)) return 2;
    }
    auto start=std::chrono::steady_clock::now();
    for(int q=1;q<rank_limit;++q) {
        if(argc>=4 && q!=std::stoi(argv[3])) continue;
        bool squarefree=true;
        for(auto p:primes) if(q%(p*p)==0) squarefree=false;
        if(!squarefree) continue;
        Mask qm=prime_mask(q);
        fiber_size=complementary_floor(qm,primes);
        for(int b=1;q*b*b<rank_limit;++b) {
            if(q*b*b<minimum_codegree || (b==1 && __builtin_popcountll(qm)==1)) continue;
            if(rank_limit<4) continue;
            if(!fixed.empty() && b!=int(fixed[0])) continue;
            Mask bm=prime_mask(b);
            auto pkt=append(PacketState{},false,qm|bm);
            prefix={static_cast<U>(b)};
            dfs(q,qm,b,bm,q*b*b-1,rank_limit-2,pkt,bm,Ternary{});
        }
        std::cerr<<"q="<<q<<" cumulative nodes "<<nodes<<" types "<<types<<"\n";
    }
    std::cerr<<"Complete rank "<<rank_limit<<", nodes "<<nodes<<", types "<<types
             <<", seconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"\n";
}
