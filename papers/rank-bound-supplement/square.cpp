// Square-dimension search in ranks at most 30. C++17 with unsigned __int128.
// Usage: square R [number-of-fours]; the optional 0,...,13 shard counts are disjoint.
// Rows go to standard output. A successful exit certifies coverage of the requested run.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <numeric>
#include <unordered_map>
#include <stdexcept>
#include <vector>
using U = std::uint64_t;
using Wide = __uint128_t;
using Mask = std::uint64_t;
using Counts = std::array<unsigned char,10>;
static const std::array<unsigned,10> primes={2,3,5,7,11,13,17,19,23,29};
static std::vector<U> prefix;
static U nodes, types;
static int four_shard;
struct PacketState {
    Counts maximum{}, count{}, packets{};
    int run=0;
};
static PacketState append(const PacketState& old,bool equal,Mask eligible) {
    PacketState st=old;
    st.run=equal?old.run+1:1;
    for(Mask a=eligible & ~Mask(1);a;a&=a-1) {
        unsigned i=__builtin_ctzll(a),m=(primes[i]-1)/2;
        ++st.count[i]; st.maximum[i]=std::max<int>(st.maximum[i],st.run);
        st.packets[i]+=(st.run%m==0);
    }
    return st;
}
// Each odd p supplies two packets of (p-1)/2 labels of constant degree.
// With just two packets the two degrees agree.
static bool possible(const PacketState& st, Mask eligible, Mask support, int h) {
    for(Mask a=support & ~Mask(1); a; a&=a-1) {
        unsigned i=__builtin_ctzll(a), m=(primes[i]-1)/2;
        unsigned current=(eligible & (Mask(1)<<i)) ? st.run : 0;
        unsigned maximum=std::max<unsigned>(st.maximum[i], h+current);
        if (st.count[i]+h < 2*m || (maximum<2*m &&
            st.packets[i]+(h+current%m)/m<3)) return false;
    }
    return true;
}
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
static std::unordered_map<U,std::vector<U>> divisor_cache;
static const std::vector<U> &nonunit_divisors(U L,Mask support) {
    auto it=divisor_cache.find(L);
    if(it!=divisor_cache.end()) return it->second;
    std::vector<U> ds={1};
    U n=L;
    for(unsigned i=0;i<primes.size();++i) if(support&(Mask(1)<<i)) {
        U p=primes[i], power=1;
        unsigned len=ds.size();
        while(n%p==0) {
            n/=p; power*=p;
            for(unsigned j=0;j<len;++j) ds.push_back(ds[j]*power);
        }
    }
    ds.erase(std::remove_if(ds.begin(),ds.end(),[](U d) {
        return __builtin_popcountll(prime_mask(d))<2;
    }),ds.end());
    std::sort(ds.begin(),ds.end());
    if(divisor_cache.size()>15000) divisor_cache.clear();
    return divisor_cache.emplace(L,std::move(ds)).first->second;
}
static void pair_tail(U L,Mask support,Wide E,
                      const PacketState &st) {
    Mask future=support|Mask(2); // Three is the only possible new odd prime.

    auto emit=[&](U k,U d,U e) {
        Wide root=Wide(k)*L;
        if(root/e < prefix.back()) return;
        if(std::gcd(k,d)!=1 || std::gcd(k,e)!=1) return;
        Wide b=root/e,c=root/d;
        Mask km; if(!sieve(k,future,km) || !(k%2) || !(k%5))return;
        Mask bm=prime_mask(b),cm=prime_mask(c),real_support=support|km;
        auto p1=append(st,b==prefix.back(),bm);
        auto p2=append(p1,c==b,cm);
        if(!possible(p2,cm,real_support,0)) return;
        std::vector<Wide> degrees={1,d,e};
        for(auto jt=prefix.rbegin();jt!=prefix.rend();++jt) {
            Wide x=root/(*jt);
            if(__builtin_popcountll(prime_mask(x))<2) return;
            degrees.push_back(x);
        }
        for(unsigned i=0;i<primes.size();++i) if(real_support&(Mask(1)<<i)){
            int free=0;
            for(unsigned j=1;j<degrees.size();++j) {
                if(degrees[j]%primes[i]) ++free;
            }
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
    };
    Wide ll=Wide(L)*L;
    U lo=std::max(prefix.back(),square_root(ll/E)+1);
    U hi=square_root((73*ll)/(36*E));
    if(lo>hi) return;
    const auto &ds=nonunit_divisors(L,support);
    if(ds.empty()) return;
    if(hi-lo <= std::max<U>(32,ds.size()*2)) {
        for(U b=lo;b<=hi;++b) {
            Mask mask;
            if(!sieve(b,future,mask)) continue;
            U g=std::gcd(L,b),k=b/g,e=L/g;
            if(__builtin_popcountll(prime_mask(e))<2 || k%2==0 || k%5==0) continue;
            // E*k^2 <= 73 L^2 k^2/(36 b^2) = (73/36)e^2.
            Wide mass=E*k*k;
            if(mass<=Wide(e)*e+1) continue;
            Wide d2=mass-Wide(e)*e-1;
            U d=square_root(d2);
            if(Wide(d)*d==d2 && d<=e && d>=6 && L%d==0 && __builtin_popcountll(prime_mask(d))>=2)
                emit(k,d,e);
        }
    } else {
        std::vector<std::pair<Wide,U>> keyed;
        for(U e:ds) keyed.emplace_back(Wide(e)*e%E,e);
        std::sort(keyed.begin(),keyed.end());
        for(U d:ds) {
            Wide a=Wide(d)*d+1, rem=(E-a%E)%E;
            auto it=std::lower_bound(keyed.begin(),keyed.end(),std::make_pair(rem,d));
            for(;it!=keyed.end() && it->first==rem;++it) {
                U e=it->second;
                Wide k2=(a+Wide(e)*e)/E;
                if(!k2) continue;
                U k=square_root(k2);
                if(Wide(k)*k==k2) emit(k,d,e);
            }
        }
    }
}
static void dfs(U L, Mask support, Wide E, int h, const PacketState& st) {
    int fours=std::upper_bound(prefix.begin(),prefix.end(),4)-prefix.begin();
    if (fours>13 || (four_shard>=0 && (fours>four_shard ||
        (prefix.back()>4 && fours!=four_shard)))) return;
    // At thirteen fours and four fives a seven contradicts the Gauss sum.
    if (fours==13 && prefix.size()>=18 && prefix[16]==5 && prefix[17]==7) return;
    if (++nodes==0) throw std::overflow_error("node counter");
    // L<=2^58 bounds every product used to calculate a denominator interval.
    Wide ll=Wide(L)*L;
    if (E>ll*(36*h+1)/36) return;
    Mask future=support;
    for(unsigned j=0;j<primes.size() && primes[j]<=unsigned(h+1);++j) future|=Mask(1)<<j;
    if(__builtin_popcountll(future)<3) return;
    U high=square_root((36*h+1)*ll/(Wide(36)*E));
    if(__builtin_popcountll(future)==3) {
        bool fives=std::binary_search(prefix.begin(),prefix.end(),5);
        if(h==2 && ( (fours>=4 && L%8) || (fives && L%25) )) return;
        if(h==3 && fours>=4 && L%8 && ((L%4) || E%2==0)) return;
        // Every three-prime nonunit can be free of at most one support prime.
        // At two and five at least three nonunits are p-free, by reducing D modulo 4 or 5.
        unsigned slots=0;
        Wide reserved=0;
        auto cost=[&](U least) {
            Wide square=Wide(least)*least;
            return (ll+square-1)/square;
        };
        for(Mask a=future;a;a&=a-1) {
            int j=__builtin_ctzll(a); U p=primes[j];
            U power=1; for(U n=L;n%p==0;n/=p) power*=p;
            power=std::max(power,p);
            if(p==2 && fours>=4) power=std::max<U>(power,8);
            if(p==5 && fives) power=std::max<U>(power,25);
            int need=(p==2 || p==5)?3:2;
            int missing4=(fours>=4), missing5=fives;
            for(U x:prefix) if(x%power==0) {
                if(need>0) --need;
                if(x%4==0) missing4=0;
                if(x%5==0) missing5=0;
            }
            need=std::max(need,missing4|missing5);
            slots+=need;
            if(slots>unsigned(h)) return;
            if (need) {
                U largest=power;
                if(missing4) largest=std::max(largest,std::lcm(power,U(4)));
                if(missing5) largest=std::max(largest,std::lcm(power,U(5)));
                reserved+=(need-1)*cost(power)+cost(largest);
            }
        }
        // These slots have combined mass at most reserved in units of 1/L^2.
        // The vacuum contributes at most one more in the same units.
        if (E>reserved+1) {
            high=std::min(high,square_root((h-slots)*ll/(E-reserved-1)));
            if (high<prefix.back()) return;
        }
    }
    if(h==2) {pair_tail(L,support,E,st); return;}
    U low=std::max(prefix.back(),square_root(ll/E)+1);
    if(low>high) return;
    auto next=[&](U b,Mask mask) {
        U g=std::gcd(b,L), k=b/g, d=L/g;
        auto state=append(st,b==prefix.back(),mask);
        if (!possible(state,mask,support|mask,h-1)) return;
        if (L > (U(1)<<58)/k) throw std::overflow_error("prefix lcm exceeds 2^58");
        Wide next_E=E*k*k-Wide(d)*d;
        prefix.push_back(b);
        dfs(L*k,support|mask,next_E,h-1,state);
        prefix.pop_back();
    };
    if(high-low<=35) {
        for(U b=low;b<=high;++b) {Mask mask; if(sieve(b,future,mask)) next(b,mask);}
    } else smooth(low,high,1,future,0,next);
}
int main(int argc,char** argv) {
    if (argc!=2 && argc!=3) return 2;
    int rank_limit=std::stoi(argv[1]);
    four_shard=argc==3 ? std::stoi(argv[2]) : -1;
    if (rank_limit<2 || rank_limit>30 || four_shard < -1 || four_shard>13) return 2;
    for (U b: {4,5}) if (b*b<unsigned(rank_limit)) {
        Mask mask=prime_mask(b);
        prefix={b};
        dfs(b,mask,b*b-1,rank_limit-2,append(PacketState{},false,mask));
    }
    std::cout.flush();
    if (!std::cout) return 3;
    std::cerr<<"rank "<<rank_limit<<", fours "<<four_shard
             <<", nodes "<<nodes<<", profiles "<<types<<"\n";
}
