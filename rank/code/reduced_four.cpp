// Enumerate positive reduced degree types with an exact two-nonunit completion.
// Compile: c++ -O3 -std=c++17 reduced_four.cpp -o reduced_search
// Run: ./reduced_search 16 > four16.jsonl. A range error aborts the run.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
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
static int rank_limit, minimum_codegree = 3;
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
static void pair_tail(int q, Mask qmask,U L,Mask support,Wide E,
                      const Counts &packet, int current_run) {
    Mask future=support|qmask|prime_mask(3);
    if(__builtin_popcountll(future)<3) return;
    // At the two-tail only the 2 and 5 heights are frozen.
    if(__builtin_popcountll(future)==3) {
        int pos=0;
        for(unsigned i=0;i<primes.size();++i) if(future&(Mask(1)<<i)) {
            U p=primes[i];
            if(p==2 || p==5) {
                int e=(q%p==0)+(p==2);
                U n=L;
                while(n%p==0){n/=p;e+=2;}
                if(e<=5-pos) return;
            }
            ++pos;
        }
    }
    const auto &ds=nonunit_divisors(L,support);
    if(ds.empty()) return;
    auto emit=[&](U k,U d,U e) {
            Wide root=Wide(k)*L;
            if(root/e < prefix.back()) return;
            if(std::gcd(k,d)!=1 || std::gcd(k,e)!=1) return;
            Wide b=root/e,c=root/d;
            int mb=b==prefix.back()?current_run+1:1;
            int mc=b==c?mb+1:1;
            Mask bm=prime_mask(b), cm=prime_mask(c);
            Mask real_support=support|prime_mask(k);
            for(unsigned i=2;i<primes.size();++i) {
                Mask bit=Mask(1)<<i;
                if(!((qmask|real_support)&bit)) continue;
                int need=(primes[i]-1)/2;
                if(packet[i]<need && (!((qmask|bm)&bit)||mb<need) &&
                   (!((qmask|cm)&bit)||mc<need)) return;
            }
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
                if ((q%2) && (root%2)) return;
                ++types;
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
    if(hi-lo <= std::max<U>(32,ds.size()*2)) {
        Mask allowed=support|qmask|prime_mask(3);
        for(U b=lo;b<=hi;++b) {
            Mask mask;
            if(!sieve(b,allowed,mask)) continue;
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
// A row of codegree four is a signed vector on three distinct
// pairwise noncoprime nonunit degrees, all proper divisors of its degree.
// Its support lies in the principal 2-block, so at even x every hit degree
// has smaller two-height than x.
static int four_tail(int q, U L, U low, int slots) {
    if(q!=1 || prefix.front()!=2) return 0;
    int m=std::count(prefix.begin(),prefix.end(),U(2));
    int union_size=3*m-m*(m-1)/2;
    if(slots>=union_size) return 0;
    std::vector<U> candidates;
    U modulus=L%4==0?4:2;
    for(U b:prefix) if(b>2 && b%modulus==0) candidates.push_back(b);
    int available=0;
    for(U b:candidates) {
        bool usable=Wide(6)*low<=Wide(13)*b;
        if(!usable) for(U c:candidates) {
            U diff = b>c?b-c:c-b;
            if(diff && Wide(6)*low*diff<=Wide(7)*b*c) {usable=true;break;}
        }
        if(usable) ++available;
    }
    if(available+slots>=union_size) return 0;
    // All three entries of at least one row must already occur: their signed
    // sum is +/- 1/M, hence M=L and their L-scaled sum is +/-1.
    for(unsigned i=0;i<candidates.size();++i)
      for(unsigned j=i+1;j<candidates.size();++j)
       for(unsigned k=j+1;k<candidates.size();++k) {
        U a=candidates[i],b=candidates[j],c=candidates[k];
        if(a==b || b==c) continue;
        U A=L/a,B=L/b,C=L/c;
        if(A+1==B+C || A==B+C+1) return 1;
       }
    return -1;
}
static void dfs(int q,Mask qmask,U L,Mask support,Wide E,int h,
                const Counts &packet,int run) {
    ++nodes;

    if(L>(U(1)<<58)) throw std::overflow_error("Egyptian search interval too large");
    Wide ll=Wide(L)*L;
    U low=std::max(prefix.back(),square_root(ll/E)+1);
    U high=square_root((36*h+1)*ll/(Wide(36)*E));
    if(low>high) return;
    int completion=four_tail(q,L,low,h);
    if(completion<0) return;
    if(h==2) {pair_tail(q,qmask,L,support,E,packet,run); return;}
    Mask allowed=qmask|support;
    for(unsigned i=0;i<primes.size() && primes[i]<=unsigned(h+1);++i)
        allowed|=Mask(1)<<i;
    auto next=[&](U b,Mask bmask) {
        if(completion && L%b) return;
        U g=std::gcd(b,L), k=b/g, d=L/g;
        if(b==0 || L>(UINT64_MAX/k)) throw std::overflow_error("prefix LCM too large");
        Wide Enew=Wide(E)*k*k-Wide(d)*d;
        if(k>1 && E > (~Wide(0))/k/k) throw std::overflow_error("residual square too large");
        
        Counts next_counts=packet;
        int m=b==prefix.back()?run+1:1;
        Mask current=support|bmask;
        for(unsigned i=2;i<primes.size();++i) {
            Mask bit=Mask(1)<<i;
            if(!((qmask|current)&bit)) continue;
            if(!(qmask&bit) && int(primes[i])-1>h) {
                int occupied=0;
                for(U den:prefix) if(den%primes[i]==0) ++occupied;
                if(occupied+(b%primes[i]==0)+(h-1)<int(primes[i])-1) return;
            }
            int need=(primes[i]-1)/2;
            bool eligible=(bit&(qmask|bmask))!=0;
            if(eligible) next_counts[i]=std::max<int>(packet[i],m);
            if(next_counts[i]<need && h-1+(eligible?m:0)<need) return;
        }
        prefix.push_back(b);
        dfs(q,qmask,L*k,current,Enew,h-1,next_counts,m);
        prefix.pop_back();
    };
    if(high-low<=35) {
        for(U b=low;b<=high;++b) {Mask mask; if(sieve(b,allowed,mask)) next(b,mask);}
    } else smooth(low,high,1,allowed,0,next);
}
int main(int argc,char** argv) {
    if(argc<2 || argc>4) return 2;
    rank_limit=std::stoi(argv[1]);
    if(argc>=3) minimum_codegree=std::stoi(argv[2]);
    if(rank_limit<2 || rank_limit>40 || minimum_codegree<3) return 2;
    for(int p=2;p<=2*rank_limit+1;++p) {
        bool prime=true; for(int i=2;i*i<=p;++i) if(p%i==0) prime=false;
        if(prime) primes.push_back(p);
    }
    auto start=std::chrono::steady_clock::now();
    for(int q=1;q<rank_limit;++q) {
        if(argc==4 && q!=std::stoi(argv[3])) continue;
        bool squarefree=true;
        for(auto p:primes) if(q%(p*p)==0) squarefree=false;
        if(!squarefree) continue;
        Mask qm=prime_mask(q);
        for(int b=1;q*b*b<rank_limit;++b) {
            if(q*b*b<minimum_codegree) continue;
            if(rank_limit<4) continue;
            Counts packet={};
            Mask bm=prime_mask(b);
            for(unsigned i=0;i<primes.size();++i) if((qm|bm)&(Mask(1)<<i)) packet[i]=1;
            prefix={static_cast<U>(b)};
            dfs(q,qm,b,bm,q*b*b-1,rank_limit-2,packet,1);
        }
        std::cerr<<"q="<<q<<" cumulative nodes "<<nodes<<" types "<<types<<"\n";
    }
    std::cerr<<"Complete rank "<<rank_limit<<", nodes "<<nodes<<", types "<<types
             <<", seconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"\n";
}
