// Exact square-dimension prefix search with fixed-width arithmetic. Requires C++17 and unsigned __int128.
// See README for the newline-delimited job, handoff, and completion formats.
#include <algorithm>
#include <cerrno>
#include <system_error>
#include <unistd.h>
#include <charconv>
#include <string_view>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <numeric>
#include <unordered_map>
#include <stdexcept>
#include <vector>
#include <set>
using U = std::uint64_t;
using Wide = __uint128_t;
using Mask = std::uint64_t;
// This rank-44 experiment retains the published rank <=43 branches unchanged. All support primes are below rank.
static constexpr unsigned native_rank_limit=44;
static constexpr std::array<unsigned,13> primes={2,3,5,7,11,13,17,19,23,29,31,37,41};
using Counts = std::array<unsigned char,primes.size()>;
static unsigned num_primes;
static Mask fixed_support_branch=0;
static std::set<std::vector<U>> wide_prefixes;
// Exact continuation is exercised in tests by compiling a smaller ceiling.
#ifndef SQ_LCM_BITS
#define SQ_LCM_BITS 58
#endif
static_assert(SQ_LCM_BITS >= 3 && SQ_LCM_BITS <= 58);
static constexpr U lcm_limit = U(1) << SQ_LCM_BITS;

static std::vector<U> prefix;
static U nodes, types;
static unsigned rank_limit;
static bool listing=false;
static bool bounded_pilot=false;
struct PilotLimit {};
static unsigned split_depth=22;
static std::size_t frontier_per_job=128;
static bool list_cost=false;
static std::vector<U> fixed;
static U handoffs=0, jobs=0;

// Flush only complete records. Writes within a pipe's atomic-write limit
// cannot publish a truncated JSON line if the producer is killed.
static std::size_t pipe_limit=512;
static std::string records;

static void write_all(const std::string& bytes) {
    std::size_t done=0;
    while(done<bytes.size()) {
        auto n=::write(STDOUT_FILENO,bytes.data()+done,bytes.size()-done);
        if(n<0 && errno==EINTR) continue;
        if(n<=0) throw std::system_error(n<0?errno:EIO,std::generic_category(),"stdout");
        done+=static_cast<std::size_t>(n);
    }
}

static void flush_records() {
    write_all(records);
    records.clear();
}

static void record(std::string line) {
    line.push_back('\n');
    if(records.size()+line.size()>pipe_limit) flush_records();
    if(line.size()>pipe_limit) write_all(line); else records+=line;
}

static void decimal(std::string& out,Wide number) {
    char buffer[40];
    if(number<=UINT64_MAX) {
        auto result=std::to_chars(buffer,buffer+sizeof(buffer),static_cast<U>(number));
        out.append(buffer,result.ptr);
    } else {
        auto ptr=buffer+sizeof(buffer);
        do {*--ptr=static_cast<char>('0'+number%10); number/=10;} while(number);
        out.append(ptr,buffer+sizeof(buffer));
    }
}
static void write_prefix(const std::vector<U>& bs) {
    std::cout << '[';
    for (std::size_t i=0;i<bs.size();++i) std::cout << (i?",":"") << bs[i];
    std::cout << ']';
}

// DFS visits each child's complete subtree contiguously, even when the
// smooth-number iterator visits the children in non-numerical order. Keep
// at most frontier_per_job unreported child counts at each open depth.
// When a subtree exceeds the threshold, no ancestor can coalesce it.
struct Group {
    U count=0;
    std::vector<std::pair<U,U>> pending;
};
static std::vector<U> planner_path;
static std::vector<Group> groups(1);

static void write_job(U weight) {
    if(list_cost) std::cout << '[';
    write_prefix(planner_path);
    if(list_cost) std::cout << ',' << weight << ']';
    std::cout << '\n';
    ++jobs;
}

static void close_group() {
    const U child=planner_path.back(), weight=groups.back().count;
    groups.pop_back();
    planner_path.pop_back();
    auto& parent=groups.back();
    const bool can_group=!planner_path.empty() && planner_path.size()>=fixed.size();
    const bool grouped_before=can_group && parent.count<=frontier_per_job;
    if (UINT64_MAX-parent.count<weight) throw std::overflow_error("planner frontier count");
    parent.count+=weight;
    if(can_group && parent.count<=frontier_per_job) {
        parent.pending.emplace_back(child,weight);
        return;
    }
    if(grouped_before) {
        for(auto [label,cost]:parent.pending) {
            planner_path.push_back(label);
            write_job(cost);
            planner_path.pop_back();
        }
        parent.pending.clear();
    }
    if(weight<=frontier_per_job && planner_path.size()+1>=fixed.size()) {
        planner_path.push_back(child);
        write_job(weight);
        planner_path.pop_back();
    }
}

static void task() {
    const auto& key=prefix.size()<fixed.size()?fixed:prefix;
    std::size_t common=0;
    while(common<planner_path.size() && common<key.size() &&
          planner_path[common]==key[common]) ++common;
    while(planner_path.size()>common) close_group();
    for(auto i=common;i<key.size();++i) {
        planner_path.push_back(key[i]);
        groups.emplace_back();
    }
    ++groups.back().count;
}

static void exact_task() {
    if (listing) {
        task();
        return;
    }
    if(fixed_support_branch && !wide_prefixes.insert(prefix).second) return;
    std::string line="{\"exact\":[";
    const auto& bs=prefix.size()<fixed.size()?fixed:prefix;
    for(std::size_t i=0;i<bs.size();++i) {
        if(i) line.push_back(',');
        decimal(line,bs[i]);
    }
    line+="]}";
    record(std::move(line));
    ++handoffs;
}

struct PacketState {
    Counts maximum{}, count{}, packets{}, pure_max{}, pure_count{};
    int run=0;
};
static PacketState append(const PacketState& old,bool equal,Mask eligible) {
    PacketState st=old;
    st.run=equal?old.run+1:1;
    for(Mask a=eligible & ~Mask(1);a;a&=a-1) {
        unsigned i=__builtin_ctzll(a),m=(primes[i]-1)/2;
        ++st.count[i]; st.maximum[i]=std::max<int>(st.maximum[i],st.run);
        if(eligible==(Mask(1)<<i)) {
            st.pure_max[i]=std::max<int>(st.pure_max[i],st.run);
            ++st.pure_count[i];
        }
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
        if (st.count[i]+unsigned(h) < 2*m || (maximum<2*m &&
            st.packets[i]+(h+current%m)/m<3)) return false;
    }
    return true;
}
// At these ranks primes >=7 are isolated, and a denominator five also
// isolates five. On a nontrivial complementary phase fiber the isolated
// prime does not divide the codegree. Thus a mixed multiple of that
// prime cannot supply any odd packet.
static Mask active_mask(Mask mask, bool has_five) {
    const Mask odd=mask&~Mask(1);
    if(odd&(odd-1)) return 0;
    const Mask isolated=(mask & ~Mask(7)) | (has_five ? mask & Mask(4) : 0);
    return isolated && (mask & (mask-1)) ? 0 : mask;
}

static U square_root(Wide x) {
    U y = x<=UINT64_MAX ? static_cast<U>(std::sqrt(static_cast<double>(static_cast<U>(x))))
                        : static_cast<U>(std::sqrt(static_cast<long double>(x)));
    while (Wide(y)*y > x) --y;
    while (y < UINT64_MAX && Wide(y+1)*(y+1) <= x) ++y;
    return y;
}
static constexpr unsigned mask_table_size=4096;
static Mask small_masks[mask_table_size]{};
static bool small_smooth[mask_table_size]{};
static void initialize_masks() {
    for (unsigned n=1;n<mask_table_size;++n) {
        unsigned rest=n;
        for (unsigned j=0;j<num_primes && rest>1;++j) {
            const unsigned p=primes[j];
            if(rest%p) continue;
            small_masks[n]|=Mask(1)<<j;
            do {rest/=p;} while(rest%p==0);
        }
        small_smooth[n]=(rest==1);
    }
}
static Mask prime_mask(Wide x) {
    if (x<mask_table_size) return small_masks[static_cast<U>(x)];
    Mask mask = 0;
    for (unsigned i=0; i<num_primes; ++i) if (x%primes[i]==0)
        mask |= Mask(1)<<i;
    return mask;
}
static bool sieve(U x, Mask allowed, Mask &mask) {
    // The DFS tests many adjacent small denominators. Precompute both
    // their exact support and whether a prime exceeded the allowed list.
    if(x<mask_table_size) {
        mask=small_masks[x];
        return small_smooth[x] && !(mask & ~allowed);
    }
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
// At four fixed primes, one capped list serves every denominator interval
// of any width in this search process. A prime denominator is known to
// occur to exponent one in M. The table entries include their support
// masks, just as in the unrestricted smooth iterator above.
struct SmoothValue {U value;Mask mask;};
static std::unordered_map<Mask,std::vector<SmoothValue>> ordered_smooth;
static const std::vector<SmoothValue>& smooth_table(Mask allowed,Mask capped) {
    const Mask key=allowed|(capped<<16);
    auto it=ordered_smooth.find(key);
    if(it!=ordered_smooth.end()) return it->second;
    std::vector<SmoothValue> out;
    auto generate=[&](auto&& self,Mask remaining,U value,Mask used)->void {
        if(!remaining) {out.push_back({value,used}); return;}
        const Mask bit=remaining&-remaining;
        const U p=primes[__builtin_ctzll(bit)];
        self(self,remaining^bit,value,used);
        while(value<=UINT64_MAX/p) {
            value*=p;
            self(self,remaining^bit,value,used|bit);
            if(capped&bit) break;
        }
    };
    generate(generate,allowed,1,0);
    std::sort(out.begin(),out.end(),[](const auto&a,const auto&b){return a.value<b.value;});
    return ordered_smooth.emplace(key,std::move(out)).first->second;
}

static std::unordered_map<U,std::vector<U>> divisor_cache;
static const std::vector<U> &nonunit_divisors(U L,Mask support) {
    auto it=divisor_cache.find(L);
    if(it!=divisor_cache.end()) return it->second;
    std::vector<U> ds={1};
    U n=L;
    for(unsigned i=0;i<num_primes;++i) if(support&(Mask(1)<<i)) {
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
    // Bound the dynamic storage per worker, including slack in each vector.
    static std::size_t cached_bytes=0;
    const auto bytes=ds.capacity()*sizeof(U)+sizeof(ds)+64;
    if(cached_bytes+bytes>4*1024*1024) {
        divisor_cache.clear();
        cached_bytes=0;
    }
    cached_bytes+=bytes;
    return divisor_cache.emplace(L,std::move(ds)).first->second;
}
#include "prefilter.hpp"

// A full row of degree M/p at p || M fixes a signed integral column.
static Wide gcd_wide(Wide x,Wide y) {
    while(y) {Wide rem=x%y; x=y; y=rem;}
    return x;
}

static bool prime_column(const std::vector<Wide>& degrees,Wide M,unsigned p,unsigned m) {
    if(!m) return true;
    if(M%(p*p)==0) return false; // Prime-denominator theorem through rank 42.
    Wide n=M/p,common=0,last=0;
    unsigned distinct=0;
    std::vector<Wide> proper;
    for(unsigned i=1;i<degrees.size();++i) {
        const Wide d=degrees[i];
        if(d==n || d%p==0) continue;
        common=gcd_wide(common,d);
        distinct+=d!=last; last=d;
        proper.push_back(d);
    }
    if(common!=1 || distinct<3) return false;
    for(unsigned j=1;j<p && (p-1)*j<=m;++j) {
        const int z=int(p*j)-int(m);
        const int norm=int(j*(p-j))+z;
        if(norm<int(proper.size()) ||
           Wide(p*std::abs(z))*n>Wide(j*(p-j))*n+(p+1)*j) continue;
        const unsigned residue=(n%p)*j%p;
        if(residue!=1 && residue!=p-1) continue;
        const int epsilon=residue==1 ? -1 : 1;
        const __int128 target=__int128(n)*z+epsilon*int(j);
        std::vector<int> choices;
        for(Wide d:proper) {
            const unsigned k=(n/d)%p;
            int a=1;
            while((a*k)%p!=1) ++a;
            choices.push_back(a);
        }
        auto solve=[&](auto&& self,unsigned i,int remaining,__int128 sum)->bool {
            if(i==choices.size()) return remaining==0 && sum==target;
            for(int a:{choices[i],choices[i]-int(p)}) {
                if(remaining-a*a<int(choices.size()-i-1)) continue;
                if(self(self,i+1,remaining-a*a,sum+__int128(a)*__int128(proper[i])))
                    return true;
            }
            return false;
        };
        if(solve(solve,0,norm,0)) return true;
    }
    return false;
}

static void pair_tail(U L,Mask support,Wide E,U low,U high,
                      const PacketState &st) {
    Mask future=support; // A new prime needs at least three distinct nonunit denominators.

    auto emit=[&](U k,U d,U e) {
        Wide root=Wide(k)*L;
        if(root/e < low || root/e > high) return;
        if(std::gcd(k,d)!=1 || std::gcd(k,e)!=1) return;
        Wide b=root/e,c=root/d;
        Mask km; if(!sieve(k,future,km) || !(k%2) || !(k%5))return;
        Mask bm=prime_mask(b),cm=prime_mask(c),real_support=support|km;
        const bool has_five=std::binary_search(prefix.begin(),prefix.end(),5) || b==5 || c==5;
        auto p1=append(st,b==prefix.back(),active_mask(bm,has_five));
        auto p2=append(p1,c==b,active_mask(cm,has_five));
        if(!possible(p2,active_mask(cm,has_five),real_support,0)) return;
        unsigned n5=std::count(prefix.begin(),prefix.end(),5)+unsigned(b==5)+unsigned(c==5);
        if(root%2 || __builtin_popcountll(real_support)!=4) return;
        if(n5 && (root%25==0 || n5!=5 ||
                  __builtin_popcountll(real_support)<4)) return;
        if(n5==5) {
            int other=0;
            Wide previous=0,positive=0,negative=0;int distinct=0;
            auto matches=[&](Wide w) {
                if(w<=5 || w%5!=0) return true;
                if(w%25!=5 && w%25!=20) return false;
                if(w%25==5) positive+=root/w; else negative+=root/w;
                ++other;distinct+=previous!=w;previous=w;
                return true;
            };
            for(U w:prefix) if(!matches(w)) return;
            if(!matches(b) || !matches(c) || other!=4 || distinct<3 || (positive>negative ? positive-negative : negative-positive)!=1) return;
        }
        std::vector<Wide> degrees={1,d,e};
        for(auto jt=prefix.rbegin();jt!=prefix.rend();++jt) {
            Wide x=root/(*jt);
            if(__builtin_popcountll(prime_mask(x))<2) return;
            degrees.push_back(x);
        }
        for(unsigned i=0;i<num_primes;++i) if(real_support&(Mask(1)<<i)){
            int free=0;
            for(unsigned j=1;j<degrees.size();++j) {
                if(degrees[j]%primes[i]) ++free;
            }
            if(free<2) return;
        }
        // All four-row hit degrees divide M/4 if it is even or occurs
        // four times. Their nonunit degrees have gcd one.
        unsigned fours=std::count(prefix.begin(),prefix.end(),4)+(b==4)+(c==4);
        if(fours && (root%8==0 || fours>=4)) {
            Wide gcd_hits=0;
            for(unsigned j=1;j<degrees.size();++j)
                if((root/degrees[j])%4==0) gcd_hits=gcd_wide(gcd_hits,degrees[j]);
            if(gcd_hits!=1) return;
        }
        for(unsigned i=0;i<num_primes;++i) if(real_support & (Mask(1)<<i)) {
            const unsigned prime=primes[i];
            Wide common=0;
            for(unsigned j=1;j<degrees.size();++j)
                if((root/degrees[j])%prime==0) common=gcd_wide(common,degrees[j]);
            if(common!=1) return;
        }
        for(unsigned i=2;i<num_primes;++i) if(real_support & (Mask(1)<<i)) {
            const unsigned p=primes[i];
            const unsigned m=std::count(prefix.begin(),prefix.end(),p)+(b==p)+(c==p);
            if(!prime_column(degrees,root,p,m)) return;
        }
        if (++types == 0) throw std::overflow_error("type count exceeds 64 bits");
        std::string line="[";
        for(std::size_t i=0;i<degrees.size();++i) {
            if(i) line.push_back(',');
            decimal(line,degrees[i]);
        }
        line.push_back(']');
        record(std::move(line));
    };
    // k is odd. Reducing E*k^2=1+d^2+e^2 modulo eight gives these five residues.
    if (((110ULL >> unsigned(E%8)) & 1)==0) return;
    Wide ll=Wide(L)*L;
    U lo=low;
    U hi=high;
    if(lo>hi || !tail_cache::possible(L,support,E)) return;
    auto from_interval=[&] {
        for(U b=lo;b<=hi;++b) {
            U g=std::gcd(L,b),k=b/g,e=L/g;
            if(k%2==0 || k%5==0) continue;
            Wide mass=E*k*k;
            if(mass<=Wide(e)*e+1) continue;
            Wide d2=mass-Wide(e)*e-1;
            // The count of residue classes cut out modulo 64 is 52/64.
            constexpr U residues=0x202021202030213ULL;
            if (!(residues & (U(1)<<unsigned(d2%64)))) continue;
            U d=square_root(d2);
            if(Wide(d)*d==d2 && d>=6 && d<=e && L%d==0 &&
               __builtin_popcountll(prime_mask(d))>=2 && __builtin_popcountll(prime_mask(e))>=2)
                emit(k,d,e);
        }
    };
    if(hi-lo<=8) {from_interval(); return;}
    const auto& ds=nonunit_divisors(L,support);
    if(ds.empty()) return;
    if(hi-lo<=ds.size()/4) {from_interval(); return;}

    // Only smooth k coprime to 10 can complete this prefix. If p divides k,
    // both degrees divide the part of L prime to p.
    U kmax=square_root((1+2*ll)/E);
    Mask eligible=future & ~Mask(5); // exclude primes 2 and 5
    std::array<U,primes.size()> limit{};
    for(Mask a=eligible;a;a&=a-1) {
        unsigned i=__builtin_ctzll(a);
        U rest=L;
        while(rest%primes[i]==0) rest/=primes[i];
        limit[i]=square_root((1+2*Wide(rest)*rest)/E);
        if (limit[i]<primes[i]) eligible &= ~(Mask(1)<<i);
    }
    auto try_k=[&](auto&& self,U k,Mask remaining,U bound)->void {
        Wide target=E*k*k-1;
        U last=std::min(square_root(target),static_cast<U>(std::min<Wide>(L,Wide(L)*k/prefix.back())));
        U first=square_root(target/2);
        if(Wide(2)*first*first<target) ++first;
        auto stop=std::upper_bound(ds.begin(),ds.end(),last);
        for(auto it=std::lower_bound(ds.begin(),stop,first);it!=stop;++it) {
            U e=*it;
            Wide d2=target-Wide(e)*e;
            constexpr U residues=0x202021202030213ULL;
            if (!(residues & (U(1)<<unsigned(d2%64)))) continue;
            U d=square_root(d2);
            if(d>=6 && Wide(d)*d==d2 && L%d==0 && __builtin_popcountll(prime_mask(d))>=2)
                emit(k,d,e);
        }
        while(remaining) {
            unsigned i=__builtin_ctzll(remaining);
            remaining &= remaining-1;
            U cap=std::min(bound,limit[i]);
            U power=k, p=primes[i];
            while(power<=cap/p) {
                power*=p;
                self(self,power,remaining,cap);
            }
        }
    };
    try_k(try_k,1,eligible,kmax);
}
// Charge the least square allowed by the component-column congruence.
static constexpr U norm_scale=U(1)<<48, norm_root=U(1)<<24;
static U upper_mass(U denominator) {
    if(denominator>=norm_root) return 1;
    return 1+(norm_scale-1)/(denominator*denominator);
}
static U inverse_mod(U a,U n) {
    U r=n,s=a%n;
    __int128 x=0,y=1;
    while(s) {
        U q=r/s,t=r%s; r=s;s=t;
        __int128 v=x-__int128(q)*y; x=y;y=v;
    }
    return x<0 ? U(x+n) : U(x);
}
struct SmallClasses {
    unsigned char fours=0, fives=0, sixes=0, sevens=0, proper_four_hits=0;
    bool isolated_five() const {return fives || fours>=9 || (sevens && fours+sixes+27>rank_limit);}
    std::array<std::array<unsigned char,8>,2> paired_heights{};
    std::array<U,primes.size()> pure_mass{},mixed_cost{};
    std::array<std::int64_t,primes.size()> closed_net{};
    U last=0;
    unsigned char runs=0, pure_index=0;
    SmallClasses child(U b) const {
        SmallClasses next=*this;
        if(b!=last) {
            if(pure_index) {
                const int p=primes[pure_index];
                const int coefficient=2*p*(runs/(p-1))-runs;
                const U lower=last>=norm_root ? 0 : norm_scale/(last*last);
                const U weight=coefficient>0 ? upper_mass(last) : lower;
                next.closed_net[pure_index]+=coefficient*std::int64_t(weight);
            }
            next.last=b;next.runs=0;next.pure_index=0;
        }
        ++next.runs;
        for(Mask seen=prime_mask(b)&~Mask(3);seen;seen&=seen-1) {
            unsigned j=__builtin_ctzll(seen),p=primes[j];
            U k=b;
            do {k/=p;} while(k%p==0);
            if(k==1) {next.pure_index=j;next.pure_mass[j]+=upper_mass(b);continue;}
            U u=b/k,a=inverse_mod(k,u);
            U v=Wide(norm_root)*std::min(a,u-a)/u;
            next.mixed_cost[j]+=v*v;
        }
        next.fours+=b==4;
        next.fives+=b==5;next.sixes+=b==6;next.sevens+=b==7;
        next.proper_four_hits+=b>4 && b%4==0;
        if(b%2==0) {
            unsigned height=__builtin_ctzll(b);
            U odd=b>>height;
            if(odd>1) for(unsigned j=0;j<2;++j) {
                U p=j ? 5 : 3, rest=odd;
                while(rest%p==0) rest/=p;
                if(rest==1) ++next.paired_heights[j][std::min(7u,height-1)];
            }
        }
        return next;
    }
};
// The pure p-power labels contribute their full mass to the norm,
// while g_p is bounded by p/(p-1) times that mass. The remaining
// p-power labels have denominators at least the next power >= last.
static bool component_room(const SmallClasses& st,Mask known,U last,int h) {
    for(Mask u=known&~Mask(3);u;u&=u-1) {
        unsigned j=__builtin_ctzll(u),p=primes[j];
        if(!st.mixed_cost[j] || (p==5 && !st.isolated_five())) continue;
        U floor=p;
        while(floor<last && floor<norm_root) floor*=p;
        U pure=st.pure_mass[j]+h*upper_mass(floor);
        if((p-1)*st.mixed_cost[j]>(p+1)*pure) return false;
    }
    return true;
}

// Test the actual bounded trace digits at the closed degree classes.
// All rounding of unassigned reciprocal mass is upwards. Exceeding a
// storage or lattice bound accepts the prefix for ordinary recursion.
static bool gauss_digits(U L,Wide E,int h,const SmallClasses& st) {
    U last=prefix.back();
    const bool twelve_fours=st.fours==12 && st.fives==5;
    const bool ten=st.fours==10 && st.fives==5 && st.sixes==6 && last>6;
    const bool coupled=st.fours && st.fours%2==0 && st.sixes<=8 && last>6;
    const bool parity_eight=last>8 && std::count(prefix.begin(),prefix.end(),8)%2;
    const bool singleton_eighteen=last>18 && std::count(prefix.begin(),prefix.end(),18)==1;
    if((last<32 || h>17) && !(twelve_fours && last>8) && !coupled && !ten && !parity_eight && !singleton_eighteen) return true;
    const Mask isolated=(~Mask(7))|(st.isolated_five() ? Mask(4):0);
    struct Digit {U denominator,step;unsigned maximum;U upper;};
    std::array<Digit,native_rank_limit> digits{};
    constexpr U accuracy=U(1)<<55;
    auto mass=[&](U b)->U {return 1+Wide(accuracy-1)/(Wide(b)*b);};
    unsigned used=0;U P=1,A=1;
    unsigned neutral=0;
    if(st.fives && last>7) {
        unsigned five=0,seven=0,shared=0,n7=0;
        for(U b:prefix) {five+=b>5 && b%5==0; seven+=b>7 && b%7==0;
                         n7+=b==7;shared+=b%35==0;}
        neutral=4-std::min(4u,five);
        if(n7==7) {
            unsigned min7=st.fours>=8 ? 6:3;
            neutral+=min7-std::min(min7,seven);
            if(!shared && neutral) --neutral;
        }
    }
    if(neutral>unsigned(h)) return false;
    U tail=2*(h-neutral)*mass(last);
    if(E) {
        Wide ll=Wide(L)*L;
        unsigned shift=0;
        while((ll>>shift)>(U(1)<<62)) shift+=8;
        U den=ll>>shift, num=1+((E-1)>>shift);
        tail=std::min<U>(tail,1+(2*Wide(accuracy)*num-1)/den);
    }
    const U Mmin=std::max<U>(L,last<=UINT64_MAX/6 ? 6*last : UINT64_MAX);
    const U error=1+(accuracy-1)/Mmin;
    std::array<bool,native_rank_limit> zero{};
    U remaining=accuracy;
    auto known_zero=[&](U b,unsigned n) {
        return b==5 || ((b==7 || b==11 || b==13) && n==b) ||
               (b==4 && (ten || twelve_fours)) || (b==6 && ten) ||
               (b==8 && n==2 && twelve_fours);
    };
    for(unsigned i=0;i<prefix.size();) {
        U b=prefix[i];unsigned j=i+1;while(j<prefix.size() && prefix[j]==b) ++j;
        if(b<last && known_zero(b,j-i)) {
            zero[i]=true;remaining-=U(Wide(j-i)*accuracy/(Wide(b)*b));
        }
        i=j;
    }
    // The mean-square norm of a nonzero p-primary root sum is at
    // least s(p-s)/(p-1), where s is its multiplicity modulo p, and
    // at least two when s=0. Reapply the triangle bound after each
    // zero class; all mass rounding is outwards.
    bool changed=true;
    while(changed) {
        changed=false;
        // The same weighted residue test uses closed divisors of 72.
        // At two the odd-sized dyadic 8-class alone survives. At three
        // the singleton 18-class alone survives, and its root is a unit.
        for(unsigned p: {2u,3u}) {
            if(p==2 ? !parity_eight : !singleton_eighteen) continue;
            U upper=remaining+error;
            for(unsigned i=0;i<prefix.size();) {
                const unsigned start=i;const U b=prefix[i];
                while(i<prefix.size() && prefix[i]==b) ++i;
                if(b==last || zero[start] || 72%b) continue;
                if(p==2 ? (b%8==0 && b!=8) : (b%9==0 && b!=18)) continue;
                upper-=U(Wide(i-start)*accuracy/(Wide(b)*b));
            }
            if(Wide(5184)*upper<accuracy) return false;
        }
        // Short dyadic projections of the complete 4-, 6-, and
        // 8-classes. The bounds use the coefficient l1 norm after
        // normalized relative trace; every remaining root costs <=1.
        if(coupled) {
            const unsigned i6=std::lower_bound(prefix.begin(),prefix.end(),6)-prefix.begin();
            const U w4=zero[0] ? 0 : U(Wide(st.fours)*accuracy/16);
            const U w6=!st.sixes || zero[i6] ? 0 : U(Wide(st.sixes)*accuracy/36);
            U upper=remaining-w4-w6+error;
            if((w4 || w6) && (72*upper<accuracy || (st.sixes<=3 && 36*upper<accuracy))) {
                if(st.sixes==1) return false;
                zero[0]=true;if(st.sixes) zero[i6]=true;
                remaining-=w4+w6;changed=true;
            } else if(w4 && st.sixes<=5 && last>8) {
                auto first8=std::lower_bound(prefix.begin(),prefix.end(),8);
                auto end8=std::upper_bound(first8,prefix.end(),8);
                const unsigned n8=end8-first8,i8=first8-prefix.begin();
                if(n8==0 || n8==2) {
                    const U w8=!n8 || zero[i8] ? 0 : U(Wide(n8)*accuracy/64);
                    if(96*(upper-w8)<accuracy) {
                        zero[0]=true;remaining-=w4;changed=true;
                    }
                }
            }
        }
        // Dyadic cancellation applies to any two complete classes.
        // Floors on the norm thresholds and an upper bound on the
        // complementary mass make both strict tests conservative.
        for(unsigned ia=0;ia<prefix.size();) {
            const unsigned astart=ia;const U a=prefix[ia];
            while(ia<prefix.size() && prefix[ia]==a) ++ia;
            const unsigned na=ia-astart;
            if(a==last || zero[astart] || (a&(a-1)) || (na&1)) continue;
            const U wa=U(Wide(na)*accuracy/(Wide(a)*a));
            const U ua=U(Wide(accuracy)/(Wide(a)*a));
            bool removed=false;
            for(unsigned ib=0;ib<prefix.size();) {
                const unsigned bstart=ib;const U b=prefix[ib];
                while(ib<prefix.size() && prefix[ib]==b) ++ib;
                const unsigned nb=ib-bstart;
                if(b==last || bstart==astart || zero[bstart]) continue;
                const Wide bb=Wide(b)*b;
                const U wb=U(Wide(nb)*accuracy/bb),ub=U(Wide(accuracy)/bb);
                const U upper=remaining-wa-wb+error;
                if(upper>=ub) continue;
                const U half=U(1+(Wide(nb)*accuracy-1)/(2*bb));
                if(ua<=half) continue;
                const U delta=ua-half;
                if(Wide(upper)*upper>=2*Wide(delta)*delta) continue;
                const Mask ps=prime_mask(b);
                const bool neutral=(ps&isolated) && (ps&(ps-1));
                if(nb==1 || neutral || (!(ps&(ps-1)) && nb%primes[__builtin_ctzll(ps)])) return false;
                zero[astart]=zero[bstart]=true;
                remaining-=wa+wb;changed=removed=true;break;
            }
            if(removed) break;
        }
        for(unsigned i=0;i<prefix.size();) {
            unsigned begin=i;U b=prefix[i];while(i<prefix.size() && prefix[i]==b) ++i;
            unsigned n=i-begin;
            if(zero[begin] || last==b) continue;
            Wide bb=Wide(b)*b;
            U weight=U(Wide(n)*accuracy/bb);
            Mask ps=prime_mask(b);
            bool neutral=(ps&isolated) && (ps&(ps-1));
            U upper=remaining-weight+error,unit=accuracy/bb;
            if(neutral && upper<weight) return false;
            bool forced_zero=upper<unit;
            if(ps && !(ps&(ps-1))) {
                unsigned p=primes[__builtin_ctzll(ps)],residue=n%p;
                unsigned numerator=residue ? residue*(p-residue):2;
                unsigned denominator=residue ? p-1:1;
                forced_zero=Wide(upper)*upper*denominator < Wide(unit)*unit*numerator;
            }
            if(forced_zero) {
                if(n==1 || neutral || (!(ps&(ps-1)) && n%primes[__builtin_ctzll(ps)])) return false;
                zero[begin]=true;remaining-=weight;changed=true;
            }
        }
    }
    auto exact=[&](U n,U b) {
        if(b>1000000) return false;
        U den=b*b,g=std::gcd(P,den),multiplier=den/g;
        if(P>(U(1)<<55)/multiplier) return false;
        A=A*multiplier-n*(P/g);P*=multiplier;
        return true;
    };
    for(unsigned i=0;i<prefix.size();) {
        U b=prefix[i]; unsigned j=i+1;
        while(j<prefix.size() && prefix[j]==b) ++j;
        const unsigned n=j-i;bool trace_zero=zero[i];i=j;
        Mask ps=prime_mask(b);
        if((ps&isolated) && (ps&(ps-1))) continue;
        if(trace_zero) {if(!exact(n,b)) return true;continue;}
        unsigned scale=1,step=1,maximum=2*n;
        if(ps&isolated) {
            unsigned p=primes[__builtin_ctzll(ps)];
            // With five absent and >=9 fours fewer than 25 labels
            // can be five-active. At p>=7 the whole rank is <p².
            if(p>=7 || st.fours>=9 || (!st.fives && st.isolated_five())) {
                if(n<p-1 && b<last) continue;
                step=p; maximum=n/(p-1);
            } else {scale=p-1;maximum=2*n*scale;}
        } else {
            if(ps&2) scale*=2;
            if(ps&4) scale*=4;
            maximum=2*n*scale;
            if(ps==2 && n<3 && b<last) {step=3;maximum=n;}
        }
        U upper=2*n*mass(b);
        tail+=upper;
        if(b<=1000000 && b<last) digits[used++]={scale*b*b,step,maximum,upper};
    }
    U H=P;
    std::array<std::int64_t,160> old{},next{};
    unsigned count=1;old[0]=A;
    for(unsigned i=0;i<used;++i) {
        auto d=digits[i];
        U multiplier=d.denominator/std::gcd(H,d.denominator);
        if(H>(U(1)<<56)/multiplier) break;
        H*=multiplier;
        U weight=(H/d.denominator)*d.step;
        tail-=d.upper;
        U bound=Wide(H)*(tail+error)/accuracy,lower=Wide(H)*error/accuracy;
        unsigned nn=0;
        for(unsigned j=0;j<count;++j) {
            const auto value=old[j]*std::int64_t(multiplier);
            if(value < -std::int64_t(lower)) continue;
            U start=value<=std::int64_t(bound) ? 0 : U(value-bound+weight-1)/weight;
            unsigned end=std::min<U>(d.maximum,U(value+std::int64_t(lower))/weight);
            for(unsigned k=start;k<=end;++k) {
                auto r=value-std::int64_t(k*weight);
                if(!r) continue;
                if(nn==next.size()) return true;
                next[nn++]=r;
            }
        }
        if(!nn) return false;
        std::sort(next.begin(),next.begin()+nn);
        count=std::unique(next.begin(),next.begin()+nn)-next.begin();
        std::copy_n(next.begin(),count,old.begin());
    }
    return true;
}

static U mixed_floor(U last,U prime) {
    U best=UINT64_MAX;
    for(U power=2*prime;power<best;) {
        U q=power;
        while(q<last) q*=2;
        best=std::min(best,q);
        if(power>best/prime) break;
        power*=prime;
    }
    return best;
}

// Each mixed p phase must contain four denominators of one dyadic
// height. Pooling all heights >=8 can only reduce the missing count.
static unsigned mixed_need(const std::array<unsigned char,8>& heights,unsigned phases) {
    unsigned full=0,residues[4]={};
    for(unsigned n:heights) {full+=n/4;++residues[n%4];}
    if(full>=phases) return 0;
    unsigned missing=4*(phases-full), left=phases-full;
    for(unsigned r=3;r && left;--r) {
        unsigned used=std::min(left,residues[r]);
        missing-=used*r;left-=used;
    }
    return missing;
}

// Primitive ternary phase masses agree. Each initial contribution with
// denominator lcm W is an integral multiple of 1/W². If the entire
// unknown difference is smaller, the prefix itself balances. In each
// completed balanced segment a unique maximal prime height forces its
// coefficient to vanish. A singleton cannot have a higher ternary phase.
static unsigned separated_three(const SmallClasses& st,U last,int h) {
    struct Item {U b,W;unsigned n;};
    std::array<Item,native_rank_limit> eligible{}; unsigned used=0;
    auto bins=st.paired_heights[0];
    U W=1;
    for(unsigned i=0;i<prefix.size();) {
        U b=prefix[i]; unsigned end=i+1;
        while(end<prefix.size() && prefix[end]==b) ++end;
        if(b%3==0) {
            U odd=b>>__builtin_ctzll(b);
            while(odd%3==0) odd/=3;
            if(odd==1) {
                U q=b/std::gcd(W,b);
                if(W>lcm_limit/q) break; // skip all later cuts; retain their possible mass below
                W*=q; eligible[used++]={b,W,end-i};
            }
        }
        i=end;
    }
    unsigned begin=0;
    for(unsigned j=0;j<used;++j) {
        const U bound=eligible[j].W;
        if(bound>=last) break;
        const Wide unit=Wide(256)*bound*bound;
        Wide tail=Wide(h)*(1+(unit-1)/(Wide(last)*last));
        for(U b:prefix) {
            if(b<=eligible[j].b || b%3) continue;
            U odd=b>>__builtin_ctzll(b);
            while(odd%3==0) odd/=3;
            if(odd==1) tail+=1+(unit-1)/(Wide(b)*b);
            if(tail>=256) break;
        }
        if(tail>=256) continue;
        for(bool changed=true;changed;) {
            changed=false;
            for(U prime:{2,3}) {
                U high=0; unsigned idx=0,count=0;
                for(unsigned k=begin;k<=j;++k) if(eligible[k].n) {
                    U part=1,n=eligible[k].b;
                    while(n%prime==0) {n/=prime;part*=prime;}
                    if(part>high) {high=part;idx=k;count=0;}
                    if(part==high) count+=eligible[k].n;
                }
                if(count==1) {
                    U b=eligible[idx].b;eligible[idx].n=0;changed=true;
                    if(!(b%2)) --bins[std::min(7,__builtin_ctzll(b)-1)];
                }
            }
        }
        begin=j+1;
    }
    return mixed_need(bins,2);
}

// High-four omission profiles constrain the dimension support, the
// admissible prime columns and the mixed ternary alternative.
struct SupportCases {
    Mask common,available;
    bool force_three_mixed=false;
    bool force_large_power=false;
};
// The integer masks encode the dimension support and the odd prime
// denominators. A positive root identifies an all-primary completion.
struct FourCase {unsigned fours; Mask support,primary; U root;};
static constexpr FourCase four_cases[]={
    {9,15,8,0},
    {9,15,12,15482880},
    {9,15,4,0},
    {9,15,12,0},
    {9,23,4,0},
    {9,23,20,0},
    {10,15,0,0},
    {10,15,8,0},
    {10,15,4,0},
    {10,15,12,0},
    {10,23,4,0},
    {11,15,0,0},
    {11,15,8,0},
    {11,15,4,0},
    {12,15,0,0},
    {12,15,8,0},
    {12,15,4,0},
};
static constexpr FourCase four_cases42[]={
    {9,15,8,0},
    {9,27,24,0},
    {9,15,12,136080},
    {9,15,12,15482880},
    {9,15,4,0},
    {9,15,12,0},
    {9,23,4,0},
    {9,23,20,0},
    {10,15,0,0},
    {10,15,8,0},
    {10,23,16,0},
    {10,15,4,0},
    {10,15,12,0},
    {10,23,4,0},
    {10,23,20,0},
    {11,15,0,0},
    {11,15,8,0},
    {11,23,16,0},
    {11,15,4,0},
    {12,15,0,0},
    {12,15,8,0},
    {12,15,4,0},
    {13,15,0,0},
};
static constexpr FourCase four_cases43[]={
    {9,15,0,0},
    {9,15,8,0},
    {9,27,24,0},
    {9,15,12,136080},
    {9,15,12,15482880},
    {9,15,4,0},
    {9,15,12,0},
    {9,23,4,0},
    {9,23,20,0},
    {9,39,4,0},
    {9,39,36,0},
    {10,15,0,0},
    {10,15,8,0},
    {10,23,16,0},
    {10,27,8,0},
    {10,27,24,0},
    {10,15,4,0},
    {10,15,12,0},
    {10,23,4,0},
    {10,23,20,0},
    {10,39,4,0},
    {11,15,0,0},
    {11,15,8,0},
    {11,23,16,0},
    {11,15,4,0},
    {11,23,4,0},
    {11,23,20,0},
    {12,15,0,0},
    {12,15,8,0},
    {12,23,16,0},
    {12,15,4,0},
    {13,15,0,0},
};
// Rank44:39 mixed signatures and the two all-primary roots, both
// closed by the exact rank44 divisor checks in the structural receipts.
static constexpr FourCase four_cases44[]={
    {9,15,0,0},
    {9,15,4,0},
    {9,15,8,0},
    {9,15,12,0},
    {9,15,12,136080},
    {9,15,12,15482880},
    {9,23,4,0},
    {9,23,16,0},
    {9,23,20,0},
    {9,27,8,0},
    {9,27,24,0},
    {9,39,4,0},
    {9,39,36,0},
    {9,43,40,0},
    {10,15,0,0},
    {10,15,4,0},
    {10,15,8,0},
    {10,15,12,0},
    {10,23,4,0},
    {10,23,16,0},
    {10,23,20,0},
    {10,27,8,0},
    {10,27,24,0},
    {10,39,4,0},
    {10,39,36,0},
    {11,15,0,0},
    {11,15,4,0},
    {11,15,8,0},
    {11,23,4,0},
    {11,23,16,0},
    {11,23,20,0},
    {11,27,8,0},
    {11,27,16,0},
    {11,27,24,0},
    {12,15,0,0},
    {12,15,4,0},
    {12,15,8,0},
    {12,23,0,0},
    {12,23,4,0},
    {12,23,16,0},
    {13,15,0,0},
};

static SupportCases admissible_support(U L,U last,Mask mask,unsigned fours,unsigned fives) {
    if(fixed_support_branch && (mask&~fixed_support_branch)) return {0,0};
    mask|=fixed_support_branch;
    if(fours<9) return {mask,fixed_support_branch ? fixed_support_branch : ~Mask(0)};
    if(rank_limit>=41) {
        SupportCases joined{~Mask(0),0,true,false};
        Mask present=0,closed=0;
        for(unsigned i=2;i<=5;++i) {
            const unsigned p=primes[i];
            if(last>p) closed|=Mask(1)<<i;
            if(last>=p && std::binary_search(prefix.begin(),prefix.end(),p)) present|=Mask(1)<<i;
        }
        auto cases=rank_limit==41 ? four_cases : rank_limit==42 ? four_cases42 : rank_limit==43 ? four_cases43 : four_cases44;
        auto end=rank_limit==41 ? std::end(four_cases) : rank_limit==42 ? std::end(four_cases42) : rank_limit==43 ? std::end(four_cases43) : std::end(four_cases44);
        for(auto it=cases;it!=end;++it) {
            const auto& c=*it;
            if(c.fours<fours || (last>4 && c.fours!=fours) || (mask & ~c.support) ||
               (present & ~c.primary) || (closed & c.primary & ~present) ||
               c.root) continue;
            joined.common&=c.support;
            joined.available|=c.support;
            joined.force_three_mixed&=c.root==0;
        }
        return joined;
    }
    const bool seven=std::binary_search(prefix.begin(),prefix.end(),7);
    constexpr Mask four_primes=15;
    if((mask&~four_primes) || (fours<=11 && !fives && last>7 && !seven) ||
       ((fours==12 || (fours>=11 && fives)) && seven)) return {0,0};
    return {four_primes,four_primes,true,last>7 && !seven};
}

// At five or seven the signed proper degrees sum to epsilon, so their
// denominators have lcm M. Enumerate four signs at five and the
// three- or six-witness shapes at seven. A fixed dimension support
// permits a larger first denominator. Incomplete families retain the
// ordinary branch.
struct ColumnRoot {U root; std::array<U,6> witnesses; Mask support;};
struct RootFamily {bool complete=true; std::vector<ColumnRoot> rows;};
static std::unordered_map<U,RootFamily> root_families;
static const RootFamily& column_roots(U first,unsigned p,Mask allowed=0,Mask capped=0,bool six=false) {
    const Mask key=(first<<28)|(allowed<<15)|(capped<<2)|(2*six)|(p==7);
    auto found=root_families.find(key);
    if(found!=root_families.end()) return found->second;
    RootFamily family;
    std::array<U,6> ws={first};
    const bool narrow=allowed!=0;
    if(!allowed) allowed=(Mask(1)<<num_primes)-1;
    unsigned work=0;
    const unsigned total=six ? 6 : 4, length=six ? 6 : p==5 ? 4 : 3;
    auto coeff=[p,six](U w) {
        constexpr int seven[7]={0,1,0,-2,2,0,-1};
        return p==7 ? (six && std::abs(seven[(w/7)%7])>1 ? 0 : seven[(w/7)%7]) : w%25==5 ? 1 : w%25==20 ? -1 : 0;
    };
    auto search=[&](auto&& self,unsigned size,U W,std::int64_t S,unsigned charge)->void {
        // For a proper initial segment of signs the signed W/w are all
        // W/p modulo p, so S is nonzero. In the three-witness shape
        // a zero partial sum would leave a single nonunit degree
        // times its coefficient equal to one up to sign.
        if(!S || !family.complete || charge+length-size>total) return;
        const U last=ws[size-1];
        if(size+1<length) {
            const U h=total-charge;
            U high=(6*h+1)*W/(6*std::abs(S));
            auto consider=[&](U z) {
                if(++work>(six ? 15000u : 100000u)) {family.complete=false;return;}
                const int a=coeff(z);
                if(!a || (p==7 && !six && z==last)) return;
                Mask mask;
                if(!sieve(z,allowed,mask)) return;
                const U k=z/std::gcd(W,z);
                // S has absolute value at most four times W/first.
                // A wider lcm accepts the family for ordinary recursion.
                if(W>(U(1)<<52)/k) {family.complete=false;return;}
                const U V=W*k;ws[size]=z;
                self(self,size+1,V,S*std::int64_t(k)+a*std::int64_t(V/z),charge+std::abs(a));
            };
            if(narrow) {
                const auto& nums=smooth_table(allowed,capped);
                auto it=std::lower_bound(nums.begin(),nums.end(),last,
                    [](const auto& x,U v){return x.value<v;});
                for(;it!=nums.end() && it->value<=high && family.complete;++it)
                    if(it->value%p==0) consider(it->value);
            } else for(U z=last;z<=high && family.complete;z+=p) consider(z);
            return;
        }
        // Put d=W/gcd(W,z) and t=z/gcd(W,z). The next root is W*t
        // and S*t+c*d=±1. Each admissible d is a divisor of W.
        auto divisors=[&](auto&& visit,unsigned i,U remainder,U d)->void {
            if(remainder==1) {
                if(d<6) return;
                for(int c:{-2,-1,1,2}) for(int target:{-1,1}) {
                    if(charge+std::abs(c)!=total) continue;
                    const std::int64_t num=target-c*std::int64_t(d);
                    if(num%S || num/S<=0) continue;
                    const U t=num/S;
                    if(std::gcd(d,t)!=1) continue;
                    Wide full=Wide(W)*t;
                    if(full>UINT64_MAX) {family.complete=false;return;}
                    U M=static_cast<U>(full), z=M/d;
                    if(z<last || z%p || coeff(z)!=c || M%2 || M%(p*p)==0) continue;
                    Mask support;
                    if(!sieve(M,allowed,support)) continue;
                    unsigned count=__builtin_popcountll(support);
                    if(count!=4) continue;
                    bool over=false;
                    for(Mask seen=capped;seen;seen&=seen-1) {
                        U q=primes[__builtin_ctzll(seen)];over|=M%(q*q)==0;
                    }
                    if(over) continue;
                    ws[length-1]=z;
                    unsigned distinct=1;
                    for(unsigned j=1;j<length;++j) distinct+=ws[j]!=ws[j-1];
                    if(distinct<3) continue;
                    bool viable=true;
                    for(unsigned j=0;j<length && viable;++j) {
                        U e=M/ws[j];
                        if(__builtin_popcountll(prime_mask(e))<2) viable=false;
                        for(unsigned k=0;k<j;++k)
                            if(std::gcd(e,M/ws[k])==1) viable=false;
                    }
                    if(viable) family.rows.push_back({M,ws,support});
                }
                return;
            }
            while(remainder%primes[i]) ++i;
            const U p=primes[i];
            U power=1;
            while(remainder%p==0) {remainder/=p;power*=p;}
            for(U v=1;;v*=p) {
                visit(visit,i+1,remainder,d*v);
                if(v==power) break;
            }
        };
        divisors(divisors,0,W,1);
    };
    // Intermediate lcms above 2^52 and final roots above 64 bits leave
    // the family incomplete. The final product W*t is evaluated wide.
    const int a=coeff(first);
    if(a) search(search,1,first,a,std::abs(a));
    return root_families.emplace(key,std::move(family)).first->second;
}

// In simultaneous five- and seven-packets a five-column root leaves
// one common proper degree. Every other seven-free proper degree has
// the form 5d, where d divides the two remaining prime powers of M.
// The signed degree sum in the seven-column is epsilon. Meet five
// signs as a pair and a triple; the alternate shape uses one sign
// and one double. A work ceiling retains the candidate.
static bool seven_completion(U M,U common,bool six_only) {
    if(M>(U(1)<<60)) return true;
    static std::unordered_map<U,std::unordered_map<U,bool>> cache;
    auto& row=cache[M];const U key=2*common+six_only;
    if(auto it=row.find(key);it!=row.end()) return it->second;
    const U n=M/7;
    if(n%7!=1 && n%7!=6) return row[key]=false;
    const std::int64_t sign=n%7==1 ? -1 : 1;
    const U n5=M/5;
    if((sign<0 && n5%5!=1) || (sign>0 && n5%5!=4)) return row[key]=false;
    constexpr int seven[7]={0,1,0,-2,2,0,-1};
    const int shared=seven[(common/7)%7];
    if(std::abs(shared)!=1) return row[key]=false;
    const std::int64_t target=sign-shared*std::int64_t(M/common);
    std::vector<std::int64_t> ones,twos;
    const Mask ps=prime_mask(n)&~Mask(4);
    const auto& values=smooth_table(ps,0);
    for(auto z:values) {
        if(z.value>n/10) break;
        const U d=5*z.value;
        if(z.value==1 || n%d) continue;
        const int a=seven[(n/d)%7];
        if(std::abs(a)==1) ones.push_back(a*std::int64_t(d));
        if(std::abs(a)==2) twos.push_back(a*std::int64_t(d));
    }
    if(!six_only) for(auto x:ones) for(auto y:twos) if(x+y==target) return row[key]=true;
    if(ones.size()>80) return row[key]=true;
    std::vector<std::int64_t> pairs;
    for(unsigned i=0;i<ones.size();++i) for(unsigned j=0;j<=i;++j) pairs.push_back(ones[i]+ones[j]);
    std::sort(pairs.begin(),pairs.end());
    for(auto x:ones) for(auto y:pairs)
        if(std::binary_search(pairs.begin(),pairs.end(),target-x-y)) return row[key]=true;
    return row[key]=false;
}

// A width transfer preserves the full subtree. Apply inexpensive slot
// obstructions on its first child as well, before spawning exact work.
static bool width_room(U parent_L,U last,int h,const PacketState& st,
                       const SmallClasses& small, Mask mask) {
    auto cases=admissible_support(parent_L,last,mask,small.fours,small.fives);
    if(!cases.available || !component_room(small,mask,last,h) || !gauss_digits(parent_L,0,h,small)) return false;
    Mask support=cases.common|mask;
    if(__builtin_popcountll(support)>4) return false;
    if(__builtin_popcountll(support)!=4 || !cases.force_three_mixed) return true;
    unsigned active=separated_three(small,last,h), inert=0;
    U two_height=std::max<U>(8,std::max(parent_L & -parent_L,last & -last));
    unsigned have_odds=0,proper5=0,odds_other=0;
    for(U b:prefix) {
        bool odd=b%two_height==0;
        have_odds+=odd;
        odds_other+=odd && b%5!=0;
        proper5+=b>5 && b%5==0;
    }
    inert=3-std::min(3u,have_odds);
    if(small.fives && last>5) {
        if(proper5>4) return false;
        inert=std::max(inert,4-proper5+2-std::min(2u,odds_other));
    }
    unsigned pure_slots=active;
    unsigned prime_inert=0;
    for(Mask u=support&~Mask(3);u;u&=u-1) {
        const unsigned p=primes[__builtin_ctzll(u)];
        pure_slots+=p-std::min<unsigned>(p,st.pure_count[__builtin_ctzll(u)]);
        const bool prime_present=std::binary_search(prefix.begin(),prefix.end(),p);
        unsigned vp=0;
        for(U x:{parent_L,last}) {
            unsigned v=0;while(x%p==0) {x/=p;++v;}
            vp=std::max(vp,v);
        }
        if(prime_present && vp>1) return false;
        Wide power=p;
        for(unsigned k=1;k<(vp+2)/2;++k) power*=p;
        bool paid=false;
        while(power<last) {
            auto eq=std::equal_range(prefix.begin(),prefix.end(),U(power));
            if(unsigned(eq.second-eq.first)>=p-1) {paid=true;break;}
            power*=p;
        }
        if(!paid) {
            unsigned count=power==last ? st.run:0;
            active+=p-1-std::min(p-1,count);
        }
        if(!prime_present && last>p) {
            Wide eventual=power>last ? Wide(last)+1 : power*power/p;
            unsigned paid_free=0;
            if(eventual<=last) for(U b:prefix) paid_free+=b%U(eventual)==0;
            unsigned needed=(p==5 ? 3 : 2);
            prime_inert=std::max(prime_inert,needed-std::min(needed,paid_free));
        }
    }
    return pure_slots+inert<=unsigned(h) && active+std::max(inert,prime_inert)<=unsigned(h);
}

static void dfs(U L, Mask support, Wide E, int h, const PacketState& st, const SmallClasses& classes) {
    const unsigned fours=classes.fours, fives=classes.fives;
    if(fours>(rank_limit>=43 ? 13u : rank_limit>=42 ? 13u : 12u) || fives>5 || !component_room(classes,support,prefix.back(),h) || !gauss_digits(L,E,h,classes)) return;
    if(fours==10 && fives==5 && std::count(prefix.begin(),prefix.end(),6)==6) {
        // T4,T5,T6 vanish. The remaining mass is 1/120. Repeatedly
        // remove a group whose sum of twists has house less than one.
        // Seven is absent by mass, so
        // the isolated-prime height theorem gives M>=8*3*5*7^3.
        const Wide scale=Wide(L)*L;
        const Wide vacuum_error=1+(scale-1)/std::max<U>(L,rank_limit<=40 ? 8*3*5*343 : 8*3*5*7);
        Wide rest=scale/120;
        for(auto it=std::upper_bound(prefix.begin(),prefix.end(),6);it!=prefix.end();) {
            U b=*it;
            auto end=std::upper_bound(it,prefix.end(),b);
            unsigned n=end-it;
            Wide weight=Wide(L/b)*(L/b);
            rest-=n*weight;
            if(rest+vacuum_error>=weight) break;
            if(n==1 || ((n&1) && !(b&(b-1)))) return;
            it=end;
        }
    }
    if(fours==12 && fives==5 && (std::binary_search(prefix.begin(),prefix.end(),6) ||
        std::count(prefix.begin(),prefix.end(),8)>2)) return;

    if(__builtin_popcountll(support)>4) return;
    // Proven restrictions on a closed class of degree M/5.
    const bool closed5 = prefix.back()>5 || h==0;
    auto cases=admissible_support(L,prefix.back(),support,fours,fives);
    if(!cases.available || __builtin_popcountll(cases.common|support)>4) return;
    if (fives && fives!=5 && (closed5 || fives+h<5)) return;
    // Once the published prefix is fixed, separate the possible
    // fourth dimension primes. Pure and omission quotas apply even
    // before that prime first divides an observed denominator.
    // Terminal supports are unique; wide prefixes are continued once
    // with the complete, unrestricted arbitrary-precision search.
    if(!listing && !fixed_support_branch && prefix.size()>=fixed.size() &&
       prefix.back()>=7 && h<=23) {
        const Mask base=cases.common|support|(Mask(1));
        if(__builtin_popcountll(base)==3) {
            for(unsigned i=1;i<num_primes;++i) {
                const Mask bit=Mask(1)<<i;
                if((base&bit) || !(cases.available&bit) || primes[i]>unsigned(h)+1) continue;
                fixed_support_branch=base|bit;
                dfs(L,support,E,h,st,classes);
            }
            fixed_support_branch=0;return;
        }
        if(__builtin_popcountll(base)==4 && __builtin_popcountll(support)<4) {
            fixed_support_branch=base;
            dfs(L,support,E,h,st,classes);
            fixed_support_branch=0;return;
        }
    }
    int witnesses_needed=0;
    U fixed_root=0,first_witness_floor=0;
    const bool cross57=fives==5 && prefix.back()>7 && std::count(prefix.begin(),prefix.end(),7)==7;
    Mask fixed_support=(cases.common|support),capped_witness=0;
    if(__builtin_popcountll(fixed_support)!=4) fixed_support=0;
    else for(Mask a=fixed_support&~Mask(3);a;a&=a-1)
        if(std::binary_search(prefix.begin(),prefix.end(),primes[__builtin_ctzll(a)])) capped_witness|=a&-a;
    unsigned minimal=0;
    U first_pure=0,least_B=0;
    if(fixed_support && cases.force_three_mixed && fives==5 && prefix.back()>7 &&
       (fixed_support&7)==7 && (fixed_support&~Mask(7))) {
        const unsigned j=__builtin_ctzll(fixed_support&~Mask(7)),p=primes[j];
        // Fours, five prime rows, four proper five-free labels, two
        // other odd labels, eight mixed labels and the unit are disjoint.
        const unsigned seats=rank_limit-fours-20;
        if(prefix.back()>p && seats<=2*(p-1) && !std::binary_search(prefix.begin(),prefix.end(),p)) {
            minimal=j;
            U B=p*p;
            while(B<=prefix.back()) {
                if(std::binary_search(prefix.begin(),prefix.end(),B)) {first_pure=B;break;}
                B*=p;
            }
            unsigned observed=0;
            if(first_pure) {
                auto eq=std::equal_range(prefix.begin(),prefix.end(),B);observed=eq.second-eq.first;
                if(observed>p || (B<prefix.back() && observed!=p)) return;
            }
            const unsigned due=p-observed;
            if(due>unsigned(h)) return;
            if(due && B<=square_root(Wide(due)*L*L/E)) {
                if(first_pure) return;
                const U bound=square_root(Wide(p)*L*L/E);
                while(B<=bound && B<(U(1)<<58)) B*=p;
            }
            least_B=B;
            for(U z:prefix) if(z%p==0) {
                U rest=z;while(rest%p==0) rest/=p;
                const U u=z/rest;
                if(u<B || (first_pure && u>Wide(B)*B/p)) return;
            }
        }
    }
    auto compatible_packet=[&](U M) {
        if(!minimal) return true;
        const U p=primes[minimal];U n=M;
        while(n%p==0) n/=p;
        const U mp=M/n;
        const U n5=M/5;
        return mp>=Wide(least_B)*least_B/p && (!first_pure || mp==Wide(first_pure)*first_pure/p) &&
               ((n5%5==1 && n%mp==1) || (n5%5==4 && n%mp==mp-1));
    };
    auto compatible_seven=[&](const ColumnRoot& row) {
        if(!cross57) return true;
        unsigned shared=0;U common=0;
        for(unsigned i=0;i<4;++i) if(U b=row.witnesses[i]; b%35==0) {
            common=b;
            U k=b/35%35;
            if(k!=17 && k!=18) return false;
            ++shared;
        }
        return shared==1 && seven_completion(row.root,common,fours>=8);
    };
    if (closed5 && fives==5) {
        std::array<U,4> observed{};
        U positive=0, negative=0, previous=0, witness_lcm=1;
        unsigned on_fours=0;
        int number=0, distinct=0;
        for(U b:prefix) if (b>5 && b%5==0) {
            if(fours==12 && (b<80 || b%4)) return;
            ++number;
            if(number>4) return;
            observed[number-1]=b;
            if(b%25==5) positive+=L/b;
            else if(b%25==20) negative+=L/b;
            else return;
            witness_lcm=std::lcm(witness_lcm,b);
            on_fours+=b%4==0;
            if(b!=previous) ++distinct;
            previous=b;
        }
        if(cross57) {
            unsigned shared=0;
            for(unsigned i=0;i<unsigned(number);++i) if(observed[i]%35==0) {
                U k=observed[i]/35%35;
                if(k!=17 && k!=18) return;
                ++shared;
            }
            if(shared>1 || (!shared && number==4)) return;
        }
        witnesses_needed=4-number;
        if(number && number<4 && observed[0]<=(fixed_support ? 25000 : 500)) {
            const auto& family=column_roots(observed[0],5,fixed_support,capped_witness);
            if(family.complete) {
                U only=0;
                bool found=false,unique=true;
                for(const auto &row:family.rows) {
                    if(row.root%L || !compatible_packet(row.root) || !compatible_seven(row) || (cases.common&~row.support) ||
                       row.witnesses[number]<prefix.back() ||
                       !std::equal(observed.begin(),observed.begin()+number,row.witnesses.begin())) continue;
                    if(found && only!=row.root) unique=false;
                    only=row.root; found=true;
                }
                if(!found) return;
                if(unique) fixed_root=only;
            }
        }
        if(!number) {
            U lower=std::max<U>(fours==12 ? 80 : 20,prefix.back());
            auto possible_first=[&](U b) {
                if(b%25!=5 && b%25!=20) return false;
                const auto& family=column_roots(b,5,fixed_support,capped_witness);
                if(!family.complete) return true;
                for(const auto &row:family.rows)
                    if(row.root%L==0 && compatible_seven(row) && (cases.common&~row.support)==0) return true;
                return false;
            };
            if(fixed_support && lower<=25000) {
                static std::unordered_map<Mask,std::unordered_map<U,std::vector<std::pair<U,U>>>> cache;
                // The seven-column test permits a double below eight
                // fours. Its cache must retain both possible shapes.
                const Mask shape=(unsigned(cross57)<<1)|unsigned(cross57 && fours>=8);
                auto& intervals=cache[(fixed_support<<16)|(capped_witness<<2)|shape][L];
                U first=lower;
                for(auto [a,b]:intervals) if(a<=lower && lower<=b) {first=b;break;}
                if(first==lower) {
                    const auto& ns=smooth_table(fixed_support,capped_witness);
                    auto it=std::lower_bound(ns.begin(),ns.end(),lower,[](const auto& x,U v){return x.value<v;});
                    for(;it!=ns.end() && it->value<=25000;++it) if(possible_first(it->value)) break;
                    first=it->value;
                    if(first>lower && intervals.size()<256) intervals.emplace_back(lower,first);
                }
                lower=first;
            } else while(lower<=200 && !possible_first(lower)) ++lower;
            first_witness_floor=lower;
        }
        if(fours>=4 && (16-fours)*(on_fours+witnesses_needed)<fours) return;
        if(witnesses_needed>h || distinct+witnesses_needed<3) return;
        U t=positive>=negative ? positive-negative : negative-positive;
        if(!witnesses_needed) {
            if(t!=1 || witness_lcm!=L) return;
            fixed_root=L;
        } else if (Wide(6)*t*prefix.back() > Wide(6*witnesses_needed+1)*L) {
            return;
        }
    }
    // The prime-denominator projection excludes 25|M throughout this range.
    if(fives && L%25==0) return;
    // For p>=7 the projection already forbids p^2|M at these ranks.
    // All active p-labels then have b=p, and at least p-1 are needed.
    for(unsigned i=3;i<num_primes && primes[i]<=prefix.back();++i) {
        const unsigned p=primes[i];
        auto range=std::equal_range(prefix.begin(),prefix.end(),p);
        unsigned count=range.second-range.first;
        if(count && (L%(p*p)==0 || (count<p-1 && (prefix.back()>p || count+h<p-1)))) return;
    }
    int seven_needed=0;
    U seven_floor=0;
    if (L%7==0 && prefix.back()>7) {
        const auto sevens=std::equal_range(prefix.begin(),prefix.end(),7);
        const unsigned n7=sevens.second-sevens.first;
        if(n7) {
            const U lower_n=std::max<U>(L/7,30); // three distinct primes in M/7
            // The full Weyl cost, including several primitive packets.
            // Since the proper squared norm is <=13, each coefficient
            // is the unique least absolute inverse of b/7 modulo seven.
            constexpr int full_coeff[7]={0,1,-3,-2,2,3,-1};
            unsigned spent=0,distinct=0;U previous=0;
            __int128 value=0;
            for(U b:prefix) if(b>7 && b%7==0) {
                const int a=full_coeff[(b/7)%7];
                spent+=a*a;value+=__int128(a)*(L/b);
                distinct+=b!=previous;previous=b;
            }
            const U next_k=(prefix.back()+6)/7;
            U dpos=UINT64_MAX,dneg=UINT64_MAX;
            for(U a=1;a<=3;++a) {
                U kp=next_k,kn=next_k;
                while((kp*a)%7!=1) ++kp;
                while((kn*a)%7!=6) ++kn;
                dpos=std::min(dpos,a*kp);dneg=std::min(dneg,a*kn);
            }
            const U Mmin=std::max<U>(L,6*prefix.back());
            bool valid=false;
            for(unsigned j=1;j<7 && 6*j<=n7;++j) {
                const int z=7*int(j)-int(n7);
                if(Wide(7*std::abs(z))*lower_n>Wide(j*(7-j))*lower_n+8*j) continue;
                const int budget=int(j*(7-j))+z;
                if(budget>13) {valid=true;continue;} // conservative guard
                if(budget<int(spent)) continue;
                const unsigned cost=budget-spent;
                if(distinct+std::min<unsigned>(cost,h)<3) continue;
                const __int128 residual=7*value-__int128(z)*L;
                const U d=residual<0 ? dpos:dneg;
                const Wide magnitude=residual<0 ? Wide(-residual):Wide(residual);
                const Wide error=(Wide(7*j)*L*d+Mmin-1)/Mmin;
                if(magnitude*d>Wide(cost)*L+error) continue;
                valid=true;
            }
            if(!valid) return;
        }
        
        if(n7==7) {
            // One seven-packet has T7=0. After removing it, if the
            // house of 9T4+4T6 (optionally plus a known T12 term) is
            // below one, reduction at two detects an impossible parity.
            const unsigned sixes=std::count(prefix.begin(),prefix.end(),6);
            const U lower=std::max<U>(L,210);
            if(fours%2 || std::binary_search(prefix.begin(),prefix.end(),12)) {
                const unsigned select=(fours%2==0);
                const int residue=7056-441*int(fours)-196*int(sixes)-144*7-49*int(select);
                if(residue<49 && Wide(49-residue)*lower>7056) return;
            }
            // The 7-free proper coefficients have squared norm six. The
            // alternatives are six signs, or two signs and a magnitude two.
            unsigned ones=0,twos=0,distinct=0;
            std::array<U,6> first{};
            U previous=0,witness_lcm=1;
            __int128 value=0;
            constexpr int coeff[7]={0,1,0,-2,2,0,-1};
            // inverse(2,7)=4 and inverse(5,7)=3 have cost at least nine.
            for(U w:prefix) if(w>7 && w%7==0) {
                if(ones+twos<6) first[ones+twos]=w;
                const int a=coeff[(w/7)%7];
                if (!a) return;
                witness_lcm=std::lcm(witness_lcm,w);
                value+=__int128(a)*(L/w);
                if(a==1 || a==-1) ++ones; else ++twos;
                distinct+=previous!=w; previous=w;
            }
            if(twos>1 || ones>(twos ? 2U : 6U)) return;
            bool three=fours<=10 && !(cross57 && fours>=8) && ones<=2 && ones+twos<=3;
            if(three && first[0] && first[0]<=(fixed_support ? 25000 : 500)) {
                const auto& family=column_roots(first[0],7,fixed_support,capped_witness);
                if(family.complete) {
                    three=false;
                    for(const auto& row:family.rows) {
                        const unsigned seen=ones+twos;
                        if(row.root%L || (cases.common&~row.support) ||
                           (seen<3 && row.witnesses[seen]<prefix.back()) ||
                           !std::equal(first.begin(),first.begin()+seen,row.witnesses.begin())) continue;
                        three=true;break;
                    }
                }
            }
            if(twos && !three) return;
            if(!three && fixed_support && first[0] && first[0]<=500) {
                const auto& family=column_roots(first[0],7,fixed_support,capped_witness,true);
                if(family.complete) {
                    bool found=false,unique=true;U only=0;
                    for(const auto& row:family.rows) {
                        if(row.root%L || !compatible_packet(row.root) || (fixed_root && fixed_root!=row.root) ||
                           (ones<6 && row.witnesses[ones]<prefix.back()) ||
                           !std::equal(first.begin(),first.begin()+ones,row.witnesses.begin())) continue;
                        bool valid=true;
                        if(cross57) {
                            unsigned common=0;
                            for(U z:row.witnesses) if(z%35==0) {
                                const U k=(z/35)%35;valid&=k==17 || k==18;++common;
                            }
                            valid&=common==1;
                        }
                        if(!valid) continue;
                        if(found && only!=row.root) unique=false;
                        found=true;only=row.root;
                    }
                    if(!found) return;
                    if(unique) fixed_root=only;
                }
            }
            unsigned required=three ? 3 : 6;
            seven_needed=std::max<int>(3-int(distinct),required-ones-twos);
            if(seven_needed>h) return;
            unsigned cost=ones+4*twos;
            seven_floor=std::max<U>(14,((prefix.back()-1)/7+1)*7);
            while(!coeff[(seven_floor/7)%7] || (!three && std::abs(coeff[(seven_floor/7)%7])!=1))
                seven_floor+=7;
            const Wide magnitude=value<0 ? -value : value;
            // The unseen coefficients have total l1 norm at most 6-cost.
            // The product (magnitude-1)*seven_floor can exceed 128 bits
            // only outside the fixed-width handoff, since magnitude<=6L.
            if(magnitude>1 && (magnitude-1)*seven_floor>Wide(6-cost)*L) return;
            if(cost==6) {
                if(magnitude!=1 || witness_lcm!=L) return;
                fixed_root=L;
            }
        }
    }
    if(fours==11 && fives==5 && std::count(prefix.begin(),prefix.end(),6)>=4) return;
    if (++nodes==0) throw std::overflow_error("node counter");
    if (bounded_pilot && nodes>400000) throw PilotLimit{};
    // L<=2^58 bounds every product used to calculate a denominator interval.
    Wide ll=Wide(L)*L;
    // For each p, at least (p-1) distinct labels must be active in a
    // simultaneous choice of representatives for all odd primes: the
    // union bound of the packet proposition gives Hall's condition.
    // Deduct an upper bound on prefix capacity prime by prime; the
    // residual requirements still satisfy Hall in the suffix.
    unsigned deficit=0;
    const U last=prefix.back();
    // At short suffixes use the exact tail and integral packet counts.
    const bool weight_packets=h>=6;
    const U scale=L<(U(1)<<48) ? 256 : 1;
    const Wide scaled_ll=scale*ll;
    const Wide last2=Wide(last)*last;
    Wide packet_mass=weight_packets ? std::min<Wide>(scale,1+(scaled_ll-1)/(36*last2)) : 0;
    // When five is isolated, a missing denominator five excludes all
    // denominators of five-adic valuation one. Once such a denominator
    // occurs, the mixed-five alternative is compulsory.
    bool require_mixed_five=false;
    if(!fives && L%5==0 && last>5)
        for(U z:prefix) require_mixed_five|=z%5==0 && z%25!=0;
    const Mask last_active=active_mask(prime_mask(last),fives);
    struct Choice {unsigned need=0; U floor=1;};
    Choice small[2][2]{};
    for (auto &pair:small) pair[1]=pair[0];
    const Mask required_support=cases.common|support;
    for (Mask known=required_support & ~Mask(1); known; known &= known-1) {
        unsigned i=__builtin_ctzll(known), prime=primes[i], step=(prime-1)/2;
        if (prime==3 || prime==5) {
            unsigned j=i-1;
            const bool must_pure=prime==5 && classes.isolated_five();
            const unsigned complete=2*step;
            U power=prime;
            while(power<last) power*=prime;
            // Only prime powers can occur in the pure case. The general
            // packet ceiling is also valid in each case separately.
            // Equal primitive-phase masses modulo p² at the greatest
            // active p-power denominator require p−1 equal labels.
            unsigned available=st.pure_max[i]>=complete ? complete :
                last_active==(Mask(1)<<i) ? std::min<unsigned>(complete,st.run) : 0;
            if(prime==5) {
                U rest=L,least=prime; unsigned height=0;
                while(rest%prime==0) {rest/=prime;++height;}
                for(unsigned k=1;k<(height+2)/2;++k) least*=prime;
                if(least>prime) {
                    bool paid=false;
                    for(U v=least;v<=last;) {
                        auto span=std::equal_range(prefix.begin(),prefix.end(),v);
                        if(unsigned(span.second-span.first)>=complete) {paid=true;break;}
                        if(v>last/prime) break;
                        v*=prime;
                    }
                    available=paid ? complete : last>=least && last_active==(Mask(1)<<i) ?
                              std::min<unsigned>(complete,st.run) : 0;
                    power=std::max(power,least);
                }
            }
            small[j][0]={complete-std::min(complete,available),power};
            small[j][1]={mixed_need(classes.paired_heights[j],complete),mixed_floor(last,prime)};
            if(prime==3 && cases.force_three_mixed && L>100000 && h<=13)
                small[j][1].need=separated_three(classes,last,h);
            if(must_pure) small[j][1]=small[j][0];
            if((prime==3 && cases.force_three_mixed) || (prime==5 && require_mixed_five))
                small[j][0]=small[j][1];
            if((must_pure && prime==5 && require_mixed_five) ||
               (prime==3 && cases.force_three_mixed && require_mixed_five)) return;
            if(prime==3 && require_mixed_five) small[j][1]=small[j][0];
            continue;
        }
        // Local Gauss valuation: the greatest active height K satisfies
        // 2K-1=v_p(M). It must already account for the observed height
        // of L. A completed smaller class cannot pay this requirement.
        unsigned height=0; U part=L,least=prime;
        while(part%prime==0) {part/=prime; ++height;}
        for(unsigned k=1;k<(height+2)/2;++k) least*=prime;
        bool paid=false;
        if(st.pure_max[i]>=2*step) for(U v=least;v<=last;) {
            auto span=std::equal_range(prefix.begin(),prefix.end(),v);
            if(unsigned(span.second-span.first)>=2*step) {paid=true;break;}
            if(v>last/prime) break;
            v*=prime;
        }
        if(paid) continue;
        const unsigned current=last>=least && last_active==(Mask(1)<<i) ? st.run : 0;
        const unsigned need=2*step-current;
        deficit+=need;
        if(deficit>unsigned(h)) return;
        if(weight_packets) {
            U multiple=((last-1)/prime+1)*prime;
            if(prime>=7) {
                multiple=std::max<U>(least,cases.force_large_power ? prime*prime : prime);
                if(minimal==i) multiple=std::max(multiple,least_B);
                while(multiple<last) multiple*=prime;
            }
            packet_mass+=need*(1+(scaled_ll-1)/(Wide(multiple)*multiple));
        }
    }
    // There are at least p pure-power denominators at an isolated
    // p>=5; the extra label need not be active. Combine only with
    // disjoint prime-power classes and the mixed ternary quota here.
    unsigned pure_slots=cases.force_three_mixed ? small[0][1].need : 0;
    for(Mask a=required_support&~Mask(3);a;a&=a-1) {
        const unsigned i=__builtin_ctzll(a), p=primes[i];
        if(p==5 && !classes.isolated_five()) continue;
        pure_slots+=p-std::min<unsigned>(p,st.pure_count[i]);
    }
    if(pure_slots>unsigned(h)) return;
    // On four support primes an odd nonunit omits at most one odd prime
    // and an even nonunit at most two. The four-row projection quotas
    // therefore require five or six proper hits, all inactive at odd primes.
    unsigned four_needed=0;
    if(fours>=11) {
        unsigned required=(L%15==0) ? 6 : 5;
        four_needed=required-std::min<unsigned>(required,classes.proper_four_hits);
    }
    // When the four-row degree is even, all odd degrees are proper hits
    // of dyadic twist. At least three nonunit degrees are odd. Count all
    // prefix multiples of the current two-part as possible odd witnesses.
    unsigned odd_needed=0;
    U odd_floor=1;
    if(fours && (L%8==0 || fours>=7)) {
        U two_part=std::max<U>(8,L & -L);
        unsigned possible=0;
        for(U z:prefix) possible+=z%two_part==0;
        odd_needed=3-std::min(3u,possible);
        odd_floor=((last-1)/two_part+1)*two_part;
    }
    unsigned joint=std::max(witnesses_needed,seven_needed);
    if(cross57) {
        bool shared=false;
        int miss5=0,miss7=0,miss_both=0;
        constexpr int coeff[7]={0,1,0,-2,2,0,-1};
        for(U b:prefix) {
            bool five=b>5 && b%5==0, seven=b>7 && b%7==0;
            shared|=five && seven;
            if(b%4==0) continue;
            if(five) ++miss5;
            if(seven) {
                int a=coeff[(b/7)%7];miss7+=a*a;
            }
            miss_both+=five && seven;
        }
        if(fours>=4 && L%8==0) {
            // Compress the two signed columns simultaneously onto the
            // four-row hits. Their norms are 4,6 and their dot is -1.
            int k=16-fours;
            int A=k*(4-miss5)-int(fours),B=k*(6-miss7)-int(fours);
            int C=k*(-1+miss_both)-int(fours);
            if(A<0 || B<0 || A*B<C*C) return;
        }
        unsigned both=unsigned(witnesses_needed+seven_needed);
        joint=std::max(joint,both-unsigned(!shared && both>0));
    }
    unsigned reserve=std::max({joint,four_needed,odd_needed});
    std::array<U,6> witness_floor{};
    std::array<U,primes.size()> height_floor{};
    constexpr U floor_cap=U(1)<<60; // exceeds the current prefix and is only a lower bound
    unsigned omission_total=odd_needed+unsigned(witnesses_needed)+unsigned(seven_needed);
    for(Mask u=required_support;u;u&=u-1) {
        const unsigned j=__builtin_ctzll(u); const U p=primes[j];
        unsigned v=0; U height=1,rest=L;
        while(rest%p==0) {rest/=p;height*=p;++v;}
        height=std::max(height,p);
        if(p==2 && fours>=7) height=std::max<U>(8,height);
        const bool isolated=p>=7 || (p==5 && classes.isolated_five());
        if(isolated && last>p && !std::binary_search(prefix.begin(),prefix.end(),p)) {
            U active=p*p;
            for(unsigned k=2;k<(v+2)/2;++k) active*=p;
            while(active<last) {
                auto span=std::equal_range(prefix.begin(),prefix.end(),active);
                if(unsigned(span.second-span.first)>=p-1) break;
                active*=p;
            }
            if(minimal==j) active=std::max(active,least_B);
            height=static_cast<U>(std::min<Wide>(Wide(active)*active/p,floor_cap));
            unsigned seen=0;
            for(U z:prefix) seen+=z%height==0;
            unsigned min_free=p==5 ? 3 : 2;
            const unsigned need=min_free-std::min(min_free,seen);
            omission_total+=need;
            reserve=std::max(reserve,need);
            U floor=((last-1)/height+1)*height;
            for(unsigned k=0;k<need;++k) witness_floor[k]=std::max(witness_floor[k],floor);
        }
        height_floor[j]=height;
    }
    // Each nonunit degree omits at most two primes on a four-prime support.
    unsigned hit_slots=0;
    // Select all possible p-free four-row hits. A proper hit has 4|b,
    // and for each p-free one the entire p-part of M divides b.
    // Count all compatible prefix labels optimistically.
    if(fours>=4 && L%8==0) {
        static constexpr unsigned quota[14][11]={
            {0,0,0,0,0,0,0,0,0,0,0},
            {1,1,1,1,1,1,1,1,1,1,1},
            {1,1,1,1,1,1,1,1,1,1,1},
            {1,1,1,1,1,1,1,1,1,1,1},
            {1,1,1,1,1,1,1,1,1,1,1},
            {1,1,1,1,1,1,1,1,1,1,1},
            {1,1,1,1,1,1,1,1,1,1,1},
            {1,1,1,1,1,1,1,1,1,1,1},
            {1,1,1,1,1,1,1,1,1,1,1},
            {2,2,2,2,2,2,2,2,2,2,2},
            {2,2,2,2,2,2,2,2,2,2,2},
            {3,3,3,3,3,3,3,3,3,3,3},
            {3,3,3,3,3,3,3,3,3,3,3},
            {4,4,4,4,4,4,4,4,4,4,4},
        };
        unsigned owed=0;
        for(Mask u=required_support&~Mask(1);u;u&=u-1) {
            unsigned j=__builtin_ctzll(u);
            U height=height_floor[j];unsigned have=0;
            for(U z:prefix) have+=z%4==0 && z%height==0;
            unsigned want=quota[fours][j-1];
            unsigned need=want-std::min(want,have);
            hit_slots=std::max(hit_slots,need);owed+=need;
            U step=std::min<Wide>(Wide(4)*height,floor_cap);
            U floor=((last-1)/step+1)*step;
            for(unsigned k=0;k<need;++k) witness_floor[k]=std::max(witness_floor[k],floor);
        }
        hit_slots=std::max(hit_slots,(owed+odd_needed+1)/2);
        reserve=std::max(reserve,hit_slots);
    }
    reserve=std::max(reserve,(omission_total+1)/2);
    // On a support of four primes, an omitted pair identifies exactly
    // the two primes of the nonunit degree. Distinct unordered pairs
    // therefore require distinct positions, even for different columns.
    if(fives || (L%7==0 && last>7)) {
        std::array<U,4> ps{},heights{};
        unsigned used=0;
        const bool five_column=fives==5 && closed5;
        const bool seven_column=last>7 && std::binary_search(prefix.begin(),prefix.end(),7);
        for(Mask a=required_support;a;a&=a-1) {
            const unsigned k=__builtin_ctzll(a);
            ps[used]=primes[k];heights[used++]=height_floor[k];
        }
        std::array<U,6> forced{};
        unsigned positions=0, for_odd=0, for_five=0, for_seven=0;
        for(unsigned i=0;i<used;++i) for(unsigned j=0;j<i;++j) {
            if(!(five_column && (ps[i]==5 || ps[j]==5)) &&
               !(seven_column && (ps[i]==7 || ps[j]==7))) continue;
            // A clipped bound exceeds every current denominator. In this
            // case use its lower floor directly instead of a residue loop.
            const bool large=Wide(heights[i])*heights[j]>floor_cap;
            const U step=large ? floor_cap : heights[i]*heights[j];
            if(std::any_of(prefix.begin(),prefix.end(),[&](U b){return b%step==0;})) continue;
            U floor=((last-1)/step+1)*step;
            if(!large && five_column && (ps[i]==5 || ps[j]==5))
                while(floor%25!=5 && floor%25!=20) floor+=step;
            if(cross57 && ps[i]*ps[j]==35) {
                while(floor/35%35!=17 && floor/35%35!=18) floor+=step;
            }
            forced[positions++]=floor;
            for_odd+=ps[i]==2 || ps[j]==2;
            for_five+=ps[i]==5 || ps[j]==5;
            for_seven+=ps[i]==7 || ps[j]==7;
        }
        std::sort(forced.rbegin(),forced.rend());
        for(unsigned i=0;i<positions;++i) witness_floor[i]=std::max(witness_floor[i],forced[i]);
        unsigned d2=odd_needed>for_odd ? odd_needed-for_odd : 0;
        unsigned d5=unsigned(witnesses_needed)>for_five ? witnesses_needed-for_five : 0;
        unsigned d7=unsigned(seven_needed)>for_seven ? seven_needed-for_seven : 0;
        reserve=std::max(reserve,positions+std::max({d2,d5,d7,(d2+d5+d7+1)/2}));
    }
    if(fives==5 && odd_needed && (fours>=7 || L%8==0)) {
        U height=std::max<U>(8,L & -L);
        unsigned seen=0;
        for(U b:prefix) seen+=b%5!=0 && b%height==0;
        unsigned outside=2-std::min(2u,seen);
        reserve=std::max(reserve,unsigned(witnesses_needed)+outside);
    }
    // Odd degrees, proper four-hits and the proper column witnesses
    // have non-primary denominators or even denominators, disjoint
    // from all pure odd-prime powers and from mixed ternary activity.
    unsigned outside_pure=std::max({hit_slots,odd_needed,four_needed,joint});
    if(fives==5 && odd_needed && (fours>=7 || L%8==0)) {
        U height=std::max<U>(8,L & -L);unsigned seen=0;
        for(U z:prefix) seen+=z%5!=0 && z%height==0;
        outside_pure=std::max(outside_pure,unsigned(witnesses_needed)+2-std::min(2u,seen));
    }
    if(pure_slots+outside_pure>unsigned(h)) return;
    // A closed pure p-power class of n labels supplies a trace deficit
    // pc/b² with (p-1)c<=n. Remove its mandatory n/b² from the
    // component norm. Open pure classes can supply at most p+1
    // units of net norm per packet; independent reservations bound
    // how many such packets fit. All comparisons round outwards.
    for(Mask a=required_support&~Mask(3);a;a&=a-1) {
        const unsigned i=__builtin_ctzll(a),p=primes[i];
        if((p==5 && !classes.isolated_five()) || std::binary_search(prefix.begin(),prefix.end(),p)) continue;
        const unsigned own=p-std::min<unsigned>(p,st.pure_count[i]);
        unsigned room=unsigned(h)-outside_pure-pure_slots+own;
        if(classes.pure_index==i) room+=classes.runs;
        U floor=p;
        while(floor<last && floor<norm_root) floor*=p;
        const auto capacity=std::int64_t((p+1)*(room/(p-1))*upper_mass(floor));
        if(std::int64_t(classes.mixed_cost[i])>classes.closed_net[i]+capacity) return;
    }
    if(first_witness_floor) for(unsigned i=0;i<unsigned(witnesses_needed);++i)
        witness_floor[i]=std::max(witness_floor[i],first_witness_floor);
    if(reserve) {
        deficit+=reserve;
        if(deficit>unsigned(h)) return;
        if(weight_packets) {
            U minimum=last;
            if(odd_needed==reserve) minimum=odd_floor;
            if(four_needed==reserve) minimum=std::max<U>(minimum,((last-1)/4+1)*4);
            if(seven_needed==int(reserve)) minimum=std::max(minimum,seven_floor);
            if(witnesses_needed==int(reserve)) {
                U v=std::max<U>({last,first_witness_floor,20});
                while((v%25!=5 && v%25!=20)) ++v;
                minimum=std::max(minimum,v);
            }
            for(unsigned i=0;i<reserve;++i) {
                const U floor=i<witness_floor.size()?std::max(minimum,witness_floor[i]):minimum;
                packet_mass+=1+(scaled_ll-1)/(Wide(floor)*floor);
            }
        }
    }
    const unsigned compulsory=std::min(small[0][0].need,small[0][1].need)+
                              std::min(small[1][0].need,small[1][1].need);
    if(deficit+compulsory>unsigned(h)) return;
    Mask future=required_support|Mask(1);
    for(unsigned j=1;j<num_primes && primes[j]<=unsigned(h)+1-deficit-compulsory;++j)
        future|=Mask(1)<<j;
    future&=cases.available;
    if(__builtin_popcountll(required_support)==4) future&=required_support;
    const int min_support=4;
    if(__builtin_popcountll(future)<min_support) return;
    // At least two odd primes occur even if two is itself available.
    // The smallest missing primes bound below any missing actual ones.
    int missing=min_support-1-__builtin_popcountll(required_support & ~Mask(1));
    for(Mask a=future & ~required_support & ~Mask(1);missing>0 && a;a&=a-1,--missing) {
        const unsigned prime=primes[__builtin_ctzll(a)], need=prime-1;
        deficit+=need;
        if(deficit>unsigned(h)) return;
        if(weight_packets) {
            U floor=std::max<U>(last,prime);
            if (prime>=7) {
                // The actual missing prime may exceed this minimal choice.
                // The next prime power above last is not monotone in its base.
                floor=UINT64_MAX;
                for(Mask opts=a;opts;opts&=opts-1) {
                    U t=primes[__builtin_ctzll(opts)];
                    U power=t;
                    while(power<last) power*=t;
                    floor=std::min(floor,power);
                }
            }
            packet_mass+=need*(1+(scaled_ll-1)/(Wide(floor)*floor));
        }
    }
    U high=square_root((36*h+1)*ll/(Wide(36)*E));
    // Bound the remaining norm using p pure labels, including the
    // extra inactive label. Proper prime-column witnesses and odd
    // degrees have already been proved disjoint from these reservations.
    if(weight_packets && pure_slots+outside_pure>0) {
        Wide charged=1+(scaled_ll-1)/(36*last2);
        auto book=[&](unsigned count,U floor) {
            if(count) charged+=count*(1+(scaled_ll-1)/(Wide(floor)*floor));
        };
        if(cases.force_three_mixed) book(small[0][1].need,small[0][1].floor);
        for(Mask a=required_support&~Mask(3);a;a&=a-1) {
            unsigned i=__builtin_ctzll(a),p=primes[i];
            if(p==5 && !classes.isolated_five()) continue;
            U floor=p;while(floor<last) floor*=p;
            if(minimal==i) floor=std::max(floor,least_B);
            book(p-std::min<unsigned>(p,st.pure_count[i]),floor);
        }
        U f5=std::max<U>({last,first_witness_floor,20});
        while(f5%25!=5 && f5%25!=20) ++f5;
        for(unsigned i=0;i<outside_pure;++i) {
            U floor=last;
            if(i<unsigned(witnesses_needed)) floor=std::max(floor,f5);
            if(i<unsigned(seven_needed)) floor=std::max(floor,seven_floor);
            if(i<odd_needed) floor=std::max(floor,odd_floor);
            if(i<four_needed) floor=std::max(floor,((last-1)/4+1)*4);
            book(1,floor);
        }
        if(scale*E>charged) high=std::min(high,square_root(
            (h-pure_slots-outside_pure)*scaled_ll/(scale*E-charged)));
    }
    U largest=0;
    for(auto a:small[0]) for(auto b:small[1]) {
        unsigned required=deficit+a.need+b.need;
        if(required>unsigned(h)) continue;
        U bound=high;
        if(weight_packets) {
            Wide mass=packet_mass;
            if(a.need) mass+=a.need*(1+(scaled_ll-1)/(Wide(a.floor)*a.floor));
            if(b.need) mass+=b.need*(1+(scaled_ll-1)/(Wide(b.floor)*b.floor));
            if(scale*E>mass) bound=std::min(bound,square_root((h-required)*scaled_ll/(scale*E-mass)));
        }
        largest=std::max(largest,bound);
    }
    high=largest;
    if(witnesses_needed) {
        U minimum=std::max<U>({last,first_witness_floor,20});
        while((minimum%25!=5 && minimum%25!=20)) ++minimum;
        // Every outstanding witness is at least minimum, and M is at
        // least six times each nonunit denominator.
        Wide denominator=36*Wide(minimum)*minimum;
        Wide reserved=1+((36*witnesses_needed+1)*scaled_ll-1)/denominator;
        if(scale*E>reserved) {
            U cap=square_root((h-witnesses_needed)*scaled_ll/(scale*E-reserved));
            high=std::min(high,cap);
        }
    }
    if(high<last) return;
    // Split every start at the same pilot depth, including long runs of fives.
    if (listing && (h==2 || prefix.size()>=split_depth)) {
        task();
        return;
    }
    U low=std::max(prefix.back(),square_root(ll/E)+1);
    if(low>high) return;
    if(h==2) {pair_tail(L,support,E,low,high,st); return;}
    Mask capped=0;
    for(Mask a=support&~Mask(3);a;a&=a-1) {
        U p=primes[__builtin_ctzll(a)];
        if(std::binary_search(prefix.begin(),prefix.end(),p)) capped|=a&-a;
    }
    auto next=[&](U b,Mask mask) {
        for(Mask a=mask&capped;a;a&=a-1) {
            U p=primes[__builtin_ctzll(a)];
            if(b%(p*p)==0) return;
        }
        // An absent isolated odd prime denominator of size >=7 forbids
        // every multiple of this prime having valuation one.
        for(Mask a=mask&~Mask(7);a;a&=a-1) {
            U p=primes[__builtin_ctzll(a)];
            if(b!=p && b%(p*p) && !std::binary_search(prefix.begin(),prefix.end(),p)) return;
        }
        if(fixed_root && fixed_root%b) return;
        if(minimal && b%primes[minimal]==0) {
            const U p=primes[minimal];U rest=b;while(rest%p==0) rest/=p;
            if(b/rest<least_B || (first_pure && b/rest>Wide(first_pure)*first_pure/p)) return;
        }
        U g=std::gcd(b,L), k=b/g, d=L/g;
        Mask eligible=active_mask(mask,fives || b==5);
        auto state=append(st,b==prefix.back(),eligible);
        if (!possible(state,eligible,support|mask,h-1)) return;
        if (L > lcm_limit/k) {
            prefix.push_back(b);
            if(width_room(L,b,h-1,state,classes.child(b),support|mask)) exact_task();
            prefix.pop_back();
            return;
        }
        Wide next_E=E*k*k-Wide(d)*d;
        prefix.push_back(b);
        dfs(L*k,support|mask,next_E,h-1,state,classes.child(b));
        prefix.pop_back();
    };
    if(prefix.size()<fixed.size()) {
        U b=fixed[prefix.size()];
        Mask mask;
        if(low<=b && b<=high && sieve(b,future,mask)) next(b,mask);
    } else if(high-low<=35) {
        for(U b=low;b<=high;++b) {Mask mask; if(sieve(b,future,mask)) next(b,mask);}
    } else if(__builtin_popcountll(future)==4 &&
              (ordered_smooth.size()<24 || ordered_smooth.count(future|(capped<<16)))) {
        const auto& table=smooth_table(future,capped);
        auto it=std::lower_bound(table.begin(),table.end(),low,
                                 [](const SmoothValue& a,U b){return a.value<b;});
        for(;it!=table.end() && it->value<=high;++it) next(it->value,it->mask);
    } else smooth(low,high,1,future,0,next);
}
static U integer(std::string_view text) {
    U value=0;
    auto end=text.data()+text.size();
    auto result=std::from_chars(text.data(),end,value);
    if (result.ec!=std::errc{} || result.ptr!=end || text.empty())
        throw std::invalid_argument("expected an unsigned decimal integer");
    return value;
}

static std::vector<U> parse_prefix(std::string_view text) {
    std::vector<U> bs;
    for (;;) {
        auto comma=text.find(',');
        bs.push_back(integer(text.substr(0,comma)));
        if(comma==std::string_view::npos) break;
        text.remove_prefix(comma+1);
    }
    if ((bs.front()<4 || bs.front()>square_root(rank_limit-2)) || bs.size()>rank_limit-3 ||
        !std::is_sorted(bs.begin(),bs.end())) throw std::invalid_argument("invalid prefix");
    return bs;
}

int main(int argc,char** argv) {
    if (argc<4 || argc>7 || std::string_view(argv[1])!="--rank") return 2;
    try {
        const auto rank=integer(argv[2]);
        if (rank<10 || rank>native_rank_limit) throw std::invalid_argument("native rank must be 10 through 44");
        // Once the first label is chosen, an unseen odd prime p needs p-1
        // of the other nonunit labels. A prime in the first denominator is smaller still.
        num_primes=std::lower_bound(primes.begin(),primes.end(),rank)-primes.begin();
        initialize_masks();
        rank_limit=static_cast<unsigned>(rank);
        argc-=2;
        argv+=2;
        auto atomic=::fpathconf(STDOUT_FILENO,_PC_PIPE_BUF);
        if(atomic>0) pipe_limit=static_cast<std::size_t>(atomic);
        std::string_view mode=argv[1];
        if (mode=="--tasks" || mode=="--tasks-cost") {
            if (argc>4) return 2;
            list_cost=(mode=="--tasks-cost");
            listing=true;
            if (argc>=3) {
                auto depth=integer(argv[2]);
                if (depth<8 || depth>32) return 2;
                split_depth=static_cast<unsigned>(depth);
            }
            if (argc==4) {
                auto group=integer(argv[3]);
                if (!group || group>1000000) return 2;
                frontier_per_job=group;
            }
        } else if(mode=="--refine-cost" && argc==5) {
            fixed=parse_prefix(argv[2]);
            const auto depth=integer(argv[3]), group=integer(argv[4]);
            if(depth<fixed.size() || depth>rank_limit-3 || !group || group>1000000) return 2;
            split_depth=static_cast<unsigned>(depth);
            frontier_per_job=group;
            listing=list_cost=bounded_pilot=true;
        } else if(mode=="--prefix" && argc==3) {
            fixed=parse_prefix(argv[2]);
        } else return 2;
        for (U b=4;b*b<=rank_limit-2;++b) if (fixed.empty() || b==fixed.front()) {
            Mask mask=prime_mask(b);
            prefix={b};
            dfs(b,mask,b*b-1,rank_limit-2,append(PacketState{},false,mask),SmallClasses{}.child(b));
        }
        if(listing) while(!planner_path.empty()) close_group();
        flush_records();
        std::cout.flush();
        if (!std::cout) return 3;
        if (listing) std::cerr<<"TASKS "<<jobs<<"\n";
        else std::cerr<<"DONE nodes="<<nodes<<" profiles="<<types<<" exact="<<handoffs<<"\n";
    } catch (const PilotLimit&) {
        std::cerr<<"PLAN_LIMIT\n";
        return 75;
    } catch (const std::exception& e) {
        std::cerr<<e.what()<<"\n";
        return 2;
    }
    return 0;
}
