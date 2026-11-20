// Nonsquare dimensions. Build with g++ -O3 -std=c++17 nonsquare.cpp -o nonsquare.
// ./nonsquare R [Q B top_multiplicity [--prefix comma,separated,denominators]]
// Fresh local nonsquare runner through rank 44; complete starts are supplied by nonsquare.py.
// A row is a candidate profile. An exact object transfers a whole subtree
// to nonsquare_exact.py before the fixed-width arithmetic could overflow.
#include <algorithm>
#include <array>
#include <limits>
#include <utility>
#include <cmath>
#include <cerrno>
#include <charconv>
#include <string_view>
#include <unistd.h>
#include <string>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

using U = std::uint64_t;
using W = __uint128_t;
using Mask = std::uint32_t;
std::vector<unsigned> primes;
std::vector<U> prefix, fixed_prefix;
std::vector<W> max_square;
int rank_limit, q, fiber, tile_prime;
unsigned effective_prime=0;
int selected_q=0, selected_b=0, selected_m=0;
Mask target, ramified;
U nodes=0, candidates=0, exact=0;
#ifndef NS_LCM_BITS
#define NS_LCM_BITS 58
#endif
#ifndef NS_EFFECTIVE_GCD
#define NS_EFFECTIVE_GCD 1
#endif
static_assert(NS_LCM_BITS >= 1 && NS_LCM_BITS <= 58);
// Complete the line before writing it to the driver's blocking pipe.
void record(std::string line) {
    line.push_back('\n');
    std::size_t done=0;
    while (done<line.size()) {
        auto n=::write(STDOUT_FILENO,line.data()+done,line.size()-done);
        if (n<0 && errno==EINTR) continue;
        if (n<=0) throw std::runtime_error("cannot write nonsquare output");
        done+=static_cast<std::size_t>(n);
    }
}

std::string list(const std::vector<U>& values) {
    std::string text="[";
    for (std::size_t i=0;i<values.size();++i) {
        if (i) text+=',';
        text+=std::to_string(values[i]);
    }
    return text+"]";
}

void handoff() {
    std::vector<U> support;
    for (unsigned i=0;i<primes.size();++i)
        if (target&(1u<<i)) support.push_back(primes[i]);
    // The handoff may be reached before the assigned worker prefix ends.
    const auto& continuation=prefix.size()<fixed_prefix.size() ? fixed_prefix : prefix;
    record("{\"exact\":"+list(continuation)+",\"primes\":"+list(support)+"}");
    if(++exact==0) throw std::overflow_error("handoff counter");
}

// Descendants only write at or beyond their parent's end index.  Passing the
// three indices by value restores the pending obligations on backtracking.
using Obligations=std::vector<std::pair<U,unsigned>>;
Obligations tile_due;
struct Tiles {
    int head=0, end=0, debt=0;
    bool close(U z,unsigned n,U next,int slots) {
        if (tile_prime && z%tile_prime) {
            unsigned held=0;
            if (head<end && tile_due[head].first==z) {
                held=tile_due[head++].second; debt-=held;
            }
            if (n<held) return false;
            unsigned remaining=n-held, need=0;
            // An odd denominator in the dyadic case comes in equal pairs.
            // At 3 the minimal number of mixed triples is remaining mod 3.
            // At 5 each mixed block uses two z and two 2z. Taking five
            // extra mixed blocks is equivalent to equal quintuples.
            if (tile_prime==2) {
                if (remaining%2) return false;
            } else if (tile_prime==3) {
                need=2*((q==3 && z==2) ? remaining : remaining%3);
            } else {
                unsigned mixed=(3*(remaining%5))%5;
                need=2*mixed;
                if (need>remaining) return false;
            }
            if (need) {
                // The native width check gives z < 2^64/6, also true at emission.
                tile_due[end++]={2*z,need}; debt+=need;
            }
        }
        return (head==end || tile_due[head].first>=next) && debt<=slots;
    }
};

Mask factors(U n) {
    Mask m=0;
    for (unsigned i=0;i<primes.size();++i) if (n%primes[i]==0) m|=1u<<i;
    return m;
}

U root(W n) {
    U s=static_cast<U>(std::sqrt(static_cast<long double>(n)));
    while (W(s)*s>n) --s;
    while (s<UINT64_MAX && W(s+1)*(s+1)<=n) ++s;
    return s;
}

template<class Visit>
void smooth(U lo,U hi,U a,Mask available,Visit& visit) {
    if (!available) {
        if (a>=lo) visit(a);
        return;
    }
    unsigned i=__builtin_ctz(available), p=primes[i];
    for (;;) {
        smooth(lo,hi,a,available^(1u<<i),visit);
        if (a>hi/p) break;
        a*=p;
    }
}

// Cache small smooth denominators in ascending order for each fixed support.
// Larger intervals still use complete, uncapped enumeration.
constexpr U cache_end=1000000;
std::vector<U> small_smooth;
std::vector<unsigned> support_primes, activity_primes;
struct Facts {
    std::array<U,5> prime_part{};
    std::array<unsigned,5> at_prime_part{}, principal_at_height{};
    std::array<U,5> at_height_positions{};
    bool pure_ramified=false;
    std::array<unsigned,4> odd_capacity{};
    std::array<unsigned,4> packets{}, doubles{};
    std::array<unsigned,16> eligible{};
    std::array<unsigned,2> effective_prime_count{}, effective_double_count{};
    bool effective_forced=false;
};
std::vector<Facts> levels;
std::vector<unsigned> group_demand;

void extend_facts(std::size_t depth,U z,unsigned run) {
    // The immutable parent slot makes a complete sibling rollback implicit.
    Facts& child=levels[depth];
    child=levels[depth-1];
    unsigned zmask=0;
    for (unsigned j=0;j<support_primes.size();++j) 
        if (z%support_primes[j]==0) zmask|=1u<<j;
    if (support_primes.size()==3 && (q==2 || q==3 || q==5 || q==7)) {
        U residue=z;
        while (residue%q==0) residue/=q;
        child.pure_ramified|=(residue==1);
    }
    for (unsigned i=0;i<support_primes.size();++i) {
        unsigned p=support_primes[i];
        if (!(zmask&(1u<<i))) continue;
        U part=p, rest=z/p;
        while (rest%p==0) {part*=p;rest/=p;}
        if (part>child.prime_part[i]) {
            child.prime_part[i]=part;
            child.at_prime_part[i]=1;
            child.principal_at_height[i]=zmask;
            child.at_height_positions[i]=U(1)<<(depth-1);
        } else if (part==child.prime_part[i]) {
            ++child.at_prime_part[i];
            child.principal_at_height[i]|=zmask;
            child.at_height_positions[i]|=U(1)<<(depth-1);
        }
    }
    if (effective_prime) {
        const unsigned unfilled=unsigned(rank_limit-1-depth);
        // At fixed support {2,3,p}, the prime-free trace theorem forbids
        // a full equal ternary block at denominator 2 or p. Each such
        // label must be the larger member of (n/2,n/2,n), consuming two
        // double-denominator partners. A larger partner for denominator
        // 2 would have denominator 1, whose codegree 3 is squarefree.
        // Every unfilled position is optimistically available until the
        // partner denominator is passed; this is safe for open runs.
        for(unsigned j=0;j<2;++j) {
            const U p=j==0 ? 2:effective_prime;
            child.effective_prime_count[j]+=z==p;
            child.effective_double_count[j]+=z==2*p;
            const unsigned possible=child.effective_double_count[j]
                +(z<=2*p ? unfilled:0);
            child.effective_forced|=2*child.effective_prime_count[j]>possible;
        }
    }
    unsigned active=0;
    for (unsigned i=0;i<activity_primes.size();++i) {
        unsigned p=activity_primes[i];
        if (z%p) continue;
        active|=1u<<i;
        unsigned step=(p-1)/2;
        if ((q==3 || q==5) && z%2) {
            unsigned block=(z%q==0) ? step : q*step;
            if(run%block==0) child.odd_capacity[i]+=block;
        }
        child.packets[i]+=run%step==0;
        child.doubles[i]+=run%(2*step)==0;
    }
    for (unsigned mask=1;mask<group_demand.size();++mask)
        if (mask&active) ++child.eligible[mask];
}

void prepare_support() {
    small_smooth.clear();
    auto add=[](U z) {small_smooth.push_back(z);};
    smooth(1,cache_end,1,target,add);
    std::sort(small_smooth.begin(),small_smooth.end());
    support_primes.clear(); activity_primes.clear();
    for (unsigned i=0;i<primes.size();++i) if (target&(Mask(1)<<i)) {
        unsigned p=primes[i];
        support_primes.push_back(p);
        if (p>2 && q%p) activity_primes.push_back(p);
    }
    effective_prime=NS_EFFECTIVE_GCD && q==3 && support_primes.size()==3
        && support_primes[0]==2 && support_primes[1]==3 ? support_primes[2]:0;
    // Refuse supports exceeding the exact fixed-size packet counters.
    if (support_primes.size()>5 || activity_primes.size()>4)
        throw std::logic_error("prime support exceeds packet counters");
    levels.assign(rank_limit,Facts{});
    group_demand.assign(1u<<activity_primes.size(),0);
    for (unsigned mask=1;mask<group_demand.size();++mask) {
        unsigned i=__builtin_ctz(mask);
        group_demand[mask]=group_demand[mask^(1u<<i)]+fiber*(activity_primes[i]-1);
    }
}

template<class Visit>
void visit_denominators(U lo,U hi,Visit& visit) {
    if (lo<=cache_end) {
        auto it=std::lower_bound(small_smooth.begin(),small_smooth.end(),lo);
        auto end=std::upper_bound(it,small_smooth.end(),hi);
        for (;it!=end;++it) visit(*it);
    }
    if (hi>cache_end) smooth(std::max(lo,cache_end+1),hi,1,target,visit);
}

// Keep integer upper bounds on the remaining reciprocal-square mass. If
// A future slots must carry multiples of v, those slots can't all use lo.
// ceil((scale*L^2)/z^2) rounds upwards in every term. The scaled products
// fit because search caps (36*h+1)*L^2 and the scale itself.
W ceiling_ratio(W n,W d) {return 1+(n-1)/d;}

bool mass_possible(W ll,W E,U lo,int h,const Facts& facts) {
    W scale=std::min<W>(W(1)<<32,(~W(0)/ll)/(h+1));
    W scaled=ll*scale, base=ceiling_ratio(scaled,W(lo)*lo);
    W unit=ceiling_ratio(scaled,W(36)*lo*lo);
    W wanted=E*scale;
    if (W(h)*base+unit<wanted) return false;
    auto affordable=[&](int missing,U threshold) {
        if (missing<=0) return true;
        if (missing>h) return false;
        W limited=ceiling_ratio(scaled,W(threshold)*threshold);
        return W(h-missing)*base+W(missing)*limited+unit>=wanted;
    };
    // A p-free degree requires a denominator divisible by the full p-part
    // of M, hence by the p-part already forced by L. There are at least
    // two such nonunit labels; mod 4 or mod 5 gives three at 2 and 5.
    for (unsigned i=0;i<support_primes.size();++i) {
        unsigned p=support_primes[i];
        U part=std::max<U>(p,facts.prime_part[i]);
        int missing=((p==2 || p==5)?3:2)-int(facts.at_prime_part[i]);
        U threshold=((lo-1)/part+1)*part;
        if (!affordable(missing,threshold)) return false;
    }
    // At exactly three primes, nonunit degrees free of different primes
    // are disjoint. Each t-free family must include a member of every
    // principal p-block, hence a denominator divisible also by p.
    if (support_primes.size()==3) {
        int places=0, missing_ram=0;
        W gain=0;
        for (unsigned j=0;j<3;++j) {
            const U t=support_primes[j];
            const U part=std::max<U>(facts.prime_part[j],t*t);
            const bool prefix_t_free=facts.prime_part[j]>=t*t;
            const unsigned seen=prefix_t_free ? facts.at_prime_part[j] : 0;
            const unsigned mask=prefix_t_free ? facts.principal_at_height[j] : 0;
            const unsigned missing_mask=(7u ^ (1u<<j)) & ~mask;
            const int required=t==2||t==5 ? 3 : 2;
            int free_count=std::max(0,required-int(seen));
            // If the height stays fixed, all seen labels are t-free.
            // If it grows, at least required future labels are t-free.
            if (seen && (t==2||t==3)) {
                const int period=t==2 ? 4 : 3;
                int extras=(required+period-(seen%period))%period;
                if (int(seen)+extras<required) extras+=period;
                free_count=std::min(extras,required);
            }
            if (missing_mask) free_count=std::max(free_count,1);
            if (int(t)==q) missing_ram=free_count;
            places+=free_count;
            if (places>h) return false;
            if (free_count) {
                const U rounded=((lo-1)/part+1)*part;
                const W unit_gain=base-ceiling_ratio(scaled,W(rounded)*rounded);
                gain+=free_count*unit_gain;
                if (missing_mask) {
                    U strong=rounded;
                    for (unsigned k=0;k<3;++k) if (missing_mask&(1u<<k)) {
                        const W mult=W(part)*support_primes[k];
                        const W next=((lo-1)/mult+1)*mult;
                        // A neighboring prime supplies a further divisor.
                        // Clip to lo on overflow for a conservative bound.
                        if(next<=UINT64_MAX) strong=std::max(strong,U(next));
                    }
                    gain+=base-ceiling_ratio(scaled,W(strong)*strong)-unit_gain;
                }
            }
        }
        if (W(h)*base+unit-gain<wanted) return false;
        if (facts.pure_ramified && missing_ram) {
            unsigned lastmask=group_demand.size()-1;
            int act=std::max(0,int(group_demand[lastmask])-int(facts.eligible[lastmask]));
            if (act+missing_ram>h) return false;
        }
    }
    // With three support primes, the sets of t-free nonunits are pairwise
    // disjoint. Decide for each prefix height whether it is the final one
    // or must rise. At a fixed height every position at that height is
    // already t-free, including its exact residue-class contribution.
    if (support_primes.size()==3 && h<=9) {
        struct Choice { U positions; unsigned needed; };
        Choice choice[3][2]{};
        for (unsigned j=0;j<3;++j) {
            const unsigned t=support_primes[j];
            const unsigned minimum=(t==2 || t==5) ? 3 : 2;
            choice[j][1]={0,minimum};
            unsigned seen=facts.at_prime_part[j];
            if (facts.prime_part[j]<U(t)*t) {
                choice[j][0]={0,unsigned(h)+1};
                continue;
            }
            unsigned need=seen<minimum ? minimum-seen : 0;
            unsigned missing=(7u ^ (1u<<j)) & ~facts.principal_at_height[j];
            if (missing) need=std::max(1U,need);
            if (t==2) need+=(3u-seen-need)&3u;
            if (t==3) need+=(5u-(seen+need)%3)%3;
            choice[j][0]={facts.at_height_positions[j],need};
        }
        bool ok=false;
        for (auto x:choice[0]) for (auto y:choice[1]) {
            if (x.positions & y.positions) continue;
            for (auto z:choice[2])
                if (!((x.positions|y.positions)&z.positions) &&
                      x.needed+y.needed+z.needed<=unsigned(h)) ok=true;
        }
        if (!ok) return false;
    }
    // The union of the active labels for any subset J of the odd
    // unramified primes costs fiber*sum(p-1) places by the joint packet bound.
    std::array<U,4> next{};
    for (unsigned j=0;j<activity_primes.size();++j) {
        unsigned p=activity_primes[j];
        next[j]=((lo-1)/p+1)*p;
    }
    for (unsigned mask=1;mask<group_demand.size();++mask) {
        int missing=int(group_demand[mask])-int(facts.eligible[mask]);
        if (missing<=0) continue;
        U threshold=UINT64_MAX;
        for (unsigned bits=mask;bits;bits&=bits-1)
            threshold=std::min(threshold,next[__builtin_ctz(bits)]);
        if (!affordable(missing,threshold)) return false;
    }
    for (unsigned i=0;i<activity_primes.size();++i) {
        const unsigned p=activity_primes[i];
        if ((q==3 || q==5) && p>=5 && facts.eligible[1u<<i]+h<unsigned(2*fiber*(p-1))) {
            unsigned step=(p-1)/2;
            unsigned block=(prefix.back()%q==0)?step:q*step;
            // Count every open copy conservatively: unlike enough_packets,
            // this function does not receive its run length.
            int available=facts.odd_capacity[i]+
                    ((prefix.back()%2 && prefix.back()%p==0)?block-1:0);
            int missing=std::max(0,fiber*int(p-1)-available);
            U mult=((lo-1)/p+1);
            if (!(mult%2)) ++mult;
            if (!affordable(missing,mult*p)) return false;
        }
    }
    return true;
}

// Necessary joint activity counts and a stronger mass test at a fixed top
// multiplicity. For any threshold t use the primes whose least eligible
// future denominator is >=t. Maximal labels at B can already be eligible.
bool seed_mass_possible(unsigned b,int m) {
    std::vector<U> thresholds={b+1};
    U common=6;
    for (unsigned p:activity_primes) thresholds.push_back(((b/p)+1)*p);
    std::sort(thresholds.begin(),thresholds.end());
    thresholds.erase(std::unique(thresholds.begin(),thresholds.end()),thresholds.end());
    for (U t:thresholds) common=std::lcm(common,t);
    W scale=W(common)*common;
    W upper=scale/36+W(m)*scale;
    int assigned=0;
    for (auto it=thresholds.rbegin();it!=thresholds.rend();++it) {
        U t=*it;
        int required=0;
        for (unsigned mask=1;mask<(1u<<activity_primes.size());++mask) {
            bool qualifies=true, maximal_eligible=false;
            int demand=0;
            for (unsigned i=0;i<activity_primes.size();++i) if (mask&(1u<<i)) {
                unsigned p=activity_primes[i];
                if ((b/p+1)*p<t) {qualifies=false;break;}
                demand+=fiber*(p-1);
                if (b%p==0) maximal_eligible=true;
            }
            if (qualifies) required=std::max(required,demand-(maximal_eligible?m:0));
        }
        if (required>rank_limit-1-m) return false;
        if (t==b+1) required=rank_limit-1-m;
        required=std::max(required,assigned);
        if (required>rank_limit-1-m) return false;
        U normalized=b*(common/t);
        upper+=W(required-assigned)*normalized*normalized;
        assigned=required;
    }
    return W(q)*b*b*scale<=upper;
}

// Closed runs contribute complete packets. Only the last run can grow.
bool enough_packets(int slots, U last, unsigned run, const Facts& facts) {
    int debt=0;
    for (unsigned i=0;i<activity_primes.size();++i) {
        unsigned p=activity_primes[i];
        int step=(p-1)/2, tail=last%p==0 ? run : 0;
        int packets=facts.packets[i], doubles=facts.doubles[i];
        debt+=std::max(0,int(fiber*(p-1))-step*packets-tail%step);
        if (debt>slots) return false;
        if((q==3 || q==5) && p>=5 && facts.eligible[1u<<i]+slots < unsigned(2*fiber*(p-1))) {
            // No p-primary order above p fits. A mixed p/2-support
            // would give all p-1 primitive roots with trivial dyadic phase
            // and again at one fixed nontrivial dyadic phase. Each needs
            // fiber labels. Thus every p-active denominator is odd.
            // At a defect-one ramified label, the mixed block would have
            // an even denominator; only equal q-tuples can be p-active.
            int block=(last%q==0)?step:q*step;
            int open=(last%2 && last%p==0)?run%block:0;
            if(int(facts.odd_capacity[i])+open+slots<fiber*(p-1)) return false;
        }
        if (q%2==0 && packets+(slots+tail%step)/step
                +doubles+(slots+tail%(2*step))/(2*step)<6) return false;
    }
    return true;
}

// The effective-row projection theorem: in D=2M^2, any prime
// denominator p>=5 through rank59 has p||M. Prefix L divides final M.
// The two optional denominators are needed at the terminal join.
bool effective_prime_heights(U M,U b=0,U c=0) {
    if(q!=2) return true;
    for(unsigned p:primes) if(p>=5 && M%(U(p)*p)==0 &&
            (std::binary_search(prefix.begin(),prefix.end(),U(p)) || b==p || c==p))
        return false;
    return true;
}

void emit(U L,U k,U d,U e,Tiles pending,unsigned run) {
    if (std::gcd(k,d)!=1 || std::gcd(k,e)!=1 || k%2==0 || k%5==0) return;
    if (L>UINT64_MAX/k) throw std::overflow_error("terminal scale");
    U M=k*L;
    if (factors(M)!=target) return;
    U b=M/e, c=M/d;
    if(!effective_prime_heights(M,b,c)) return;
    // Consume the last two denominators, keeping the current run open.
    U z=prefix.back();
    if (b!=z) {
        if (!pending.close(z,run,b,2)) return;
        z=b; run=0;
    }
    ++run;
    const auto depth=prefix.size();
    extend_facts(depth+1,b,run);
    if (c!=z) {
        if (!pending.close(z,run,c,1)) return;
        z=c; run=0;
    }
    ++run;
    if (!pending.close(z,run,UINT64_MAX,0)) return;
    extend_facts(depth+2,c,run);
    const Facts& full=levels[depth+2];
    if(full.effective_forced) return;
    if (!enough_packets(0,c,run,full)) return;
    for (unsigned mask=1;mask<group_demand.size();++mask)
        if (full.eligible[mask]<group_demand[mask]) return;
    std::vector<U> degrees={1,d,e};
    for (auto it=prefix.rbegin();it!=prefix.rend();++it) degrees.push_back(M / *it);
    auto m=std::count(degrees.begin()+1,degrees.end(),degrees.back());
    if (selected_m && m!=selected_m) return;
    for (unsigned j=1;j<degrees.size();++j)
        if (__builtin_popcount(factors(degrees[j]))<2) return;
    for (unsigned i=0;i<primes.size();++i) if (target&(1u<<i)) {
        unsigned p=primes[i], free=0;
        for (unsigned j=1;j<degrees.size();++j) free+=(degrees[j]%p!=0);
        if (free<((p==2 || p==5)?3u:2u)) return;
    }
    record(list(degrees));
    if(++candidates==0) throw std::overflow_error("profile counter");
}

void search(U L,W E,int h,unsigned run,const Tiles& pending) {
    if(++nodes==0) throw std::overflow_error("node counter");
    if(!effective_prime_heights(L)) return;
    if(levels[prefix.size()].effective_forced) return;
    if (!enough_packets(h,prefix.back(),run,levels[prefix.size()])) return;
    if (prefix.back()==prefix.front()) {
        U f=q*prefix.front()*prefix.front();
        for (unsigned p:primes) if (p%2 && f%p==0 && f%(p*p)
                                      && prefix.size()>f-1-(p-1)/2) return;
    }
    if (L>(U(1)<<NS_LCM_BITS)) {handoff(); return;}
    W ll=W(L)*L;
    if (ll>max_square[h]) {handoff(); return;}
    if (E>(36*h+1)*ll/36) return;
    U lo=std::max(prefix.back(),root(ll/E)+1);
    U hi=root((36*h+1)*ll/(36*E));
    if (selected_m && prefix.front()==prefix.back()) {
        if (prefix.size()<static_cast<unsigned>(selected_m)) hi=std::min(hi,prefix.front());
        else lo=std::max(lo,prefix.front()+1);
    }
    if (pending.head!=pending.end) {
        auto [z,need]=tile_due[pending.head];
        if (z>prefix.back()) {
            hi=std::min(hi,z);
            if (pending.debt==h) lo=std::max(lo,z);
        } else if (z==prefix.back() && run<need) {
            hi=std::min(hi,z);
            if (need-run>unsigned(h)) return;
        }
    }
    if (lo>hi) return;
    // A degree at f=18 must be prime to three.
    if(q==2 && prefix.back()==3 && prefix.size()>2 && std::count(prefix.begin(),prefix.end(),3)>(prefix.front()==2 ? 6 : 10)) return;
    if(q==2 && L%9==0 && std::binary_search(prefix.begin(),prefix.end(),3)) return;
    if (!mass_possible(ll,E,lo,h,levels[prefix.size()])) return;
    // At two remaining places k is odd, hence E = 1+d*d+e*e mod 8.
    if (h==2 && (E%8==0 || E%8==4 || E%8==7)) return;
    // At the terminal pair M=kL, k<=hi. Hand off before producing any
    // native result if this bound could exceed the representation of M.
    if (h==2 && hi>UINT64_MAX/L) {handoff(); return;}
    auto visit=[&](U b) {
        U g=std::gcd(b,L), k=b/g, e=L/g;
        // E*k*k <= E*b*b <= (h+1/36)L^2, strictly below 2^128.
        W rest=E*k*k-W(e)*e;
        if (h==2) {
            if (e<6 || rest<37 || rest>W(e)*e+1 || k%2==0 || k%5==0) return;
            U d=root(rest-1);
            if (W(d)*d==rest-1 && L%d==0) emit(L,k,d,e,pending,run);
        } else {
            bool same=(b==prefix.back());
            Tiles next=pending;
            if (!same && !next.close(prefix.back(),run,b,h)) return;
            prefix.push_back(b);
            extend_facts(prefix.size(),b,same?run+1:1);
            if (L>UINT64_MAX/k) handoff();
            else search(L*k,rest,h-1,same?run+1:1,next);
            prefix.pop_back();
        }
    };
    if (prefix.size()<fixed_prefix.size()) {
        U b=fixed_prefix[prefix.size()];
        if (lo<=b && b<=hi) visit(b);  // its support was checked at the root
    } else visit_denominators(lo,hi,visit);
}

// The packet bound pays p-1 for each odd unramified prime and zero for 2.
// Traverse only feasible supports, in increasing mask order.
void supports(int i,Mask support,int budget,int minimum,unsigned b,int f) {
    if (i<0) {
        if (__builtin_popcount(support)>=minimum) {
            target=support;
            for (U z:fixed_prefix) {
                U rem=z;
                for (unsigned j=0;j<primes.size();++j) if (target&(Mask(1)<<j))
                    while (rem%primes[j]==0) rem/=primes[j];
                if (rem!=1) return;
            }
            prepare_support();
            if (selected_m && !seed_mass_possible(b,selected_m)) return;
            prefix={b};
            extend_facts(1,b,1);
            search(b,f-1,rank_limit-2,1,Tiles{});
        }
        return;
    }
    supports(i-1,support,budget,minimum,b,f);
    Mask bit=Mask(1)<<i;
    int cost=i ? primes[i]-1 : 0;
    if (!(support&bit) && cost<=budget)
        supports(i-1,support|bit,budget-cost,minimum,b,f);
}

int integer(std::string_view token) {
    int value=0;
    auto end=token.data()+token.size();
    auto parsed=std::from_chars(token.data(),end,value);
    if (parsed.ec!=std::errc{} || parsed.ptr!=end || value<0)
        throw std::invalid_argument("expected a nonnegative decimal integer");
    return value;
}

int main(int argc,char** argv) try {
    if (argc!=5 && argc!=7) return 2;
    rank_limit=integer(argv[1]);
    if (rank_limit<1 || rank_limit>44) return 2;
    for (int p=2;p<=rank_limit;++p) {
        if (std::none_of(primes.begin(),primes.end(),[p](unsigned d){return p%d==0;})) {
            if (primes.size()==std::numeric_limits<Mask>::digits) {
                std::cerr<<"prime support exceeds native mask, use exact replay\n";
                return 75;
            }
            primes.push_back(p);
        }
    }
    tile_due.resize(rank_limit);
    max_square.resize(rank_limit);
    for (int h=0;h<rank_limit;++h) max_square[h]=~W(0)/(36*h+1);
    if (argc>=5) {
        selected_q=integer(argv[2]); selected_b=integer(argv[3]);
        selected_m=integer(argv[4]);
        if (selected_q<2 || selected_q>=rank_limit || selected_b<1 ||
            selected_m<0 || selected_m>=rank_limit) return 2;
    }
    if (argc==7) {
        if (std::string_view(argv[5])!="--prefix") return 2;
        std::string_view rest=argv[6];
        for (;;) {
            auto comma=rest.find(',');
            auto token=rest.substr(0,comma);
            U z=0;
            auto last=token.data()+token.size();
            auto parsed=std::from_chars(token.data(),last,z);
            if (parsed.ec!=std::errc{} || parsed.ptr!=last || z==0) return 2;
            fixed_prefix.push_back(z);
            if (comma==std::string_view::npos) break;
            rest.remove_prefix(comma+1);
        }
        if (fixed_prefix.front()!=U(selected_b) ||
            fixed_prefix.size()>unsigned(rank_limit-3) ||
            !std::is_sorted(fixed_prefix.begin(),fixed_prefix.end())) return 2;
    }
    // Work partitions cover every possible qB^2, independently of the
    // rank-34 five-start reduction. The odd-dimension rank theorem forces 2|M.
    q=selected_q;
    unsigned b=selected_b;
    const int f=q*b*b;
    if(f>rank_limit-2 || b<2 || f==8)
        throw std::invalid_argument("inadmissible starting codegree (codegree eight requires rank at least 46)");
    for(unsigned p:primes) if(q%(p*p)==0)
        throw std::invalid_argument("the dimension kernel must be squarefree");
    ramified=factors(q);
    int product=1;
    for(unsigned p:primes) if(p>2 && q%p==0) product*=(p+1)/2;
    fiber=(q%2 ? 1:2)*product;
    for(unsigned p:primes) if(p>2 && q%p==0) {
        int factor=(p+1)/2;
        fiber=std::max(fiber,(product/factor)*(factor+1));
    }
    // Applying one ramified tiling suffices for prefix pruning; the Python
    // terminal test applies all available tilings for composite q as well.
    tile_prime=q%2==0 ? 2 : q%3==0 ? 3 : q%5==0 ? 5 : 0;
    const int inactive=q==3 ? std::max(4,fiber) : fiber;
    int budget=(rank_limit-inactive)/fiber;
    Mask must=ramified|factors(b)|Mask(1);
    for(unsigned i=1;i<primes.size();++i)
        if((must&~ramified)&(Mask(1)<<i)) budget-=primes[i]-1;
    if((f==18 && selected_m>10) || (f==12 && selected_m>7)) budget=-1;
    if(budget>=0) supports(int(primes.size())-1,must,budget,(f==8 || f==18) ? 4 : 3,b,f);
    std::cerr<<"rank "<<rank_limit;
    if (selected_q) std::cerr<<", q "<<selected_q<<", B "<<selected_b<<", multiplicity "<<selected_m;
    std::cerr<<", nodes "<<nodes<<", candidates "<<candidates<<", exact "<<exact<<"\n";
} catch (const std::exception& e) {
    std::cerr<<"nonsquare generator failed: "<<e.what()<<'\n';
    return 1;
}
