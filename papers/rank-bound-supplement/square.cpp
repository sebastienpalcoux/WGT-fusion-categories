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
using U = std::uint64_t;
using Wide = __uint128_t;
using Mask = std::uint64_t;
// The 64-bit position mask covers ranks through 64. At most 18 primes
// occur below 64; unfinished branches at larger ranks use Python integers.
static constexpr unsigned native_rank_limit=64;
static constexpr std::array<unsigned,18> primes={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61};
using Counts = std::array<unsigned char,primes.size()>;
static unsigned num_primes;
static unsigned max_fours;
static std::vector<U> prefix;
static U nodes, types;
static unsigned rank_limit;
static bool listing=false;
static unsigned split_depth=22;
static std::size_t frontier_per_job=512;
static bool list_cost=false;
static std::vector<std::vector<U>> frontier;
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

static void task() {
    frontier.push_back(prefix);
}

// Each range has a common prefix of length depth. Grouping complete
// descendant sets gives disjoint jobs even when branches terminate early.
static void emit_groups(std::size_t first, std::size_t end, unsigned depth) {
    if (depth && end-first<=frontier_per_job) {
        auto& example=frontier[first];
        if(list_cost) std::cout << '[';
        write_prefix(std::vector<U>(example.begin(),example.begin()+depth));
        if(list_cost) std::cout << ',' << end-first << ']';
        std::cout << '\n';
        ++jobs;
        return;
    }
    for(std::size_t i=first;i<end;) {
        auto j=i+1;
        while(j<end && frontier[j][depth]==frontier[i][depth]) ++j;
        emit_groups(i,j,depth+1);
        i=j;
    }
}

static void exact_task() {
    if (listing) {
        task();
        return;
    }
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
        if (st.count[i]+unsigned(h) < 2*m || (maximum<2*m &&
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
    for (unsigned i=0; i<num_primes; ++i) if (x%primes[i]==0)
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
        for(unsigned i=0;i<num_primes;++i) if(real_support&(Mask(1)<<i)){
            int free=0;
            for(unsigned j=1;j<degrees.size();++j) {
                if(degrees[j]%primes[i]) ++free;
            }
            if(free<2) return;
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
    U lo=std::max(prefix.back(),square_root(ll/E)+1);
    U hi=square_root((73*ll)/(36*E));
    if(lo>hi) return;
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
    if(hi-lo<=32) {from_interval(); return;}
    const auto& ds=nonunit_divisors(L,support);
    if(ds.empty()) return;
    if(hi-lo<=2*ds.size()) {from_interval(); return;}

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
// Highest power attained in the prefix and counts of the entries at that power.
struct Heights {
    std::array<U,primes.size()> power{};
    Counts count{}, flags{};
    std::array<Mask,primes.size()> positions{};
    bool has_five=false;
    Heights child(U b, Mask mask, unsigned at) const {
        Heights next=*this;
        next.has_five|=b==5;
        const unsigned witnessed=(b%4==0 ? 1U : 0U) | (b%25==0 ? 2U : 0U);
        while(mask) {
            unsigned i=__builtin_ctzll(mask); mask&=mask-1;
            U prime=primes[i], power=1, rest=b;
            do { power*=prime; rest/=prime; } while(rest%prime==0);
            if(power<next.power[i]) continue;
            if(power>next.power[i]) {
                next.power[i]=power; next.count[i]=0; next.flags[i]=0; next.positions[i]=0;
            }
            ++next.count[i]; next.flags[i]|=witnessed; next.positions[i]|=Mask(1)<<at;
        }
        return next;
    }
};
static void dfs(U L, Mask support, Wide E, int h, const PacketState& st, const Heights& heights) {
    int fours=std::upper_bound(prefix.begin(),prefix.end(),4)-prefix.begin();
    if (unsigned(fours)>max_fours) return;
    // The mixed-degree Gauss-sum exclusions in the coverage notes.
    if ((fours==13 && prefix.size()>=17 && (prefix[16]==5 ||
            (prefix[15]==5 && prefix[16]==6))) ||
        (fours==12 && prefix.size()>=18 && prefix[17]==5) ||
        (fours==10 && prefix.size()>=19 && prefix[18]==5)) return;
    if (++nodes==0) throw std::overflow_error("node counter");
    // L<=2^58 bounds every product used to calculate a denominator interval.
    Wide ll=Wide(L)*L;
    if (E>ll*(36*h+1)/36) return;
    unsigned deficit=0;
    for (Mask known=support & ~Mask(1); known; known &= known-1) {
        unsigned i=__builtin_ctzll(known), step=(primes[i]-1)/2;
        unsigned cap=step*st.packets[i];
        if (prefix.back()%primes[i]==0) cap+=st.run%step;
        if (cap<2*step) deficit+=2*step-cap;
    }
    if (deficit>unsigned(h)) return;
    Mask future=support|Mask(1);
    for(unsigned j=1;j<num_primes && primes[j]<=unsigned(h)+1-deficit;++j)
        future|=Mask(1)<<j;
    if(__builtin_popcountll(future)<3) return;
    U high=square_root((36*h+1)*ll/(Wide(36)*E));
    if(__builtin_popcountll(future)==3) {
        bool fives=heights.has_five;
        if(h==2 && ( (fours>=4 && L%8) || (fives && L%125) )) return;
        if(h==3 && fours>=4 && L%8 && ((L%4) || E%2==0)) return;
        // Every three-prime nonunit can be free of at most one support prime.
        // At two and five at least three nonunits are p-free, by reducing D modulo 4 or 5.
        struct Option { unsigned need; Wide cost; Mask known; };
        Option options[3][2]{};
        unsigned index=0;
        for(Mask a=future;a;a&=a-1,++index) {
            unsigned j=__builtin_ctzll(a); U p=primes[j];
            U seen=heights.power[j], floor=p;
            if (p==2 && fours>=4) floor=8;
            if (p==5 && fives) floor=125;
            for(int grown=0;grown<=1;++grown) {
                auto &o=options[index][grown];
                if(!grown && seen<floor) {o.need=h+1;continue;}
                unsigned prefix_free=grown?0:heights.count[j];
                unsigned flags=grown?0:heights.flags[j];
                o.known=grown?0:heights.positions[j];
                unsigned need_min=(p==2 || p==5)?3:2;
                o.need=need_min>prefix_free?need_min-prefix_free:0;
                const bool missing4=fours>=4 && !(flags&1);
                const bool missing5=fives && !(flags&2);
                if(missing4 || missing5) o.need=std::max(1U,o.need);
                if(p==2) o.need+=(3U-prefix_free-o.need)&3U;
                if(p==3) o.need+=(5U-(prefix_free+o.need)%3)%3;
                if(o.need>unsigned(h) || !o.need) continue;
                U power=grown?std::max(floor,seen ? seen*p : p):seen;
                Wide strongest=power;
                if(missing4) strongest=std::max(strongest,Wide(power/std::gcd(power,U(4)))*4);
                if(missing5) strongest=std::max(strongest,Wide(power/std::gcd(power,U(25)))*25);
                auto cost=[&](Wide div) -> Wide {
                    if(div>=L) return 1;
                    const U step=static_cast<U>(div);
                    const Wide b=Wide((prefix.back()-1)/step+1)*step;
                    if(b>=L) return 1;
                    return 1+(ll-1)/(b*b);
                };
                o.cost=(o.need-1)*cost(power)+cost(strongest);
            }
        }
        U branch_high=0;
        for(auto& a:options[0]) for(auto& b:options[1]) {
            if(a.known & b.known || a.need+b.need>unsigned(h)) continue;
            for(auto& c:options[2]) {
                unsigned slots=a.need+b.need+c.need;
                if(slots>unsigned(h) || ((a.known|b.known)&c.known)) continue;
                Wide reserved=1+a.cost+b.cost+c.cost;
                U bound=high;
                if(E>reserved) bound=std::min(bound,square_root((h-slots)*ll/(E-reserved)));
                branch_high=std::max(branch_high,bound);
            }
        }
        high=branch_high;
        if(high<prefix.back()) return;
    }
    // Split only after the initial run of fours. Its low-multiplicity
    // branches are short enough to serve directly as portable jobs.
    if (listing && (h==2 || (prefix.back()>4 && (fours<10 || prefix.size()>=split_depth)))) {
        task();
        return;
    }
    if(h==2) {pair_tail(L,support,E,st); return;}
    U low=std::max(prefix.back(),square_root(ll/E)+1);
    if(low>high) return;
    auto next=[&](U b,Mask mask) {
        U g=std::gcd(b,L), k=b/g, d=L/g;
        auto state=append(st,b==prefix.back(),mask);
        if (!possible(state,mask,support|mask,h-1)) return;
        if (L > (U(1)<<58)/k) {
            prefix.push_back(b);
            exact_task();
            prefix.pop_back();
            return;
        }
        Wide next_E=E*k*k-Wide(d)*d;
        prefix.push_back(b);
        dfs(L*k,support|mask,next_E,h-1,state,heights.child(b,mask,prefix.size()-1));
        prefix.pop_back();
    };
    if(prefix.size()<fixed.size()) {
        U b=fixed[prefix.size()];
        Mask mask;
        if(low<=b && b<=high && sieve(b,future,mask)) next(b,mask);
    } else if(high-low<=35) {
        for(U b=low;b<=high;++b) {Mask mask; if(sieve(b,future,mask)) next(b,mask);}
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
    if ((bs.front()<4 || bs.front()>square_root(rank_limit-1)) || bs.size()>rank_limit-3 ||
        !std::is_sorted(bs.begin(),bs.end())) throw std::invalid_argument("invalid prefix");
    return bs;
}

int main(int argc,char** argv) {
    if (argc<4 || argc>6 || std::string_view(argv[1])!="--rank") return 2;
    try {
        const auto rank=integer(argv[2]);
        if (rank<9 || rank>native_rank_limit) throw std::invalid_argument("native rank must be between 9 and 64");
        // Once the first label is chosen, an unseen odd prime p needs p-1
        // of the other nonunit labels. A prime in the first denominator is smaller still.
        num_primes=std::lower_bound(primes.begin(),primes.end(),rank)-primes.begin();
        max_fours=rank<=31 ? 13 : 14;
        rank_limit=static_cast<unsigned>(rank);
        argc-=2;
        argv+=2;
        auto atomic=::fpathconf(STDOUT_FILENO,_PC_PIPE_BUF);
        if(atomic>0) pipe_limit=static_cast<std::size_t>(atomic);
        std::string_view mode=argv[1];
        if (mode=="--tasks" || mode=="--tasks-cost") {
            list_cost=(mode=="--tasks-cost");
            listing=true;
            if (argc>=3) {
                auto depth=integer(argv[2]);
                if (depth<14 || depth>22) return 2;
                split_depth=static_cast<unsigned>(depth);
            }
            if (argc==4) {
                auto group=integer(argv[3]);
                if (!group || group>1000000) return 2;
                frontier_per_job=group;
            }
        } else if(mode=="--prefix" && argc==3) {
            fixed=parse_prefix(argv[2]);
        } else return 2;
        for (U b=4;b*b<rank_limit;++b) if (fixed.empty() || b==fixed.front()) {
            Mask mask=prime_mask(b);
            prefix={b};
            dfs(b,mask,b*b-1,rank_limit-2,append(PacketState{},false,mask),Heights{}.child(b,mask,0));
        }
        if(listing) {
            std::sort(frontier.begin(),frontier.end());
            emit_groups(0,frontier.size(),0);
        }
        flush_records();
        std::cout.flush();
        if (!std::cout) return 3;
        if (listing) std::cerr<<"TASKS "<<jobs<<"\n";
        else std::cerr<<"DONE nodes="<<nodes<<" profiles="<<types<<" exact="<<handoffs<<"\n";
    } catch (const std::exception& e) {
        std::cerr<<e.what()<<"\n";
        return 2;
    }
    return 0;
}
