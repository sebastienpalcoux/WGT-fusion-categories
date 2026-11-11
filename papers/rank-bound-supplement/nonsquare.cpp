// Nonsquare dimensions. Build with g++ -O3 -std=c++17 nonsquare.cpp -o nonsquare.
// ./nonsquare R [Q B initial_multiplicity] (zero multiplicity selects the full Q,B root).
// ./nonsquare R writes the candidate degree lists at rank R in JSON Lines format.
#include <algorithm>
#include <limits>
#include <utility>
#include <cmath>
#include <cerrno>
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
std::vector<U> prefix;
std::vector<W> max_square;
int rank_limit, q, fiber;
int selected_q=0, selected_b=0, selected_m=0;
Mask target, ramified;
U nodes=0, candidates=0, exact=0;
#ifndef NS_LCM_BITS
#define NS_LCM_BITS 58
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
    record("{\"exact\":"+list(prefix)+",\"primes\":"+list(support)+"}");
    ++exact;
}

// Descendants only write at or beyond their parent's end index.  Passing the
// three indices by value restores the pending obligations on backtracking.
std::vector<std::pair<U,unsigned>> ternary_due;
struct Ternary {
    int head=0, end=0, debt=0;
    bool close(U z,unsigned n,U next,int slots) {
        if (z%3) {
            unsigned held=0;
            if (head<end && ternary_due[head].first==z) {
                held=ternary_due[head++].second; debt-=held;
            }
            if (n<held) return false;
            unsigned need=2*((n-held)%3);
            if (need) {
                // The native width check gives z < 2^64/6, also true at emission.
                ternary_due[end++]={2*z,need}; debt+=need;
            }
        }
        return (head==end || ternary_due[head].first>=next) && debt<=slots;
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

// A closed run contributes complete square-Galois packets. The last run
// can continue into the h remaining positions.
bool enough_packets(int h) {
    int debt=0;
    for (unsigned i=1;i<primes.size();++i) if ((target&~ramified)&(1u<<i)) {
        int p=primes[i], step=(p-1)/2, packets=0, doubles=0, tail=0;
        for (unsigned j=0;j<prefix.size();) {
            unsigned k=j+1;
            while (k<prefix.size() && prefix[k]==prefix[j]) ++k;
            if (prefix[j]%p==0) {
                int count=k-j;
                packets+=count/step;
                doubles+=count/(2*step);
                if (k==prefix.size()) tail=count;
            }
            j=k;
        }
        debt+=std::max(0,fiber*(p-1)-step*packets-tail%step);
        if (debt>h) return false;
        if (q%2==0 && packets+(h+tail%step)/step
                           +doubles+(h+tail%(2*step))/(2*step)<6) return false;
    }
    return true;
}

void emit(U L,U k,U d,U e,Ternary pending,unsigned run) {
    if (std::gcd(k,d)!=1 || std::gcd(k,e)!=1 || k%2==0 || k%5==0) return;
    if (L>UINT64_MAX/k) throw std::overflow_error("terminal scale");
    U M=k*L;
    if (factors(M)!=target) return;
    U b=M/e, c=M/d;
    if (q==3) {
        // Consume the two last denominators and close the last run.
        U z=prefix.back();
        if (b!=z) {
            if (!pending.close(z,run,b,2)) return;
            z=b; run=0;
        }
        ++run;
        if (c!=z) {
            if (!pending.close(z,run,c,1)) return;
            z=c; run=0;
        }
        if (!pending.close(z,run+1,UINT64_MAX,0)) return;
    }
    std::vector<U> degrees={1,d,e};
    for (auto it=prefix.rbegin();it!=prefix.rend();++it) degrees.push_back(M / *it);
    auto m=std::count(degrees.begin()+1,degrees.end(),degrees.back());
    if (selected_m && m!=selected_m) return;
    if (q==2 && prefix.front()==2 &&
        ((m>4 && __builtin_popcount(target)<=4) || (m>2 && M%4==0))) return;
    for (unsigned j=1;j<degrees.size();++j)
        if (__builtin_popcount(factors(degrees[j]))<2) return;
    for (unsigned i=0;i<primes.size();++i) if (target&(1u<<i)) {
        unsigned p=primes[i], free=0;
        for (unsigned j=1;j<degrees.size();++j) free+=(degrees[j]%p!=0);
        if (free<2) return;
    }
    prefix.push_back(M/e);
    prefix.push_back(M/d);
    bool possible=enough_packets(0);
    prefix.pop_back();
    prefix.pop_back();
    if (!possible) return;
    record(list(degrees));
    ++candidates;
}

void search(U L,W E,int h,unsigned maxima,unsigned run,const Ternary& pending) {
    ++nodes;
    bool eight=q==2 && prefix.front()==2;
    if (eight && ((maxima>4 && __builtin_popcount(target)<=4) || (maxima>2 && L%4==0))) return;
    if (!enough_packets(h)) return;
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
    if (lo>hi) return;
    // At the terminal pair M=kL, k<=hi. Hand off before producing any
    // native result if this bound could exceed the representation of M.
    if (h==2 && hi>UINT64_MAX/L) {handoff(); return;}
    auto visit=[&](U b) {
        if (eight && maxima>2 && b%4==0) return;
        U g=std::gcd(b,L), k=b/g, e=L/g;
        // E*k*k <= E*b*b <= (h+1/36)L^2, strictly below 2^128.
        W rest=E*k*k-W(e)*e;
        if (h==2) {
            if (e<6 || rest<37 || rest>W(e)*e+1 || k%2==0 || k%5==0) return;
            U d=root(rest-1);
            if (W(d)*d==rest-1 && L%d==0) emit(L,k,d,e,pending,run);
        } else {
            bool same=(b==prefix.back());
            Ternary next=pending;
            if (q==3 && !same && !next.close(prefix.back(),run,b,h)) return;
            prefix.push_back(b);
            if (L>UINT64_MAX/k) handoff();
            else search(L*k,rest,h-1,maxima+(b==prefix.front()),same?run+1:1,next);
            prefix.pop_back();
        }
    };
    smooth(lo,hi,1,target,visit);
}

// The packet bound pays p-1 for each odd unramified prime and zero for 2.
// Traverse only feasible supports, in increasing mask order.
void supports(int i,Mask support,int budget,int minimum,unsigned b,int f) {
    if (i<0) {
        if (__builtin_popcount(support)>=minimum) {
            target=support;
            prefix={b};
            search(b,f-1,rank_limit-2,1,1,Ternary{});
        }
        return;
    }
    supports(i-1,support,budget,minimum,b,f);
    Mask bit=Mask(1)<<i;
    int cost=i ? primes[i]-1 : 0;
    if (!(support&bit) && cost<=budget)
        supports(i-1,support|bit,budget-cost,minimum,b,f);
}

int main(int argc,char** argv) try {
    if (argc!=2 && argc!=5) return 2;
    rank_limit=std::stoi(argv[1]);
    if (rank_limit<2) return 2;
    for (int p=2;p<=rank_limit;++p) {
        if (std::none_of(primes.begin(),primes.end(),[p](unsigned d){return p%d==0;})) {
            if (primes.size()==std::numeric_limits<Mask>::digits) {
                std::cerr<<"prime support exceeds native mask, use exact replay\n";
                return 75;
            }
            primes.push_back(p);
        }
    }
    ternary_due.resize(rank_limit);
    max_square.resize(rank_limit);
    for (int h=0;h<rank_limit;++h) max_square[h]=~W(0)/(36*h+1);
    if (argc==5) {
        selected_q=std::stoi(argv[2]); selected_b=std::stoi(argv[3]);
        selected_m=std::stoi(argv[4]);
        if (selected_q<2 || selected_q>=rank_limit || selected_b<1 ||
            selected_m<0 || selected_m>=rank_limit) return 2;
    }
    for (q=2;q<rank_limit;++q) {
        if (selected_q && selected_q!=q) continue;
        bool squarefree=true;
        for (unsigned p:primes) if (q%(p*p)==0) squarefree=false;
        if (!squarefree) continue;
        ramified=factors(q);
        int product=1;
        for (unsigned i=1;i<primes.size();++i)
            if (ramified&(1u<<i)) product*=(primes[i]+1)/2;
        fiber=(q%2 ? 1 : 2)*product;
        for (unsigned i=1;i<primes.size();++i) if (ramified&(1u<<i)) {
            int a=(primes[i]+1)/2;
            fiber=std::max(fiber,product/a*(a+1));
        }
        for (int b=1;q*b*b<rank_limit;++b) {
            if (selected_b && selected_b!=b) continue;
            int f=q*b*b;
            if (f<8 || (b==1 && __builtin_popcount(ramified)==1)) continue;
            Mask must=ramified|factors(b);
            int budget=rank_limit/fiber-1;
            for (unsigned i=1;i<primes.size();++i)
                if ((must&~ramified)&(Mask(1)<<i)) budget-=primes[i]-1;
            if (budget>=0) supports(int(primes.size())-1,must,budget,f==8?4:3,b,f);
        }
    }
    std::cerr<<"rank "<<rank_limit;
    if (selected_q) std::cerr<<", q "<<selected_q<<", B "<<selected_b<<", multiplicity "<<selected_m;
    std::cerr<<", nodes "<<nodes<<", candidates "<<candidates<<", exact "<<exact<<"\n";
} catch (const std::exception& e) {
    std::cerr<<"nonsquare generator failed: "<<e.what()<<'\n';
    return 1;
}
