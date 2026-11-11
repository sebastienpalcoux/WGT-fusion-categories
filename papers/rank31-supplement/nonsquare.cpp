// Nonsquare dimensions. Build with g++ -O3 -std=c++17 nonsquare.cpp -o nonsquare.
// ./nonsquare R writes the candidate degree lists at rank R in JSON Lines format.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

using U = std::uint64_t;
using W = __uint128_t;
using Mask = unsigned;
const std::vector<unsigned> primes = {2,3,5,7,11,13,17,19,23,29};
std::vector<U> prefix;
int rank_limit, q, fiber;
Mask target, ramified;
U nodes=0, candidates=0;

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

void emit(U L,U k,U d,U e) {
    if (std::gcd(k,d)!=1 || std::gcd(k,e)!=1 || k%2==0 || k%5==0) return;
    if (L>UINT64_MAX/k) throw std::overflow_error("terminal scale");
    U M=k*L;
    if (factors(M)!=target) return;
    std::vector<U> degrees={1,d,e};
    for (auto it=prefix.rbegin();it!=prefix.rend();++it) degrees.push_back(M / *it);
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
    std::cout<<"[";
    for (unsigned j=0;j<degrees.size();++j) std::cout<<(j?",":"")<<degrees[j];
    std::cout<<"]\n";
    ++candidates;
}

void search(U L,W E,int h) {
    ++nodes;
    if (!enough_packets(h)) return;
    if (prefix.back()==prefix.front()) {
        U f=q*prefix.front()*prefix.front();
        for (unsigned p:primes) if (p%2 && f%p==0 && f%(p*p)
                                      && prefix.size()>f-1-(p-1)/2) return;
    }
    if (L>(U(1)<<58)) throw std::overflow_error("prefix scale");
    W ll=W(L)*L;
    if (E>(36*h+1)*ll/36) return;
    U lo=std::max(prefix.back(),root(ll/E)+1);
    U hi=root((36*h+1)*ll/(36*E));
    if (lo>hi) return;
    auto visit=[&](U b) {
        U g=std::gcd(b,L), k=b/g, e=L/g;
        if (E>~W(0)/k/k) throw std::overflow_error("residual");
        W rest=E*k*k-W(e)*e;
        if (h==2) {
            if (e<6 || rest<37 || rest>W(e)*e+1 || k%2==0 || k%5==0) return;
            U d=root(rest-1);
            if (W(d)*d==rest-1 && L%d==0) emit(L,k,d,e);
        } else {
            if (L>UINT64_MAX/k) throw std::overflow_error("scale");
            prefix.push_back(b);
            search(L*k,rest,h-1);
            prefix.pop_back();
        }
    };
    smooth(lo,hi,1,target,visit);
}

int main(int argc,char** argv) {
    if (argc!=2) return 2;
    rank_limit=std::stoi(argv[1]);
    if (rank_limit<6 || rank_limit>30) return 2;
    for (q=2;q<rank_limit;++q) {
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
            int f=q*b*b;
            if (f<8 || (b==1 && __builtin_popcount(ramified)==1)) continue;
            Mask must=ramified|factors(b);
            for (target=must;target<(1u<<primes.size());++target) {
                if ((target&must)!=must || __builtin_popcount(target)<(f==8?4:3)) continue;
                int phases=1;
                for (unsigned i=1;i<primes.size();++i)
                    if ((target&~ramified)&(1u<<i)) phases+=primes[i]-1;
                if (fiber*phases>rank_limit) continue;
                prefix={U(b)};
                search(b,f-1,rank_limit-2);
            }
        }
    }
    std::cerr<<"rank "<<rank_limit<<", nodes "<<nodes<<", candidates "<<candidates<<"\n";
}
