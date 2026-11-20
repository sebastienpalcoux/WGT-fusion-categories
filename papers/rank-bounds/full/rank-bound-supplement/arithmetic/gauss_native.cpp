// Exact native continuation of the high-denominator odd Gauss recursion.
// All number-theoretic quantities use GMP integers; an interrupted run is unresolved.
#include <gmpxx.h>
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>
using Z=mpz_class;
using State=std::array<int,4>;
using Balance=std::pair<Z,Z>;
using States=std::vector<State>;
using Balances=std::vector<Balance>;
Z floorq(const Z&a,const Z&b){Z q;mpz_fdiv_q(q.get_mpz_t(),a.get_mpz_t(),b.get_mpz_t());return q;}
Z ceilq(const Z&a,const Z&b){Z q;mpz_cdiv_q(q.get_mpz_t(),a.get_mpz_t(),b.get_mpz_t());return q;}
Z az(const Z&a){return a<0?-a:a;}
Z mn(const Z&a,const Z&b){return a<b?a:b;}
Z mx(const Z&a,const Z&b){return a>b?a:b;}
Z lcmz(const Z&a,const Z&b){Z c;mpz_lcm(c.get_mpz_t(),a.get_mpz_t(),b.get_mpz_t());return c;}
Z rootz(const Z&a){Z b;mpz_sqrt(b.get_mpz_t(),a.get_mpz_t());return b;}
Z units(const Z&scale,const Z&b){return ceilq(scale,b*b);}
bool divides(const Z&a,const Z&b){return mpz_divisible_p(a.get_mpz_t(),b.get_mpz_t());}
unsigned mod(const Z&a,unsigned b){return mpz_fdiv_ui(a.get_mpz_t(),b);}
Z pure_next(const Z&b,unsigned p){Z v=p*p;while(v<b)v*=p;return v;}


bool zero_moments(int p,const Z&V,const Z&ZZ,const Z&lo,const Z&hi,int step){
    if(V<=ZZ)return false;
    int q=p-1;Z J=mx(lo,mn(hi,floorq((p*p-1)*V-q*q*ZZ,2)));
    if(V*((p*p-1)*J-q*q*ZZ)<J*J+q*q*ZZ*J)return false;
    Z T=V-ZZ;std::set<Z> hs={floorq(q*V-lo,q*step*T),floorq(q*V-hi,q*step*T),-1,1};
    for(const Z&h:hs){if(h==0)continue;J=q+step*(2*h+1)>=0?hi:lo;
        if(q*(V+J)-step*(2*h+1)*(q*V-J)+q*step*step*h*(h+1)*T<0)return false;}
    return true;
}

struct Packet {int a,f,charge;std::vector<std::pair<int,int>> phases;};
std::map<std::tuple<int,int,int>,std::vector<Packet>> packet_cache;
const std::vector<Packet>& packets(int g,int slots,int order){
    auto key=std::make_tuple(g,slots,order);auto old=packet_cache.find(key);
    if(old!=packet_cache.end())return old->second;
    std::array<std::vector<std::array<int,3>>,2> local;
    for(int side=0;side<2;++side){int p=side?5:3;
        if(g%p){local[side].push_back({0,0,0});continue;}
        for(int j=0;j<=slots/(p*(p-1));++j)
          for(int i=0;i<=slots/(p-1)-p*j;++i){
            if(g==15 && ((p==3&&order==45)?i!=0:j!=0))continue;
            local[side].push_back({(p-1)*(i+p*j),p*i+p*(p-1)*j,i});
          }
    }
    std::map<std::tuple<int,int,int>,std::vector<std::pair<int,int>>> groups;
    for(auto a:local[0])for(auto f:local[1])if(a[0]+f[0]<=slots)
        groups[{a[0],f[0],a[1]+f[1]}].push_back({a[2],f[2]});
    std::vector<Packet> out;
    for(auto &entry:groups){auto[a,f,c]=entry.first;out.push_back({a,f,c,entry.second});}
    return packet_cache.emplace(key,std::move(out)).first->second;
}
std::map<unsigned,std::vector<Z>> level_cache;
const std::vector<Z>& levels(const Z&cap){
    unsigned bits=mpz_sizeinbase(cap.get_mpz_t(),2);auto old=level_cache.find(bits);
    if(old!=level_cache.end())return old->second;
    Z high=Z(1)<<bits;std::set<Z> all;
    for(Z a=1;a<=high;a*=3)for(Z b=a;b<=high;b*=5)if(b>=15)all.insert(b);
    for(Z b=49;b<=high;b*=7)all.insert(b);
    return level_cache.emplace(bits,std::vector<Z>(all.begin(),all.end())).first->second;
}
Z mixed_next(const Z&b){
    Z bound=3*b;const auto&all=levels(bound);
    for(auto it=std::lower_bound(all.begin(),all.end(),b);it!=all.end();++it)
        if(mod(*it,15)==0)return *it;
    throw std::runtime_error("mixed semigroup bound failed");
}
Z part(const Z&L,unsigned p){Z r=L,v=1;while(mod(r,p)==0){r/=p;v*=p;}return v;}
int support_size,wq,initial_o3,initial_o5,total_neutral;
std::map<std::pair<Z,bool>,std::vector<Z>> neutral_cache;
const std::vector<Z>& neutral_floors(const Z&L,bool nine){
    auto key=std::make_pair(L,nine);auto old=neutral_cache.find(key);
    if(old!=neutral_cache.end())return old->second;
    Z p3=part(L,3),p5=part(L,5),p7=part(L,7);
    std::vector<Z> floors;
    if(support_size==3){
        if(nine)floors.assign(total_neutral/2-1,9*p7);
        else {floors.assign(wq/2,9*p7);floors.insert(floors.end(),(total_neutral-2-wq)/2,p7);}
        return neutral_cache.emplace(key,std::move(floors)).first->second;
    }
    floors={11*p3,11*p5,11*p7,p7,p7,p7};
    if(nine)for(auto&b:floors)b=lcmz(9,b);
    return neutral_cache.emplace(key,std::move(floors)).first->second;
}
Z neutral_mass(const Z&L,const Z&scale,bool nine=false){
    Z total=0;for(const auto&b:neutral_floors(L,nine))total+=2*units(scale,b);return total;
}

int nines,z,big_z,initial_qcap;
bool has_small_caps;
std::set<std::pair<int,int>> small_caps_three,small_caps_five;
bool qexact,small_five_rigid,minimal_five_allowed;
int joint_num,packet_size,joint_order;
std::uint64_t nodes=0,terminals=0,node_limit=0;
struct Limit:std::exception{};

bool terminal(const Z&M,const Z&L,int k,const Z&mass,const Z&scale,const Z&first,int o3,int o5){
    ++terminals;
    std::cout<<"T "<<M<<' '<<L<<' '<<k<<' '<<mass<<' '<<scale<<' '<<first<<' '<<o3<<' '<<o5<<std::endl;
    int answer;if(!(std::cin>>answer)||!(answer==0||answer==1))throw std::runtime_error("invalid terminal response");
    return answer==1;
}

bool visit(const Z&last,int k,const Z&U,const Z&mass,const Z&L,const Z&S,
           States states,int zleft,int qleft,int qmax,const Z&qfirst,int o3,int o5,Balances balances,
           const Z&eligible,const Z&cosine,const Z&cosine_square,
           const Z&mixed_three,const Z&mixed_five,const Z&mixed_joint,int small_three,int small_five,int count_three,int count_five,int three45,int five45){
    if(node_limit && ++nodes>node_limit)throw Limit();if(!node_limit)++nodes;
    Z scale=8*S*S,bound=scale/L,charge_room=U+bound;
    if(states.empty()||charge_room<0||U==0)return false;
    int outside=qleft+zleft;
    Z next_mix=units(scale,mixed_next(last+1));
    Z future_three=z==48?charge_room:ceilq(2*charge_room,3);
    Z future_five=ceilq(4*charge_room,5);
    States kept;
    for(auto[sx,sy,shx,shy]:states){
        if(k<std::max(sx,o3+o5)+sy+outside)continue;
        int fslots=k-sx-sy-outside;
        Z highmix=units(scale,mixed_next(mx(last+1,45)));
        Z Au=(sx-shx)*next_mix+shx*highmix,Bu=(sy-shy)*next_mix+shy*highmix;
        Z pfloor=mn(pure_next(last+1,3),mn(pure_next(last+1,5),pure_next(last+1,7)));
        Z upper=(z!=48?ceilq(Au,2):Z(0))+ceilq(Bu,4)+ceilq(fslots*units(scale,pfloor),2)+ceilq(qleft*units(scale,pure_next(last+1,7)),6);
        if(U-bound+mass-scale>upper)continue;
        Z V=mixed_three+mn(sx*next_mix,future_three);
        if(24*V<mixed_joint)continue;
        if(z!=48 && 8*(mixed_five+mn(sy*next_mix,future_five))<mixed_joint)continue;
        Z lo=mixed_joint,hi=lo+zleft*highmix,ZZ=small_three*(scale/225);int lattice=5;
        if(z==32&&big_z==16){
            int tx=sx+count_three,ty=sy+count_five;
            if(tx<=20){Z odd=(tx-16)*units(scale,45);
                if(22*V<3*lo-50*odd||18*mixed_three>7*hi+50*odd)continue;
                ZZ+=std::max(0,three45-(tx-16))*(scale/2025);if(tx==16)lattice=10;}
            Z W=mixed_five+mn(sy*next_mix,future_five),odd=(ty-32)*units(scale,45);
            if(7*W<2*lo-9*odd||5*mixed_five>4*hi+9*odd)continue;
            Z Z5=small_five*(scale/225);if(ty<=44)Z5+=std::max(0,five45-(ty-32))*(scale/2025);
            if(!zero_moments(3,W,Z5,lo,hi,3))continue;
            if(tx<=20){Z odd3=(tx-16)*units(scale,45);
                Z excess=86*lo-308*mixed_three-154*mixed_five-700*odd3-198*odd;
                Z neutral=neutral_mass(L,scale)+floorq(scale,L*L)+(!(qfirst!=0&&qexact)?2*units(scale,qfirst==0?pure_next(last+1,7):qfirst):Z(0));
                Z room=mass+charge_room+neutral+fslots*units(scale,15)-scale;
                if(616*room<excess)continue;}
        }
        if(z==32||z==48||(z==64&&big_z))if(!zero_moments(5,V,ZZ,lo,hi,lattice))continue;
        kept.push_back({sx,sy,shx,shy});
    }
    states=std::move(kept);if(states.empty())return false;
    int x=k+1,y=k+1,hx=k+1,hy=k+1,maxx=0,maxy=0,maxsum=0;
    for(auto[sx,sy,shx,shy]:states){x=std::min(x,sx);y=std::min(y,sy);hx=std::min(hx,shx);hy=std::min(hy,shy);maxx=std::max(maxx,sx);maxy=std::max(maxy,sy);maxsum=std::max(maxsum,sx+sy);}
    int left=x,free=k-left-y-outside;if(free<0||U < -bound)return false;
    if(z==32&&big_z==0&&small_five_rigid&&last>=15){
        bool good=false;Z totalZ=small_five*(scale/225),maximum=mixed_five+mn(maxy*next_mix,future_five);
        for(int saturated=0;saturated<3;++saturated){
            Z sat=20*saturated*(scale/225);if(totalZ<sat||mixed_five<sat)continue;
            Z J=(32-16*saturated)*(scale/225),ZZ=totalZ-sat,V=mn(maximum-sat,5*J/4);
            if(J==0){if(minimal_five_allowed&&small_five==40){good=true;break;}continue;}
            if(zero_moments(3,V,ZZ,J,J,3)){good=true;break;}
        }
        if(!good)return false;
    }
    Z pure_three=std::max(0,k-x-y-outside)*units(scale,pure_next(last+1,3));
    Z imbalance_room=9*(scale-mass)+9*mixed_joint-16*mixed_three+16*pure_three;
    Balances retained;
    for(const auto&v:balances)if(z==48?2*v.first<=pure_three:50*v.first<=imbalance_room)retained.push_back(v);
    balances=std::move(retained);if(balances.empty())return false;
    Z tail=neutral_mass(L,scale,true)+floorq(scale,L*L)
          +mn((k-outside-y-x)*(scale/729)+x*mn(scale/2025,floorq(scale,last*last)),(k-outside-y-o5)*(scale/729)+o5*(scale/50625));
    Z H=nines*(scale/81)+eligible+tail;
    Z u=mn(4*tail,mx(-2*tail,H-4*cosine));
    if((4*cosine+u)*(4*cosine+u)>2*H*(8*cosine_square+4*tail+u))return false;
    if(std::find(states.begin(),states.end(),State{0,0,0,0})!=states.end()&&!outside
       &&std::find(balances.begin(),balances.end(),Balance{0,0})!=balances.end()
       &&divides(scale,az(U))){Z M=scale/az(U);if(divides(M,L)&&terminal(M,L,k,mass,scale,qfirst,o3,o5))return true;}
    Z base_book=mass+neutral_mass(L,scale)+floorq(scale,L*L);
    Z book=base_book+(!(qfirst!=0&&qexact)?2*units(scale,qfirst==0?pure_next(last+1,7):qfirst):Z(0));
    if(book+mx(0,U+bound)+(8-joint_num)*zleft*(scale/(8*2025))+free*(scale/225)<scale)return false;
    if(3*mass+2*(U-bound)>3*scale)return false;
    Z p7=part(L,7),factor=mn(L/(part(L,3)*part(L,5)),7*(L/p7));
    Z cap=floorq(scale+rootz(scale*scale+6*k*factor*factor*az(U)*scale),2*factor*az(U));
    const auto&all=levels(cap);
    for(auto it=std::upper_bound(all.begin(),all.end(),last);it!=all.end()&&*it<=cap;++it){
        const Z b=*it,b2=b*b;
        if(2*(U-bound)*b2>3*k*scale)break;
        bool prime=mod(b,7)==0,mixed=mod(b,15)==0;
        Z multiple=lcmz(L,prime?Z(b2/7):b),S1=lcmz(S,multiple),ratio=S1/S,up=ratio*ratio,scale1=scale*up,unit=scale1/b2;
        Z book1=base_book*up+(!(qfirst!=0&&qexact)?2*units(scale1,qfirst==0?pure_next(b,7):qfirst):Z(0));
        Z mix=units(scale1,mixed_next(b)),highmix=units(scale1,mixed_next(mx(b,45)));
        Z isolated_mass=qleft*units(scale1,pure_next(b,7));
        Z three_mass=(left-hx)*mix+hx*highmix,five_mass=(y-hy)*mix+hy*highmix;
        Z room=book1+free*(scale1/225)+three_mass+five_mass+zleft*highmix+isolated_mass-scale1;
        if(room<0)break;
        Z required=(z==48?24:36)*three_mass+30*five_mass+28*isolated_mass+3*joint_num*zleft*highmix;
        if(675*24*room<2*(b2-225)*(24*(U-bound)*up-required))break;
        if(support_size==3){
            Z best=0;Z overlap3=units(scale1,mixed_next(mx(b,135))),overlap5=units(scale1,mixed_next(mx(b,225)));
            // These overlaps require divisibility, not just the numerical floor.
            auto next_divisible=[&](int divisor)->Z{Z target=mx(b,divisor),cap=3*target;
                while(true){const auto& ls=levels(cap);for(auto it=std::lower_bound(ls.begin(),ls.end(),target);it!=ls.end();++it)if(mod(*it,divisor)==0)return *it;cap*=3;}};
            overlap3=units(scale1,next_divisible(135));overlap5=units(scale1,next_divisible(225));
            for(int i=0;i<=std::min(left,o5);i+=2)for(int j=0;j<=std::min(left-i,o3);j+=2){
                int excess=std::max(0,hx-i-j),ordinary=left-i-j-excess,available=k-outside-y-left-o3-o5+i+j;
                if(available<0)continue;
                Z upper=book1+y*mix+isolated_mass+zleft*highmix+excess*highmix+ordinary*mix+j*overlap3+i*overlap5+(o3-j)*(scale1/729)+(o5-i)*(scale1/50625)+available*(scale1/225);
                best=mx(best,upper);}
            if(best<scale1)break;
        }else if(book1+o3*(scale1/729)+o5*(scale1/50625)+(k-o3-o5-y-outside)*(scale1/225)+five_mass+isolated_mass+zleft*highmix<scale1)break;
        Z max3=floorq((k-y-outside-o3-o5)*unit+o3*mn(unit,scale1/729)+o5*mn(unit,scale1/50625),2),max5=floorq((k-outside-std::max(x,o3+o5))*unit,4);
        Balances feasible;for(auto[d,e]:balances)if(d*up<=max3&&e*up<=max5)feasible.push_back({d*up,e*up});
        if(feasible.empty())break;
        Z cap_charge=U*up+scale1/multiple;
        auto choose=[&](int a,int f,int j,int q,int ncharge,const std::vector<std::pair<int,int>>&phases)->bool{
            if(mass*up+(a+f+j+q)*unit>=scale1)return false;
            Z charge=ncharge*(unit/8);if(charge>cap_charge)return false;
            Z next_three=mixed_three*up+(mixed?a*unit:Z(0));
            Z next_five=mixed_five*up+(mixed?f*unit:Z(0));
            Z next_joint=mixed_joint*up+j*unit,joint_upper=next_joint+(zleft-j)*unit;
            if(16*next_three>9*joint_upper || (z==48?2*next_five>joint_upper:4*next_five>5*joint_upper))return false;
            int k1=k-a-f-j-q;States next_states;
            for(auto[sx,sy,shx,shy]:states){
                if(last==9 && has_small_caps){int a15=b==15?a:0,f15=b==15?f:0;
                    if(!small_caps_three.count({a15,sx})||!small_caps_five.count({f15,sy}))continue;}
                if(mixed&&(a>sx||f>sy))continue;
                int nx=mixed?sx-a:sx,ny=mixed?sy-f:sy;
                int nhx=mixed&&b>15?std::max(0,shx-a):shx;
                int nhy=mixed&&b>15?std::max(0,shy-f):shy;
                if(nhx<=nx&&nhy<=ny)next_states.push_back({nx,ny,nx,ny});
            }
            std::sort(next_states.begin(),next_states.end());next_states.erase(std::unique(next_states.begin(),next_states.end()),next_states.end());
            if(next_states.empty())return false;
            int x1=k+1,y1=k+1;for(auto[sx,sy,shx,shy]:next_states){x1=std::min(x1,sx);y1=std::min(y1,sy);}
            int q1=std::max(0,qleft-q),other=q1+zleft-j;
            int o3new=mod(b,27)==0?std::max(0,o3-a):o3,o5new=mod(b,225)==0?std::max(0,o5-a):o5;
            if(k1<std::max(x1,o3new+o5new)+y1+other)return false;
            Z m3=floorq((k1-y1-other-o3new-o5new)*unit+o3new*mn(unit,scale1/729)+o5new*mn(unit,scale1/50625),2),m5=floorq((k1-other-std::max(x1,o3new+o5new))*unit,4);
            std::set<Balance> next_set;
            for(auto[i,t]:phases)for(const auto&[d,e]:feasible)
                for(int v=-i;v<=i;v+=2){Z dd=az(d+v*unit);if(dd>m3)continue;
                    for(int w=-t;w<=t;w+=2){Z ee=az(e+w*unit);if(ee<=m5)next_set.insert({dd,ee});}}
            if(next_set.empty())return false;
            int active=mod(b,9)==0?a:0,first=active?phases[0].first:0;
            return visit(b,k1,U*up-charge,mass*up+(a+f+j+q+(prime&&qexact?2:0))*unit,multiple,S1,
              std::move(next_states),zleft-j,q1,qmax-q,(prime&&qfirst==0)?b:qfirst,o3new,o5new,Balances(next_set.begin(),next_set.end()),
              eligible*up+active*unit,cosine*up-first*unit,cosine_square*up+(active-first)*(unit/2),next_three,next_five,next_joint,small_three+(b==15?a:0),small_five+(b==15?f:0),count_three+(mixed?a:0),count_five+(mixed?f:0),three45+(b==45?a:0),five45+(b==45?f:0));
        };
        if(prime){for(int q=12;q<=std::min(k,qmax);q+=12)if(choose(0,0,0,q,56*(q/6),{{0,0}}))return true;continue;}
        int g=std::gcd(static_cast<int>(mod(b,15)),15);
        for(int j=0;j<=std::min(k,zleft);j+=packet_size){
            if(j && !(mixed&&b>15))break;
            if(z==32&&big_z==16&&j&&j!=zleft)continue;
            int room_slots=mixed?std::min(k-j,maxsum):k-j;
            for(const auto&p:packets(g,room_slots,joint_order)){
                if(!j&&!p.a&&!p.f)continue;
                if(mixed&&(p.a>maxx||p.f>maxy))continue;
                if(choose(p.a,p.f,j,0,8*p.charge+joint_num*j,p.phases))return true;
            }
        }
    }
    return false;
}

int main(int argc,char**argv){
    if(argc>1)node_limit=std::stoull(argv[1]);
    try{
        int slots,number;Z root,S,initial_mass,initial_charge;
        if(!(std::cin>>support_size>>wq>>initial_o3>>initial_o5>>nines>>z>>big_z>>total_neutral>>initial_qcap>>slots>>root>>S>>initial_mass>>initial_charge>>number))return 2;
        if((support_size!=3&&support_size!=4) || (support_size==4&&total_neutral!=14) || !(z==32||z==48||z==64))throw std::runtime_error("unsupported native context");
        qexact=initial_qcap==12;
        joint_num=z==48?8:7;packet_size=z==48?12:4;joint_order=z==48?45:15;
        States states;for(int i=0;i<number;++i){int x,y,hx,hy;std::cin>>x>>y>>hx>>hy;states.push_back({x,y,hx,hy});}
        small_five_rigid=true;minimal_five_allowed=false;
        for(auto st:states){if(st[1]>=72)small_five_rigid=false;if(st[1]==40)minimal_five_allowed=true;}
        std::cin>>has_small_caps;
        if(has_small_caps){int count;std::cin>>count;for(int i=0;i<count;++i){int a,x;std::cin>>a>>x;small_caps_three.insert({a,x});}
            std::cin>>count;for(int i=0;i<count;++i){int a,x;std::cin>>a>>x;small_caps_five.insert({a,x});}}
        int charges;std::cin>>charges;Z scale=8*S*S;
        std::vector<std::pair<int,Balances>> starts;
        for(int t=0;t<charges;++t){int c,n;std::cin>>c>>n;Balances balances;for(int i=0;i<n;++i){int d;std::cin>>d;balances.push_back({d*(scale/81),0});}starts.push_back({c,std::move(balances)});}
        for(auto&[c,balances]:starts){
            if(visit(9,slots,initial_charge-c*(scale/27),initial_mass,root,S,states,big_z,12,initial_qcap,0,initial_o3,initial_o5,balances,0,(nines-3*c)*(scale/81),0,0,0,(z-big_z)*(scale/225),0,0,0,0,0,0)){
                std::cout<<"DONE 1 "<<nodes<<' '<<terminals<<std::endl;return 0;}
        }
        std::cout<<"DONE 0 "<<nodes<<' '<<terminals<<std::endl;
    }catch(const Limit&){std::cout<<"UNRESOLVED nodes "<<nodes<<std::endl;return 3;}
    catch(const std::exception&e){std::cout<<"UNRESOLVED error "<<e.what()<<std::endl;return 2;}
}
