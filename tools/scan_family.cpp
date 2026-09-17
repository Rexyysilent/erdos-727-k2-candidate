// Exact finite scanner for Erdős 727, k=2. Not an asymptotic proof.
// g++ -std=c++17 -O3 tools/scan_family.cpp -o scan_family
// ./scan_family START STOP [--rows]  (STOP is exclusive)
// Memory is approximately 24*STOP bytes; STOP is capped at 20,000,000.
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
using u64=std::uint64_t;
using u128=__uint128_t;
u64 parse(const char* s) {
    std::string v(s);
    if(v.empty() || v.find_first_not_of("0123456789")!=std::string::npos)
        throw std::invalid_argument("arguments must be nonnegative integers");
    return std::stoull(v);
}
u64 valuation(u64 n,u64 p){u64 total=0;while(n){n/=p;total+=n;}return total;}
void factor(u64 n,const std::vector<std::uint32_t>& spf,std::map<u64,unsigned>& factors){
    while(n>1){u64 p=spf[n]?spf[n]:n;do{++factors[p];n/=p;}while(n%p==0);}
}
int main(int argc,char**argv){
    try{
        if(argc<3 || argc>4)throw std::invalid_argument("use START STOP [--rows]");
        const u64 start=parse(argv[1]),stop=parse(argv[2]);
        if(start>=stop || stop>20000000)throw std::invalid_argument("require 0 <= START < STOP <= 20000000");
        bool rows=argc==4;
        if(rows && std::string(argv[3])!="--rows")throw std::invalid_argument("unknown option");
        if(rows && stop-start>10000)throw std::invalid_argument("--rows limited to 10000 parameters");
        const u64 limit=6*stop-1;
        std::vector<std::uint32_t> spf(limit+1,0);
        for(u64 p=2;p*p<=limit;++p)if(!spf[p])
            for(u64 i=p*p;i<=limit;i+=p)if(!spf[i])spf[i]=p;
        u64 good=0,large=0,fixed=0,overlap=0,ones=0,excess=0,remainder=0;
        bool dyadic=start>0 && stop==2*start;
        std::vector<u64> first;
        if(rows)std::cout<<"[";
        for(u64 t=start;t<stop;++t){
            u64 n=6*t*t+11*t+3,m=6*t+5;
            std::map<u64,unsigned> fs;
            for(u64 v:{2*t+1,3*t+4,t+1,m})factor(v,spf,fs);
            bool pass=true,H=false,F=false;u64 Hprime=0,W=0;
            for(auto [p,e]:fs){
                u64 cp=valuation(2*n,p)-2*valuation(n,p);
                bool ok=cp>=2*e;pass=pass&&ok;
                if(!ok && p<=5)F=true;
                if(p>=7 && cp==1){++W;if(e!=1)throw std::logic_error("initial carry identity failed");}
                if(dyadic && p>=7 && m%p==0 && p*p>12*start+3){
                    if(Hprime)throw std::logic_error("large primes were not unique");
                    Hprime=p;u128 q=u128(p)*p*p;
                    H=2*(u128(n)%q)<q;
                    if(H!=(cp==1) || 2*(n%(p*p))>=p*p || 2*u128(n)>=q*p)
                        throw std::logic_error("large-prime identity failed");
                }
            }
            good+=pass;large+=H;fixed+=F;overlap+=(H&&F);ones+=W;
            excess+=W>0?W-1:0;remainder+=(!pass && W==0);
            if(H && pass)throw std::logic_error("large obstruction admitted good value");
            if(pass && first.size()<8)first.push_back(n);
            if(rows){
                if(t>start)std::cout<<",";
                std::cout<<"{\"t\":"<<t<<",\"n\":"<<n<<",\"passes\":"<<(pass?"true":"false")
                    <<",\"simple_obstructions\":"<<W<<",\"large_obstruction\":"<<(H?"true":"false")<<"}";
            }
        }
        if(rows){std::cout<<"]\n";return 0;}
        if(good!=(stop-start)-ones+excess-remainder)throw std::logic_error("exact accounting failed");
        std::cout<<std::setprecision(12)<<"{\"start_inclusive\":"<<start<<",\"stop_exclusive\":"<<stop
            <<",\"parameters_checked\":"<<(stop-start)<<",\"valid_parameters\":"<<good
            <<",\"observed_fraction\":"<<double(good)/(stop-start)
            <<",\"simple_obstruction_incidences\":"<<ones<<",\"excess_simple_obstructions\":"<<excess
            <<",\"failed_without_simple_obstruction\":"<<remainder
            <<",\"large_obstruction_count\":"<<(dyadic?std::to_string(large):"null")
            <<",\"large_obstruction_fraction\":";
        if(dyadic)std::cout<<double(large)/(stop-start);else std::cout<<"null";
        std::cout<<",\"fixed_2_3_5_failure_count\":"<<fixed
            <<",\"large_and_fixed_overlap\":"<<(dyadic?std::to_string(overlap):"null")
            <<",\"exact_count_identity_passed\":true,\"first_n_values\":[";
        for(unsigned i=0;i<first.size();++i){if(i)std::cout<<",";std::cout<<first[i];}
        std::cout<<"],\"scope\":\"finite exact arithmetic only, not asymptotic verification\"}\n";
    }catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 2;}
}
