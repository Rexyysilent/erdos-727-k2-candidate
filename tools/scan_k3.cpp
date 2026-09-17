// Exact scanner for n=210*t^2+391*t+179, k=3; finite evidence, not an analytic proof.
// Requires GCC/Clang's unsigned 128-bit integer extension.
// g++ -std=c++17 -O3 tools/scan_k3.cpp -o scan_k3
// ./scan_k3 START STOP [--rows]
// Largest sieve index <=50 million; at most about 200 MB for the sieve.
#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
using U64=std::uint64_t;
using I64=std::int64_t;
using U128=__uint128_t;
U64 number(const char* arg){
    std::string s(arg);
    if(s.empty()||s.find_first_not_of("0123456789")!=std::string::npos)
        throw std::invalid_argument("require unsigned decimal integers");
    return std::stoull(s);
}
U64 valuation(U64 n,U64 p){U64 total=0;while(n){n/=p;total+=n;}return total;}
int main(int argc,char**argv){
 try{
    if(argc!=3&&argc!=4)throw std::invalid_argument("use START STOP [--rows]");
    U64 u=14,r=195,start=number(argv[1]),stop=number(argv[2]);
    bool rows=argc==4;
    if(rows&&std::string(argv[3])!="--rows")throw std::invalid_argument("unknown option");
    if(start>=stop||stop>25000000)throw std::invalid_argument("require START<STOP<=25000000, subject to the sieve cap");
    U64 D=u*(u+1);
    if(r>=D||r*(r+1)%D)throw std::invalid_argument("R must satisfy 0<=R<D and D divides R(R+1)");
    if(rows&&stop-start>5000)throw std::invalid_argument("--rows limited to 5000 values");
    I64 k=I64(r*(r+1)/D)-3;
    std::vector<std::pair<I64,I64>> lines;
    U64 limit=1;
    for(I64 c:{I64(0),I64(1),-I64(u),I64(u+1),I64(-20),I64(21)}){
        I64 g=std::gcd(I64(D),I64(r)+c),A=I64(D)/g,B=(I64(r)+c)/g;
        I64 first=A*I64(start)+B,last=A*I64(stop-1)+B;
        if(first<1||last>50000000)throw std::invalid_argument("require positive linear factors and largest sieve index <=50000000");
        lines.emplace_back(A,B);limit=std::max(limit,U64(last));
    }
    std::vector<std::uint32_t> spf(limit+1,0);
    for(U64 p=2;p*p<=limit;++p)if(!spf[p])
        for(U64 j=p*p;j<=limit;j+=p)if(!spf[j])spf[j]=p;
    U64 good=0;std::vector<U64>firsts;
    if(rows)std::cout<<"[";
    for(U64 t=start;t<stop;++t){
        I64 ni=I64(D*t*t+(2*r+1)*t)+k;
        if(ni<0)throw std::invalid_argument("n must be nonnegative");
        U64 n=U64(ni);std::map<U64,unsigned>fs;
        std::vector<U64> values;
        for(auto [A,B]:lines){U64 v=U64(A*I64(t)+B);values.push_back(v);
            while(v>1){U64 p=spf[v]?spf[v]:v;do{++fs[p];v/=p;}while(v%p==0);}}
        if(U128(values[0])*values[1]!=n+3 || U128(values[2])*values[3]!=n+2 || U128(values[4])*values[5]!=n+1)throw std::logic_error("factorization identity failed");
        bool pass=true;
        for(auto [p,e]:fs)pass=pass&&(valuation(2*n,p)-2*valuation(n,p)>=2*e);
        good+=pass;if(pass&&firsts.size()<8)firsts.push_back(n);
        if(rows){if(t>start)std::cout<<",";std::cout<<"{\"t\":"<<t<<",\"n\":"<<n<<",\"passes\":"<<(pass?"true":"false")<<"}";}
    }
    if(rows){std::cout<<"]\n";return 0;}
    std::cout<<std::setprecision(12)<<"{\"k\":3,\"u\":"<<u<<",\"D\":"<<D<<",\"r\":"<<r
        <<",\"start_inclusive\":"<<start<<",\"stop_exclusive\":"<<stop
        <<",\"parameters_checked\":"<<(stop-start)<<",\"valid_parameters\":"<<good
        <<",\"observed_fraction\":"<<double(good)/(stop-start)<<",\"first_n_values\":[";
    for(unsigned i=0;i<firsts.size();++i){if(i)std::cout<<",";std::cout<<firsts[i];}
    std::cout<<"],\"scope\":\"finite exact arithmetic only; not asymptotic verification\"}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 2;}
}
