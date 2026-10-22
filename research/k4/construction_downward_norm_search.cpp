// Search actual square roots, including either sign of M=4PQ.
// Delta=2521 corresponds to n=m(m+1)/210-4 (new shift m^2+m-630).
// Delta=-839 corresponds to n=m(m+1)/210-3 (new shift m^2+m+210).
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <vector>
using I=long long;
I root(I x){if(x<0)return -1;I s=std::sqrt((long double)x);while((__int128)s*s<x)++s;while((__int128)s*s>x)--s;return (__int128)s*s==x?s:-1;}
struct Row{I p,q,m,n,c,mask;std::array<I,6>r;};
int main(int argc,char**argv){
 I limit=argc>1?std::atoll(argv[1]):3000,delta=argc>2?std::atoll(argv[2]):2521;
 bool allsix=argc>3&&std::string(argv[3])=="allsix";
 std::vector<Row>rows;
 #pragma omp parallel for schedule(dynamic)
 for(I ap=1;ap<=limit;++ap){
  std::vector<Row>local;
  for(I q=1;q<=limit;++q){if(std::gcd(ap,q)!=1)continue;
   for(I sign:std::array<I,2>{1,-1}){I p=ap*sign,m=4*p*q,n=-delta*q*q-p*p-2*p*q,c=p*p-delta*q*q;
    std::array<I,6>r={root(-n),root(-n+14*m),root(-n-21*m),root(-n-m),root(-n-15*m),root(-n+20*m)};
    I mask=0;for(int k=0;k<(allsix?6:3);++k)if(r[k]>=0)mask|=1<<k;
    if(__builtin_popcount((unsigned)mask)>=(allsix?3:2))local.push_back(Row{p,q,m,n,c,mask,r});
   }
  }
  #pragma omp critical
  rows.insert(rows.end(),local.begin(),local.end());
 }
 std::sort(rows.begin(),rows.end(),[](const Row&a,const Row&b){return std::abs(a.m)<std::abs(b.m);});
 for(auto&r:rows){std::cout<<r.p<<' '<<r.q<<' '<<r.m<<' '<<r.n<<' '<<r.c<<' '<<r.mask;for(int k=0;k<(allsix?6:3);++k)std::cout<<' '<<r.r[k];std::cout<<'\n';}
 std::cerr<<"Delta "<<delta<<"; abs(P),Q <= "<<limit<<"; "<<(allsix?"at least three of all six splits ":"at least two selected splits ")<<rows.size()<<'\n';
}
