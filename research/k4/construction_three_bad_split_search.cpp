// Complete bounded search within the norm parametrization for simultaneous
// splitting of m, m-14 and m+21 as rational linear factors in y.
// m=(y^2+N)/M, N=839Q^2-P^2-2PQ, M=4PQ, gcd(P,Q)=1.
// Required: -N, -N+14M, -N-21M are squares. The latter forces P/Q>91.
// Parametrize -N=U^2 by (P+Q-U)(P+Q+U)=840Q^2.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>
using I=long long;
bool sq(I v){
 if(v<0)return false;
 I s=std::sqrt((long double)v);
 while(s*s<v)++s;while(s*s>v)--s;return s*s==v;
}
int main(int argc,char**argv){
 I bound=argc>1?std::atoll(argv[1]):10000000;
 I maxQ=bound/91;
 std::vector<int> spf(maxQ+1);
 for(I i=2;i<=maxQ;++i)if(!spf[i])for(I j=i;j<=maxQ;j+=i)if(!spf[j])spf[j]=i;
 I candidates=0,two=0,hits=0;
 for(I Q=1;Q<=maxQ;++Q){
  std::vector<std::pair<I,int>> fs={{2,3},{3,1},{5,1},{7,1}};
  I z=Q;
  while(z>1){I p=spf[z];int e=0;while(z%p==0){z/=p;++e;}
   auto it=std::find_if(fs.begin(),fs.end(),[p](auto x){return x.first==p;});
   if(it==fs.end())fs.emplace_back(p,2*e);else it->second+=2*e;
  }
  std::vector<I> ds={1};
  for(auto [p,e]:fs){size_t n=ds.size();I mul=1;for(int j=1;j<=e;++j){mul*=p;for(size_t h=0;h<n;++h)ds.push_back(ds[h]*mul);}}
  I D=840*Q*Q;
  for(I d:ds){
   I e=D/d;if(d>e||(d+e)%2)continue;
   I P=(d+e)/2-Q;if(P<1||P>bound||std::gcd(P,Q)!=1)continue;
   I U=(e-d)/2,M=4*P*Q;
   I v2=U*U-21*M;if(v2<0)continue;
   ++candidates;
   if(!sq(U*U+14*M))continue;
   ++two;
   if(!sq(v2))continue;
   ++hits;std::cout<<P<<' '<<Q<<' '<<U<<' '<<M<<'\n';
  }
 }
 std::cerr<<"P bound "<<bound<<"; Q bound "<<maxQ<<"; candidates "<<candidates
          <<"; second-square passes "<<two<<"; three-square hits "<<hits<<'\n';
}
