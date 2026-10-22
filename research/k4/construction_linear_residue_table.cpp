// Exact second/third/fourth carry residue table for the four linear factors
// in m=(y^2-1681)/552, n=m(m+1)/210-3. Fixed y modulo D.
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
using I=long long;
constexpr I D=210*552*552;
I mod(I x){x%=D;return x<0?x+D:x;}
I mul(I a,I b){return a*b%D;}
bool bad(I s,int sign){
 if(s==0)return sign>0;
 if(2*s==D)return sign<0;
 return 2*s<D;
}
int main(int argc,char**argv){
 I Y=argc>1?std::atoll(argv[1]):23971;
 std::array<I,4> roots={-97,-41,41,97};
 std::array<std::array<I,3>,4> K;
 std::array<std::array<I,5>,4> cnt{};
 for(int j=0;j<4;++j){
  I r=roots[j],a=mod(Y-r),a2=mul(a,a),a3=mul(a2,a);
  K[j]={mul(a2,mod(Y*Y+2*r*Y+3*r*r-2810)),
        mul(a3,mod(Y+3*r)),mul(a2,a2)};
 }
 I total=0;
 for(I u=1;u<D;++u){
  if(u%2==0||u%3==0||u%5==0||u%7==0||u%23==0)continue;
  ++total;I u2=mul(u,u),u3=mul(u2,u),u4=mul(u2,u2);
  for(int j=0;j<4;++j){
   int sign=roots[j]>0?1:-1;
   bool b2=bad(mul(K[j][0],u2),sign);
   bool b3=bad(mul(K[j][1],u3),1);
   bool b4=bad(mul(K[j][2],u4),sign);
   cnt[j][0]+=b2;cnt[j][1]+=b3;cnt[j][2]+=b4;
   cnt[j][3]+=b2&&b3;cnt[j][4]+=b2&&b3&&b4;
  }
 }
 std::cout<<"{\n  \"Y\": "<<Y<<",\n  \"modulus\": "<<D
          <<",\n  \"phi\": "<<total<<",\n  \"rows\": [\n";
 for(int j=0;j<4;++j){
  std::cout<<"    {\"root\": "<<roots[j]<<", \"K2\": "<<K[j][0]
           <<", \"K3\": "<<K[j][1]<<", \"K4\": "<<K[j][2]<<", \"bad_counts\": [";
  for(int k=0;k<5;++k)std::cout<<(k?", ":"")<<cnt[j][k];
  std::cout<<"]}"<<(j<3?",":"")<<'\n';
 }
 std::cout<<"  ],\n  \"bad_count_order\": [\"level2\",\"level3\",\"level4\",\"levels23\",\"levels234\"]\n}\n";
}
