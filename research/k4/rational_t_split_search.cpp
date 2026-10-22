// Search two explicit rational one-parameter quartic factorization families.
// A: u=(1,3,t,5), L=15,
// B=(-t^6+336t^2+225)/(2t^2(t^2-17)).
// B: u=(7,5,1,t), L=35,
// B=(-t^6-3024t^2+1225)/(2t^2(t^2+23)).
// Here t=p/q in lowest terms. The output counts the eight rational quadratic
// factors whose positive discriminants are rational squares. Hits do not
// certify integer values or the factorial-divisibility property.
// Build: g++ -O3 -std=c++17 rational_t_split_search.cpp -o search-t
// Run: ./search-t 5000
// At bound 5000: 15,200,915 reduced fractions; the only >=4-split hit is
// family B, t=2, with four split quadratics. No >=5-split hit was found.
// Integer arithmetic is exact; floating sqrt only seeds a corrected exact
// integer-square test. Bounds above 10000 are rejected to prevent overflow.
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
using I=__int128_t;
bool square(I z){if(z<=0)return false;I r=(I)sqrtl((long double)z);while(r*r<z)++r;while(r*r>z)--r;return r*r==z;}
int main(int argc,char**argv){int n=argc>1?std::stoi(argv[1]):5000;if(n>10000)return 2;long long tested=0,hits=0;
for(int p=1;p<=n;p++)for(int q=1;q<=n;q++){
 if(std::gcd(p,q)>1)continue;I P=I(p)*p,Q=I(q)*q,p4=P*P,q4=Q*Q,p6=p4*P,q6=q4*Q;tested++;
 for(int f=0;f<2;f++){
  if(q==1&&(p==1||p==5||p==(f?7:3)))continue;
  I m=P+(f?23:-17)*Q;if(m==0)continue;
  I N=p6+(f?3024:-336)*P*q4+(f?-1225:-225)*q6,den=P*Q*m;
  std::array<int,5>off=f?std::array<int,5>{0,12,32,60,140}:std::array<int,5>{0,12,32,60,0};
  std::array<int,5>weight=f?std::array<int,5>{1,1,2,1,1}:std::array<int,5>{2,2,1,1,0};
  int count=0,mask=0;
  for(int j=0;j<(f?5:4);j++){if(square((N+(off[j]-(f?71:31))*den)*m)){count+=weight[j];mask|=1<<j;}}
  if(count<2)continue;
  int extra=0;for(int s:{-1,1})if(square((N-p4*m+s*(f?70:30)*I(p)*q*Q*m)*m))extra++;
  if(count+extra>=4){hits++;std::cout<<"family="<<(f?'B':'A')<<" p="<<p<<" q="<<q<<" split="<<count+extra<<" oldmask="<<mask<<" extra="<<extra<<'\n';}
 }
}
std::cerr<<"fractions="<<tested<<" hits="<<hits<<'\n';
}
