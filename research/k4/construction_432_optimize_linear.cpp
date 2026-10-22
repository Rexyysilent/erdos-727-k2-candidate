// Exhaustive exact optimization over the 1728 refinements of the four
// y mod 90720 classes with minimal quadratic large-prime cost 1/3.
// Pairing u with -u makes the third-carry bad event have probability 1/2
// conditional on u^2, except for the explicitly handled endpoint classes.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <set>
#include <vector>
using I=long long;
constexpr I D=39191040,L=90720;
I mod(I x){x%=D;return x<0?x+D:x;}
I mul(I a,I b){return a*b%D;}
I inv(I a,I m){I b=m,u=1,v=0;while(b){I t=a/b;a-=t*b;std::swap(a,b);u-=t*v;std::swap(u,v);}u%=m;return u<0?u+m:u;}
bool bad(I s,int sign){if(s==0)return sign>0;if(2*s==D)return sign<0;return 2*s<D;}
struct Row {I Y;std::array<I,6>b2{},b23twice{},b234twice{};double score;};
int main(){
 std::vector<I>squares={0};
 for(I q:std::array<I,4>{512,2187,5,7}){
  std::set<I>sq;for(I u=1;u<q;++u)if(std::gcd(u,q)==1)sq.insert(u*u%q);
  I basis=(D/q)*inv(D/q%q,q);std::vector<I>next;
  for(I a:squares)for(I b:sq)next.push_back((a+mul(b,basis))%D);
  squares=std::move(next);
 }
 std::vector<I>fourths;for(I s:squares)fourths.push_back(mul(s,s));
 std::vector<Row> rows;
 for(I Y:std::array<I,4>{16273,29087,61633,74447})
  for(I k=0;k<432;++k)rows.push_back(Row{Y+k*L});
 double c2=(1.0/16)*log(7.0/6)+(1.0/32)*log(8.0/7);
 double c23=log(9.0/8)/4,c234=log(16.0/15)/2+log(5.0/4);
 double small=0;for(int J=8;J<80;++J)small+=6*pow(2.0,1-J)*log((J+1.0)/J);
 #pragma omp parallel for schedule(dynamic)
 for(int j=0;j<(int)rows.size();++j){
  Row &row=rows[j];I Y=row.Y;
  std::array<I,6> roots={-143,-89,-73,73,89,143};
  for(int k=0;k<6;++k){
   I r=roots[k],a=mod(Y-r),a2=mul(a,a),a3=mul(a2,a);
   I K2=mul(a2,mod(Y*Y+2*r*Y+3*r*r-28370));
   I K3=mul(a3,mod(Y+3*r));I K4=mul(a2,a2);
   I E1=4*r*r*r-56740*r;int sign2=E1>0?1:-1,sign4=r>0?1:-1;
   // Number of bad members of the pair {u,-u} at carry level 3.
   int pairbad=K3==0?2:(2*K3==D?0:1);
   I b2=0,b24=0;
   for(int h=0;h<(int)squares.size();++h){
    bool b=bad(mul(K2,squares[h]),sign2);b2+=b;
    if(b&&bad(mul(K4,fourths[h]),sign4))++b24;
   }
   row.b2[k]=b2;row.b23twice[k]=pairbad*b2;row.b234twice[k]=pairbad*b24;
  }
  I B2=std::accumulate(row.b2.begin(),row.b2.end(),0LL);
  I B23=std::accumulate(row.b23twice.begin(),row.b23twice.end(),0LL);
  I B234=std::accumulate(row.b234twice.begin(),row.b234twice.end(),0LL);
  row.score=small+(c2*B2+(c23*B23+c234*B234)/2)/squares.size();
 }
 std::sort(rows.begin(),rows.end(),[](const Row&a,const Row&b){return a.score<b.score||(a.score==b.score&&a.Y<b.Y);});
 std::cout<<"{\n  \"modulus\": "<<D<<",\n  \"square_classes\": "<<squares.size()<<",\n  \"candidates\": "<<rows.size()<<",\n  \"root_order\": [-143,-89,-73,73,89,143],\n  \"rows\": [\n";
 std::cout.precision(17);
 for(int j=0;j<(int)rows.size();++j){
  auto &row=rows[j];std::cout<<"    {\"Y\": "<<row.Y<<", \"linear_cost\": "<<row.score;
  for(auto pair:{std::pair<const char*,std::array<I,6>>("b2",row.b2),{"b23_twice",row.b23twice},{"b234_twice",row.b234twice}}){
   std::cout<<", \""<<pair.first<<"\": [";for(int k=0;k<6;++k)std::cout<<(k?", ":"")<<pair.second[k];std::cout<<']';
  }
  std::cout<<'}'<<(j+1<(int)rows.size()?",":"")<<'\n';
 }
 std::cout<<"  ]\n}\n";
 std::cerr<<"Exact candidates "<<rows.size()<<"; square classes "<<squares.size()<<"; best Y "<<rows.front().Y<<"; best linear cost "<<rows.front().score<<'\n';
}
