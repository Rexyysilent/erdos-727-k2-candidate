// Exact residue optimization for either of two non-even quartics.
// Default D=483840; argument 40320 selects the eight-linear-factor family.
// A second argument "reverse" analyzes F(-y), including endpoint signs.
// norm120,norm208,norm432,norm152,norm688 analyze coarse representatives
// for the five even norm families with three split old quadratics.
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
#include <set>
#include <vector>
using I=long long;
I D=483840;
I mod(I x){x%=D;return x<0?x+D:x;}
I mul(I a,I b){return a*b%D;}
I inv(I a,I m){I b=m,u=1,v=0;while(b){I t=a/b;a-=t*b;std::swap(a,b);u-=t*v;std::swap(u,v);}u%=m;return u<0?u+m:u;}
bool bad(I s,int sign){if(s==0)return sign>0;if(2*s==D)return sign<0;return 2*s<D;}
struct Row {I Y;std::array<I,6> high{};std::array<I,8>b2{},b23twice{},b234twice{};I H;double score;};
int main(int argc,char**argv){
 bool eight=argc>1&&std::string(argv[1])=="40320";
 I B=eight?-810:-3010,C=eight?2240:7680,E=eight?-41751:3969;
 if(eight)D=40320;
 std::vector<I>roots=eight?std::vector<I>{-29,-27,-9,1,9,13,19,27}:std::vector<I>{-53,-51,21,27};
 std::vector<std::pair<I,I>>quad=eight?std::vector<std::pair<I,I>>{{28,-53},{-20,-149},{4,-677},{-8,-233}}:std::vector<std::pair<I,I>>{{40,-801},{-40,-609},{-32,-873},{-24,-1057},{8,-1953},{-8,-993}};
 std::vector<int>qsign;for(auto q:quad)qsign.push_back(q.first>0?-1:1);
 bool norm=argc>1&&std::string(argv[1]).rfind("norm",0)==0;
 I Ylimit=D;
 if(norm){
  I M=std::stoll(std::string(argv[1]).substr(4)),P=M/4,delta=(M==152||M==688)?-839:2521;
  I N=-delta-P*P-2*P,Cnorm=P*P-delta,offset=delta==2521?4:3;
  D=210*M*M;B=2*N+M;C=0;E=N*(N+M)-offset*D;Ylimit=210*M;
  roots.clear();quad.clear();qsign.clear();
  for(I c:std::array<I,6>{0,1,-14,15,-20,21}){
   I z=N+c*M,r=z<=0?(I)std::sqrt((long double)-z):-1;
   if(r>=0&&r*r==-z){roots.push_back(-r);roots.push_back(r);}
   else{quad.emplace_back(0,z);qsign.push_back((1-2*c)*M>0?1:-1);}
  }
  quad.emplace_back(-2*P,Cnorm);qsign.push_back(1);
  quad.emplace_back(2*P,Cnorm);qsign.push_back(-1);
 }
 bool reverse=argc>2&&std::string(argv[2])=="reverse";
 if(reverse){C=-C;for(I&r:roots)r=-r;for(int k=0;k<(int)quad.size();++k)if(quad[k].first){quad[k].first=-quad[k].first;qsign[k]=-qsign[k];}}
 std::vector<I>squares={0};
 std::vector<I>primepowers;I unfactored=D;
 for(I p=2;p*p<=unfactored;++p)if(unfactored%p==0){I pe=1;while(unfactored%p==0){unfactored/=p;pe*=p;}primepowers.push_back(pe);}if(unfactored>1)primepowers.push_back(unfactored);
 for(I q:primepowers){
  std::set<I>sq;for(I u=1;u<q;++u)if(std::gcd(u,q)==1)sq.insert(u*u%q);
  I basis=(D/q)*inv(D/q%q,q);std::vector<I>next;
  for(I a:squares)for(I b:sq)next.push_back((a+mul(b,basis))%D);
  squares=std::move(next);
 }
 std::vector<I>fourths;for(I s:squares)fourths.push_back(mul(s,s));
 std::vector<Row>rows;
 for(I Y=0;Y<Ylimit;++Y){I Y2=mul(Y,Y);if(mod(mul(Y2,Y2)+B*Y2+C*Y+E)==0)rows.push_back(Row{Y});}
 double c2=(1.0/16)*log(7.0/6)+(1.0/32)*log(8.0/7);
 double c23=log(9.0/8)/4,c234=log(16.0/15)/2+log(5.0/4);
 double small=0;for(int J=8;J<80;++J)small+=roots.size()*pow(2.0,1-J)*log((J+1.0)/J);
 #pragma omp parallel for schedule(dynamic)
 for(int j=0;j<(int)rows.size();++j){
  Row&row=rows[j];I Y=row.Y,Y2=mul(Y,Y);
  for(int k=0;k<(int)quad.size();++k){auto [u,v]=quad[k];I q=mod(Y2+u*Y+v),K=mul(q,q);int sign=qsign[k];
   for(I s:squares)row.high[k]+=bad(mul(K,s),sign);
  }
  row.H=std::accumulate(row.high.begin(),row.high.end(),0LL);
  for(int k=0;k<(int)roots.size();++k){I r=roots[k],a=mod(Y-r),a2=mul(a,a),a3=mul(a2,a);
   I K2=mul(a2,mod(Y*Y+2*r*Y+3*r*r+B));I K3=mul(a3,mod(Y+3*r));I K4=mul(a2,a2);
   I E1=4*r*r*r+2*B*r+C,E2=6*r*r+B;
   int sign2=E1>0?1:-1,sign3=E2>0?1:-1,sign4=r>0?1:-1;
   int pairbad=K3==0?(sign3>0?2:0):(2*K3==D?(sign3<0?2:0):1);
   I b2=0,b24=0;
   for(int h=0;h<(int)squares.size();++h){bool b=bad(mul(K2,squares[h]),sign2);b2+=b;if(b&&bad(mul(K4,fourths[h]),sign4))++b24;}
   row.b2[k]=b2;row.b23twice[k]=pairbad*b2;row.b234twice[k]=pairbad*b24;
  }
  I B2=std::accumulate(row.b2.begin(),row.b2.end(),0LL),B23=std::accumulate(row.b23twice.begin(),row.b23twice.end(),0LL),B234=std::accumulate(row.b234twice.begin(),row.b234twice.end(),0LL);
  row.score=small+(c2*B2+(c23*B23+c234*B234)/2)/squares.size();
 }
 std::sort(rows.begin(),rows.end(),[](const Row&a,const Row&b){return a.H<b.H||(a.H==b.H&&(a.score<b.score||(a.score==b.score&&a.Y<b.Y)));});
 std::cout.precision(17);
 std::cout<<"{\n  \"modulus\": "<<D<<",\n  \"direction\": "<<(reverse?-1:1)<<",\n  \"residue_search_bound\": "<<Ylimit<<",\n  \"coarse_representatives_only\": "<<(norm?"true":"false")<<",\n  \"square_classes\": "<<squares.size()<<",\n  \"candidates\": "<<rows.size()<<",\n  \"root_order\": [";
 for(int k=0;k<(int)roots.size();++k)std::cout<<(k?", ":"")<<roots[k];
 std::cout<<"],\n  \"quadratic_order_uv\": [";for(int k=0;k<(int)quad.size();++k)std::cout<<(k?", ":"")<<'['<<quad[k].first<<','<<quad[k].second<<']';std::cout<<"],\n  \"rows\": [\n";
 for(int j=0;j<(int)rows.size();++j){auto&row=rows[j];std::cout<<"    {\"Y\": "<<row.Y<<", \"high_count\": "<<row.H<<", \"linear_cost\": "<<row.score<<", \"high_counts\": [";
  for(int k=0;k<(int)quad.size();++k)std::cout<<(k?", ":"")<<row.high[k];std::cout<<']';
  for(auto pair:{std::pair<const char*,std::array<I,8>>("b2",row.b2),{"b23_twice",row.b23twice},{"b234_twice",row.b234twice}}){std::cout<<", \""<<pair.first<<"\": [";for(int k=0;k<(int)roots.size();++k)std::cout<<(k?", ":"")<<pair.second[k];std::cout<<']';}
  std::cout<<'}'<<(j+1<(int)rows.size()?",":"")<<'\n';
 }
 std::cout<<"  ]\n}\n";
 std::cerr<<"Candidates "<<rows.size()<<"; square classes "<<squares.size()<<"; best Y "<<rows.front().Y<<"; high coefficient "<<(double)rows.front().H/squares.size()<<"; linear cost "<<rows.front().score<<'\n';
}
