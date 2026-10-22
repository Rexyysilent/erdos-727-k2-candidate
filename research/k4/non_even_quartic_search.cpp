// Exact rational-factorization search for centered quartics with a linear term.
// Build: g++ -O3 -std=c++17 non_even_quartic_search.cpp -o search
// Run: ./search 60
// A hit supplies four ordered u values and rational C such that
// ((u^2+B)^2-C^2/u^2)/4 is an arithmetic progression for some rational B.
// This proves only rational factorization. Integrality, positive n, carry
// conditions, and infinitude must be checked separately.
// The default bound 60 tested 5,851,620 orderings and found 5,038 hits.
// Bound 120 tested 98,574,840 orderings and found 25,196 hits.
// Signed 128-bit arithmetic is ample through bound 300. Larger inputs are
// rejected rather than silently risking integer overflow.
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <string>
using I=__int128_t;using U=__uint128_t;
I gcd(I a,I b){if(a<0)a=-a;if(b<0)b=-b;while(b){I c=a%b;a=b;b=c;}return a;}
std::string out(I v){if(!v)return"0";bool s=v<0;if(s)v=-v;std::string r;while(v){r.push_back('0'+v%10);v/=10;}if(s)r.push_back('-');std::reverse(r.begin(),r.end());return r;}
bool square(I z,I&r){if(z<0)return false;r=(I)sqrtl((long double)z);while(r*r<z)++r;while(r*r>z)--r;return r*r==z;}
int main(int argc,char**argv){int bound=argc>1?std::stoi(argv[1]):60;if(bound<4||bound>300){std::cerr<<"bound must be between 4 and 300\n";return 2;}long long tested=0,hit=0;
for(int x=1;x<=bound;x++)for(int y=x+1;y<=bound;y++)for(int z=y+1;z<=bound;z++)for(int w=z+1;w<=bound;w++){
 std::array<int,4>u{x,y,z,w};do{if(u[0]>u[3])continue;tested++;I a=I(u[0])*u[0],b=I(u[1])*u[1],c=I(u[2])*u[2],d=I(u[3])*u[3];
 I L1=a-2*b+c,L2=b-2*c+d,M1=a*a-2*b*b+c*c,M2=b*b-2*c*c+d*d;
 I N1=b*c-2*a*c+a*b,N2=c*d-2*b*d+b*c;
 I num=(M1*L2-M2*L1)*a*b*c*d,den=N1*d*L2-N2*a*L1;
 if(!num||!den||(num>0)!=(den>0))continue;if(den<0){num=-num;den=-den;}I g=gcd(num,den);num/=g;den/=g;I sn,sd;
 if(square(num,sn)&&square(den,sd)){hit++;std::cout<<u[0]<<' '<<u[1]<<' '<<u[2]<<' '<<u[3]<<" C="<<out(sn)<<'/'<<out(sd)<<'\n';}
 }while(std::next_permutation(u.begin(),u.end()));
}
std::cerr<<"tested="<<tested<<" hits="<<hit<<'\n';
}
