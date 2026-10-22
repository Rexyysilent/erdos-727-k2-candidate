// Search monic depressed cubics with four completely split translates in AP.
// A split x^3-A*x+B has roots r,s,-r-s, hence A=r^2+r*s+s^2,
// B=r*s*(r+s). Symmetry gives -B as well. Scaling by the AP difference
// would turn any four-term AP of such B values into consecutive translates.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>
using Pair=std::pair<long long,long long>;
int main(int argc,char** argv){
    long long maxA=argc>1?std::atoll(argv[1]):20000000;
    std::vector<Pair> values;
    for(long long r=1;r*r<=maxA;++r)
      for(long long s=0;s<=r;++s){
        long long A=r*r+r*s+s*s;
        if(A>maxA)break;
        values.emplace_back(A,r*s*(r+s));
      }
    std::sort(values.begin(),values.end());
    std::cerr<<"pairs "<<values.size()<<"; max A "<<maxA<<'\n';
    size_t begin=0,groups=0;long long hits=0;
    while(begin<values.size()){
      size_t end=begin+1;
      while(end<values.size()&&values[end].first==values[begin].first)++end;
      if(end-begin>=2){
        ++groups;
        std::vector<long long> bs;
        for(size_t j=begin;j<end;++j){
          bs.push_back(values[j].second);bs.push_back(-values[j].second);
        }
        std::sort(bs.begin(),bs.end());
        bs.erase(std::unique(bs.begin(),bs.end()),bs.end());
        for(size_t i=0;i<bs.size();++i)
          for(size_t j=i+1;j<bs.size();++j){
            long long a=bs[i],d=bs[j]-a;
            if(a+3*d>bs.back())break;
            if(std::binary_search(bs.begin(),bs.end(),a+2*d)&&
               std::binary_search(bs.begin(),bs.end(),a+3*d)){
              ++hits;
              std::cout<<values[begin].first<<' '<<a<<' '<<d<<'\n';
            }
          }
      }
      begin=end;
    }
    std::cerr<<"multiple-representation groups "<<groups<<"; AP hits "<<hits<<'\n';
}
