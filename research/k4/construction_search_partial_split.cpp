#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <cstdlib>
using I = __int128_t;
using U = __uint128_t;
bool square(U z) {
    U s = (U) std::sqrt((long double) z);
    while (s*s < z) ++s;
    while (s*s > z) --s;
    return s*s == z;
}
bool same_square_class(I a, I b) {
    return a*b>0 && square((U)(a*b));
}
int main(int argc, char** argv) {
    int limit = argc>1 ? std::atoi(argv[1]) : 2000;
    int found=0;
    for (long long p=1;p<=limit;++p) {
        for (long long q=1;q<=limit;++q) {
            if (std::gcd(p,q)!=1) continue;
            I num=(I)839*q*q-(I)p*p-(I)2*p*q;
            I den=(I)4*p*q;
            I v[3]={num,num-14*den,num+21*den};
            int mask=0;
            if(same_square_class(v[0],v[1])) mask|=1;
            if(same_square_class(v[0],v[2])) mask|=2;
            if(same_square_class(v[1],v[2])) mask|=4;
            if(mask) {
                std::cout << p << ' ' << q << ' ' << (long long)num
                          << ' ' << (long long)den << ' ' <<mask<< '\n';
                ++found;
            }
        }
    }
    std::cerr << "norm parameters bounded by " << limit << "; matches " <<found<<'\n';
}
