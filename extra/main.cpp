//
// Created by zhenyuzhang on 2/4/22.
//

#include <vector>

#include "../src/AbstractZigZag.h"
#include "../src/ZigZag.h"
#include "../src/NoUTurn.h"
#include "../src/UniformGenerator.h"
using namespace std;
int main(int argCount, char** args) {
    int dimension = 4;
    std::vector<double> mask = {1, 1, 1, 1};
    std::vector<double> observed = {1, 1, 1, 1};
    std::vector<double> parameterSign = {1, 1, 1, 1};

    std::vector<double> lb = {0.0, 0.0, 0.0, 0.0};
    std::vector<double> ub = {INFINITY, INFINITY, INFINITY, INFINITY};
    if (INFINITY == INFINITY){
        cerr << " inf = inf " ;
    }
    std::vector<double> momentum = {1.138124890, -0.485223884 , 0.009229834 ,-0.551786161};
    long flags = 128L;
    long info = 1L;
    long seed = 666L;

    std::vector<double> position={0.1, 0.1, 0.1, 0.1};
//    std::fill (position.begin(),position.end(),0.1);

    std::vector<double> mean={0, 0, 0, 0};
    std::vector<double> precision={1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};

    auto zz = zz::make_unique<zz::ZigZag<zz::DoubleSseTypeInfo>>(
            dimension, mask.data(), lb.data(), ub.data(), flags, info, seed);
    std::shared_ptr<zz::ZigZag<zz::DoubleSseTypeInfo>> shared = std::move(zz);
    shared->setMean(mean);
    shared->setPrecision(precision);
    std::unique_ptr<nuts::NoUTurn> nuts = nuts::dispatchNuts(100, 100, 666, true, 0.1, shared);

    //' Set mean for MTN
//'
//' @param sexp pointer to zigzag object
//' @param mean a numerica vector containing the MTN mean
//' @export
// [[Rcpp::export(setMean)]]

    for (int n = 0; n < 1; ++n) {
        auto res = shared -> operate(zz::DblSpan(position), zz::DblSpan(momentum), 1.414214);
        for (int i = 0; i < position.size(); ++i) {
            cerr << position[i] << " ";
        }
        cerr << "\n";
    }
    return 0;
}
