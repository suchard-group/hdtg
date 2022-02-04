//
// Created by zhenyuzhang on 2/4/22.
//

#include <vector>

#include "../src/AbstractZigZag.hpp"
#include "../src/ZigZag.hpp"
#include "../src/NoUTurn.hpp"
using namespace std;
int main(int argCount, char** args) {

    int dimension = 2;
    std::vector<double> mask(dimension);
    std::vector<double> observed(dimension);
    std::vector<double> parameterSign(dimension);

    long flags = 0L;
    long info = 0L;
    long seed = 666L;

    std::vector<double> mean(dimension);
    std::vector<double> precision(dimension * dimension);

    cout << "hello" << endl;
//    zz::DblSpan tmp(mean);


    auto zz =
            zz::dispatch(dimension, mask.data(), observed.data(), parameterSign.data(), flags, info, seed,
                         zz::DblSpan(mean), zz::DblSpan(precision));

    return 1;
}
