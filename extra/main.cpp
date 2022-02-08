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
    double time = 100;
    std::vector<double> mask = {1,1};
    std::vector<double> observed = {1,1};
    std::vector<double> parameterSign = {1,1};

    std::vector<double> position = {0.1,0.1};
    std::vector<double> velocity = {-1,1};
    std::vector<double> action = {-1,1};
    std::vector<double> gradient = {-0.1,-0.1};
    std::vector<double> momentum = {-0.3980153048394741, 2.2483327571774736};

    long flags = 128L;
    long info = 1L;
    long seed = 666L;

    std::vector<double> mean={0, 0};
    std::vector<double> precision={1,0,0,1};

//    std::unique_ptr<zz::ZigZag<zz::DoubleSseTypeInfo>> zz = zz::dispatch(dimension, mask.data(), observed.data(), parameterSign.data(), flags, info, seed,
//                         zz::DblSpan(mean), zz::DblSpan(precision));

    auto zz = zz::make_unique<zz::ZigZag<zz::DoubleSseTypeInfo>>(
            dimension, mask.data(), observed.data(), parameterSign.data(), flags, info, seed, mean, precision);
    std::shared_ptr<zz::ZigZag<zz::DoubleSseTypeInfo>> shared = std::move(zz);

    std::unique_ptr<nuts::NoUTurn> nuts = nuts::dispatchNuts(100, 100, 10, 666, 15.182, shared);
    cout << "nuts one step starts:" << endl;

//    nuts->testOneStep(position, momentum, gradient);

    nuts->takeOneStep(position, momentum, gradient);
//    shared->operate(zz::DblSpan(position),
//                zz::DblSpan(velocity),
//                zz::DblSpan(action),
//                zz::DblSpan(gradient),
//                zz::DblSpan(momentum),
//                time);
    cout << "position after:" << endl;
    for (int i = 0; i < position.size(); ++i) {
        cout << position[i] << endl;
    }

    return 666;
}
