//
// Created by Zhenyu Zhang on 2/1/22.
//
#include <random>

#ifndef NOUTURN_HPP_UNIFORMGENERATOR_H
#define NOUTURN_HPP_UNIFORMGENERATOR_H


class UniformGenerator {
public:
    UniformGenerator(int seed, bool randomFlg) : randomFlg(randomFlg) {
        count = 1;
        generator = std::mt19937(seed);
        distribution = std::uniform_real_distribution<double>(0, 1);
//        std::cerr << "uniform generator constructed" << '\n' << std::endl;
    }

    double getUniform() {
        double tmp;
        if (randomFlg) {
            tmp = distribution(generator);
        } else {
            if (count % 10 == 0) {
                ++count;
            }
            tmp = count % 10 / 10.;
            ++count;
        }
        //std::cerr << tmp << std::endl;
        return tmp;
    }

    std::mt19937 generator;
    std::uniform_real_distribution<double> distribution;
    bool randomFlg;
    int count;
};


#endif //NOUTURN_HPP_UNIFORMGENERATOR_H
