//
// Created by Zhenyu Zhang on 2/1/22.
//

#ifndef NOUTURN_HPP_UNIFORMGENERATOR_H
#define NOUTURN_HPP_UNIFORMGENERATOR_H


class UniformGenerator{
public:
    UniformGenerator(int seed){
        generator = std::mt19937(seed);
        distribution = std::uniform_real_distribution<double>(0, 1);
        //std::cerr << "uniform generator constructed" << '\n' << std::endl;
    }
    double getUniform(){
        return distribution(generator);
    }
    std::mt19937 generator;
    std::uniform_real_distribution<double> distribution;
};


#endif //NOUTURN_HPP_UNIFORMGENERATOR_H
