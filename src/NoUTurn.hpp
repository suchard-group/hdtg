//
// Created by Marc Suchard on 2019-12-03.
//

#ifndef NUTS_HPP
#define NUTS_HPP

#pragma clang diagnostic push
#pragma ide diagnostic ignored "OCUnusedMacroInspection" // Turn off warning for TBB_PREVIEW_GLOBAL_CONTROL

#include <vector>
#include <cmath>

#define TBB_PREVIEW_GLOBAL_CONTROL 1

#define TIMING

#ifdef TIMING

#include <map>
#include <iomanip>
#include "Timing.h"

#endif // TIMING

#include "threefry.h"
//#include "dr_evomodel_operators_NativeZigZag.h"
#include "MemoryManagement.hpp"
#include "Simd.hpp"
#include "ZigZag.hpp"
#include "NutsTreeState.hpp"
#include "UniformGenerator.h"

namespace nuts {
    using DblSpan = tcb::span<double>;
    using UniPtrTreeState = std::unique_ptr<TreeState>;
    using SharedPtrTreeState = std::shared_ptr<TreeState>;

    class NoUTurn {
    public:
        NoUTurn(double logProbErrorTol,
                int findMax,
                int maxHeight,
                int seed,
                double stepSize,
                std::shared_ptr<zz::ZigZag<zz::DoubleSseTypeInfo>> zigzag) : logProbErrorTol(logProbErrorTol),
                                                                             findMax(findMax),
                                                                             maxHeight(maxHeight),
                                                                             uniSeed(seed),
                                                                             uniGenerator(UniformGenerator(seed)),
                                                                             stepSize(stepSize),
                                                                             zzEngine(*zigzag) {
            std::cerr << "nuts constructed" << '\n' << std::endl;
        }

        ~NoUTurn() = default;

        void printDblSpan(zz::DblSpan &span) {
            for (auto e: span) std::cout << e << ' ';
            std::cout << "\n";
        }

        DblSpan takeOneStep(DblSpan initialPosition, DblSpan initialMomentum, DblSpan gradient) {
            DblSpan endPosition = initialPosition;

            const double initialJointDensity = zzEngine.getJointProbability(initialPosition,
                                                                            initialMomentum);
            //double tmp = log(uniGenerator.getUniform());
            //double tmp = -1.3230210224799723;
            double logSliceU = log(uniGenerator.getUniform()) + initialJointDensity;
            //std::cerr << "logSlice u is " << logSliceU << std::endl;
            //printDblSpan(initialMomentum);

            TreeState *newState = new TreeState(initialPosition, initialMomentum, gradient, 1, true,
                                                0, 0, ++uniSeed);
            SharedPtrTreeState trajectoryTree = std::move(zz::make_unique<TreeState>(*newState));

//        TreeState trajectoryTree = TreeState(initialPosition, initialMomentum, gradient, 1, true,
//                                             0, 0, ++uniSeed); //todo to check

            int height = 0;

            while (trajectoryTree->flagContinue) {
                //std::cerr << "***************appear once" << std::endl;
                DblSpan tmp = updateTrajectoryTree(trajectoryTree, height, logSliceU, initialJointDensity);
                if (!tmp.empty()) {
                    endPosition = tmp;
                    printDblSpan(tmp);
                }

                height++;

                if (height > maxHeight) {
                    trajectoryTree->flagContinue = false;
                }
            }
            return endPosition;
        }

        DblSpan updateTrajectoryTree(SharedPtrTreeState trajectoryTree,
                                     int depth,
                                     double logSliceU,
                                     double initialJointDensity) {
            DblSpan endPosition;
            //zz::DblSpan tmpspan = trajectoryTree->getMomentum(-1);
            //printDblSpan(tmpspan);
            int direction = (uniGenerator.getUniform() < 0.5) ? -1 : 1;
//        int direction = 1;
//        std::cout << "fixed direction:" << direction << std::endl;

            UniPtrTreeState nextTrajectoryTree = buildTree(
                    trajectoryTree->getPosition(direction), trajectoryTree->getMomentum(direction),
                    trajectoryTree->getGradient(direction),
                    direction, logSliceU, depth, stepSize, initialJointDensity);

            if ((*nextTrajectoryTree).flagContinue) {

                const double acceptProb = (double) (*nextTrajectoryTree).numNodes / (double) trajectoryTree->numNodes;
                if (uniGenerator.getUniform() < acceptProb) {
                    endPosition = (*nextTrajectoryTree).getSample();
                }
            }

            trajectoryTree->mergeNextTree((*nextTrajectoryTree), direction);

            return endPosition;
        }

        UniPtrTreeState buildTree(DblSpan position, DblSpan momentum, DblSpan gradient, int direction,
                                  double logSliceU, int height, double stepSize, double initialJointDensity) {
            //std::cerr << "height is" << height << std::endl;
            if (height == 0) {
                return buildBaseCase(position, momentum, gradient, direction, logSliceU, stepSize, initialJointDensity);
            } else {
                return buildRecursiveCase(position, momentum, gradient, direction, logSliceU, height, stepSize,
                                          initialJointDensity);
            }
        }

        UniPtrTreeState buildBaseCase(DblSpan inPosition, DblSpan inMomentum, DblSpan inGradient, int direction,
                                      double logSliceU, double stepSize, double initialJointDensity) {
            // Make deep copy of position and momentum
            std::vector<double> positionVec;
            std::vector<double> momentumVec;
            std::vector<double> gradientVec;

            positionVec.assign(inPosition.begin(), inPosition.end());
            momentumVec.assign(inMomentum.begin(), inMomentum.end());
            gradientVec.assign(inGradient.begin(), inGradient.end());

            DblSpan position{positionVec};
            DblSpan momentum{momentumVec};
            DblSpan gradient{gradientVec};

            // "one reversibleHMC integral
//            std::cout << "position before reversible update:" << std::endl;
//            printDblSpan(inPosition);

            zzEngine.reversiblePositionMomentumUpdate(position, momentum, gradient, direction, stepSize);

//            std::cout << "position after reversible update:" << std::endl;
//            printDblSpan(position);
            double logJointProbAfter = zzEngine.getJointProbability(position, momentum);

            const int numNodes = (logSliceU <= logJointProbAfter ? 1 : 0);

            const bool flagContinue = (logSliceU < logProbErrorTol + logJointProbAfter);

            // Values for dual-averaging
            const double acceptProb = std::min(1.0, exp(logJointProbAfter - initialJointDensity));
            const int numAcceptProbStates = 1;

            //hmcProvider.setParameter(inPosition);
            TreeState *newState = new TreeState(position, momentum, gradient, numNodes,
                                                flagContinue, acceptProb,
                                                numAcceptProbStates, ++uniSeed);
            return zz::make_unique<TreeState>(*newState);
        }


        UniPtrTreeState buildRecursiveCase(DblSpan inPosition, DblSpan inMomentum, DblSpan gradient, int direction,
                                           double logSliceU, int height, double stepSize, double initialJointDensity) {

            UniPtrTreeState subtree = buildTree(inPosition, inMomentum, gradient, direction, logSliceU,
                                                height - 1, // Recursion
                                                stepSize, initialJointDensity);

            if ((*subtree).flagContinue) {

                UniPtrTreeState nextSubtree = buildTree((*subtree).getPosition(direction),
                                                        (*subtree).getMomentum(direction),
                                                        (*subtree).getGradient(direction), direction,
                                                        logSliceU, height - 1, stepSize,
                                                        initialJointDensity);

                (*subtree).mergeNextTree((*nextSubtree), direction);
            }
            return subtree;
        }

        double logProbErrorTol = 100.0;
        const int findMax = 100;
        const int maxHeight = 10;
        double stepSize;
        zz::ZigZag<zz::DoubleSseTypeInfo> zzEngine;
        int uniSeed;
        UniformGenerator uniGenerator;
    };

    std::unique_ptr<nuts::NoUTurn> dispatchNuts(
            double logProbErrorTol,
            int findMax,
            int maxHeight,
            int seed,
            double stepSize,
            std::shared_ptr<zz::ZigZag<zz::DoubleSseTypeInfo>> ptr) {
        std::cerr << "Factory: SSE" << std::endl;
        //NoUTurn *newNuts = new NoUTurn(logProbErrorTol, findMax, maxHeight, seed, stepSize, ptr);
        return zz::make_unique<nuts::NoUTurn>(logProbErrorTol, findMax, maxHeight, seed, stepSize, ptr);
    }
}

#pragma clang diagnostic pop

#endif //ZIG_ZAG_ZIGZAG_HPP
