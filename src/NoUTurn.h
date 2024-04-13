//
// Created by Marc Suchard on 2019-12-03.
//

#ifndef NUTS_HPP
#define NUTS_HPP

//#pragma clang diagnostic push
//#pragma ide diagnostic ignored "OCUnusedMacroInspection" // Turn off warning for TBB_PREVIEW_GLOBAL_CONTROL

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
#include "MemoryManagement.h"
#include "Simd.h"
#include "ZigZag.h"
#include "NutsTreeState.h"
#include "UniformGenerator.h"

namespace nuts {
    using DblSpan = tcb::span<double>;
    using UniPtrTreeState = std::unique_ptr<TreeState>;
    using SharedPtrTreeState = std::shared_ptr<TreeState>;

    class NoUTurn {
    public:
        NoUTurn(double logProbErrorTol,
                int maxHeight,
                int seed,
                bool randomFlg,
                double stepSize,
                std::shared_ptr<zz::ZigZag<zz::DoubleSseTypeInfo>> zigzag) : logProbErrorTol(logProbErrorTol),
                                                                             maxHeight(maxHeight),
                                                                             stepSize(stepSize),
                                                                             zzEngine(*zigzag),
                                                                             uniGenerator(UniformGenerator(seed,
                                                                                                           randomFlg)) {
        }

        ~NoUTurn() = default;

        // Utility function for debugging purpose
        template<typename T>
        void printDblSpan(T &span) {
           for (auto e: span) std::cout << e << ' ';
        }


        std::vector<double> takeOneStep(DblSpan initialPosition, DblSpan initialMomentum) {

            
            const double initialJointDensity = zzEngine.getLogPDFnoDet(initialPosition, initialMomentum);
            double logSliceU = log(uniGenerator.getUniform()) + initialJointDensity;

            std::unique_ptr<Eigen::VectorXd> gPtr = zzEngine.getLogdGradient(initialPosition);
            DblSpan gradient(*gPtr);

            TreeState *newState = new TreeState(initialPosition, initialMomentum, gradient, 1, true,
                                                0, 0, uniGenerator);
            SharedPtrTreeState trajectoryTree = std::shared_ptr<TreeState>(newState);

            int height = 0;

            while (trajectoryTree->flagContinue) {
                updateTrajectoryTree(trajectoryTree, height, logSliceU, initialJointDensity);
                height++;
                if (height > maxHeight) {
                    trajectoryTree->flagContinue = false;
                }
            }
            DblSpan sampleSpan = trajectoryTree->getSample();
            std::vector<double> sample(sampleSpan.begin(), sampleSpan.end());
            return sample;
        }

        void updateTrajectoryTree(SharedPtrTreeState trajectoryTree,
                                     int depth,
                                     double logSliceU,
                                     double initialJointDensity) {
            int direction = (uniGenerator.getUniform() < 0.5) ? -1 : 1;
            UniPtrTreeState nextTrajectoryTree = buildTree(
                    trajectoryTree->getPosition(direction), trajectoryTree->getMomentum(direction),
                    trajectoryTree->getGradient(direction),
                    direction, logSliceU, depth, stepSize, initialJointDensity);

            if ((*nextTrajectoryTree).flagContinue) {
                bool swapSampling = true;
                trajectoryTree->mergeNextTree((*nextTrajectoryTree), direction, swapSampling);
            } else {
                trajectoryTree->flagContinue = false;
            }
        }

        UniPtrTreeState buildTree(DblSpan position, DblSpan momentum, DblSpan gradient, int direction,
                                  double logSliceU, int height, double stepSize, double initialJointDensity) {
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

            zzEngine.reversiblePositionMomentumUpdate(position, momentum, gradient, direction, stepSize);
            
            double logJointProbAfter = zzEngine.getLogPDFnoDet(position, momentum);

            const int numNodes = (logSliceU <= logJointProbAfter ? 1 : 0);

            const bool flagContinue = (logSliceU < logProbErrorTol + logJointProbAfter);

            // Values for dual-averaging
            const double acceptProb = std::min(1.0, exp(logJointProbAfter - initialJointDensity));
            const int numAcceptProbStates = 1;

            //hmcProvider.setParameter(inPosition);
            TreeState *newState = new TreeState(position, momentum, gradient, numNodes,
                                                flagContinue, acceptProb,
                                                numAcceptProbStates, uniGenerator);
            return UniPtrTreeState(newState);
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

                if ((*nextSubtree).flagContinue) {
                    bool swapSampling = false;
                    (*subtree).mergeNextTree((*nextSubtree), direction, swapSampling);
                } else {
                    (*subtree).flagContinue = false;
                }
            }
            return subtree;
        }

        double logProbErrorTol = 100.0;
        const int maxHeight = 10;
        double stepSize;
        zz::ZigZag<zz::DoubleSseTypeInfo> zzEngine;
        UniformGenerator uniGenerator;
    };

    std::unique_ptr<nuts::NoUTurn> dispatchNuts(
            double logProbErrorTol,
            int maxHeight,
            int seed,
            bool randomFlg,
            double stepSize,
            std::shared_ptr<zz::ZigZag<zz::DoubleSseTypeInfo>> ptr) {
        return zz::make_unique<nuts::NoUTurn>(logProbErrorTol, maxHeight, seed, randomFlg, stepSize, ptr);
    }
}

//#pragma clang diagnostic pop

#endif //ZIG_ZAG_ZIGZAG_HPP
