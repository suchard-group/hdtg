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

    class NoUTurn {
    public:
        NoUTurn(double logProbErrorTol,
                int findMax,
                int maxHeight,
                int seed) : logProbErrorTol(logProbErrorTol),
                                 findMax(findMax),
                                 maxHeight(maxHeight),
                                 uniSeed(seed),
                                 uniGenerator(UniformGenerator(seed)){
            std::cerr << "nuts constructed" << '\n' << std::endl;
        }

        ~NoUTurn() = default;

        DblSpan takeOneStep(DblSpan initialPosition, DblSpan initialMomentum, DblSpan gradient) {

            DblSpan endPosition = initialPosition;

            const double initialJointDensity = zigzag.getJointProbability(initialPosition, initialMomentum); //todo to check

            double logSliceU = log(uniGenerator.getUniform()) + initialJointDensity;

            TreeState trajectoryTree = TreeState(initialPosition, initialMomentum, gradient, 1, true,
                                                 0, 0, ++uniSeed); //todo to check

            int height = 0;

            while (trajectoryTree.flagContinue) {

                DblSpan tmp = updateTrajectoryTree(trajectoryTree, height, logSliceU, initialJointDensity);
                if (!tmp.empty()) {
                    endPosition = tmp;
                }

                height++;

                if (height > options.maxHeight) {
                    trajectoryTree.flagContinue = false;
                }
            }
            return endPosition;
        }

        DblSpan updateTrajectoryTree(TreeState trajectoryTree,
                                     int depth,
                                     double logSliceU,
                                     double initialJointDensity) {
            DblSpan endPosition;

            int direction = (uniGenerator.getUniform() < 0.5) ? -1 : 1;
            TreeState nextTrajectoryTree = buildTree(
                    trajectoryTree.getPosition(direction), trajectoryTree.getMomentum(direction),
                    trajectoryTree.getGradient(direction),
                    direction, logSliceU, depth, stepSizeInformation.getStepSize(), initialJointDensity);

            if (nextTrajectoryTree.flagContinue) {

                const double acceptProb = (double) nextTrajectoryTree.numNodes / (double) trajectoryTree.numNodes;
                if (uniGenerator.getUniform() < acceptProb) {
                    endPosition = nextTrajectoryTree.getSample();
                }
            }

            trajectoryTree.mergeNextTree(nextTrajectoryTree, direction);

            return endPosition;
        }

        TreeState buildTree(DblSpan position, DblSpan momentum, DblSpan gradient, int direction,
                            double logSliceU, int height, double stepSize, double initialJointDensity) {

            if (height == 0) {
                return buildBaseCase(position, momentum, gradient, direction, logSliceU, stepSize, initialJointDensity);
            } else {
                return buildRecursiveCase(position, momentum, gradient, direction, logSliceU, height, stepSize,
                                          initialJointDensity);
            }
        }

        static TreeState buildBaseCase(DblSpan inPosition, DblSpan inMomentum, DblSpan inGradient, int direction,
                                double logSliceU, double stepSize, double initialJointDensity) {
            // Make deep copy of position and momentum
            std::vector<double> positionVec;
            std::vector<double> momentumVec;
            std::vector<double> gradientVec;

            position.assign(inPosition.begin(), inPosition.end());
            momentum.assign(inMomentum.begin(), inMomentum.end());
            gradient.assign(inGradient.begin(), inGradient.end());

            DblSpan position{positionVec};
            DblSpan momentum{momentumVec};
            DblSpan gradient{gradientVec};
            //hmcProvider.setParameter(position.getBuffer());
            // "one reversibleHMC integral
            ZigZag.reversiblePositionMomentumUpdate(position, momentum, gradient, direction, stepSize);

            double logJointProbAfter = ZigZag.getJointProbability(momentum);

            const int numNodes = (logSliceU <= logJointProbAfter ? 1 : 0);

            const bool flagContinue = (logSliceU < options.logProbErrorTol + logJointProbAfter);

            // Values for dual-averaging
            const double acceptProb = Math.min(1.0, Math.exp(logJointProbAfter - initialJointDensity));
            const int numAcceptProbStates = 1;

            //hmcProvider.setParameter(inPosition);

            return new TreeState(position, momentum, gradient, numNodes,
                                 flagContinue, acceptProb,
                                 numAcceptProbStates);
        }

        TreeState buildRecursiveCase(DblSpan inPosition, DblSpan inMomentum, DblSpan gradient, int direction,
                                     double logSliceU, int height, double stepSize, double initialJointDensity) {

            TreeState subtree = buildTree(inPosition, inMomentum, gradient, direction, logSliceU,
                                          height - 1, // Recursion
                                          stepSize, initialJointDensity);

            if (subtree.flagContinue) {

                TreeState nextSubtree = buildTree(subtree.getPosition(direction), subtree.getMomentum(direction),
                                                  subtree.getGradient(direction), direction,
                                                  logSliceU, height - 1, stepSizeInformation.getStepSize(),
                                                  initialJointDensity);

                subtree.mergeNextTree(nextSubtree, direction);
            }
            return subtree;
        }

        double logProbErrorTol = 100.0;
        const int findMax = 100;
        const int maxHeight = 100;
        zz::ZigZag<zz::DoubleSseTypeInfo> zigzag;
        int uniSeed;
        UniformGenerator uniGenerator;
    };

}

#pragma clang diagnostic pop

#endif //ZIG_ZAG_ZIGZAG_HPP
